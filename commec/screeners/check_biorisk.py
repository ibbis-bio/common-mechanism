#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks the output from hmmscan and prints to screen the results

Usage:
 python check_biorisk.py -i INPUT.biorisk.hmmscan -d databases/biorisk_db/
"""

import logging
import os
import sys
import argparse
import pandas as pd
from commec.config.constants import BIORISK_E_VALUE_THRESHOLD
from commec.tools.hmmer import (
    readhmmer,
    remove_overlaps,
    HmmerHandler,
)
from commec.config.result import (
    ScreenResult,
    HitResult,
    ScreenStep,
    Recommendation,
    HitRecommendationContainer,
    MatchRange,
    compare
)


def _guess_domain(search_string : str) -> str:
    """
    Given a string description, try to determine
    which domain of life this has come from. Temporary work around
    until we can retrieve this data directly from biorisk outputs.
    """
    def contains(search_string : str, search_terms):
        for token in search_terms:
            if search_string.find(token) == -1:
                continue
            return True
        return False

    search_token = search_string.lower()
    if contains(search_token, ["vir", "capsid", "RNA Polymerase"]):
        return "Virus"
    if contains(search_token, ["cillus","bact","coccus","phila","ella","cocci","coli"]):
        return "Bacteria"
    if contains(search_token, ["eukary","nucleus","sona","odium","myces"]):
        return "Eukaryote"
    return "not assigned"

def parse_biorisk_hits(search_handler : HmmerHandler, result : ScreenResult):
    """
    Parse the outputs of a biorisk search and update screen results.

    Inputs:
        search_handler: The handler which has performed a biorisk search.
        result : ScreenResult to be updated with information from biorisk search.
    """

    # Read annotations CSV
    annotations_csv = os.path.join(search_handler.db_directory,"biorisk_annotations.csv")
    if not os.path.exists(annotations_csv):
        logging.error("\t...biorisk_annotations.csv does not exist\n %s", annotations_csv)
        return

    # Confirm that search succeeded
    if not search_handler.check_output():
        logging.error("\t...biorisk screening output doesn't exist:\n %s", search_handler.out_file)
        result.set_recommendation(ScreenStep.BIORISK, Recommendation.ERROR)
        return
    if search_handler.is_empty(search_handler.out_file):
        logging.error("\t...ERROR: biorisk search results empty\n")
        result.set_recommendation(ScreenStep.BIORISK, Recommendation.ERROR)
        return

    result.set_recommendation(ScreenStep.BIORISK, Recommendation.PASS)

    # If there are no hits, everything passes, and we can just return
    if not search_handler.has_hits(search_handler.out_file):
        return

    # Read and parse outputs
    hmmer : pd.DataFrame = readhmmer(search_handler.out_file)
    keep1 = [i for i, x in enumerate(hmmer['E-value']) if x < BIORISK_E_VALUE_THRESHOLD]
    hmmer = hmmer.iloc[keep1,:]

    hmmer['nt len'] = _get_nt_len_from_query(hmmer, result)
    hmmer = remove_overlaps(hmmer)

    # Read annotations
    lookup : pd.DataFrame = pd.read_csv(annotations_csv)
    lookup.fillna(False, inplace=True)

    # Append description, and must_flag columns from annotations:
    hmmer['description'] = ''
    hmmer['Must flag'] = False
    hmmer = hmmer.reset_index(drop=True)
    for model in range(hmmer.shape[0]):
        name_index = [i for i, x in enumerate([lookup['ID'] == hmmer['target name'][model]][0]) if x]
        hmmer.loc[model, 'description'] = lookup.iloc[name_index[0], 1]
        hmmer.loc[model, 'Must flag'] = lookup.iloc[name_index[0], 2]

    # Update the results to capture the biorisk outputs
    queries_with_hits = hmmer['query name'].unique()
    for query_name in queries_with_hits:
        biorisk_overall : Recommendation = Recommendation.PASS

        query_result = result.get_query(query_name)
        if not query_result:
            logging.error("Could not find query matching query name in hmmscan output: [%s]", query_name)
            continue

        # Grab a list of unique queries, and targets for iteration.
        unique_query_data : pd.DataFrame = hmmer[hmmer['query name'] == query_name]
        unique_targets = unique_query_data['target name'].unique()

        for target_name in unique_targets:
            unique_target_data : pd.DataFrame = unique_query_data[unique_query_data['target name'] == target_name]
            target_description = ", ".join(set(unique_target_data['description'])) # First should be unique.
            must_flag = unique_target_data['Must flag'].iloc[0] # First should be unique.
            match_ranges = []
            for _, region in unique_target_data.iterrows():
                match_range = MatchRange(
                    float(region['E-value']),
                    int(region['hmm from']), int(region['hmm to']),
                    int(region['ali from']), int(region['ali to'])
                )
                match_ranges.append(match_range)

            target_recommendation = Recommendation.FLAG if must_flag > 0 else Recommendation.WARN

            biorisk_overall = compare(target_recommendation, biorisk_overall)

            hit_data : HitResult = query_result.get_hit(target_name)
            if hit_data:
                hit_data.ranges.extend(match_ranges)
                continue

            regulation_str = "Regulated Gene" if must_flag else "Virulence Factor"
            domain = _guess_domain(""+str(target_name)+target_description)

            new_hit : HitResult = HitResult(
                HitRecommendationContainer(
                    target_recommendation,
                    ScreenStep.BIORISK
                ),
                target_name,
                target_description,
                match_ranges,
                {"domain" : [domain],"regulated":[regulation_str]},
            )
            query_result.hits[target_name] = new_hit

        query_result.set_recommendation(ScreenStep.BIORISK, biorisk_overall)

def _get_nt_len_from_query(hmmer: pd.DataFrame, result: ScreenResult) -> pd.Series:
    """
    Get a Series of nucleotide lengths of original (i.e. untranslated) queries. If no matching
    query is found, use qlen * 3 (may lead to off-by-1 or off-by-2 errors depending on frames)
    """
    matched_queries = hmmer['query name'].apply(result.get_query)
    # Get query length where matches were found
    lengths = matched_queries.apply(lambda x: getattr(x, 'query_length', None))
    # If no match, use qlen * 3
    return lengths.where(lengths.notna(), hmmer['qlen'] * 3)

def check_biorisk(hmmscan_input_file : str, biorisk_annotations_directory : str):
    """
    Checks an HMM scan output, and parses it for biorisks, according to those found in the biorisk_annotations.csv.
    INPUTS:
        - hmmscan_input_file - the file output from hmmscan, containing information about potential hits.
        - hmm_folder - the directory containing biorisk_annotations.csv
    """

    # check input files
    hmm_folder_csv = biorisk_annotations_directory + "/biorisk_annotations.csv"
    if not os.path.exists(hmmscan_input_file):
        logging.error("\t...input file does not exist\n")
        return 1
    if not os.path.exists(hmm_folder_csv):
        logging.error("\t...biorisk_annotations.csv does not exist\n" + hmm_folder_csv)
        return 1

    # Specify input file and read in database file
    lookup = pd.read_csv(hmm_folder_csv)
    lookup.fillna(False, inplace=True)

    # read in HMMER output and check for valid hits
    if HmmerHandler.is_empty(hmmscan_input_file):
        logging.info("\t...ERROR: biorisk search results empty\n")
        return 1

    if not HmmerHandler.has_hits(hmmscan_input_file):
        logging.info("\t\t --> Biorisks: no hits detected, PASS\n")
        return 1

    hmmer = readhmmer(hmmscan_input_file)

    keep1 = [i for i, x in enumerate(hmmer["E-value"]) if x < 1e-20]
    hmmer = hmmer.iloc[keep1, :]

    # Recalculate hit ranges into query based nucleotide coordinates
    hmmer = remove_overlaps(hmmer)

    hmmer["description"] = ""
    hmmer["Must flag"] = False
    hmmer = hmmer.reset_index(drop=True)

    for model in range(hmmer.shape[0]):
        name_index = [
            i
            for i, x in enumerate([lookup["ID"] == hmmer["target name"][model]][0])
            if x
        ]
        hmmer.loc[model, "description"] = lookup.iloc[name_index[0], 1]
        hmmer.loc[model, "Must flag"] = lookup.iloc[name_index[0], 2]

    if hmmer.shape[0] == 0:
        logging.info("\t\t --> Biorisks: no significant hits detected, PASS\n")
        return

    if sum(hmmer["Must flag"]) > 0:
        for region in hmmer.index[hmmer["Must flag"] != 0]:
            if hmmer["ali from"][region] > hmmer["qlen"][region]:
                hmmer["ali from"][region] = divmod(
                    hmmer["ali from"][region], hmmer["qlen"][region]
                )[0]
                hmmer["ali to"][region] = divmod(
                    hmmer["ali to"][region], hmmer["qlen"][region]
                )[0]
            logging.info(
                "\t\t --> Biorisks: Regulated gene in bases "
                + str(hmmer["q. start"][region])
                + " to "
                + str(hmmer["q. end"][region])
                + ": FLAG\n\t\t     Gene: "
                + ", ".join(set(hmmer["description"][hmmer["Must flag"] == True]))
                + "\n"
            )

    else:
        logging.info("\t\t --> Biorisks: Regulated genes not found, PASS\n")
        return 0

    if sum(hmmer["Must flag"]) != hmmer.shape[0]:
        for region in hmmer.index[hmmer["Must flag"] == 0]:
            logging.info(
                "\t\t --> Virulence factor found in bases "
                + str(hmmer["q. start"][region])
                + " to "
                + str(hmmer["q. end"][region])
                + ", WARNING\n\t\t     Gene: "
                + ", ".join(set(hmmer["description"][hmmer["Must flag"] == False]))
                + "\n"
            )

    return 0


def main():
    """
    Wrapper for parsing arguments direction to check_biorisk if called as main.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="in_file",
        required=True,
        help="Input file - hmmscan output file",
    )
    parser.add_argument(
        "-d",
        "--database",
        dest="db",
        required=True,
        help="HMM folder (must contain biorisk_annotations.csv)",
    )
    parser.add_argument(
        "-o", "--out", dest="output_json", required=True, help="output_json_filepath"
    )
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    return_value = check_biorisk(args.in_file, args.db)
    return return_value


if __name__ == "__main__":
    main()
