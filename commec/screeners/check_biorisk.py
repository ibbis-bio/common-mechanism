#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks the output from hmmscan and prints to screen the results

Usage:
 python check_biorisk.py -i INPUT.biorisk.hmmscan -d databases/biorisk_db/ 
"""
import logging
import os
import pandas as pd
from commec.tools.hmmer import readhmmer, remove_overlaps, recalculate_hmmer_query_coordinates, append_nt_querylength_info, HmmerHandler
from commec.config.query import Query
from commec.config.result import (
    ScreenResult,
    HitResult,
    ScreenStep,
    ScreenStatus,
    HitScreenStatus,
    MatchRange,
    compare
)

logger = logging.getLogger(__name__)

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

def update_biorisk_data_from_database(search_handle : HmmerHandler, data : ScreenResult, queries : dict[str, Query]):
    """
    Takes an input database, reads its outputs, and updates the input data to contain
    biorisk hits from the database. Also requires passing of the biorisk annotations CSV file.
    Inputs:
        search : search_handle - The handler which has performed a search on a database.
        biorisk_annotations_csv_file : str - directory/filename of the biorisk annotations provided by Commec.
        data : ScreenResult - The ScreenResult to be updated with information from database, interpeted as Biorisks.
    """
    # Check for annocations.csv, as well as whether the 
    logging.debug("Directory: %s", search_handle.db_directory)
    logging.debug("Directory/file: %s", search_handle.db_file)
    #logging.debug("Directory/file: %s", search_handle.db_file)
    hmm_folder_csv = os.path.join(search_handle.db_directory,"biorisk_annotations.csv")
    if not os.path.exists(hmm_folder_csv):
        logger.error("\t...biorisk_annotations.csv does not exist\n %s", hmm_folder_csv)
        return 1
    if not search_handle.check_output():
        logger.error("\t...database output file does not exist\n %s", search_handle.out_file)
        return 1
    if search_handle.is_empty(search_handle.out_file):
        logger.error("\t...ERROR: biorisk search results empty\n")
        return 1

    for query in data.queries.values():
        query.recommendation.biorisk_status = ScreenStatus.PASS

    if not search_handle.has_hits(search_handle.out_file):
        return 0

    # Read in Output, and parse.
    hmmer : pd.DataFrame = readhmmer(search_handle.out_file)
    keep1 = [i for i, x in enumerate(hmmer['E-value']) if x < 1e-20]
    hmmer = hmmer.iloc[keep1,:]
    
    append_nt_querylength_info(hmmer, queries)
    recalculate_hmmer_query_coordinates(hmmer)
    hmmer = remove_overlaps(hmmer)

    # Read in annotations.
    lookup : pd.DataFrame = pd.read_csv(hmm_folder_csv)
    lookup.fillna(False, inplace=True)

    # Append description, and must_flag columns from annotations:
    hmmer['description'] = ''
    hmmer['Must flag'] = False
    hmmer = hmmer.reset_index(drop=True)
    for model in range(hmmer.shape[0]):
        name_index = [i for i, x in enumerate([lookup['ID'] == hmmer['target name'][model]][0]) if x]
        hmmer.loc[model, 'description'] = lookup.iloc[name_index[0], 1]
        hmmer.loc[model, 'Must flag'] = lookup.iloc[name_index[0], 2]

    # Update the data state to capture the outputs from biorisk search:
    unique_queries = hmmer['query name'].unique()
    for affected_query in unique_queries:

        biorisk_overall : ScreenStatus = ScreenStatus.PASS

        query_data = data.get_query(affected_query)
        if not query_data:
            logging.error("Query during hmmscan could not be found! [%s]", affected_query)
            continue

        # Grab a list of unique queries, and targets for iteration.
        unique_query_data : pd.DataFrame = hmmer[hmmer['query name'] == affected_query]
        unique_targets = unique_query_data['target name'].unique()

        for affected_target in unique_targets:
            unique_target_data : pd.DataFrame = unique_query_data[unique_query_data['target name'] == affected_target]
            target_description = ", ".join(set(unique_target_data['description'])) # First should be unique.
            must_flag = unique_target_data['Must flag'].iloc[0] # First should be unique.
            match_ranges = []
            for _, region in unique_target_data.iterrows():
                match_range = MatchRange(
                    float(region['E-value']),
                    int(region['hmm from']), int(region['hmm to']),
                    int(region['q. start']), int(region['q. end'])
                )
                match_ranges.append(match_range)

            target_recommendation : ScreenStatus = ScreenStatus.FLAG if must_flag > 0 else ScreenStatus.WARN

            biorisk_overall = compare(target_recommendation, biorisk_overall)

            hit_data : HitResult = query_data.get_hit(affected_target)
            if hit_data:
                hit_data.ranges.extend(match_ranges)
                continue

            regulation_str : str = "Regulated Gene" if must_flag else "Virulance Factor"
            
            domain : str = _guess_domain(""+str(affected_target)+target_description)
            
            new_hit : HitResult = HitResult(
                HitScreenStatus(
                    target_recommendation,
                    ScreenStep.BIORISK
                ),
                affected_target,
                target_description,
                match_ranges,
                {"domain" : [domain],"regulated":[regulation_str]},
            )
            query_data.hits[affected_target] = new_hit

        # Update the recommendation for this query for biorisk.
        query_data.recommendation.biorisk_status = biorisk_overall
    return 0

def check_biorisk(hmmscan_input_file : str, biorisk_annotations_directory : str, queries : dict[str,Query]):
    """
    LEGACY .screen output content.

    Checks an HMM scan output, and parses it for biorisks, according to those found in the biorisk_annotations.csv.
    INPUTS:
        - hmmscan_input_file - the file output from hmmscan, containing information about potential hits.
        - hmm_folder - the directory containing biorisk_annotations.csv
    RETURNS:
        0 or 1 depending on whether execution was successful
    """

    # check input files
    hmm_folder_csv = biorisk_annotations_directory + "/biorisk_annotations.csv"
    if not os.path.exists(hmmscan_input_file):
        logger.error(f"\t...hmmscan file does not exist: {hmmscan_input_file}")
        return 1
    if not os.path.exists(hmm_folder_csv):
        logger.error(f"\t...Biorisk annotations file does not exist: {hmm_folder_csv}")
        return 1

    # Specify input file and read in database file
    lookup = pd.read_csv(hmm_folder_csv)
    lookup.fillna(False, inplace=True)

    # read in HMMER output and check for valid hits
    if HmmerHandler.is_empty(hmmscan_input_file):
        logger.error("\t... biorisk search results empty\n")
        return 1

    if not HmmerHandler.has_hits(hmmscan_input_file):
        logger.info("\t\t --> Biorisks: no hits detected, PASS\n")
        return 0

    hmmer = readhmmer(hmmscan_input_file)

    keep1 = [i for i, x in enumerate(hmmer["E-value"]) if x < 1e-20]
    hmmer = hmmer.iloc[keep1, :]

    # Recalculate hit ranges into query based nucleotide coordinates, and trim overlaps.
    append_nt_querylength_info(hmmer, queries)
    recalculate_hmmer_query_coordinates(hmmer)
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

