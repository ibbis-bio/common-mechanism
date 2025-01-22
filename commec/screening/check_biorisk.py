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
from commec.tools.hmmer import readhmmer, trimhmmer, HmmerHandler
from commec.config.json_io import (
    ScreenData,
    HitDescription,
    CommecScreenStep,
    CommecRecommendation,
    CommecScreenStepRecommendation,
    MatchRange,
    guess_domain,
    compare
)

def update_biorisk_data_from_database(search_handle : HmmerHandler, data : ScreenData):
    """
    Takes an input database, reads its outputs, and updates the input data to contain
    biorisk hits from the database. Also requires passing of the biorisk annotations CSV file.
    Inputs:
        search : search_handle - The handler which has performed a search on a database.
        biorisk_annotations_csv_file : str - directory/filename of the biorisk annotations provided by Commec.
        data : ScreenData - The ScreenData to be updated with information from database, interpeted as Biorisks.
    """
    # Check for annocations.csv, as well as whether the 
    logging.debug("Directory: %s", search_handle.db_directory)
    logging.debug("Directory/file: %s", search_handle.db_file)
    #logging.debug("Directory/file: %s", search_handle.db_file)
    hmm_folder_csv = os.path.join(search_handle.db_directory,"biorisk_annotations.csv")
    if not os.path.exists(hmm_folder_csv):
        logging.error("\t...biorisk_annotations.csv does not exist\n %s", hmm_folder_csv)
        return
    if not search_handle.check_output():
        logging.error("\t...database output file does not exist\n %s", search_handle.out_file)
        return
    if search_handle.is_empty(search_handle.out_file):
        logging.error("\t...ERROR: biorisk search results empty\n")
        return

    for query in data.queries:
        query.recommendation.biorisk_screen = CommecRecommendation.PASS

    if not search_handle.has_hits(search_handle.out_file):
        return 0

    # Read in Output, and parse.
    hmmer : pd.DataFrame = readhmmer(search_handle.out_file)
    keep1 = [i for i, x in enumerate(hmmer['E-value']) if x < 1e-20]
    hmmer = hmmer.iloc[keep1,:]
    hmmer = trimhmmer(hmmer)

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

        biorisk_overall : CommecRecommendation = CommecRecommendation.PASS

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
                    int(region['ali from']), int(region['ali to'])
                )
                match_ranges.append(match_range)

            target_recommendation : CommecRecommendation = CommecRecommendation.FLAG if must_flag > 0 else CommecRecommendation.WARN

            biorisk_overall = compare(target_recommendation, biorisk_overall)

            hit_data : HitDescription = query_data.get_hit(affected_target)
            if hit_data:
                hit_data.ranges.extend(match_ranges)
                continue

            regulation_str : str = "Regulated Gene" if must_flag else "Virulance Factor"
            
            domain : str = guess_domain(""+str(affected_target)+target_description)
            
            new_hit : HitDescription = HitDescription(
                CommecScreenStepRecommendation(
                    target_recommendation,
                    CommecScreenStep.BIORISK
                ),
                affected_target,
                target_description,
                match_ranges,
                {"domain" : [domain],"regulated":[regulation_str]},
            )
            query_data.hits.append(new_hit)

        # Update the recommendation for this query for biorisk.
        query_data.recommendation.biorisk_screen = biorisk_overall

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

    #Specify input file and read in database file
    lookup = pd.read_csv(hmm_folder_csv)
    lookup.fillna(False, inplace=True)

    # read in HMMER output and check for valid hits
    if HmmerHandler.is_empty(hmmscan_input_file):
        logging.info("\t...ERROR: biorisk search results empty\n")
        return 0

    if not HmmerHandler.has_hits(hmmscan_input_file):
        logging.info("\t\t --> Biorisks: no hits detected, PASS\n")
        return 0

    hmmer = readhmmer(hmmscan_input_file)
    keep1 = [i for i, x in enumerate(hmmer['E-value']) if x < 1e-20]
    hmmer = hmmer.iloc[keep1,:]
    hmmer = trimhmmer(hmmer)
    hmmer['description'] = ''
    hmmer['Must flag'] = False
    hmmer = hmmer.reset_index(drop=True)

    for model in range(hmmer.shape[0]):
        name_index = [i for i, x in enumerate([lookup['ID'] == hmmer['target name'][model]][0]) if x]
        hmmer.loc[model, 'description'] = lookup.iloc[name_index[0], 1]
        hmmer.loc[model, 'Must flag'] = lookup.iloc[name_index[0], 2]

    if hmmer.shape[0] == 0:
        logging.info("\t\t --> Biorisks: no significant hits detected, PASS\n")
        return 0
    
    if sum(hmmer['Must flag']) == 0:
        logging.info("\t\t --> Biorisks: Regulated genes not found, PASS\n")
        return 0

    if sum(hmmer['Must flag']) > 0:
        for region in hmmer.index[hmmer['Must flag'] != 0]:
            if hmmer['ali from'][region] > hmmer['qlen'][region]:
                hmmer['ali from'][region] = divmod(hmmer['ali from'][region], hmmer['qlen'][region])[0]
                hmmer['ali to'][region] = divmod(hmmer['ali to'][region], hmmer['qlen'][region])[0]

            logging.info("\t\t --> Biorisks: Regulated gene in bases " + str(hmmer['ali from'][region]) +
                            " to " + str(hmmer['ali to'][region]) + 
                            ": FLAG\n\t\t     Gene: " + 
                            ", ".join(set(hmmer['description'][hmmer['Must flag'] == True])) + "\n")

    if sum(hmmer['Must flag']) != hmmer.shape[0]:
        for region in hmmer.index[hmmer['Must flag'] == 0]:
            logging.info("\t\t --> Virulence factor found in bases " + str(hmmer['ali from'][region]) +
                                " to " + str(hmmer['ali to'][region]) +
                                ", WARNING\n\t\t     Gene: " +
                                ", ".join(set(hmmer['description'][hmmer['Must flag'] == False])) + "\n")

    return 0

def main():
    '''
    Wrapper for parsing arguments direction to check_biorisk if called as main.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", dest="in_file",
        required=True, help="Input file - hmmscan output file")
    parser.add_argument("-d","--database", dest="db",
        required=True, help="HMM folder (must contain biorisk_annotations.csv)")
    parser.add_argument("-o","--out", dest="output_json",
        required=True,help="output_json_filepath")
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
