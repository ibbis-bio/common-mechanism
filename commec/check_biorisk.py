#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks the output from hmmscan and prints to screen the results

Usage:
 python check_biorisk.py -i INPUT.biorisk.hmmscan -d databases/biorisk_db/ 
"""
import os
import sys
import argparse
import pandas as pd
import logging

from commec.utils import (
    has_hits,
    readhmmer,
    trimhmmer,)

from commec.file_tools import FileTools
from commec.json_io import (
    ScreenData,
    BioRisk,
    MatchRange,
    encode_screen_data_to_json,
    get_screen_data_from_json
)

def check_biorisk(hmmscan_input_file : str, biorisk_annotations_directory : str, output_json : str):
    '''
    Checks an HMM scan output, and parses it for biorisks, according to those found in the biorisk_annotations.csv.
    INPUTS:
        - hmmscan_input_file - the file output from hmmscan, containing information about potential hits.
        - hmm_folder - the directory containing biorisk_annotations.csv
    '''

    #check input files
    hmm_folder_csv  = biorisk_annotations_directory + "/biorisk_annotations.csv"
    if not os.path.exists(hmmscan_input_file):
        logging.error("\t...input file does not exist\n")
        return 1
    if not os.path.exists(hmm_folder_csv):
        logging.error("\t...biorisk_annotations.csv does not exist\n" + hmm_folder_csv)
        return 1

    output_data : ScreenData = get_screen_data_from_json(output_json)
    
    #Specify input file and read in database file
    lookup = pd.read_csv(hmm_folder_csv)
    lookup.fillna(False, inplace=True)

    # read in HMMER output and check for valid hits
    if FileTools.is_empty(hmmscan_input_file):
        logging.info("\t...ERROR: biorisk search results empty\n")
        return

    if not has_hits(hmmscan_input_file):
        logging.info("\t\t --> Biorisks: no hits detected, PASS\n")
        return

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
        return

    if sum(hmmer['Must flag']) > 0:
        for region in hmmer.index[hmmer['Must flag'] != 0]:
            if hmmer['ali from'][region] > hmmer['qlen'][region]:
                hmmer['ali from'][region] = divmod(hmmer['ali from'][region], hmmer['qlen'][region])[0]
                hmmer['ali to'][region] = divmod(hmmer['ali to'][region], hmmer['qlen'][region])[0]
            logging.info("\t\t --> Biorisks: Regulated gene in bases " + str(hmmer['ali from'][region]) +
                                " to " + str(hmmer['ali to'][region]) + 
                                ": FLAG\n\t\t     Gene: " + 
                                ", ".join(set(hmmer['description'][hmmer['Must flag'] == True])) + "\n")
            
            new_match = BioRisk(
                ", ".join(set(hmmer['description'][hmmer['Must flag'] == True])),
                MatchRange(
                    0,0,
                    int(hmmer['ali from'][region]),
                    int(hmmer['ali to'][region]),
                    0
                )
            )
            output_data.biorisks.regulated_genes.append(new_match)
    else:
        logging.info("\t\t --> Biorisks: Regulated genes not found, PASS\n")

    if sum(hmmer['Must flag']) != hmmer.shape[0]:
        for region in hmmer.index[hmmer['Must flag'] == 0]:
            logging.info("\t\t --> Virulence factor found in bases " + str(hmmer['ali from'][region]) +
                                " to " + str(hmmer['ali to'][region]) +
                                ", WARNING\n\t\t     Gene: " +
                                ", ".join(set(hmmer['description'][hmmer['Must flag'] == False])) + "\n")

            new_match = BioRisk(
                ", ".join(set(hmmer['description'][hmmer['Must flag'] == False])),
                MatchRange(
                    0,0,
                    int(hmmer['ali from'][region]),
                    int(hmmer['ali to'][region]),
                    0
                )
            )
            output_data.biorisks.virulance_factors.append(new_match)
    encode_screen_data_to_json(output_data, output_json)
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

    return_value = check_biorisk(args.in_file, args.db, args.output_json)
    return return_value
if __name__ == "__main__":
    main()
