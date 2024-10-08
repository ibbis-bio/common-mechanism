#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks the output from hmmscan and prints to screen the results

Usage:
    python check_benign.py -i INPUT -s SEQUENCE -d DATABASE FOLDER
      -i, --input = input sample name (will check for sample.benign.hmmscan file)
      -s, --sequence = input sequence file
      -d, --database = database folder location/path (will check for benign_annotations.tsv)
"""
import logging
import argparse
import os
import sys
import pandas as pd

from commec.tools.blastn import BlastNHandler # For has_hits.
from commec.tools.blast_tools import tophits, readblast, trimblast
from commec.tools.hmmer import readhmmer, HmmerHandler
from commec.tools.cmscan import readcmscan

from commec.config.json_io import (
    ScreenData,
    HitDescription,
    CommecScreenStep,
    CommecRecomendation,
    CommecScreenStepRecommendation,
    MatchRange,
    LifeDomainFlag,
    RegulationFlag,
    guess_domain,
    compare
)

logger = logging.getLogger(__name__)

def update_benign_data_from_database(search_handle : HmmerHandler, data : ScreenData):

    if not search_handle.check_output():
        logging.info("\t...no housekeeping protein data\n")
    if not search_handle.has_hits(search_handle.out_file):
        logging.info("\t...no housekeeping protein hits\n")
    else:
        hmmer = readhmmer(search_handle.out_file)
        hmmer = hmmer[hmmer["E-value"] < 1e-20]
        # print(hmmer)

def check_for_benign(query, coords, benign_desc):
    '''
    Checks a query against taxonomy 
    '''
    cleared = [0] * coords.shape[0]

    # PROTEIN HITS
    # for each set of hits, need to pull out the coordinates covered by benign entries
    hmmscan = query + ".benign.hmmscan"
    if not BlastNHandler.has_hits(hmmscan):
        logger.info("\t...no housekeeping protein hits\n")
    else:
        hmmer = readhmmer(hmmscan)
        hmmer = hmmer[hmmer["E-value"] < 1e-20]
        for region in range(0, coords.shape[0]):  # for each regulated pathogen region
            # look at only the hmmer hits that overlap with it
            htrim = hmmer[
                ~(
                    (hmmer["ali from"] > coords["q. end"][region])
                    & (hmmer["ali to"] > coords["q. end"][region])
                )
                & ~(
                    (hmmer["ali from"] < coords["q. start"][region])
                    & (hmmer["ali to"] < coords["q. start"][region])
                )
            ]
            if htrim.shape[0] > 0:
                htrim = htrim.assign(coverage=abs(htrim["ali to"] - htrim["ali from"]))
                if any(htrim["qlen"] - htrim["coverage"] < 50):
                    htrim = htrim[htrim["coverage"] > 0.80]
                    htrim = htrim.reset_index(drop=True)
                    descriptions = []
                    # for row in range(htrim.shape[0]):
                    for row in [0]:  # just print the top hit
                        hit = htrim["target name"][row]
                        hit_msg = (
                            hit
                            + ": "
                            + str(*benign_desc["Description"][benign_desc["ID"] == hit])
                            + f" (E-value: {htrim['E-value'][row]:.3g}"
                        )
                        descriptions.append(hit_msg + "\n")
                    annot_string = "\n".join(str(v) for v in descriptions)
                    logger.info(
                        "\t\t -->Housekeeping proteins covering "
                        + str(coords["q. start"][region])
                        + " to "
                        + str(coords["q. end"][region])
                        + " = PASS\n"
                    )
                    logger.info("\t\t   " + annot_string)
                    cleared[region] = 1
                else:
                    logger.info(
                        "\t\t -->Housekeeping proteins - not enough coverage = FAIL\n"
                    )

    # RNA HITS
    # for each set of hits, need to pull out the coordinates covered by benign entries
    cmscan = query + ".benign.cmscan"
    if not BlastNHandler.has_hits(cmscan):
        logger.info("\t...no benign RNA hits\n")
    else:
        cmscan = readcmscan(cmscan)
        for region in range(0, coords.shape[0]):  # for each regulated pathogen region
            # look at only the cmscan hits that overlap with it
            qlen = abs(coords["q. start"][region] - coords["q. end"][region])
            # filter hits for ones that overlap with the regulated region
            htrim = cmscan[
                ~(
                    (cmscan["seq from"] < coords["q. start"][region])
                    & (cmscan["seq to"] < coords["q. start"][region])
                )
                & ~(
                    (cmscan["seq from"] > coords["q. end"][region])
                    & (cmscan["seq to"] > coords["q. end"][region])
                )
            ]
            if htrim.shape[0] > 0:
                # bases unaccounted for based method
                htrim = htrim.assign(
                    coverage=qlen - abs(htrim["seq to"] - htrim["seq from"])
                )
                if any(htrim["coverage"] < 50):
                    htrim = htrim[htrim["coverage"] < 50]
                    htrim = htrim.reset_index(drop=True)
                    descriptions = []
                    for row in range(htrim.shape[0]):
                        hit = htrim["target name"][row]
                        descriptions.append(hit)
                    annot_string = "\n\t...".join(str(v) for v in descriptions)
                    logger.info(
                        "\t\t -->Housekeeping RNAs - <50 bases unaccounted for: PASS\n"
                    )
                    logger.info("\t\t   RNA family: " + annot_string + "\n")
                    cleared[region] = 1
                else:
                    logger.info(
                        "\t\t -->Housekeeping RNAs - >50 bases unaccounted for = FAIL\n"
                    )

    # SYNBIO HITS
    # annotate and clear benign nucleotide sequences
    blast = query + ".benign.blastn"
    if not BlastNHandler.has_hits(blast):
        logger.info("\t...no Synbio sequence hits\n")
    else:
        blastn = read_blast(blast)  # synbio parts
        blastn = get_top_hits(blastn)
        for region in range(0, coords.shape[0]):  # for each regulated pathogen region
            htrim = blastn[
                ~(
                    (blastn["q. start"] > coords["q. end"][region])
                    & (blastn["q. end"] > coords["q. end"][region])
                )
                & ~(
                    (blastn["q. start"] < coords["q. start"][region])
                    & (blastn["q. end"] < coords["q. start"][region])
                )
            ]
            if any(htrim["q. coverage"] > 0.80):
                htrim = htrim[htrim["q. coverage"] > 0.80]
                htrim = htrim.reset_index(drop=True)
                descriptions = []
                for row in range(htrim.shape[0]):
                    hit = htrim["subject title"][row]
                    descriptions.append(hit)
                annot_string = "\n\t\t   ".join(str(v) for v in descriptions)
                logger.info(
                    "\t\t -->Synbio sequences - >80% coverage achieved = PASS\n"
                )
                logger.info("\t\t   Synbio parts: " + annot_string + "\n")
                cleared[region] = 1
            else:
                logger.info(
                    "\t\t -->Synbio sequences - <80% coverage achieved = FAIL\n"
                )

    for region in range(0, coords.shape[0]):
        if cleared[region] == 0:
            logger.info(
                "\t\t -->Regulated region at bases "
                + str(int(coords.iloc[region, 0]))
                + " to "
                + str(int(coords.iloc[region, 1]))
                + " failed to clear: FLAG\n"
            )
    if sum(cleared) == len(cleared):
        logger.info("\n\t\t -->all regulated regions cleared: PASS\n")

    return 0


def main():
    """
    Alternative legacy entry point wrapper for calling check_benign as main,
    as isolated python script. No longer used in commec screen.
    """
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

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i", "--input", dest="sample_name", required=True, help="Sample name"
    )
    parser.add_argument(
        "-d",
        "--database",
        dest="db",
        required=True,
        help="Benign HMM database folder (must contain benign_annotations.tsv)",
    )
    parser.add_argument(
        "-o", "--out", dest="output_json", required=True, help="output_json_filepath"
    )

    args = parser.parse_args()

    if not os.path.exists(args.db + "/benign_annotations.tsv"):
        sys.stderr.write("\t...benign_annotations.tsv does not exist\n")
        return 1

    # Check database file
    pd.set_option("max_colwidth", 200)
    benign_desc = pd.read_csv(args.db + "/benign_annotations.tsv", sep="\t")

    if not os.path.exists(args.sample_name + ".reg_path_coords.csv"):
        logger.info("\t...no regulated regions to clear\n")
        return 0

    # Read database file
    coords = pd.read_csv(args.sample_name + ".reg_path_coords.csv")

    if coords.shape[0] == 0:
        logger.info("\t...no regulated regions to clear\n")
        return 0

    coords.sort_values(by=["q. start"], inplace=True)
    coords.reset_index(drop=True, inplace=True)
    rv = check_for_benign(args.sample_name, coords, benign_desc)
    sys.exit(rv)


if __name__ == "__main__":
    main()
