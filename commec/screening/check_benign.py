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
from commec.tools.blastn import BlastNHandler  # For has_hits.
from commec.tools.hmmer import HmmerHandler
from commec.tools.blast_tools import get_top_hits, read_blast
from commec.tools.hmmer import readhmmer
from commec.tools.cmscan import CmscanHandler, readcmscan
from commec.tools.search_handler import SearchHandler

from commec.config.json_io import (
    ScreenData,
    HitDescription,
    QueryData,
    CommecScreenStep,
    CommecRecommendation,
    CommecScreenStepRecommendation,
    MatchRange,
    guess_domain,
    compare
)

# Constants determining Commec's sensitivity for benign screen.
BENIGN_PROTEIN_EVALUE_CUTOFF : float = 1e-20
MINIMUM_PEPTIDE_COVERAGE : int = 50 # Consider the implications of whether coverage in this context represents the translated aa, or bp count.
MINIMUM_QUERY_COVERAGE_FRACTION : float = 0.80
MINIMUM_RNA_BASEPAIR_COVERAGE : int = 50
MINIMUM_SYNBIO_COVERAGE_FRACTION : float = 0.80

def _update_benign_data_for_query(query : QueryData,
                                  benign_protein : pd.DataFrame,
                                  benign_rna : pd.DataFrame,
                                  benign_synbio : pd.DataFrame,
                                  benign_descriptions : pd.DataFrame):
    """
    For a single query, look at all three benign database outputs, and update the 
    single queries hit descriptions to record all the benign hits, as well as clear
    any overlapping WARN or FLAG hits.
    """
    # We only care about the benign data for this query.
    benign_protein_for_query = benign_protein[benign_protein["query name"] == query.query]
    benign_rna_for_query = benign_rna[benign_rna["query name"] == query.query]
    benign_synbio_for_query = benign_synbio[benign_synbio["query acc."] == query.query]

    # Check every region, of every hit that is a FLAG or WARN, against the Benign screen outcomes.
    for hit in query.hits:
        if hit.recommendation.outcome not in {
            CommecRecommendation.FLAG,
            CommecRecommendation.WARN
            }:
            continue

        for region in hit.ranges:
            benign_protein_for_query = _trim_to_region(benign_protein_for_query, region)

            # Filter benign proteins for relevance...
            benign_protein_for_query = benign_protein_for_query.assign(
                coverage=abs(
                    benign_protein_for_query["q. start"] - benign_protein_for_query["q. end"]
                    ))

            benign_protein_for_query = benign_protein_for_query[
                benign_protein_for_query["query length"] - benign_protein_for_query["coverage"]
                < MINIMUM_PEPTIDE_COVERAGE
                ]

            benign_protein_for_query = benign_protein_for_query[
                benign_protein_for_query["coverage"] > MINIMUM_QUERY_COVERAGE_FRACTION]
            benign_protein_for_query = benign_protein_for_query.reset_index(drop=True)

            # Filter benign RNA for relevance...
            benign_rna_for_query = _trim_to_region(benign_rna_for_query, region)
            benign_rna_for_query = benign_rna_for_query.assign(
                        coverage=region.query_length - abs(benign_rna_for_query["q. end"] - benign_rna_for_query["q. start"])
                    )
            
            benign_rna_for_query = benign_rna_for_query[
                benign_rna_for_query["coverage"] < MINIMUM_RNA_BASEPAIR_COVERAGE]
            benign_rna_for_query = benign_rna_for_query.reset_index(drop=True)

            # Filter benign SynBio for relevance... 
            benign_synbio_for_query = _trim_to_region(benign_synbio_for_query, region)
            benign_synbio_for_query = benign_synbio_for_query[
                benign_synbio_for_query["q. coverage"] > MINIMUM_SYNBIO_COVERAGE_FRACTION]
            benign_synbio_for_query = benign_synbio_for_query.reset_index(drop=True)

            # Report top hit for 
            if not benign_protein_for_query.empty:
                benign_hit = benign_protein_for_query["subject title"][0]
                benign_hit_description = str(*benign_descriptions["Description"][benign_descriptions["ID"] == benign_hit])
                match_ranges = [
                    MatchRange(
                    float(benign_protein_for_query['evalue'][0]),
                    int(benign_protein_for_query['s. start'][0]), int(benign_protein_for_query['s. end'][0]),
                    int(benign_protein_for_query['q. start'][0]), int(benign_protein_for_query['q. end'][0])
                    )
                ]
                benign_hit_outcome = HitDescription(
                        CommecScreenStepRecommendation(
                            CommecRecommendation.PASS,
                            CommecScreenStep.BENIGN_PROTEIN
                        ),
                        benign_hit,
                        benign_hit_description,
                        match_ranges,
                    )
                query.add_new_hit_information(benign_hit_outcome)
                hit.recommendation.outcome = hit.recommendation.outcome.clear()

            if not benign_rna_for_query.empty:
                benign_hit = benign_rna_for_query["subject title"][0]
                benign_hit_description =  benign_rna_for_query["description of target"][0]
                match_ranges = [
                    MatchRange(
                    float(benign_rna_for_query['evalue'][0]),
                    int(benign_rna_for_query['s. start'][0]), int(benign_rna_for_query['s. end'][0]),
                    int(benign_rna_for_query['q. start'][0]), int(benign_rna_for_query['q. end'][0])
                    )
                ]
                benign_hit_outcome = HitDescription(
                        CommecScreenStepRecommendation(
                            CommecRecommendation.PASS,
                            CommecScreenStep.BENIGN_RNA
                        ),
                        benign_hit,
                        benign_hit_description,
                        match_ranges,
                    )
                query.add_new_hit_information(benign_hit_outcome)
                hit.recommendation.outcome = hit.recommendation.outcome.clear()

            if not benign_synbio_for_query.empty:
                benign_hit = benign_synbio_for_query["subject title"][0]
                benign_hit_description =  benign_synbio_for_query["subject title"][0]
                match_ranges = [
                    MatchRange(
                    float(benign_synbio_for_query['evalue'][0]),
                    int(benign_synbio_for_query['s. start'][0]), int(benign_synbio_for_query['s. end'][0]),
                    int(benign_synbio_for_query['q. start'][0]), int(benign_synbio_for_query['q. end'][0])
                    )
                ]
                benign_hit_outcome = HitDescription(
                        CommecScreenStepRecommendation(
                            CommecRecommendation.PASS,
                            CommecScreenStep.BENIGN_SYNBIO
                        ),
                        benign_hit,
                        benign_hit_description,
                        match_ranges,
                    )
                query.add_new_hit_information(benign_hit_outcome)
                hit.recommendation.outcome = hit.recommendation.outcome.clear()


def update_benign_data_from_database(benign_protein_handle : HmmerHandler,
                                     benign_rna_handle : CmscanHandler,
                                     benign_synbio_handle : BlastNHandler,
                                     data : ScreenData, 
                                     benign_desc : pd.DataFrame):
    """
    Parse the outputs from the protein, rna, and synbio database searches, and populate
    the benign hits into a Screen dataset. Marks those hits that are cleared for benign
    as cleared if benign screen passes them.
    """
    # Reading empty outcomes should result in empty DataFrames, not errors.
    benign_protein_screen_data = benign_protein_handle.read_output()
    benign_rna_screen_data = benign_rna_handle.read_output()
    benign_synbio_screen_data = benign_synbio_handle.read_output()

    for query in data.queries:
        _update_benign_data_for_query(query,
                                      benign_protein_screen_data,
                                      benign_rna_screen_data,
                                      benign_synbio_screen_data,
                                      benign_desc)

        # Calculate the Benign Screen outcomes for each query.
        query.recommendation.benign_screen = CommecRecommendation.PASS
        # If any hits are still warnings, or flags, propagate that.
        for flagged_hit in query.get_flagged_hits():
            query.recommendation.benign_screen = compare(
                flagged_hit.recommendation.outcome,
                query.recommendation.benign_screen
                )

def _trim_to_coords(data : pd.DataFrame, coords, region):
    datatrim = data[
        ~(
            (data["q. start"] > coords["q. end"][region])
            & (data["q. end"] > coords["q. end"][region])
        )
        & ~(
            (data["q. start"] < coords["q. start"][region])
            & (data["q. end"] < coords["q. start"][region])
        )
    ]
    return datatrim

def _trim_to_region(data : pd.DataFrame, region : MatchRange):
    datatrim = data[
        ~(
            (data["q. start"] > region.query_end)
            & (data["q. end"] > region.query_end)
        )
        & ~(
            (data["q. start"] < region.query_start)
            & (data["q. end"] < region.query_end)
        )
    ]
    return datatrim

def check_for_benign(query, coords, benign_desc):
    """
    Checks a query against taxonomy
    """
    cleared = [0] * coords.shape[0]

    # PROTEIN HITS
    # for each set of hits, need to pull out the coordinates covered by benign entries
    hmmscan = query + ".benign.hmmscan"
    if not HmmerHandler.has_hits(hmmscan):
        logging.info("\t...no housekeeping protein hits\n")
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
                    logging.info(
                        "\t\t -->Housekeeping proteins covering "
                        + str(coords["q. start"][region])
                        + " to "
                        + str(coords["q. end"][region])
                        + " = PASS\n"
                    )
                    logging.info("\t\t   " + annot_string)
                    cleared[region] = 1
                else:
                    logging.info(
                        "\t\t -->Housekeeping proteins - not enough coverage = FAIL\n"
                    )

    # RNA HITS
    # for each set of hits, need to pull out the coordinates covered by benign entries
    cmscan = query + ".benign.cmscan"
    if not BlastNHandler.has_hits(cmscan):
        logging.info("\t...no benign RNA hits\n")
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
                    logging.info(
                        "\t\t -->Housekeeping RNAs - <50 bases unaccounted for: PASS\n"
                    )
                    logging.info("\t\t   RNA family: " + annot_string + "\n")
                    cleared[region] = 1
                else:
                    logging.info(
                        "\t\t -->Housekeeping RNAs - >50 bases unaccounted for = FAIL\n"
                    )

    # SYNBIO HITS
    # annotate and clear benign nucleotide sequences
    blast = query + ".benign.blastn"
    if not BlastNHandler.has_hits(blast):
        logging.info("\t...no Synbio sequence hits\n")
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
                logging.info(
                    "\t\t -->Synbio sequences - >80% coverage achieved = PASS\n"
                )
                logging.info("\t\t   Synbio parts: " + annot_string + "\n")
                cleared[region] = 1
            else:
                logging.info(
                    "\t\t -->Synbio sequences - <80% coverage achieved = FAIL\n"
                )

    for region in range(0, coords.shape[0]):
        if cleared[region] == 0:
            logging.info(
                "\t\t -->Regulated region at bases "
                + str(int(coords.iloc[region, 0]))
                + " to "
                + str(int(coords.iloc[region, 1]))
                + " failed to clear: FLAG\n"
            )
    if sum(cleared) == len(cleared):
        logging.info("\n\t\t -->all regulated regions cleared: PASS\n")

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
        logging.info("\t...no regulated regions to clear\n")
        return 0

    # Read database file
    coords = pd.read_csv(args.sample_name + ".reg_path_coords.csv")

    if coords.shape[0] == 0:
        logging.info("\t...no regulated regions to clear\n")
        return 0

    coords.sort_values(by=["q. start"], inplace=True)
    coords.reset_index(drop=True, inplace=True)
    rv = check_for_benign(args.sample_name, coords, benign_desc)
    sys.exit(rv)

if __name__ == "__main__":
    main()
