#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks results for regulated pathogen and prints any matched coordinates. Ignores any
synthetic constructs.

Usage:
  python check_reg_path.py -i INPUT -d database_folder -t threads

"""
import argparse
import logging
import os
import re
import sys
import textwrap
import pandas as pd
#from commec.config.screen_tools import ScreenIO
from commec.tools.blast_tools import read_blast, get_taxonomic_labels, get_top_hits
from commec.tools.blastn import BlastNHandler
from commec.tools.search_handler import SearchHandler
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

pd.set_option("display.max_colwidth", 10000)

logger = logging.getLogger(__name__)

def _check_inputs(
        search_handle : SearchHandler,
        benign_taxid_path : str | os.PathLike,
        biorisk_taxid_path : str | os.PathLike,
        taxonomy_directory : str | os.PathLike
        ):
    """ 
    Simply check for the existance of files, 
    returns True if it is safe to continue. 
    """
    # check input files
    if not search_handle.check_output():
        logger.info("\t...ERROR: Taxonomic search results empty\n %s", search_handle.out_file)
        return False

    if not os.path.exists(benign_taxid_path):
        logger.error("\t...benign db file %s does not exist\n", benign_taxid_path)
        return False

    if not os.path.exists(biorisk_taxid_path):
        logger.error("\t...biorisk db file %s does not exist\n", biorisk_taxid_path)
        return False
    
    if not os.path.exists(taxonomy_directory):
        logger.error("\t...taxonomy directory %s does not exist\n", taxonomy_directory)
        return False

    if search_handle.is_empty(search_handle.out_file):
        logger.info("\tERROR: Homology search has failed\n")
        return False
    
    return True

def update_taxonomic_data_from_database(
        search_handle : SearchHandler,
        benign_taxid_path : str | os.PathLike,
        biorisk_taxid_path : str | os.PathLike,
        taxonomy_directory : str | os.PathLike,
        data : ScreenResult,
        queries : dict[str, Query],
        step : ScreenStep,
        n_threads : int
        ):
    """
    Given a Taxonomic database screen output, update the screen data appropriately.
        search_handle : The handle of the search tool used to screen taxonomic data.
        benign_taxid_path : Path to benign taxid csv.
        biorisk_taxid_path : Path to regulated taxid csv.
        taxonomy_directory : The location of taxonom directory.
        data : the Screen data object, to be updated.
        step : Which taxonomic step this is (Nucleotide, Protein, etc)
        n_threads : maximum number of available threads for allocation.
    """
    logger.debug("Acquiring Taxonomic Data for JSON output:")

    if not _check_inputs(search_handle, benign_taxid_path, 
                         biorisk_taxid_path, taxonomy_directory):
        return 1

    # The default is to pass, its up to the data to over-write this.
    # Some Queries may already be set to skip, which we ignore.
    if step == ScreenStep.TAXONOMY_AA:
        for query in data.queries.values():
            if query.recommendation.protein_taxonomy_status != ScreenStatus.SKIP:
                query.recommendation.protein_taxonomy_status = ScreenStatus.PASS
    if step == ScreenStep.TAXONOMY_NT:
        for query in data.queries.values():
            if query.recommendation.nucleotide_taxonomy_status != ScreenStatus.SKIP:
                query.recommendation.nucleotide_taxonomy_status = ScreenStatus.PASS

    if not search_handle.has_hits(search_handle.out_file):
        logger.info("\t...no hits\n")
        return 0

    # We delay non-debug logging to sort messages via query.
    log_container = {key : [] for key in data.queries.keys()}

    # Read in lists of regulated and benign tax ids
    vax_taxids = pd.read_csv(benign_taxid_path, header=None).squeeze().astype(str).tolist()
    reg_taxids = pd.read_csv(biorisk_taxid_path, header=None).squeeze().astype(str).tolist()

    blast = read_blast(search_handle.out_file)
    logger.debug("%s Blast Import: shape: %s preview:\n%s", step, blast.shape, blast.head())

    blast = get_taxonomic_labels(blast, reg_taxids, vax_taxids, taxonomy_directory, n_threads)
    logger.debug("%s TaxLabels: shape: %s preview:\n%s", step, blast.shape, blast.head())

    blast = blast[blast["species"] != ""]  # ignore submissions made above the species level
    logger.debug("%s RemoveSpecies: shape: %s preview:\n%s", step, blast.shape, blast.head())

    # label each base with the top matching hit, but include different taxids attributed to same hit
    top_hits = get_top_hits(blast)
    logger.debug("%s Top Hits: shape: %s preview:\n%s", step, top_hits.shape, top_hits.head())

    if top_hits["regulated"].sum() == 0:
        logger.info("\t...no regulated hits\n")
        return 0
    
    # if ANY of the trimmed hits are regulated
    with pd.option_context('display.max_rows', None,
                    'display.max_columns', None,
                    'display.precision', 3,
                    ):

        unique_queries = top_hits['query acc.'].unique()
        logger.debug("%s Unique Queries: shape: %s preview:\n%s", step, unique_queries.shape, unique_queries)
        for query in unique_queries:
            logger.debug("\tProcessing query: %s", query)
            query_write = data.get_query(query)
            if not query_write:
                logger.error("Query during %s could not be found! [%s]", str(step), query)
                continue

            unique_query_data : pd.DataFrame = top_hits[top_hits['query acc.'] == query]
            unique_query_data.dropna(subset = ['species'])
            regulated_only_data = unique_query_data[unique_query_data["regulated"] == True]
            regulated_hits = regulated_only_data['subject acc.'].unique()
            logger.debug("\t%s Regulated hits: shape: %s preview:\n%s", step, regulated_hits.shape, regulated_hits)

            for hit in regulated_hits:
                logger.debug("\t\tProcessing Hit: %s", hit)
                regulated_hit_data : pd.DataFrame = regulated_only_data[regulated_only_data["subject acc."] == hit]
                logger.debug("%s Regulated Hit Data: shape: %s preview:\n%s", step, regulated_hit_data.shape, regulated_hit_data.head())
                hit_description = regulated_hit_data['subject title'].values[0]

                n_regulated_bacteria = 0
                n_regulated_virus = 0
                n_regulated_eukaryote = 0
                
                reg_taxids = [] # Regulated Taxonomy IDS
                non_reg_taxids = [] # Non-regulated Taxonomy IDS.
                reg_species = [] # List of species
                domains = [] # List of domains.
                match_ranges = [] # Ranges where hit matches query.

                for _, region in regulated_hit_data.iterrows():
                    match_range = MatchRange(
                        float(region['evalue']),
                        int(region['s. start']), int(region['s. end']),
                        int(region['q. start']), int(region['q. end'])
                    )
                    logger.debug("Processing region from hit: %s", region)
                    # Convert from non-coding to nt query coordinates if we're doing a NT taxonomy step.
                    if step == ScreenStep.TAXONOMY_NT:
                        match_range.query_start = queries[query].nc_to_nt_query_coords(match_range.query_start)
                        match_range.query_end = queries[query].nc_to_nt_query_coords(match_range.query_end)

                    match_ranges.append(match_range)

                    # Filter shared_site based on 'q. start' or 'q. end' (Previously only shared starts were used)
                    shared_site = top_hits[
                        (top_hits['q. start'] == region['q. start']) |
                        (top_hits['q. end'] == region['q. end'])
                        ]

                    # Filter for regulated and non-regulated entries
                    regulated = shared_site[shared_site["regulated"] == True]
                    non_regulated = shared_site[shared_site["regulated"] == False]

                    # Count domain information.
                    domain = region['superkingdom']
                    if domain == "Viruses":
                        n_regulated_virus += 1
                        logger.debug("\t\t\tAdded Virus.")
                    if domain == "Bacteria":
                        n_regulated_bacteria +=1
                        logger.debug("\t\t\tAdded Bacteria.")
                    if domain == "Eukaryota":
                        n_regulated_eukaryote+=1
                        logger.debug("\t\t\tAdded Eukaryote.")
                    domains.append(domain)

                    # Collect unique species from both regulated and non-regulated
                    reg_species.extend(regulated["species"])
                    # JSON serialization requires int, not np.int64, hence the map()
                    reg_taxids.extend(map(str, regulated["subject tax ids"]))
                    non_reg_taxids.extend(map(str, non_regulated["subject tax ids"]))

                    # These are the old values, now we simply count the size of the regulated, and non_regulated taxid arrays.
                    #n_reg += (top_hits["regulated"][top_hits['q. start'] == region['q. start']] != False).sum()
                    #n_total += len(top_hits["regulated"][top_hits['q. start'] == region['q. start']])

                # Uniquefy.
                reg_species = list(set(reg_species))
                reg_taxids = list(set(reg_taxids))
                non_reg_taxids = list(set(non_reg_taxids))
                match_ranges = list(set(match_ranges))

                reg_species_text = ", ".join(reg_species)
                reg_taxids_text = ", ".join(reg_taxids)
                non_reg_taxids_text = ", ".join(non_reg_taxids)
                match_ranges_text = ", ".join(map(str,match_ranges))
                domains_text = ", ".join(set(domains))

                reg_species.sort()
                reg_taxids.sort()
                non_reg_taxids.sort()

                logger.debug("\t\tRegulated Species: %s", reg_species)
                logger.debug("\t\tRegulated Taxids: %s", reg_taxids)
                logger.debug("\t\tNon Regulated Taxids: %s", non_reg_taxids)
                logger.debug("\t\tRanges: %s", match_ranges)

                recommendation : ScreenStatus = ScreenStatus.FLAG

                # TODO: Currently, we recapitulate old behaviour,
                # # " no top hit exclusive to a regulated pathogen: PASS"
                #  however in the future:
                # if all hits are in the same genus n_reg > 0, and n_total > n_reg, WARN, or other logic.
                # the point is, this is where you do it.

                logger.debug("Checking number of non regulated taxids: %i", len(non_reg_taxids))
                if len(non_reg_taxids) > 0:
                    logger.debug("Non-regulated taxids present, treating as MIXED result.")
                    recommendation = ScreenStatus.PASS

                # Update the query level recommendation of this step.
                if step == ScreenStep.TAXONOMY_AA:
                    query_write.recommendation.protein_taxonomy_status = compare(
                        query_write.recommendation.protein_taxonomy_status,
                        recommendation)
                if step == ScreenStep.TAXONOMY_NT:
                    query_write.recommendation.nucleotide_taxonomy_status = compare(
                        query_write.recommendation.nucleotide_taxonomy_status,
                        recommendation)

                regulation_dict = {"number_of_regulated_taxids" : str(len(reg_taxids)),
                                   "number_of_unregulated_taxids" : str(len(non_reg_taxids)),
                                   "regulated_eukaryotes": str(n_regulated_eukaryote),
                                   "regulated_bacteria": str(n_regulated_bacteria),
                                   "regulated_viruses": str(n_regulated_virus),
                                   "regulated_taxids": reg_taxids,
                                   "non_regulated_taxids" : non_reg_taxids,
                                   "regulated_species" : reg_species}

                # Logging logic.
                alt_text = "only " if recommendation == ScreenStatus.FLAG else "both regulated and non-"
                s = "" if len(reg_taxids) == 1 else "'s"
                log_message = (
                    f"\t --> {recommendation} at bases ({match_ranges_text}) found in {alt_text}regulated {domains_text}.\n"
                    f"\t   (Regulated Species: {reg_species_text}. Regulated TaxID{s}: {reg_taxids_text})"
                )
                logger.debug(log_message)
                log_container[query].append(log_message)

                # Append our hit information to Screen data.
                new_hit = HitResult(
                    HitScreenStatus(
                        recommendation,
                        step
                    ),
                    hit,
                    hit_description,
                    match_ranges,
                    {"domain" : [domain],"regulated_taxonomy":[regulation_dict]},
                )

                logger.debug("Hit information summary: %s", new_hit)

                if query_write.add_new_hit_information(new_hit):
                    write_hit = query_write.get_hit(hit)
                    if write_hit:
                        write_hit.annotations["domain"] = domains
                        write_hit.annotations["regulated_taxonomy"].append(regulation_dict)
                        write_hit.recommendation.status = compare(write_hit.recommendation.status, recommendation)
                        write_hit.description += ","+hit_description

    # Do all non-verbose logging in order of query:
    for query_name, log_list in log_container.items():
        if len(log_list) == 0:
            continue
        s = "" if len(log_list) == 1 else "s"
        taxtype = "protein" if step == ScreenStep.TAXONOMY_AA else "nucleotide"
        logger.info(f" Regulated {taxtype}{s} in {query_name}:")
        for log_text in log_list:
            logger.info(log_text)
    
    return 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="in_file",
        required=True,
        help="Input query file (e.g. QUERY.nr.dmnd)",
    )
    parser.add_argument(
        "-d",
        "--database",
        dest="db",
        required=True,
        help="top-level database folder (assumes /taxonomy, /benign_db, /biorisk_db dirs",
    )
    parser.add_argument("-t", "--threads", dest="threads", required=True, help="number of threads")
    args = parser.parse_args()

    # Legacy - assume hardcoded database locations; to adjust, call via screen.py
    input_database_dir = args.db
    taxonomy_db_path = f"{input_database_dir}/taxonomy/"
    benign_taxid_path = f"{input_database_dir}/benign_db/vax_taxids.txt"
    biorisk_taxid_path = f"{input_database_dir}/biorisk_db/reg_taxids.txt"

    exit_code = check_for_regulated_pathogens(
        args.in_file, taxonomy_db_path, benign_taxid_path, biorisk_taxid_path, args.threads
    )
    sys.exit(exit_code)


def check_for_regulated_pathogens(
        input_file: str | os.PathLike,
        taxonomy_path: str | os.PathLike,
        benign_taxid_path: str | os.PathLike,
        regulated_taxid_path: str | os.PathLike,
        threads: int
    ):
    """
    Check an input file (output from a database search) for regulated pathogens, from the benign and
    biorisk database taxids.
    """
    logger = logging.getLogger(__name__)

    # Check input file
    if not os.path.exists(input_file):
        logger.error("\t...input query file %s does not exist\n", input_file)
        return 1
    sample_name = re.sub(r"\.nr.*|\.nt\.blastn", "", input_file)

    # Read in lists of regulated and benign tax ids
    if not os.path.exists(benign_taxid_path):
        logger.error("\t...List of benign taxids %s does not exist\n", benign_taxid_path)
        return 1
    vax_taxids = pd.read_csv(benign_taxid_path, header=None).squeeze().astype(str).tolist()

    if not os.path.exists(regulated_taxid_path):
        logger.error("\t...List of regulated taxids %s does not exist\n", regulated_taxid_path)
        return 1
    reg_taxids = pd.read_csv(regulated_taxid_path, header=None).squeeze().astype(str).tolist()

    # if there are already regulated regions written to file for this query, add to them
    hits1 = None
    if os.path.exists(sample_name + ".reg_path_coords.csv"):
        hits1 = pd.read_csv(sample_name + ".reg_path_coords.csv")

    if BlastNHandler.is_empty(input_file):
        logger.error(
            "Cannot check for regulated pathogens in empty or non-existent file: %s\n",
            input_file,
        )
        return 1

    if not BlastNHandler.has_hits(input_file):
        logger.info("\t... Skipping regulated pathogens check, no hits in: %s\n", input_file)
        return 0

    blast = read_blast(input_file)
    blast = get_taxonomic_labels(blast, reg_taxids, vax_taxids, taxonomy_path, threads)
    blast = blast[blast["species"] != ""]  # ignore submissions made above the species level

    # label each base with the top matching hit, but include different taxids attributed to same hit
    blast2 = get_top_hits(blast)

    reg_bac = 0
    reg_vir = 0
    reg_fung = 0

    # if this is the nucleotide screen, check if any weak protein flags can be negated with strong non-regulated nt ones
    if re.findall(".nt.blastn", input_file):
        if hits1 is not None:
            for region in range(0, hits1.shape[0]):  # for each regulated pathogen region
                # look at only the hits that overlap it
                htrim = blast2[
                    ~(
                        (blast2["q. start"] > hits1["q. end"][region])
                        & (blast2["q. end"] > hits1["q. end"][region])
                    )
                    & ~(
                        (blast2["q. start"] < hits1["q. start"][region])
                        & (blast2["q. end"] < hits1["q. start"][region])
                    )
                ]
                species_list = textwrap.fill(", ".join(set(htrim["species"])), 100).replace(
                    "\n", "\n\t\t     "
                )
                taxid_list = textwrap.fill(
                    ", ".join(map(str, set(htrim["subject tax ids"]))), 100
                ).replace("\n", "\n\t\t     ")
                percent_ids = " ".join(map(str, set(htrim["% identity"])))
                if htrim.shape[0] > 0:
                    if any(htrim["q. coverage"] > 0.90):
                        htrim = htrim[htrim["q. coverage"] > 0.90]
                        htrim = htrim.reset_index(drop=True)
                        descriptions = []
                        for row in range(htrim.shape[0]):
                            hit = htrim["subject title"][row]
                            descriptions.append(hit)
                        annot_string = "\n\t...".join(str(v) for v in descriptions)
                        logger.info(
                            "\t...Regulated protein region at bases "
                            + str(int(hits1["q. start"][region]))
                            + " to "
                            + str(int(hits1["q. end"][region]))
                            + " overlapped with a nucleotide hit\n"
                        )
                        logger.info(
                            "\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n",
                            species_list,
                            taxid_list,
                            percent_ids,
                        )

    if blast2["regulated"].sum():  # if ANY of the trimmed hits are regulated
        hits = pd.DataFrame(columns=["q. start", "q. end"])
        with pd.option_context(
            "display.max_rows",
            None,
            "display.max_columns",
            None,
            "display.precision",
            3,
        ):
            # for each hit (subject acc) linked with at least one regulated taxid
            for site in set(blast2["q. start"][blast2["regulated"] != False]):
                subset = blast2[(blast2["q. start"] == site)]
                subset = subset.sort_values(by=["regulated"], ascending=False)
                subset = subset.reset_index(drop=True)
                org = ""

                blast2 = blast2.dropna(subset=["species"])
                n_reg = (blast2["regulated"][blast2["q. start"] == site] != False).sum()
                n_total = len(blast2["regulated"][blast2["q. start"] == site])
                gene_names = ", ".join(set(subset["subject acc."]))
                end = blast2["q. end"][blast2["q. start"] == site].max()
                coordinates = str(int(site)) + " - " + str(int(end))

                species_list = textwrap.fill(
                    ", ".join(set(blast2["species"][blast2["q. start"] == site])), 100
                ).replace("\n", "\n\t\t     ")
                desc = blast2["subject title"][blast2["q. start"] == site].values[0]
                taxid_list = textwrap.fill(
                    ", ".join(
                        map(
                            str,
                            set(blast2["subject tax ids"][blast2["q. start"] == site]),
                        )
                    ),
                    100,
                ).replace("\n", "\n\t\t     ")
                percent_ids = " ".join(
                    map(str, set(blast2["% identity"][blast2["q. start"] == site]))
                )
                reg_ids = " ".join(
                    map(
                        str,
                        set(
                            blast2["regulated"][
                                (blast2["q. start"] == site) & (blast2["regulated"] != False)
                            ]
                        ),
                    )
                )

                # if some of the organisms with this sequence aren't regulated, say so
                if n_reg < n_total:
                    logger.info(
                        "\t\t --> Best match to sequence(s) %s at bases %s found in both regulated and non-regulated organisms\n"
                        % (gene_names, coordinates)
                    )
                    logger.info(
                        "\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n"
                        % (species_list, taxid_list, percent_ids)
                    )
                    logger.info("\t\t     Description: %s\n" % (desc))
                    # could explicitly list which are and aren't regulated?
                # otherwise, raise a flag and say which superkingdom the flag belongs to
                elif n_reg == n_total:
                    if subset["superkingdom"][0] == "Viruses":
                        reg_vir = 1
                        org = "virus"
                    elif subset["superkingdom"][0] == "Bacteria":
                        reg_bac = 1
                        org = "bacteria"
                    elif "superkingdom" in subset:
                        if subset["superkingdom"][0] == "Eukaryota":
                            org = "eukaryote"
                            reg_fung = 1

                    new_hits = subset[["q. start", "q. end"]].dropna()
                    if not new_hits.empty and not hits.empty:
                        hits = pd.concat([hits, new_hits], ignore_index=True)
                    elif not new_hits.empty:
                        hits = new_hits.copy()
                    logger.info(
                        "\t\t --> Best match to sequence(s) %s at bases %s found in only regulated organisms: FLAG (%s)\n"
                        % (gene_names, coordinates, org)
                    )
                    logger.info(
                        "\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n"
                        % (species_list, taxid_list, percent_ids)
                    )
                    logger.info("\t\t     Description: %s\n" % (desc))
                else:  # something is wrong, n_reg > n_total
                    logger.info("\t...gene: %s\n" % gene_names)
                    logger.info("%s\n" % (blast["regulated"][blast["subject acc."] == gene_names]))
        hits = hits.drop_duplicates()
        # Create output file
        if hits1 is not None:
            hits = pd.concat([hits1, hits])
        hits.to_csv(sample_name + ".reg_path_coords.csv", index=False)

    if reg_vir == 0 and reg_bac == 0 and reg_fung == 0 and reg_fung == 0:
        logger.info("\t\t --> no top hit exclusive to a regulated pathogen: PASS\n")

    return 0


if __name__ == "__main__":
    main()