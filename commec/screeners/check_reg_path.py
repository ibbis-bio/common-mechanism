#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks results for regulated pathogen and prints any matched coordinates. Ignores any
synthetic constructs.

Usage:
  python check_reg_path.py -i INPUT -d database_folder -t threads

"""
import logging
import os
import pandas as pd
from commec.tools.search_handler import SearchHandler
from commec.config.query import Query
from commec.tools.blast_tools import (
    read_blast,
    get_taxonomic_labels,
    get_top_hits
)
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
        search_handler : SearchHandler,
        low_concern_taxid_path : str | os.PathLike,
        biorisk_taxid_path : str | os.PathLike,
        taxonomy_directory : str | os.PathLike
        ):
    """ 
    Simply check for the existance of files, 
    returns True if it is safe to continue. 
    """
    # check input files
    if not search_handler.validate_output():
        logger.info("\t...ERROR: Taxonomic search results empty\n %s", search_handler.out_file)
        return False

    if not os.path.exists(low_concern_taxid_path):
        logger.error("\t...low-concern database file %s does not exist\n", low_concern_taxid_path)
        return False

    if not os.path.exists(biorisk_taxid_path):
        logger.error("\t...biorisk database file %s does not exist\n", biorisk_taxid_path)
        return False

    if not os.path.exists(taxonomy_directory):
        logger.error("\t...taxonomy directory %s does not exist\n", taxonomy_directory)
        return False

    return True

def parse_taxonomy_hits(
        search_handler : SearchHandler,
        low_concern_taxid_path : str | os.PathLike,
        biorisk_taxid_path : str | os.PathLike,
        taxonomy_directory : str | os.PathLike,
        data : ScreenResult,
        queries : dict[str, Query],
        step : ScreenStep,
        n_threads : int
        ):
    """
    Given a Taxonomic database screen output, update the screen data appropriately.
        search_handler : The handle of the search tool used to screen taxonomic data.
        low_concern_taxid_path : Path to low-concern taxid csv.
        biorisk_taxid_path : Path to regulated taxid csv.
        taxonomy_directory : The location of taxonomy directory.
        data : the Screen data object, to be updated.
        step : Which taxonomic step this is (Nucleotide, Protein, etc)
        n_threads : maximum number of available threads for allocation.
    """
    logger.debug("Acquiring Taxonomic Data for JSON output:")

    if not _check_inputs(search_handler, low_concern_taxid_path,
                         biorisk_taxid_path, taxonomy_directory):
        return 1

    # The default is to pass, its up to the data to over-write this.
    for query in data.queries.values():
        query.status.set_step_status(step, ScreenStatus.PASS)

    if not search_handler.has_hits():
        logger.info("\t...no hits\n")
        return 0

    # We delay non-debug logging to sort messages via query.
    log_container = {key : [] for key in data.queries.keys()}

    # Read in lists of regulated and low_concern tax ids
    vax_taxids = pd.read_csv(low_concern_taxid_path, header=None).squeeze().astype(str).tolist()
    reg_taxids = pd.read_csv(biorisk_taxid_path, header=None).squeeze().astype(str).tolist()

    blast = read_blast(search_handler.out_file)
    logger.debug("%s Blast Import: shape: %s preview:\n%s", step, blast.shape, blast.head())

    # Initial check for query to be identified as anything known.
    unique_queries = blast['query acc.'].unique()
    for query_acc in unique_queries:
        query_obj = queries.get(query_acc)
        if query_obj:
            logger.debug("Confirming hits for query %s.", query_acc)
            query_obj.confirm_has_hits()
        else:
            logger.error("Could not mark query %s for confirmation of hit, "
                            "query not found in input queries.", query_acc)

    # Add taxonomic labels, and filter synthetic constructs
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

                screen_status : ScreenStatus = ScreenStatus.FLAG

                # TODO: Currently, we recapitulate old behaviour,
                # # " no top hit exclusive to a regulated pathogen: PASS"
                #  however in the future:
                # if all hits are in the same genus n_reg > 0, and n_total > n_reg, WARN, or other logic.
                # the point is, this is where you do it.

                logger.debug("Checking number of non regulated taxids: %i", len(non_reg_taxids))
                if len(non_reg_taxids) > 0:
                    logger.debug("Non-regulated taxids present, treating as MIXED result.")
                    screen_status = ScreenStatus.PASS

                # Update the query level recommendation of this step.
                query_write.status.update_step_status(step, screen_status)

                regulation_dict = {"number_of_regulated_taxids" : str(len(reg_taxids)),
                                   "number_of_unregulated_taxids" : str(len(non_reg_taxids)),
                                   "regulated_eukaryotes": str(n_regulated_eukaryote),
                                   "regulated_bacteria": str(n_regulated_bacteria),
                                   "regulated_viruses": str(n_regulated_virus),
                                   "regulated_taxids": reg_taxids,
                                   "non_regulated_taxids" : non_reg_taxids,
                                   "regulated_species" : reg_species}

                # Logging logic.
                alt_text = "only " if screen_status == ScreenStatus.FLAG else "both regulated and non-"
                s = "" if len(reg_taxids) == 1 else "'s"
                log_message = (
                    f"\t --> {screen_status} at bases ({match_ranges_text}) found in {alt_text}regulated {domains_text}.\n"
                    f"\t   (Regulated Species: {reg_species_text}. Regulated TaxID{s}: {reg_taxids_text})"
                )
                logger.debug(log_message)
                log_container[query].append(log_message)

                # Append our hit information to Screen data.
                new_hit = HitResult(
                    HitScreenStatus(
                        screen_status,
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
                        write_hit.recommendation.status = compare(write_hit.recommendation.status, screen_status)
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
