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
from dataclasses import asdict

import pandas as pd

from commec.tools.search_handler import SearchHandler
from commec.config.query import Query
from commec.tools.blast_tools import (
    read_blast,
    get_taxonomic_labels,
    get_controlled_labels,
    get_top_hits
)
from commec.config.result import (
    ScreenResult,
    TaxonomyAnnotation,
    HitResult,
    ScreenStep,
    ScreenStatus,
    HitScreenStatus,
    MatchRange,
    compare
)

from commec.control_list import (
    get_control_lists,
    get_regulation,
    ListMode,
    ControlList,
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

    def unique_annotation_set(df) -> set[TaxonomyAnnotation]:
        """
        Convenience function to extract taxonomy annotations into a 
        common structure, from BLAST results with the appropriate column headings.
        """
        return {
            TaxonomyAnnotation(*row)
            for row in df[[
                "evalue",
                "subject tax ids",
                #"species",
                #"genus",
                #"superkingdom",
                "subject acc.",
                "subject title",
            ]].itertuples(index=False, name=None)
        }
    
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
            query_obj.mark_as_hit()
        else:
            logger.error("Could not mark query %s for confirmation of hit, "
                            "query not found in input queries.", query_acc)

    # Add taxonomic labels, and filter synthetic constructs
    #blast = get_taxonomic_labels(blast, reg_taxids, vax_taxids, taxonomy_directory, n_threads)
    #logger.debug("%s TaxLabels: shape: %s preview:\n%s", step, blast.shape, blast.head())

    #blast = blast[blast["species"] != ""]  # ignore submissions made above the species level
    #logger.debug("%s RemoveSpecies: shape: %s preview:\n%s", step, blast.shape, blast.head())

    blast = get_controlled_labels(blast)

    # label each base with the top matching hit, but include different taxids attributed to same hit
    top_hits = get_top_hits(blast)
    logger.debug("%s Top Hits: shape: %s preview:\n%s", step, top_hits.shape, top_hits.head())

    if top_hits["regulated"].sum() == 0:
        logger.info("\t...no regulated hits\n")
        return 0

    unique_queries = top_hits['query acc.'].unique()
    logger.debug("%s Unique Queries: shape: %s preview:\n%s", step, unique_queries.shape, unique_queries)
    for query in unique_queries:
        logger.debug("\tProcessing query: %s", query)
        query_write = data.get_query(query)
        if not query_write:
            logger.error("Query during %s could not be found! [%s]", str(step), query)
            continue

        unique_query_data : pd.DataFrame = top_hits[top_hits['query acc.'] == query]
        #unique_query_data.dropna(subset = ['species'])
        regulated_only_data = unique_query_data[unique_query_data["regulated"] == True]
        regulated_hits = regulated_only_data['subject acc.'].unique()
        logger.debug("\t%s Regulated hits: shape: %s preview:\n%s",
                     step, regulated_hits.shape, regulated_hits)

        for hit in regulated_hits:
            logger.debug("\t\tProcessing Hit: %s", hit)
            regulated_hit_data : pd.DataFrame = regulated_only_data[regulated_only_data["subject acc."] == hit]
            logger.debug("%s Regulated Hit Data: shape: %s preview:\n%s",
                         step, regulated_hit_data.shape, regulated_hit_data.head())

            n_regulated_bacteria = 0
            n_regulated_virus = 0
            n_regulated_eukaryote = 0
            
            reg_taxids = [] # Regulated Taxonomy IDS
            non_reg_taxids = [] # Non-regulated Taxonomy IDS.
            reg_species = [] # List of species
            non_reg_species = [] # List of species
            domains = [] # List of domains.
            match_ranges = [] # Ranges where hit matches query.

            regulated_annotations = set()
            non_regulated_annotations = set()

            for _, region in regulated_hit_data.iterrows():
                # Record region information:
                match_range = MatchRange(
                    float(region['evalue']),
                    int(region['s. start']), int(region['s. end']),
                    int(region['q. start']), int(region['q. end'])
                )

                # Convert from non-coding to nt query coordinates if we're doing a NT taxonomy step.
                if step == ScreenStep.TAXONOMY_NT:
                    match_range.query_start = queries[query].nc_to_nt_query_coords(match_range.query_start)
                    match_range.query_end = queries[query].nc_to_nt_query_coords(match_range.query_end)
                match_ranges.append(match_range)
                logger.debug("Processing region from hit: %s", region)

                # Filter shared_site based on 'q. start' or 'q. end'
                # (Previously only shared starts were used)
                shared_site = unique_query_data[
                    (unique_query_data['q. start'] == region['q. start']) |
                    (unique_query_data['q. end'] == region['q. end'])
                    ]

                # Filter for regulated and non-regulated entries
                regulated_for_region = shared_site[shared_site["regulated"] == True]
                non_regulated_for_region = (
                    shared_site[shared_site["regulated"] == False]
                    .sort_values(by="evalue", ascending=True)
                    .head(10) # we only care for max 10 non-regulated.
                )

                # Optimise for conciseness by unique TaxID
                regulated_for_region = regulated_for_region.drop_duplicates(subset=["subject tax ids"], keep="first")
                non_regulated_for_region = non_regulated_for_region.drop_duplicates(subset=["subject tax ids"], keep="first")

                # Update what is a regulated region based on imported control lists.
                regulated_for_region, non_regulated_for_region = update_using_control_lists(
                    regulated_for_region, non_regulated_for_region
                )

                # Collect unique species from both regulated and non-regulated - legacy logging
                #reg_species.extend(regulated_for_region["species"])
                #non_reg_species.extend(non_regulated_for_region["species"])

                # JSON serialization requires int, not np.int64, hence the map()
                reg_taxids.extend(map(str, regulated_for_region["subject tax ids"]))
                non_reg_taxids.extend(map(str, non_regulated_for_region["subject tax ids"]))

                regulated_annotations = regulated_annotations | unique_annotation_set(regulated_for_region)
                non_regulated_annotations = non_regulated_annotations | unique_annotation_set(non_regulated_for_region)

            # Convert our structures to a dictionary for JSON export, sorted by taxid.
            regulated_annotation_list = [asdict(t) for t in regulated_annotations]
            non_regulated_annotation_list = [asdict(t) for t in non_regulated_annotations]
            regulated_annotation_list = sorted(
                regulated_annotation_list,
                key=lambda d: d["taxid"]
            )
            non_regulated_annotation_list = sorted(
                non_regulated_annotation_list,
                key=lambda d: d["taxid"]
            )

            for reg_annotation in regulated_annotation_list:
                control_info, _context_info = get_regulation(reg_annotation["taxid"])
                reg_annotation["control_list"] = control_info

                for control_output_info in control_info:
                    # Record domain information.
                    control_output_info.list = str(get_control_lists(control_output_info.list))
                    domain = control_output_info.category
                    if domain == "Viruses":
                        n_regulated_virus += 1
                        logger.debug("\t\t\tAdded Virus.")
                        break
                    if domain == "Bacteria":
                        n_regulated_bacteria +=1
                        logger.debug("\t\t\tAdded Bacteria.")
                        break
                    if domain == "Eukaryota":
                        n_regulated_eukaryote+=1
                        logger.debug("\t\t\tAdded Eukaryote.")
                        break
                    domains.append(domain)

            # Useful for when a single conditional control list compliance occured.
            for nonreg_annotation in non_regulated_annotation_list:
                control_info, _context_info = get_regulation(nonreg_annotation["taxid"])
                for control_output_info in control_info:
                    # Record domain information.
                    control_output_info.list = str(get_control_lists(control_output_info.list))
                if len(control_info) > 0:
                    nonreg_annotation["control_list"] = control_info



            # If the category in the control list isn't set, we default to the non-descript, "sequences".
            domains_text = ", ".join(set(domains)) or "sequence"
            # Set the default hit description, this is changed if result is mixed etc.
            hit_description = f"Regulated {domains_text} - {regulated_hit_data['subject title'].values[0]}"

            # TODO: Currently, we recapitulate old behaviour,
            # # " no top hit exclusive to a regulated pathogen: PASS"
            #  however in the future:
            # if all hits are in the same genus n_reg > 0, and n_total > n_reg, WARN, or other logic.
            # the point is, this is where you do it.

            screen_status : ScreenStatus = ScreenStatus.FLAG # Default is to flag.

            logger.debug("Checking number of non regulated taxids: %i", len(non_regulated_annotation_list))
            if len(non_regulated_annotation_list) > 0:
                logger.debug("Non-regulated taxids present, treating as MIXED result.")
                screen_status = ScreenStatus.PASS
                hit_description = (f"Mix of {len(regulated_annotation_list)} controlled {domains_text}"
                f" and {len(non_regulated_annotation_list)} non-regulated {domains_text}")

            # We might have 0 regulated annotations, due to removal based on regional context:
            if len(regulated_annotation_list) == 0:
                screen_status = ScreenStatus.PASS
                logger.debug("Only non-controlled entities due to control list compliance regional context")
                hit_description = (f"Externally controlled {domains_text}")
            
            # Update the query level recommendation of this step.
            query_write.status.update_step_status(step, screen_status)

            regulation_dict = {"number_of_regulated_taxids" : str(len(regulated_annotation_list)),
                                "number_of_unregulated_taxids" : str(len(non_regulated_annotation_list)),
                                "regulated_eukaryotes": str(n_regulated_eukaryote),
                                "regulated_bacteria": str(n_regulated_bacteria),
                                "regulated_viruses": str(n_regulated_virus),
                                "regulated_taxa": regulated_annotation_list,
                                "non_regulated_taxa" : non_regulated_annotation_list}
            
            # Ensure each match range is unique.
            match_ranges = list(set(match_ranges))

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

            if query_write.add_new_hit_information(new_hit):
                write_hit = query_write.get_hit(hit)
                if write_hit:
                    write_hit.ranges.extend(match_ranges)
                    write_hit.annotations["domain"] = domains
                    write_hit.annotations["regulated_taxonomy"].append(regulation_dict)
                    write_hit.recommendation.status = compare(write_hit.recommendation.status, screen_status)
                    write_hit.description += ","+hit_description

            logger.debug("Hit information summary: %s", new_hit)

            # Logging logic - somewhat convolutedly placed but for the intention of recreating 
            # legacy-like .screen.log file logging experience.
            reg_species = list(set(reg_species))
            reg_taxids = list(set(reg_taxids))
            non_reg_taxids = list(set(non_reg_taxids))

            reg_species.sort()
            reg_taxids.sort()
            non_reg_taxids.sort()
            
            reg_species_text = ", ".join(reg_species)
            non_reg_species_text = ", ".join(non_reg_species)
            reg_taxids_text = ", ".join(reg_taxids)
            non_reg_taxids_text = ", ".join(non_reg_taxids)
            match_ranges_text = ", ".join(map(str,match_ranges))

            alt_text = "only " if screen_status == ScreenStatus.FLAG else "both regulated and non-" if len(reg_taxids) > 0 else "externally "
            s = "" if len(reg_taxids) == 1 else "'s"
            ss = "" if len(non_reg_taxids) == 1 else "'s"
            log_message = (
                f"\t --> {screen_status} at bases ({match_ranges_text}) found in {alt_text}regulated {domains_text}.\n"
                f"\t   (Regulated Species: {reg_species_text or non_reg_species_text}.\n\t    Regulated TaxID{s}: {reg_taxids_text}\n"
                f"\t   Non-Regulated TaxID{ss}: {non_reg_taxids_text})"
            )
            if len(non_reg_taxids) > 0:
                log_message += f"\n\t    Non-Regulated TaxID{ss}: {non_reg_taxids_text}"
            log_message += ")"
            logger.debug(log_message)
            logger.debug("\t\tRegulated Species: %s", reg_species)
            logger.debug("\t\tRegulated Taxids: %s", regulated_annotation_list)
            logger.debug("\t\tNon Regulated Taxids: %s", non_regulated_annotation_list)
            logger.debug("\t\tRanges: %s", match_ranges)

            log_container[query].append(log_message)


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


def update_using_control_lists(regulated : pd.DataFrame, non_regulated : pd.DataFrame):
    """
    Takes dataframes for both of regulated, and non_regulated candidates.
    For each regulated candidate - check its list compliance, and only keep it
    if it matches full compliance, or conditional compliance.

    If not, push it to the non-regulated lists.

    If a single CONDITIONAL is found, it is moved to non-regulated.

    """

    CONDITIONAL_NUMBER_TO_ALLOW = 2 # i.e. if something appears in 2 or more control lists.
    indexes_to_move = []

    for index, row in regulated.iterrows():
        control_data, _context_data = get_regulation(row["subject tax ids"])
        for info in control_data:
            clist : ControlList = get_control_lists(info.list)
            regulated.at[index, "list_acronym"] = clist.acronym
            regulated.at[index, "category"] = info.category
            if clist.status == ListMode.COMPLIANCE:
                continue
            if clist.status == ListMode.COMPLIANCE_WARN:
                continue
            if clist.status == ListMode.CONDITIONAL_NUM:
                indexes_to_move.append(index)

    # If more than 1 index is removed, then treat as flag.
    if len(indexes_to_move) >= CONDITIONAL_NUMBER_TO_ALLOW:
        indexes_to_move = []

    # Move the out of context rows to the non-regulated DataFrame
    if indexes_to_move:
        rows_to_move = regulated.loc[indexes_to_move]
        non_regulated = pd.concat([non_regulated, rows_to_move], ignore_index=True)
        regulated = regulated.drop(index=indexes_to_move)

    return regulated, non_regulated

