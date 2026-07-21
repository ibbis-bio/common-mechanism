#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for parsing taxonomic BLAST screen output and annotating regulated pathogen hits.

Provides ``parse_taxonomy_hits`` as the public entry point, which validates inputs,
loads and labels BLAST results, applies control-list compliance rules, and updates
a ``ScreenResult`` in-place with ``HitResult`` objects for each regulated hit found.
"""

import logging

import pandas as pd

from commec.tools.search_handler import SearchHandler
from commec.config.query import Query
from commec.tools.blast_tools import (
    read_blast,
    get_controlled_labels,
    find_clusters,
    get_top_hits,
)
from commec.config.result import (
    ScreenResult,
    HitResult,
    ScreenStep,
    ScreenStatus,
    HitScreenStatus,
    MatchRange,
)
from commec.control_list import (
    get_control_lists,
    get_regulation,
    is_regulated,
    ListMode,
)

from commec.config.constants import N_NON_REGIONAL_HITS_TO_WARN, BAD_ACCESSIONS

pd.set_option("display.max_colwidth", 10000)

logger = logging.getLogger(__name__)


# ── Public entry point ─────────────────────────────────────────────────────────


def parse_taxonomy_hits(
    search_handler: SearchHandler,
    data: ScreenResult,
    queries: dict[str, Query],
    step: ScreenStep,
    n_threads: int,
) -> int:
    """
    Parse a taxonomic BLAST screen output and update the ScreenResult in-place.

    Process:
        1. Load in the dataframe from provided SearchHandler.db_directory
        2. Annotate the dataframe with control list annotations based on Accession (TaxID)
            Includes removing synthetic or vaccine or other ingored taxa.
            Parsing multi-species records via label, or multiple TaxID entries.
            Provides a cluster ID for each hit based on its Control List connections
        3. Split the dataframe into per Query processing.
         per Query:
        4. Filter to Top Hits and Trim Edges.
        4. Split into Controlled vs non-controlled.
        5. Creates clusters of controlled non-overlapping regions.
        6. Remove any identically described hits based on subject title.
        7. Generate a HitResults result object containing
            > Longest display_name for controlled hit.
            > ScreenResult
            > Taxonomy annotations containing:
                > Summary statistics
                > List of Hits Controlled
                > List of Hits non-regulated (Top 10 non-regulated hits).
    Args:
        search_handler: Handle of the search tool used for the taxonomic screen.
        data: ScreenResult to be updated in-place.
        queries: Mapping from query accession to Query objects.
        step: Taxonomy step (TAXONOMY_NT or TAXONOMY_AA).
        n_threads: Maximum threads available (reserved for future use).

    Returns:
        0 on success, 1 on unrecoverable input error.
    """
    logger.debug("Acquiring Taxonomic Data for JSON output:")

    if not search_handler.validate_output():
        logger.info(
            "\t...ERROR: Taxonomic search results empty\n %s", search_handler.out_file
        )
        return 1

    # Default all queries to PASS; downstream hits will override where needed.
    for query in data.queries.values():
        query.status.set_step_status(step, ScreenStatus.PASS)

    if not search_handler.has_hits():
        logger.info("\t...no hits\n")
        return 0

    # Data load and preparation.
    blast = read_blast(search_handler.out_file)
    logger.debug("%s Blast Import: shape %s\n%s", step, blast.shape, blast.head())

    # Remove very bad Accessions
    blast = blast[~blast["subject acc."].isin(BAD_ACCESSIONS)]

    # Label data with control information.
    logger.info("    Identifying controlled taxa ...")
    blast = get_controlled_labels(blast)

    # Early exit if no regulated hits.
    if blast["regulated"].sum() == 0:
        logger.info("\t...no controlled hits\n")
        return 0

    logger.debug(
        "%s Controlled Labels applied: shape %s\n%s", step, blast.shape, blast.head()
    )

    # Step 2: Per-query analysis
    unique_query_accs = blast["query acc."].unique()
    logger.debug(
        "%s: %d unique queries with controlled hits", step, len(unique_query_accs)
    )

    # Each Query should be threaded:
    for query_acc in unique_query_accs:
        query_write, query_name = data.get_query(query_acc)
        if not query_write:
            logger.error(
                "Query '%s' not found in ScreenResult during %s.", query_acc, step
            )
            continue

        query_info = queries[query_name]
        logger.info("    Processing query: %s", query_name)

        # Filter to just top hits, correct for NT coords, trim the edges, sort and Clean up
        unique_query_data = blast[blast["query acc."] == query_acc]

        if unique_query_data.empty:
            logger.info("      --> No hits.")
            continue

        # If NT Taxonomy, we have some additional parsing to consider...
        if step == ScreenStep.TAXONOMY_NT:
            nc_id = int(query_acc.split("_")[-1])
            unique_query_data["q. start"] = unique_query_data["q. start"].apply(
                query_info.nc_to_nt_query_coords, index=nc_id
            )
            unique_query_data["q. end"] = unique_query_data["q. end"].apply(
                query_info.nc_to_nt_query_coords, index=nc_id
            )
            start, stop = query_info.non_coding_regions[nc_id]
            logger.info("    for the non-coding region: %s-%s", start, stop)

        unique_query_data = unique_query_data.sort_values(
            by=["% identity"], ascending=False
        )
        unique_query_data = unique_query_data.reset_index(drop=True)

        hit_results_for_query, logs = _get_hit_result_from_data(unique_query_data, step)

        # After thread is finished:
        for new_hit in hit_results_for_query:
            logger.debug("Adding hit for query %s : %s", query_acc, new_hit)
            query_write.add_new_hit_information(new_hit)

        for log in logs:
            logger.info(log)

    return 0


def _get_hit_result_from_data(
    unique_query_data: pd.DataFrame, step: ScreenStep
) -> list[HitResult]:
    """
    Given a dataframe that has been annotated with control list information,
    derive a list of HitResults that represent the outcome of this data.
    Data submitted is taken as a whole, and is assumed to represent a single Query.
    """
    output_hit_results: list[HitResult] = []

    query_acc = unique_query_data["query acc."].to_list()[0]

    # Filter data to the top hits only.
    top_hits = get_top_hits(unique_query_data)
    # top_hits = get_controlled_labels(top_hits)

    regulated_only = top_hits[top_hits["regulated"]]
    non_regulated_only = top_hits[~top_hits["regulated"]]
    unique_hash_taxids = regulated_only["control_hash"].unique()

    log_messages = []
    for hash_taxid in unique_hash_taxids:
        __regulated = regulated_only[regulated_only["control_hash"] == hash_taxid]
        if __regulated.empty:
            logger.error(
                "Empty dataframe being sent to find_clusters... for hash taxid %s",
                hash_taxid,
            )
        cluster_data, clusters = find_clusters(
            __regulated
        )  # [regulated_only["control_hash"] == hash_taxid]
        logger.debug(
            "%s: query %s has %i regulated clusters.", step, query_acc, len(clusters)
        )
        for i, cluster in enumerate(clusters):
            logger.debug("Processing cluster[%i]: %s", i, cluster)
            new_hit_result, log_message = _create_hit_result_for_cluster(
                cluster_data[cluster_data["cluster"] == i],
                non_regulated_only,
                cluster[0],
                cluster[1],
                hash_taxid,
                step,
            )

            if new_hit_result:
                output_hit_results.append(new_hit_result)
                log_messages.append(log_message)

    return output_hit_results, log_messages


def _create_hit_result_for_cluster(
    cluster_data: pd.DataFrame,
    non_regulated_only: pd.DataFrame,
    cluster_start: int,
    cluster_end: int,
    controlled_cluster_taxid: str,
    step: ScreenStep,
) -> tuple[HitResult | None, str]:
    """
    Given a cluster of controlled hits, and the common start and end point.
    Find overlapping non-regulated hits, and generate a HitResult, to describe
    the data for commec screen output.
    """
    logger.debug("Retrieving Control Information for Cluster: ")
    cluster_regulation = get_regulation(controlled_cluster_taxid)

    if not cluster_regulation:
        logger.error(
            "Tried to create a cluster from %s, but returned no control information. raw data: %s",
            controlled_cluster_taxid,
            cluster_data["subject tax ids"].unique(),
        )
        if check := is_regulated(controlled_cluster_taxid):
            logger.error(
                "Note: Status from is_regulated() is correctly identified as %s", check
            )
        return HitResult(
            HitScreenStatus(ScreenStatus.ERROR, step),
            f"Bad_Cluster_{controlled_cluster_taxid}",
        ), "Error creating hit from cluster ID"

    # species = cluster_regulation[0].species
    # genus = cluster_regulation[0].genus
    display_name = cluster_regulation[0].name

    # We get rid of the rare, but sometimes present duplicate subject titles within a cluster.
    cluster_data = cluster_data.drop_duplicates(subset=["subject title"], keep="first")

    best_evalue = min(cluster_data["evalue"])
    match_range = MatchRange(
        best_evalue,
        int(cluster_start),
        int(cluster_end),
    )

    regulated_for_region = cluster_data

    # Get the top 10 non-controlled hits for reporting.
    if non_regulated_only.empty:
        non_regulated_for_region = non_regulated_only
    else:
        shared_site = non_regulated_only[
            (non_regulated_only["q. start"] == cluster_start)
            | (non_regulated_only["q. end"] == cluster_end)
        ]
        non_regulated_for_region = (
            shared_site[~shared_site["regulated"]]
            .drop_duplicates(subset="subject tax ids")
            .drop_duplicates(subset="subject acc.")
            .sort_values(by="evalue", ascending=True)
            .head(10)
        )

    regulated_annotations = []
    regionally_regulated_annotations = []
    warn_regulated_annotations = []
    non_regulated_annotations = []

    # Gather control list information
    compliances = []
    conditional_compliances = []
    warn_compliances = []

    # For each uniquely controlled taxid, create the annotations according to list mode.
    for candidate_taxid in regulated_for_region["subject tax ids"].unique():
        control_output = get_regulation(candidate_taxid)
        data = regulated_for_region[
            regulated_for_region["subject tax ids"] == candidate_taxid
        ]

        compliances.extend(
            [
                output
                for output in control_output
                if get_control_lists(output.list).status == ListMode.COMPLIANCE
            ]
        )
        conditional_compliances.extend(
            [
                output
                for output in control_output
                if get_control_lists(output.list).status == ListMode.CONDITIONAL_NUM
            ]
        )
        warn_compliances.extend(
            [
                output
                for output in control_output
                if get_control_lists(output.list).status == ListMode.COMPLIANCE_WARN
            ]
        )

        logger.debug(
            "Parsed %s: [%i Control outputs, %i Compliances, %i Conditional, %i Warns]",
            candidate_taxid,
            len(control_output),
            len(compliances),
            len(conditional_compliances),
            len(warn_compliances),
        )

        if len(compliances) > 0:
            regulated_annotations.extend(
                [_create_hit_info(row, compliances) for i, row in data.iterrows()]
            )
            continue

        if len(warn_compliances) > 0:
            warn_regulated_annotations.extend(
                [_create_hit_info(row, warn_compliances) for i, row in data.iterrows()]
            )
            continue

        if len(conditional_compliances) > 0:
            regionally_regulated_annotations.extend(
                [
                    _create_hit_info(row, conditional_compliances)
                    for i, row in data.iterrows()
                ]
            )

    # Create the list of non-controlled taxids.
    non_regulated_annotations = [
        _create_hit_info(row) for _, row in non_regulated_for_region.iterrows()
    ]

    # update controlled and non-controlled taxids list depending on control list mode annotations.
    is_controlled = len(compliances) > 0
    is_warning_controlled = len(warn_compliances) > 0
    is_conditionally_controlled = (
        len(conditional_compliances) >= N_NON_REGIONAL_HITS_TO_WARN
    )

    if is_controlled:
        non_regulated_annotations = (
            non_regulated_annotations
            + warn_regulated_annotations
            + regionally_regulated_annotations
        )
        return _create_hit_result_from_annotations(
            display_name,
            regulated_annotations,
            non_regulated_annotations,
            compliances,
            ListMode.COMPLIANCE,
            step,
            match_range,
        )

    if is_warning_controlled:
        return _create_hit_result_from_annotations(
            display_name,
            warn_regulated_annotations,
            non_regulated_annotations,
            warn_compliances,
            ListMode.COMPLIANCE_WARN,
            step,
            match_range,
        )

    if is_conditionally_controlled:
        return _create_hit_result_from_annotations(
            display_name,
            regionally_regulated_annotations,
            non_regulated_annotations,
            conditional_compliances,
            ListMode.CONDITIONAL_NUM,
            step,
            match_range,
        )

    logger.debug(
        "No controlled entities following parsing of cluster %s",
        controlled_cluster_taxid,
    )

    return None, ""


def _create_hit_info(row: pd.Series, control_output=None) -> dict:
    """
    Return the a formated taxonomy annotation output. Optionally, if
    it is an item of a control list, we add that info.
    """
    output_dict = {
        "evalue": row["evalue"],
        "percent_identity": row["% identity"],
        "query_start": row["q. start"],
        "query_end": row["q. end"],
        "match_start": row["s. start"],
        "match_end": row["s. end"],
        "target_hit": row["subject acc."],
        "target_description": row["subject title"],
        "taxid": row["subject tax ids"],
    }
    if control_output:
        output_dict["genus"] = control_output[0].genus
        output_dict["species"] = control_output[0].species
        output_dict["controlled_by_lists"] = [
            {"list": get_control_lists(co.list).display_name, "source": co.source_text}
            for co in control_output
        ]

    return output_dict


def _create_hit_result_from_annotations(
    controlled_cluster_label,
    regulated_annotations,
    non_regulated_annotations,
    compliances,
    mode,
    step,
    match_range,
):

    status = ScreenStatus.ERROR
    control_text = "controlled"

    if mode == ListMode.COMPLIANCE:
        status = ScreenStatus.FLAG

    if mode == ListMode.COMPLIANCE_WARN:
        status = ScreenStatus.WARN
        control_text = "observed"

    if mode == ListMode.CONDITIONAL_NUM:
        status = ScreenStatus.PASS
        control_text = "regionally exempt"

    display_names = [compliance.name for compliance in compliances]
    display_name = max(display_names, key=len)
    categories = list(set([compliance.category for compliance in compliances]))

    domains_text = ", ".join(categories) or "sequence"

    hit_description = f"{control_text} {domains_text} {display_name}"

    if len(non_regulated_annotations) > 0:
        status = ScreenStatus.PASS
        logger.debug(
            "Mixed result: %d regulated, %d non-regulated annotations.",
            len(regulated_annotations),
            len(non_regulated_annotations),
        )
        hit_description = (
            f"Mix of {len(regulated_annotations)} {control_text} {domains_text}"
            f" and {len(non_regulated_annotations)} non-controlled {domains_text}"
        )

    reg_species = list(set([ra["species"] for ra in regulated_annotations]))
    reg_taxids = list(set([str(ra["taxid"]) for ra in regulated_annotations]))
    non_reg_taxids = list(set([str(ra["taxid"]) for ra in non_regulated_annotations]))

    # Compile per-domain counts for regulation_dict
    n_virus = sum(1 for d in categories if d == "Viruses")
    n_bacteria = sum(1 for d in categories if d == "Bacteria")
    n_fungi = sum(1 for d in categories if d == "Fungi")
    # Note Other is non-virus/bacteria/fungi/human parasite. Which typically means a non-human parasite
    n_parasite = sum(1 for d in categories if (d == "Human parasite" or d == "Other"))

    # Generate the Taxonomy Annotation dict for adding to HitResult
    taxonomy_annotations = {
        "controlled_taxa": regulated_annotations,
        "non-controlled_taxa": non_regulated_annotations,
        "statistics": {
            "number_of_controlled_taxids": len(reg_taxids),
            "number_of_non-controlled_taxids": len(non_reg_taxids),
            "controlled_bacteria": n_bacteria,
            "controlled_viruses": n_virus,
            "controlled_parasites": n_parasite,
            "controlled_fungi": n_fungi,
        },
    }

    log_message = _build_log_message(
        status,
        domains_text,
        [match_range],
        hit_description,
        reg_taxids,
        non_reg_taxids,
        reg_species,
    )

    new_hit = HitResult(
        HitScreenStatus(status, step),
        display_name,
        hit_description,
        match_range,
        {"category": categories, "controlled_taxonomy": taxonomy_annotations},
    )
    return new_hit, log_message


def _build_log_message(
    screen_status: ScreenStatus,
    domains_text: str,
    match_ranges: list[MatchRange],
    hit_description: str,
    reg_taxids: list[str],
    non_reg_taxids: list[str],
    reg_species: list[str],
) -> str:
    """
    Construct a legacy-style log line for a single regulated hit.

    Describes the outcome status, affected query coordinates, domain type,
    and associated tax IDs — mirroring the original .screen.log format.
    """
    match_ranges_text = ", ".join(map(str, match_ranges))
    reg_species_text = ", ".join(reg_species)
    reg_taxids_text = ", ".join(reg_taxids)
    non_reg_taxids_text = ", ".join(non_reg_taxids)

    s_reg = "" if len(reg_taxids) == 1 else "s"
    s_non_reg = "" if len(non_reg_taxids) == 1 else "s"

    msg = (
        f"\t --> {screen_status} at bases ({match_ranges_text})"
        f" {hit_description}.\n"
        f"\t   (Regulated Species: {reg_species_text}.\n"
        f"\t    Regulated TaxID{s_reg}: {reg_taxids_text}"
    )
    if non_reg_taxids:
        msg += f"\n\t    Non-Regulated TaxID{s_non_reg}: {non_reg_taxids_text}"
    msg += ")"
    return msg
