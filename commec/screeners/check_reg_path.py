#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for parsing taxonomic BLAST screen output and annotating regulated pathogen hits.

Provides ``parse_taxonomy_hits`` as the public entry point, which validates inputs,
loads and labels BLAST results, applies control-list compliance rules, and updates
a ``ScreenResult`` in-place with ``HitResult`` objects for each regulated hit found.
"""
import logging
import os
from dataclasses import asdict

import pandas as pd

from commec.tools.search_handler import SearchHandler
from commec.config.query import Query
from commec.tools.blast_tools import (
    read_blast,
    get_controlled_labels,
)
from commec.config.result import (
    ScreenResult,
    TaxonomyAnnotation,
    HitResult,
    ScreenStep,
    ScreenStatus,
    HitScreenStatus,
    MatchRange,
    compare,
)
from commec.control_list import (
    get_control_lists,
    get_regulation,
    is_regulated,
    ListMode,
    ControlList,
)

from commec.tools.data_functions import (
    find_clusters,
    get_top_hits,
    remove_synthetic_and_vaccine_taxids,
)

from commec.config.constants import N_NON_REGIONAL_HITS_TO_WARN

pd.set_option("display.max_colwidth", 10000)

logger = logging.getLogger(__name__)


# ── Input validation ───────────────────────────────────────────────────────────


def _check_inputs(
    search_handler: SearchHandler,
    low_concern_taxid_path: str | os.PathLike,
    biorisk_taxid_path: str | os.PathLike,
    taxonomy_directory: str | os.PathLike,
) -> bool:
    """
    Check existence of required input files and directories.

    Returns True if all checks pass and processing may continue.
    """
    if not search_handler.validate_output():
        logger.info("\t...ERROR: Taxonomic search results empty\n %s", search_handler.out_file)
        return False

    if not os.path.exists(low_concern_taxid_path):
        logger.error("\t...low-concern database file %s does not exist\n", low_concern_taxid_path)
        return False

    if not os.path.exists(taxonomy_directory):
        logger.error("\t...taxonomy directory %s does not exist\n", taxonomy_directory)
        return False

    return True



# ── Region-level annotation collection ────────────────────────────────────────


def _collect_region_annotations(
    q_start : int, q_end : int,
    unique_query_data: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    For a single BLAST region, find all co-localised hits (shared query start OR end),
    split them into regulated and non-regulated DataFrames, and apply control-list rules.

    Non-regulated hits are capped at the 10 lowest e-values and deduplicated by tax ID.
    Control-list compliance may reclassify regulated rows into the non-regulated set.

    Returns:
        regulated_for_region: Regulated hits at this position.
        non_regulated_for_region: Non-regulated hits at this position.
    """
    shared_site = unique_query_data[
        (unique_query_data["q. start"] == q_start)
        | (unique_query_data["q. end"] == q_end)
    ]

    regulated_for_region = shared_site[shared_site["regulated"]]
    non_regulated_for_region = (
        shared_site[~shared_site["regulated"]]
        .sort_values(by="evalue", ascending=True)
        .head(10)
    )

    # Deduplicate by taxID for conciseness
    regulated_for_region = regulated_for_region.drop_duplicates(
        subset=["subject tax ids"], keep="first"
    )
    non_regulated_for_region = non_regulated_for_region.drop_duplicates(
        subset=["subject tax ids"], keep="first"
    )

    # Reclassify entries based on control-list compliance status
    regulated_for_region, non_regulated_for_region = _update_using_control_lists(
        regulated_for_region, non_regulated_for_region
    )

    return regulated_for_region, non_regulated_for_region


def _update_using_control_lists(
    regulated: pd.DataFrame,
    non_regulated: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reclassify regulated entries based on their control-list compliance status.

    For each row in ``regulated``:
      * ``COMPLIANCE`` or ``COMPLIANCE_WARN`` — retain as regulated.
      * ``CONDITIONAL_NUM`` — mark for movement to the non-regulated set.

    If the number of entries flagged for movement reaches
    ``CONDITIONAL_NUMBER_TO_ALLOW``, all are retained as regulated (the
    conditional threshold is met and the hit remains flagged).

    Side-effects: populates ``list_acronym``, ``category``, ``species``, and
    ``genus`` columns on ``regulated`` rows during iteration.

    Args:
        regulated: DataFrame of regulated BLAST hits for a single position.
        non_regulated: DataFrame of non-regulated BLAST hits for the same position.

    Returns:
        Updated (regulated, non_regulated) tuple.
    """
    CONDITIONAL_NUMBER_TO_ALLOW = 2
    indexes_to_move: list = []

    for index, row in regulated.iterrows():
        control_data, _context = get_regulation(row["subject tax ids"])
        for info in control_data:
            clist: ControlList = get_control_lists(info.list)
            regulated.at[index, "list_acronym"] = clist.acronym
            regulated.at[index, "category"] = info.category
            regulated.at[index, "species"] = info.species
            regulated.at[index, "genus"] = info.genus
            if clist.status in (ListMode.COMPLIANCE, ListMode.COMPLIANCE_WARN):
                continue
            if clist.status == ListMode.CONDITIONAL_NUM:
                indexes_to_move.append(index)

    # If enough conditional entries are present, treat all as regulated
    if len(indexes_to_move) >= CONDITIONAL_NUMBER_TO_ALLOW:
        indexes_to_move = []

    if indexes_to_move:
        rows_to_move = regulated.loc[indexes_to_move]
        non_regulated = pd.concat([non_regulated, rows_to_move], ignore_index=True)
        regulated = regulated.drop(index=indexes_to_move)

    return regulated, non_regulated

def _extract_unique_annotations(df: pd.DataFrame) -> set[TaxonomyAnnotation]:
    """
    Build a set of TaxonomyAnnotation objects from a BLAST-result DataFrame.
    Expected columns: 
    """
    return {
        TaxonomyAnnotation(*row)
        for row in df[
            ["evalue",
             "subject tax ids",
             "subject acc.",
             "subject title"]
        ].itertuples(index=False, name=None)
    }


# ── Annotation enrichment ──────────────────────────────────────────────────────


def _enrich_regulated_annotations(
    annotation_list: list[dict],
) -> tuple[list[dict], list[str]]:
    """
    Populate each regulated annotation dict with its control-list info and
    collect domain categories encountered across all annotations.

    For each annotation the inner loop runs until the first Viruses / Bacteria /
    Eukaryota domain is found, then stops — preserving the original domain-counting
    behaviour where each annotation contributes at most one primary domain count.

    Returns:
        annotation_list: Mutated in-place; each dict gains a 'control_list' key.
        domains: Flat list of domain strings accumulated across all annotations.
    """
    domains: list[str] = []

    for annotation in annotation_list:
        control_info, _context = get_regulation(annotation["taxid"])
        annotation["control_list"] = control_info

        for control_entry in control_info:
            control_entry.list = str(get_control_lists(control_entry.list))
            domain = control_entry.category
            domains.append(domain)
            if domain in ("Viruses", "Bacteria", "Eukaryota"):
                logger.debug("\t\t\tAdded domain: %s", domain)
                break

    return annotation_list, domains


def _enrich_non_regulated_annotations(annotation_list: list[dict]) -> list[dict]:
    """
    Attach control-list info to non-regulated annotation dicts where available.

    Useful for surfacing conditional-compliance context in the JSON output.
    Mutates and returns annotation_list.
    """
    for annotation in annotation_list:
        control_info, _context = get_regulation(annotation["taxid"])
        for control_entry in control_info:
            control_entry.list = str(get_control_lists(control_entry.list))
        if control_info:
            annotation["control_list"] = control_info
    return annotation_list


# ── Screen outcome determination ───────────────────────────────────────────────


def _determine_screen_outcome(
    reg_annotations: list[dict],
    regionally_reg_annotations: list[dict],
    warn_reg_annotations: list[dict],
    non_reg_annotations: list[dict],
    domains: list[str],
    hit_title: str,
) -> tuple[ScreenStatus, str, str]:
    """
    Determine the ScreenStatus and human-readable description for a single hit.

    Three outcomes (in priority order):
      * Regionally controlled — all regulated annotations removed by compliance → PASS
      * Mixed result — both regulated and non-regulated annotations present → PASS
      * Fully regulated — only regulated annotations remain → FLAG

    Returns:
        screen_status: ScreenStatus to apply to this hit.
        hit_description: Short text description of the outcome.
        domains_text: Comma-separated unique domain categories (default "sequence").
    """
    unique_domains = list(set(domains))
    domains_text = ", ".join(unique_domains) or "sequence"

    if reg_annotations:
        logger.debug("Regulated annotations.")
        return (
            ScreenStatus.FLAG,
            f"Controlled {domains_text} - {hit_title}",
            domains_text,
        )

    if not reg_annotations:
        logger.debug("No regulated annotations remain; regionally controlled.")
        return (
            ScreenStatus.PASS,
            f"Regionally controlled {domains_text}",
            domains_text,
        )

    if non_reg_annotations:
        logger.debug(
            "Mixed result: %d regulated, %d non-regulated annotations.",
            len(reg_annotations), len(non_reg_annotations),
        )
        return (
            ScreenStatus.PASS,
            (
                f"Mix of {len(reg_annotations)} controlled {domains_text}"
                f" and {len(non_reg_annotations)} non-controlled {domains_text}"
            ),
            domains_text,
        )




# ── Per-hit orchestration ──────────────────────────────────────────────────────
def get_hit_result_from_data(unique_query_data : pd.DataFrame, step : ScreenStep) -> list[HitResult]:
    """
    Given a dataframe that has been annotated with control list information,
    derive a list of HitResults that represent the outcome of this data.
    Data submitted is taken as a whole, and is assumed to represent a single Query.
    """
    output_hit_results : list[HitResult] = []

    query_acc = unique_query_data["query acc."].to_list()[0]

    # Filter data to the top hits only.
    top_hits = get_top_hits(unique_query_data)

    # Merge top_hits based on only unique "subject title"
    top_hits = top_hits.drop_duplicates(subset=["subject title"], keep="first")

    top_hits = get_controlled_labels(top_hits)
    regulated_only = top_hits[top_hits["regulated"]]
    non_regulated_only = top_hits[~top_hits["regulated"]]
    unique_controlled_parents = regulated_only["subject tax ids"].unique()

    logger.debug("%s: %s query has %i unique control hits:\n%s",
                step, query_acc, len(unique_controlled_parents), unique_controlled_parents)

    for controlled_parent in unique_controlled_parents:
        cluster_data, clusters = find_clusters(regulated_only[regulated_only["subject tax ids"] == controlled_parent])
        logger.debug("%s: query %s has %i regulated clusters", step, query_acc, len(clusters))
        for i, cluster in enumerate(clusters):
            logger.debug("Processing cluster[%i]: %s",i, cluster)
            new_hit_result = create_hit_result_for_cluster(cluster_data[cluster_data["cluster"] == i],
                                                            non_regulated_only, cluster[0], cluster[1], controlled_parent, step)
            if new_hit_result:
                output_hit_results.append(new_hit_result)
    return output_hit_results

def create_hit_info(row : pd.Series, control_output = None) -> dict:
    """
    Return the formated TaxonomyAnnotation Output. Optionally, if 
    it is an item of a control list, we can add that info.
    """
    output_dict =  {
        "evalue" : row["evalue"],
        "percent_identity" : row["% identity"],
        "taxid" : row["subject tax ids"],
        "start" : row["q. start"],
        "end" : row["q. end"],
        "target_hit" : row["subject acc."],
        "target_description" : row["subject title"],
    }
    if control_output:
        output_dict["genus"] = control_output[0].genus
        output_dict["species"] = control_output[0].species
        output_dict["controlled_by_lists"] = [
            {
                "list" : get_control_lists(co.list).display_name, 
                "source" : co.source_text
            } for co in control_output]

    return output_dict

def create_hit_result_for_cluster(
    cluster_data : pd.DataFrame,
    non_regulated_only : pd.DataFrame,
    cluster_start : int,
    cluster_end : int,
    controlled_cluster_taxid : str,
    step : ScreenStep) -> HitResult | None:
    """
    Given a cluster of controlled hits, and the common start and end point.
    Find overlapping non-regulated hits, and generate a HitResult, to describe
    the data for commec screen output.

    """
    cluster_regulation = get_regulation(controlled_cluster_taxid)
    species = cluster_regulation[0].species
    genus = cluster_regulation[0].genus
    display_name = cluster_regulation[0].name
    hit_name = f"{controlled_cluster_taxid}_{cluster_start}_{cluster_end}"

    best_evalue = min(cluster_data["evalue"])
    match_range = MatchRange(
            best_evalue,
            0, 0,
            int(cluster_start), int(cluster_end),
        )

    # Get the top 10 non-controlled hits for reporting.
    shared_site = non_regulated_only[
        (non_regulated_only["q. start"] == cluster_start)
        | (non_regulated_only["q. end"] == cluster_end)
    ]
    
    regulated_for_region = cluster_data
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
        data = regulated_for_region[regulated_for_region["subject tax ids"] == candidate_taxid]

        compliances.extend([output for output in control_output 
            if get_control_lists(output.list).status == ListMode.COMPLIANCE])
        conditional_compliances.extend([output for output in control_output 
            if get_control_lists(output.list).status == ListMode.CONDITIONAL_NUM])
        warn_compliances.extend([output for output in control_output 
            if get_control_lists(output.list).status == ListMode.COMPLIANCE_WARN])

        if len(compliances) > 0:
            regulated_annotations.extend(
                [create_hit_info(row, compliances) for i, row in data.iterrows()])
            continue

        if len(warn_compliances) > 0:
            warn_regulated_annotations.extend(
                [create_hit_info(row, warn_compliances) for i, row in data.iterrows()])
            continue

        if len(conditional_compliances) > 0:
            regionally_regulated_annotations.extend(
                [create_hit_info(row, conditional_compliances) for i, row in data.iterrows()])

    # Create the list of non-controlled taxids.
    for i, row in non_regulated_for_region.iterrows():
        non_regulated_annotations.extend(
            [create_hit_info(row) for i, row in data.iterrows()])

    # update controlled and non-controlled taxids list depending on control list mode annotations.
    is_controlled = len(compliances) > 0
    is_warning_controlled = len(warn_compliances) > 0
    is_conditionally_controlled = len(conditional_compliances) >= N_NON_REGIONAL_HITS_TO_WARN

    if is_controlled:
        non_regulated_annotations = non_regulated_annotations + warn_regulated_annotations + regionally_regulated_annotations
        return create_hit_result_from_annotations(
            regulated_annotations,
            non_regulated_annotations,
            compliances, ListMode.COMPLIANCE, step, match_range)

    if is_warning_controlled:
        return create_hit_result_from_annotations(
            warn_regulated_annotations,
            non_regulated_annotations,
            warn_compliances, ListMode.COMPLIANCE_WARN, step, match_range)

    if is_conditionally_controlled:
        return create_hit_result_from_annotations(
            regionally_regulated_annotations,
            non_regulated_annotations,
            conditional_compliances, ListMode.CONDITIONAL_NUM, step, match_range)

    logger.debug(" ... ... No controlled hits after parsing control list rules.")
    return None


def create_hit_result_from_annotations(
    regulated_annotations,
    non_regulated_annotations,
    compliances,
    mode, step, match_range):

    status = ScreenStatus.ERROR
    control_text = "controlled"

    if mode == ListMode.COMPLIANCE:
        status = ScreenStatus.FLAG

    if mode == ListMode.COMPLIANCE_WARN:
        status = ScreenStatus.WARN
        control_text = "observed"

    if mode == ListMode.CONDITIONAL_NUM:
        status = ScreenStatus.PASS
        control_text = "non-regionally controlled"


    display_names = [compliance.name for compliance in compliances]
    display_name = max(display_names, key=len)
    categories = list(set([compliance.category for compliance in compliances]))

    domains_text = ", ".join(categories) or "sequence"

    hit_description = f"{control_text} {domains_text} - {display_name}"

    if len(non_regulated_annotations) > 0:
        status = ScreenStatus.WARN
        logger.debug(
            "Mixed result: %d regulated, %d non-regulated annotations.",
            len(reg_annotations), len(non_reg_annotations),
        )
        hit_description =  (
                f"Mix of {len(reg_annotations)} {control_text} {domains_text}"
                f" and {len(non_reg_annotations)} non-controlled {domains_text}"
            )

    reg_accessions = list(set([ra["query acc."]] for ra in regulated_annotations))
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
        "controlled_taxa" : regulated_annotations,
        "non-controlled_taxa" : non_regulated_annotations,
        "regulated_accessions": reg_accessions,
        "statistics" : {
            "number_of_controlled_taxids" : len(reg_taxids),
            "number_of_non-controlled_taxids" : len(non_reg_taxids),
            "controlled_bacteria" : n_bacteria,
            "controlled_viruses" : n_virus,
            "controlled_parasites" : n_parasite,
            "controlled_fungi" : n_fungi
        }
    }

    log_message = _build_log_message(
        status, domains_text, [match_range],
        reg_taxids, non_reg_taxids, reg_species,
    )
    logger.info(log_message)

    new_hit = HitResult(
        HitScreenStatus(status, step),
        display_name,
        hit_description,
        match_range,
        {   "category": categories,
            "controlled_taxonomy": taxonomy_annotations},
    )
    return new_hit

def _process_regulated_hit(
    hit_acc: str,
    hit_data: pd.DataFrame,
    unique_query_data: pd.DataFrame,
    query_acc: str,
    queries: dict[str, Query],
    step: ScreenStep,
) -> tuple[HitResult, dict, list[str], ScreenStatus, list[MatchRange], str]:
    """
    Build a HitResult for a single regulated BLAST hit accession.

    Iterates over all BLAST regions for this hit, collecting regulated and
    non-regulated annotations at each shared query position, then enriches
    those annotations with control-list information and determines the outcome.

    Returns:
        new_hit: Constructed HitResult.
        regulation_dict: Annotation dict for JSON output.
        unique_domains: Deduplicated domain category strings.
        screen_status: Determined ScreenStatus for this hit.
        match_ranges: Deduplicated list of MatchRange objects.
        log_message: Pre-formatted legacy-style log string (for deferred INFO emit).
    """
    hit_name = hit_data["subject title"].values[0]
    logger.debug("Processing hit accession %s: %s", hit_acc, hit_name)
    logger.debug("Regulated Hit Data shape: %s\n%s", hit_data.shape, hit_data.head())

    match_ranges: list[MatchRange] = []
    regulated_annotations: set[TaxonomyAnnotation] = set()
    non_regulated_annotations: set[TaxonomyAnnotation] = set()
    reg_species: list[str] = []
    reg_taxids: list[str] = []
    non_reg_taxids: list[str] = []

    for _, region in hit_data.iterrows():
        match_range = MatchRange(
            float(region["evalue"]),
            int(region["s. start"]), int(region["s. end"]),
            int(region["q. start"]), int(region["q. end"]),
        )
        if step == ScreenStep.TAXONOMY_NT:
            match_range.query_start = queries[query_acc].nc_to_nt_query_coords(
                match_range.query_start
            )
            match_range.query_end = queries[query_acc].nc_to_nt_query_coords(
                match_range.query_end
            )
        match_ranges.append(match_range)
        logger.debug(
            "Processing region: evalue=%.2e  q=%d-%d  s=%d-%d",
            float(region["evalue"]),
            int(region["q. start"]), int(region["q. end"]),
            int(region["s. start"]), int(region["s. end"]),
        )

        reg_df, non_reg_df = _collect_region_annotations(region, unique_query_data)

        reg_species.extend(reg_df["species"])
        reg_taxids.extend(map(str, reg_df["subject tax ids"]))
        non_reg_taxids.extend(map(str, non_reg_df["subject tax ids"]))

        regulated_annotations |= _extract_unique_annotations(reg_df)
        non_regulated_annotations |= _extract_unique_annotations(non_reg_df)

    # Sort annotation dicts by taxid for deterministic JSON output
    reg_annotation_list = sorted(
        [asdict(t) for t in regulated_annotations], key=lambda d: d["taxid"]
    )
    non_reg_annotation_list = sorted(
        [asdict(t) for t in non_regulated_annotations], key=lambda d: d["taxid"]
    )

    # Enrich with control-list metadata
    reg_annotation_list, domains = _enrich_regulated_annotations(reg_annotation_list)
    non_reg_annotation_list = _enrich_non_regulated_annotations(non_reg_annotation_list)

    # Determine screen outcome
    screen_status, hit_description, domains_text = _determine_screen_outcome(
        reg_annotation_list, non_reg_annotation_list, domains, hit_name
    )

    logger.debug(
        "Hit %s: status=%s  regulated=%d  non_regulated=%d  domains=%s",
        hit_acc, screen_status, len(reg_annotation_list), len(non_reg_annotation_list), domains,
    )

    # Compile per-domain counts for regulation_dict
    n_virus = sum(1 for d in domains if d == "Viruses")
    n_bacteria = sum(1 for d in domains if d == "Bacteria")
    n_eukaryote = sum(1 for d in domains if d == "Eukaryota")

    regulation_dict = {
        "number_of_regulated_taxids": str(len(reg_annotation_list)),
        "number_of_unregulated_taxids": str(len(non_reg_annotation_list)),
        "regulated_eukaryotes": str(n_eukaryote),
        "regulated_bacteria": str(n_bacteria),
        "regulated_viruses": str(n_virus),
        "regulated_taxa": reg_annotation_list,
        "non_regulated_taxa": non_reg_annotation_list,
    }

    unique_domains = list(set(domains))
    match_ranges = list(set(match_ranges))

    new_hit = HitResult(
        HitScreenStatus(screen_status, step),
        hit_name,
        hit_description,
        match_ranges,
        {"taxonomy": unique_domains, "regulated_taxonomy": [regulation_dict]},
    )

    # Build legacy-format log message (emitted later at INFO level, grouped by query)
    reg_taxids = sorted(set(reg_taxids))
    non_reg_taxids = sorted(set(non_reg_taxids))
    reg_species = sorted(set(reg_species))
    log_message = _build_log_message(
        screen_status, domains_text, match_ranges,
        reg_taxids, non_reg_taxids, reg_species,
    )
    logger.debug(log_message)
    logger.debug("\t\tRegulated Taxids: %s", reg_annotation_list)
    logger.debug("\t\tNon-Regulated Taxids: %s", non_reg_annotation_list)
    logger.debug("\t\tRanges: %s", match_ranges)

    return new_hit, regulation_dict, unique_domains, screen_status, match_ranges, log_message

# ── Logging helpers ────────────────────────────────────────────────────────────

def _build_log_message(
    screen_status: ScreenStatus,
    domains_text: str,
    match_ranges: list[MatchRange],
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

    alt_text = (
        "only " if screen_status == ScreenStatus.FLAG
        else "both regulated and non-" if reg_taxids
        else "externally "
    )
    s_reg = "" if len(reg_taxids) == 1 else "s"
    s_non_reg = "" if len(non_reg_taxids) == 1 else "s"

    msg = (
        f"\t --> {screen_status} at bases ({match_ranges_text})"
        f" found in {alt_text}regulated {domains_text}.\n"
        f"\t   (Regulated Species: {reg_species_text}.\n"
        f"\t    Regulated TaxID{s_reg}: {reg_taxids_text}"
    )
    if non_reg_taxids:
        msg += f"\n\t    Non-Regulated TaxID{s_non_reg}: {non_reg_taxids_text}"
    msg += ")"
    return msg


def _emit_query_logs(log_container: dict[str, list[str]], step: ScreenStep) -> None:
    """
    Emit deferred per-query log messages at INFO level, grouped by query name.

    Messages are collected during hit processing and emitted here in query order
    to produce a legacy-compatible .screen.log layout.
    """
    taxtype = "protein" if step == ScreenStep.TAXONOMY_AA else "nucleotide"
    for query_name, log_list in log_container.items():
        if not log_list:
            continue
        s = "" if len(log_list) == 1 else "s"
        logger.info(" Regulated %s%s in %s:", taxtype, s, query_name)
        for log_text in log_list:
            logger.info(log_text)


# ── Public entry point ─────────────────────────────────────────────────────────


def parse_taxonomy_hits(
    search_handler: SearchHandler,
    low_concern_taxid_path: str | os.PathLike,
    biorisk_taxid_path: str | os.PathLike,
    taxonomy_directory: str | os.PathLike,
    data: ScreenResult,
    queries: dict[str, Query],
    step: ScreenStep,
    n_threads: int,
) -> int:
    """
    Parse a taxonomic BLAST screen output and update the ScreenResult in-place.

    Process:
        1. Load in the dataframe from provided SearchHandler.db_directory
        2. Remove Synthetic and Vaccine hits.
        2. Annotation the dataframe with control list annotations based on Accession (TaxID)
        3. Split the dataframe into per Query processing.
         per Query: (Point of Threading)
        4. Filter to Top Hits and Trim Edges.
        4. Split into Controlled vs non-controlled.
        6. Identify unique list of controlled parents.
         per UCP:
        5. Find regional clusters.
        6. Identify top 10 non-regulated hits + clustered hits.
        7. Generate a HitResults containing
            > ID by "hash_parent"_"start"_"stop"
            > Longest description for control hit.
            > ScreenResult
            > List of Control Lists Hit
            > List of Hits regulated (Those inside cluster)
            > List of Hits non-regulated (Top 10 non-regulated hits).
    Args:
        search_handler: Handle of the search tool used for the taxonomic screen.
        low_concern_taxid_path: Path to the low-concern taxid CSV.
        biorisk_taxid_path: Path to the regulated taxid CSV (reserved for future use).
        taxonomy_directory: Location of the NCBI taxonomy directory.
        data: ScreenResult to be updated in-place.
        queries: Mapping from query accession to Query objects.
        step: Taxonomy step (TAXONOMY_NT or TAXONOMY_AA).
        n_threads: Maximum threads available (reserved for future use).

    Returns:
        0 on success, 1 on unrecoverable input error.
    """
    logger.debug("Acquiring Taxonomic Data for JSON output:")

    if not _check_inputs(
        search_handler, low_concern_taxid_path,
        biorisk_taxid_path, taxonomy_directory
    ):
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

    # Label data with parent hash taxids.
    logger.info("    Identifying controlled taxa ...")

    # Build a mapping {taxid: truthiness}
    taxid_to_regulated = {taxid: is_regulated(taxid) for taxid in blast["subject tax ids"].unique()}
    # Map back to the dataframe
    blast["regulated"] = blast["subject tax ids"].map(taxid_to_regulated)

    # Early exit if no regulated hits.
    #if sum(1 for ch in blast["control_hash"] if ch is not None) == 0:
    if blast["regulated"].sum() == 0:
        logger.info("\t...no regulated hits\n")
        return 0

    logger.debug("%s Controlled Labels applied: shape %s\n%s", step, blast.shape, blast.head())

    # Step 2: Per-query analysis
    unique_query_accs = blast["query acc."].unique()
    logger.debug("%s: %d unique queries with regulated hits", step, len(unique_query_accs))

    for query_acc in unique_query_accs:
        # From here should be Threaded:
        logger.info("    Processing query: %s", query_acc)
        query_write = data.get_query(query_acc)
        if not query_write:
            logger.error("Query '%s' not found in ScreenResult during %s.", query_acc, step)
            continue

        # Filter to just top hits, trim the edges, sort and Clean up
        unique_query_data = blast[blast["query acc."] == query_acc]
        unique_query_data = unique_query_data.sort_values(by=["% identity"], ascending=False)
        unique_query_data = unique_query_data.reset_index(drop=True)
        hit_results_for_query = get_hit_result_from_data(unique_query_data, step)

        for new_hit in hit_results_for_query:
            logger.debug("Adding hit for query %s : %s", query_acc, new_hit)
            query_write.add_new_hit_information(new_hit)

    for query_acc in unique_query_accs:
        query_write = data.get_query(query_acc)
        if not query_write:
            logger.error("Query '%s' not found in ScreenResult during %s.", query_acc, step)
            continue
        # Collect all threads.

    # Step 3: Emit grouped INFO logs

    return 0

