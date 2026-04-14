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
    get_top_hits,
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
    ListMode,
    ControlList,
)

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


# ── Data preparation ───────────────────────────────────────────────────────────


def _load_blast_data(
    search_handler: SearchHandler,
    step: ScreenStep,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read raw BLAST output, apply regulated/non-regulated classification,
    and trim to per-base top hits.

    Returns:
        blast: Full labelled DataFrame (for initial hit marking).
        top_hits: Trimmed top-hit DataFrame with 'regulated' column (for analysis).
    """
    blast = read_blast(search_handler.out_file)
    logger.debug("%s Blast Import: shape %s\n%s", step, blast.shape, blast.head())

    blast = get_controlled_labels(blast)
    logger.debug("%s Controlled Labels applied: shape %s\n%s", step, blast.shape, blast.head())

    top_hits = get_top_hits(blast)
    logger.debug("%s Top Hits derived: shape %s\n%s", step, top_hits.shape, top_hits.head())

    return blast, top_hits


def _mark_initial_hits(
    blast: pd.DataFrame,
    queries: dict[str, Query],
) -> None:
    """
    Mark each Query object as having at least one BLAST hit.

    Logs an error for any query accession present in the BLAST output
    that cannot be located in the supplied queries map.
    """
    for query_acc in blast["query acc."].unique():
        query_obj = queries.get(query_acc)
        if query_obj:
            logger.debug("Confirming hit for query %s.", query_acc)
            query_obj.mark_as_hit()
        else:
            logger.error(
                "Could not mark query %s as hit: query not found in input queries.",
                query_acc,
            )


# ── Region-level annotation collection ────────────────────────────────────────


def _collect_region_annotations(
    region: pd.Series,
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
        (unique_query_data["q. start"] == region["q. start"])
        | (unique_query_data["q. end"] == region["q. end"])
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

    Expected columns: evalue, subject tax ids, subject acc., subject title.
    """
    return {
        TaxonomyAnnotation(*row)
        for row in df[
            ["evalue", "subject tax ids", "subject acc.", "subject title"]
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

    logger.debug("Only regulated annotations present — FLAG.")
    return (
        ScreenStatus.FLAG,
        f"Controlled {domains_text} - {hit_title}",
        domains_text,
    )


# ── Per-hit orchestration ──────────────────────────────────────────────────────


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
        {"domain": unique_domains, "regulated_taxonomy": [regulation_dict]},
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


# ── ScreenResult update ────────────────────────────────────────────────────────


def _apply_hit_to_screen_result(
    query_write,
    hit_acc: str,
    new_hit: HitResult,
    match_ranges: list[MatchRange],
    domains: list[str],
    regulation_dict: dict,
    screen_status: ScreenStatus,
    step: ScreenStep,
) -> None:
    """
    Update the QueryResult in the ScreenResult with the new HitResult.

    Updates the query-level step status and, if the hit accession was already
    recorded (duplicate accession across different regions), merges the new
    match ranges, domain info, and taxonomy annotations into the existing entry.
    """
    query_write.status.update_step_status(step, screen_status)

    if query_write.add_new_hit_information(new_hit):
        existing_hit = query_write.get_hit(hit_acc)
        if existing_hit:
            existing_hit.ranges.extend(match_ranges)
            existing_hit.annotations["domain"] = domains
            existing_hit.annotations["regulated_taxonomy"].append(regulation_dict)
            existing_hit.recommendation.status = compare(
                existing_hit.recommendation.status, screen_status
            )
            existing_hit.description += "," + new_hit.description


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

    Processing pipeline:
        1. Validate inputs.
        2. Load BLAST data, apply regulated labels, and derive top hits.
        3. Mark any query that produced a BLAST hit.
        4. For each query × regulated hit accession × BLAST region:
           a. Collect co-localised regulated and non-regulated hits.
           b. Apply control-list compliance rules.
           c. Determine screen status and build a HitResult.
           d. Update the ScreenResult.
        5. Emit grouped INFO-level log messages.

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
        search_handler, low_concern_taxid_path, biorisk_taxid_path, taxonomy_directory
    ):
        return 1

    # Default all queries to PASS; downstream hits will override where needed.
    for query in data.queries.values():
        query.status.set_step_status(step, ScreenStatus.PASS)

    if not search_handler.has_hits():
        logger.info("\t...no hits\n")
        return 0

    # ── Phase 1: Data preparation ─────────────────────────────────────────────
    blast, top_hits = _load_blast_data(search_handler, step)

    # ── Phase 2: Mark queries with any BLAST hit ──────────────────────────────
    _mark_initial_hits(blast, queries)

    if top_hits["regulated"].sum() == 0:
        logger.info("\t...no regulated hits\n")
        return 0

    # ── Phase 3: Per-query, per-hit analysis ──────────────────────────────────
    log_container: dict[str, list[str]] = {key: [] for key in data.queries}
    unique_query_accs = top_hits["query acc."].unique()
    logger.debug("%s: %d unique queries with regulated hits", step, len(unique_query_accs))

    for query_acc in unique_query_accs:
        logger.debug("Processing query: %s", query_acc)

        query_write = data.get_query(query_acc)
        if not query_write:
            logger.error("Query '%s' not found in ScreenResult during %s.", query_acc, step)
            continue

        unique_query_data = top_hits[top_hits["query acc."] == query_acc]
        regulated_only = unique_query_data[unique_query_data["regulated"]]
        regulated_hit_accs = regulated_only["subject acc."].unique()
        logger.debug(
            "%s: query %s has %d regulated hit accessions",
            step, query_acc, len(regulated_hit_accs),
        )

        for hit_acc in regulated_hit_accs:
            logger.debug("Processing hit accession: %s", hit_acc)
            hit_data = regulated_only[regulated_only["subject acc."] == hit_acc]

            new_hit, regulation_dict, unique_domains, screen_status, match_ranges, log_msg = (
                _process_regulated_hit(
                    hit_acc, hit_data, unique_query_data, query_acc, queries, step
                )
            )

            _apply_hit_to_screen_result(
                query_write, hit_acc, new_hit,
                match_ranges, unique_domains, regulation_dict, screen_status, step,
            )

            log_container[query_acc].append(log_msg)

    # ── Phase 4: Emit grouped INFO logs ───────────────────────────────────────
    _emit_query_logs(log_container, step)

    return 0

