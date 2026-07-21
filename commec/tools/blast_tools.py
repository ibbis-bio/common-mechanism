#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for Blast related tools, a library for dealing with general blast file parsing tasks.
Useful for reading any blast related outputs, for example from Blastx or Blastn.

Also contains the abstract base class for blastX/N database search handlers.
"""

import os
import logging
from typing import BinaryIO, TextIO
import pandas as pd
import numpy as np
import glob
import commec.control_list as cl
from commec.tools.search_handler import SearchHandler, DatabaseValidationError
from commec.config.constants import NON_CODING_REGION_PERCENT_IDENTITY_THRESHOLD

logger = logging.getLogger(__name__)


class BlastHandler(SearchHandler):
    """
    A Database handler specifically for use with Blast.
    Inherit from this, and implement screen()
    """

    def read_output(self):
        output_dataframe = pd.DataFrame()
        if self.has_hits():
            output_dataframe = read_blast(self.out_file)
        return output_dataframe

    def validate_output(self) -> bool:
        """
        BLAST tabular output (outfmt 6) contains no header lines, so a valid
        run that simply found no hits is legitimately an empty file. Only require that
        the output file exists, rather than that it is non-empty.
        """
        return os.path.isfile(self.out_file)

    def _validate_db(self):
        """
        Blast expects a set of files with shared prefix, rather than a single file.
        Here we validate such directory structures for Blast related search handlers.
        """
        if not os.path.isdir(self.db_directory):
            raise DatabaseValidationError(
                f"No screening database directory found at: {self.db_directory}."
                " Directory path can be set via --databases option or --config yaml."
            )

        # Search for files of provided prefix.
        filename, extension = os.path.splitext(self.db_file)
        search_file = os.path.join(
            self.db_directory, "*" + os.path.basename(filename) + "*" + extension
        )
        files = glob.glob(search_file)
        if len(files) == 0:
            raise DatabaseValidationError(
                f"Screening files not found. No files matched search: {filename}*"
                " File location can be set via --databases option or --config yaml."
            )

        # Search for files of provided prefix.
        filename, extension = os.path.splitext(self.db_file)
        search_file = os.path.join(
            self.db_directory, "*" + os.path.basename(filename) + "*" + extension
        )
        files = glob.glob(search_file)
        if len(files) == 0:
            raise DatabaseValidationError(
                f"Mandatory screening files with {filename}* not found."
            )


def read_blast(blast_file: str | os.PathLike | BinaryIO | TextIO) -> pd.DataFrame:
    """
    Read in BLAST files and pre-format the data frame with essential info.

    BLAST is run with tabular output (outfmt 6) which contains no comment
    lines, so the file can be handed to pandas directly.
    """
    blast = pd.read_csv(blast_file, sep="\t", header=None)
    columns = [
        "query acc.",
        "subject title",
        "subject acc.",
        "subject tax ids",
        "evalue",
        "bit score",
        "% identity",
        "query length",
        "q. start",
        "q. end",
        "subject length",
        "s. start",
        "s. end",
    ]

    blast.columns = columns
    blast = blast.sort_values(by=["% identity"], ascending=False)
    blast["log evalue"] = -np.log10(pd.to_numeric(blast["evalue"]) + 1e-300)
    blast["q. coverage"] = (
        abs(blast["q. end"] - blast["q. start"]) / blast["query length"].max()
    )
    blast["s. coverage"] = (
        abs(blast["s. end"] - blast["s. start"]) / blast["subject length"]
    )
    blast["query acc."] = blast["query acc."].astype(str)

    blast = blast[blast["subject tax ids"].notna()]
    blast = blast.reset_index(drop=True)

    return blast


def get_top_hits(blast: pd.DataFrame):
    """
    Trim BLAST results down to the top hit for each base.
    """

    assert isinstance(blast, pd.DataFrame), (
        "get_top_hits expects a pandas dataframe object."
    )
    if "query acc." in blast.columns:
        assert len(blast["query acc."].unique() == 1)

    if blast.empty:
        logger.debug("Empty dataframe passed to Get Top Hits.")
        return blast

    top_hits = _trim_overlapping(blast)
    top_hits = top_hits.sort_values("% identity", ascending=False)

    # only keep coordinates of each hit that are not already covered by a better hit
    df = top_hits

    rerun = 1
    while (
        rerun == 1
    ):  # edges of hits can be moved within a higher scoring hit in the first pass
        df, rerun = _trim_edges(df)

    for j in df.index:
        top_hits.loc[j, "subject length"] = max(
            [df.loc[j, "q. start"], df.loc[j, "q. end"]]
        ) - min([df.loc[j, "q. start"], df.loc[j, "q. end"]])
        top_hits.loc[j, "q. start"] = df.loc[j, "q. start"]
        top_hits.loc[j, "q. end"] = df.loc[j, "q. end"]

    top_hits = top_hits.sort_values("q. start")
    top_hits = top_hits[top_hits["q. start"] != 0]

    # only keep annotations covering 50 bases or more
    top_hits = top_hits[top_hits["subject length"] >= 50]
    top_hits = top_hits.reset_index(drop=True)
    return top_hits


def get_high_identity_hits(
    blast_output_file, threshold=NON_CODING_REGION_PERCENT_IDENTITY_THRESHOLD
):
    """
    Read all hits with high sequence identity above a given threshold from a
    BLAST results file. This is current used to quickly grab candidates from the
    Protein Taxonomy search to identify non-coding regions of a query.
    """
    hits = read_blast(blast_output_file)
    hits = pd.concat(
        _trim_overlapping(group) for _, group in hits.groupby("query acc.")
    )
    return hits[hits["% identity"] >= threshold]


def get_controlled_labels(
    blast: pd.DataFrame, taxids_column_name="subject tax ids"
) -> pd.DataFrame:
    """
    Uses the Control Lists to label each taxid/accession supplied from imported
    blast outputs to label as regulated True/False, depending on its existance
    in the regulated lists.
    Ignores synthetic constructs.
    Returns the mutated blast dataframe to contain the following headings:
    cluster_hash - the TaxID of the "most parent" in the context of a controlled taxid hierarchy.
    species - the species of the subject taxid
    genus - genus of subject taxid.
    control_list - list of all controlled outputs.
    """
    logger.debug("Dealing with multi-taxid entries... ")

    # Deal with multiple taxids entries, including removing ignored taxids.
    blast = split_by_tax_id(blast, taxids_column_name)

    # Build a mapping {taxid: ignoriness} and Map back to the dataframe
    taxids_to_ignore = {
        taxid: cl.should_ignore(taxid) for taxid in blast[taxids_column_name].unique()
    }
    logger.debug(
        "Removing %s ignored taxids: \n%s", len(taxids_to_ignore), taxids_to_ignore
    )
    ignored_mask = blast[taxids_column_name].map(taxids_to_ignore)
    blast = blast[~ignored_mask]

    # Build a mapping {taxid: truthiness} and Map back to the dataframe
    taxid_to_regulated = {
        taxid: cl.is_regulated(taxid) for taxid in blast[taxids_column_name].unique()
    }
    blast["regulated"] = blast[taxids_column_name].map(taxid_to_regulated)

    # Get unique taxids
    unique_taxids = blast[taxids_column_name].dropna().unique()
    logger.debug("Checking %s unique taxids: \n%s", len(unique_taxids), unique_taxids)
    taxid_to_controlhash = {
        taxid: cl.get_cluster_hash(taxid) for taxid in unique_taxids
    }

    # Map control status hash to the dataframe
    blast["control_hash"] = blast[taxids_column_name].map(taxid_to_controlhash)

    return blast


def split_by_tax_id(blast: pd.DataFrame, taxids_col_name="subject tax ids"):
    """
    Some results will have multiple tax ids listed in a semicolon-separated list.
    These are handled depending on context:
    * If any of the taxids should be ignored, we ignore that entire row.
    * If the row is labelled as MULTISPECIES, we substitute the taxid with that information.
    * If the row has more than 100 taxids, we treat as MULTISPECIES.
    * The row is split into many rows for any remaining valid taxids.
    """
    # Create a list to hold all rows, including split ones
    new_rows = []
    for _, row in blast.iterrows():
        tax_ids = str(row[taxids_col_name]).split(";")

        # We don't care for any multiple taxid entry that hits an ignored taxid.
        if any(cl.should_ignore(taxid) for taxid in tax_ids):
            continue

        if str(row["subject title"]).startswith("MULTISPECIES:"):
            logger.debug(
                "     Treating labelled MULTISPECIES:, %i multi-taxid entry, as MULTISPECIES",
                len(tax_ids),
            )
            new_row = row.copy()
            new_row[taxids_col_name] = "MULTISPECIES"
            new_rows.append(new_row)
            continue

        if len(tax_ids) > 100:
            logger.debug(
                "     Treating %i multi-taxid entry as MULTISPECIES", len(tax_ids)
            )
            new_row = row.copy()
            new_row[taxids_col_name] = "MULTISPECIES"
            new_rows.append(new_row)
            continue

        if len(tax_ids) > 1:
            logger.debug("     Splitting %i multi-taxid entry", len(tax_ids))
            # If there are multiple tax IDs, create a new row for each controlled one
            for tax_id in tax_ids:
                new_row = row.copy()
                new_row[taxids_col_name] = tax_id
                new_rows.append(new_row)
            continue

        # If there's only one tax ID, keep the original row
        new_rows.append(row)

    # Create a new DataFrame from the list of rows
    split = pd.DataFrame(new_rows)
    return split


def find_clusters(
    input_data: pd.DataFrame,
    start_heading: str = "q. start",
    end_heading: str = "q. end",
    output_label_heading: str = "cluster",
) -> tuple[pd.DataFrame, list[tuple]]:
    """
    Groups rows into clusters where any set of rows with collectively overlapping
    ranges (transitively) belong to the same cluster.

    Returns a copy of input_data with a new integer cluster-ID column, and a list
    of (min_start, max_end) tuples — one per cluster.
    """
    assert not input_data.empty, "Empty dataframe given to find_clusters()"
    assert start_heading in input_data.columns
    assert end_heading in input_data.columns

    sorted_df = input_data.sort_values(by=start_heading)

    clusters: list[tuple] = []
    row_cluster_ids: dict = {}
    current_cluster = -1
    current_max_end = None

    for idx, row in sorted_df.iterrows():
        start = row[start_heading] - 1
        end = row[end_heading] + 1

        if current_cluster == -1 or start >= current_max_end:
            current_cluster += 1
            clusters.append((row[start_heading], row[end_heading]))
            current_max_end = end
        else:
            current_max_end = max(current_max_end, end)
            clusters[current_cluster] = (clusters[current_cluster][0], current_max_end)

        row_cluster_ids[idx] = current_cluster

    output_data = input_data.copy()
    output_data[output_label_heading] = output_data.index.map(row_cluster_ids)

    return output_data, clusters


def _trim_overlapping(blast: pd.DataFrame):
    """
    Remove any hits that are completely overlapped by another, higher-quality hit.
    """
    # set start to the lowest coordinate and end to the highest to avoid confusion
    blast = _shift_hits_pos_strand(blast)

    # if any multispecies hits contain regulated pathogens, put the regulated up top
    if "regulated" in blast:
        blast = blast.sort_values(by=["regulated"], ascending=False)

    # rank hits by percent identity, then bit score
    blast = blast.sort_values(by=["% identity", "bit score"], ascending=False)
    blast = blast.reset_index(drop=True)

    blast2 = blast

    # only keep top-ranked hits that don't overlap
    for i in blast.index:  # run through each hit from the top
        for j in blast.index[(i + 1) :]:  # compare to each below
            if j in blast2.index:
                i_envelopes_j = (
                    blast.loc[i, "q. start"] <= blast.loc[j, "q. start"]
                    and blast.loc[i, "q. end"] >= blast.loc[j, "q. end"]
                )
                not_shared_coords_or_identity = (
                    blast.loc[i, "q. start"] < blast.loc[j, "q. start"]
                    or blast.loc[i, "q. end"] > blast.loc[j, "q. end"]
                    or blast.loc[i, "% identity"] > blast.loc[j, "% identity"]
                )

                if i_envelopes_j and not_shared_coords_or_identity:
                    blast2 = blast2.drop([j])

    blast2 = blast2.reset_index(drop=True)

    return blast2


def _shift_hits_pos_strand(blast):
    for j in blast.index:
        if blast.loc[j, "q. start"] > blast.loc[j, "q. end"]:
            start = blast.loc[j, "q. end"]
            end = blast.loc[j, "q. start"]
            blast.loc[j, "q. start"] = start
            blast.loc[j, "q. end"] = end
    return blast


def _trim_edges(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """
    Function for filtering a Blast derived dataframe, removes weaker hits
    (based on % identity) that have extents within that of stronger hits.
    Also trims weaker hits to not overlap with stronger hits.

    input:
        - df :: pd.Dataframe, containing 'q. start', 'q. end', and '% identity' information.

    output:
        - The altered dataframe.s
        - Integer flag, 0 or 1, where 1 indicates a need to rerun _trim_edges()
    """

    assert "q. start" in df.columns, (
        'Expected column "q. start" does not exist for _trim_edges().\n'
        f"Existing columns: {', '.join(df.columns)}"
    )

    assert "q. end" in df.columns, (
        'Expected column "q. end" does not exist for _trim_edges().\n'
        f"Existing columns: {', '.join(df.columns)}"
    )

    assert "% identity" in df.columns, (
        'Expected column "query length" does not exist for _trim_edges().\n'
        f"Existing columns: {', '.join(df.columns)}"
    )

    for top, i in enumerate(df.index):  # run through each hit from the top
        for _, j in enumerate(
            df.index[top + 1 :], start=top + 1
        ):  # compare to each below
            i_start = df.loc[i, "q. start"]
            i_end = df.loc[i, "q. end"]
            j_start = df.loc[j, "q. start"]
            j_end = df.loc[j, "q. end"]

            # if the beginning of a weaker hit is inside a stronger hit, alter its start to the next base after that hit
            if j_start >= i_start and j_start <= i_end:
                # keep equivalent hits
                if df.loc[j, "% identity"] == df.loc[i, "% identity"]:
                    pass
                # if the hit extends past the end of the earlier one
                elif i_end + 1 < j_end:
                    df.loc[j, "q. start"] = i_end + 1
                elif (
                    i_end == j_end
                    and df.loc[j, "% identity"] == df.loc[i, "% identity"]
                ):
                    pass
                # remove if the hit is contained in the earlier one
                else:
                    df.loc[j, "q. start"] = 0
                    df.loc[j, "q. end"] = 0

            # if the end of a weaker hit is inside a stronger hit, alter the end to just before that hit
            if j_end >= i_start and j_end <= i_end:
                # keep equivalent hits
                if df.loc[j, "% identity"] == df.loc[i, "% identity"]:
                    pass
                elif i_start - 1 > j_start:
                    df.loc[j, "q. end"] = i_start - 1
                elif (
                    i_start == j_start
                    and df.loc[j, "% identity"] == df.loc[i, "% identity"]
                ):
                    pass
                else:
                    df.loc[j, "q. start"] = 0
                    df.loc[j, "q. end"] = 0

    rerun = 0
    mix_starts = 0
    for start in set(df["q. start"]):
        if (
            len(
                set(
                    zip(
                        df["q. start"][df["q. start"] == start],
                        df["q. end"][df["q. start"] == start],
                    )
                )
            )
            > 1
        ):
            if (
                len(set(df["% identity"][df["q. start"] == start])) > 1
            ):  # if there are overlapping annotations with different % identities, re-run
                rerun = 1
                mix_starts = mix_starts + 1
    return df, rerun
