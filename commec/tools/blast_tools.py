#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for Blast related tools, a library for dealing with general blast file parsing tasks.
Useful for reading any blast related outputs, for example from Blastx, Blastn, or diamond.

Also contains the abstract base class for blastX/N/Diamond database search handlers.
"""
import os
import logging
import glob
from typing import List, Union, BinaryIO, TextIO
import pytaxonkit
import pandas as pd
import numpy as np
from commec.tools.search_handler import SearchHandler, DatabaseValidationError

TAXID_SYNTHETIC_CONSTRUCTS = 32630
TAXID_VECTORS = 29278


class BlastHandler(SearchHandler):
    """
    A Database handler specifically for use with Blast.
    Inherit from this, and implement screen()
    """

    def read_output(self):
        output_dataframe = read_blast(self.out_file)
        return output_dataframe

    def _validate_db(self):
        """
        Blast expects a set of files with shared prefix, rather than a single file.
        Here we validate such directory structures for Blast related search handlers.
        """
        if not os.path.isdir(self.db_directory):
            raise DatabaseValidationError(
                f"Mandatory screening database directory not found: {self.db_directory}."
                " Screening directory can be set via --databases option or --config yaml."
            )

        # Search for files of provided prefix.
        filename, extension = os.path.splitext(self.db_file)
        search_file = os.path.join(
            self.db_directory, "*" + os.path.basename(filename) + "*" + extension
        )
        files = glob.glob(search_file)
        if len(files) == 0:
            raise DatabaseValidationError(
                f"Mandatory screening files not found. No files matched search: {filename}*"
                " File location can be set via --databases option or --config yaml."
            )

        # Search for files of provided prefix.
        filename, extension = os.path.splitext(self.db_file)
        search_file = os.path.join(self.db_directory, "*" + os.path.basename(filename) + "*" + extension)
        files = glob.glob(search_file)
        if len(files) == 0:
            raise DatabaseValidationError(f"Mandatory screening files with {filename}* not found.")

def _split_by_tax_id(blast: pd.DataFrame, taxids_col_name="subject tax ids"):
    """
    Some results will have multiple tax ids listed in a semicolon-separated list; split these into
    multiple rows, each with their own taxon id.
    """
    # Create a list to hold all rows, including split ones
    new_rows = []

    for _, row in blast.iterrows():
        tax_ids = str(row[taxids_col_name]).split(";")
        if len(tax_ids) > 1:
            # If there are multiple tax IDs, create a new row for each
            for tax_id in tax_ids:
                new_row = row.copy()
                new_row[taxids_col_name] = tax_id
                new_rows.append(new_row)
        else:
            # If there's only one tax ID, keep the original row
            new_rows.append(row)

    # Create a new DataFrame from the list of rows
    split = pd.DataFrame(new_rows)
    split[taxids_col_name] = split[taxids_col_name].astype("int")
    return split


def _get_lineages(taxids, db_path: str | os.PathLike, threads: int):
    """
    Get the full lineage for each unique taxid. This is needed to determine whether it belongs to
    any regulated pathogen, since pathogens might be regulated at various points in the lineage
    (i.e. not just species or genus).
    """
    lin = pytaxonkit.lineage(taxids.drop_duplicates(), data_dir=db_path, threads=threads)

    # Warn about error codes from lineage search
    taxids_not_found = lin[lin["Code"] == -1]["TaxID"]
    if not taxids_not_found.empty:
        logging.warning(
            "No information about the following taxID(s) was found in the taxonomy database: %s",
            ", ".join(taxids_not_found.astype(str).tolist()),
        )

    taxids_deleted = lin[lin["Code"] == 0]["TaxID"]
    if not taxids_deleted.empty:
        logging.warning(
            "The following taxID(s) have been deleted (in delnodes.dmp): %s",
            ", ".join(taxids_deleted.astype(str).tolist()),
        )

    # Remove non-success codes from the list
    return lin[(lin["Code"] != -1) & (lin["Code"] != 0)]


def get_taxonomic_labels(
    blast: pd.DataFrame,
    regulated_taxids: List[str],
    vaccine_taxids: List[str],
    db_path: str | os.PathLike,
    threads: int,
):
    """
    Fetch the full lineage for each taxonomy id returned in a similarity search, check if any
    taxonomy id in the lineage is regulated (filtering out synthetic constructs), and return a new
    dataframe with taxonomy information.
    """
    TAXIDS_COL = "subject tax ids"

    # prevent truncation of taxonomy results
    pd.set_option("display.max_colwidth", None)

    blast = _split_by_tax_id(blast, TAXIDS_COL)

    # Add new columns, which will later be used to classify the hits as regulated or non-regulated
    blast["regulated"] = False
    blast["superkingdom"] = ""
    blast["phylum"] = ""
    blast["genus"] = ""
    blast["species"] = ""

    blast = blast[blast[TAXIDS_COL] != TAXID_SYNTHETIC_CONSTRUCTS]
    blast = blast[blast[TAXIDS_COL] != TAXID_VECTORS]
    blast = blast.reset_index(drop=True)

    lin = _get_lineages(blast[TAXIDS_COL], db_path, threads)

    # Check if any rows will be removed due to not finding a valid lineage for them
    rows_to_remove = blast[~blast[TAXIDS_COL].isin(lin["TaxID"])]
    if not rows_to_remove.empty:
        logging.warning(
            "Removing %i rows from BLAST results due to invalid taxID(s): %s"
            " - check that taxonomy and protein databases are up to date!",
            len(rows_to_remove),
            ", ".join(map(str, rows_to_remove[TAXIDS_COL].unique())),
        )
    # Filter to only those rows which have a matching taxonomic lineage
    blast = blast[blast[TAXIDS_COL].isin(lin["TaxID"])]

    # Process each hit
    rows_to_drop = []
    for index, row in blast.iterrows():
        row_lin = lin[lin["TaxID"] == row[TAXIDS_COL]].iloc[0]
        full_lineage = pd.DataFrame(
            {
                "Lineage": row_lin["FullLineage"].split(";"),
                "TaxID": row_lin["FullLineageTaxIDs"].split(";"),
                "Rank": row_lin["FullLineageRanks"].split(";"),
            }
        )
        full_lineage.set_index("Rank", inplace=True)
        full_lineage_taxids = list(map(str, full_lineage["TaxID"]))

        # If any organism in the lineage is synthetic, drop the row
        if any(
            taxid in [str(TAXID_SYNTHETIC_CONSTRUCTS), str(TAXID_VECTORS)]
            for taxid in full_lineage_taxids
        ):
            rows_to_drop.append(index)
            continue

        # If any organism in the lineage is regulated, set this hit as regulated
        if any(taxid in regulated_taxids for taxid in full_lineage_taxids):
            blast.at[index, "regulated"] = True

        # Unless we're dealing with a known vaccine strain
        if any(taxid in vaccine_taxids for taxid in full_lineage_taxids):
            blast.at[index, "regulated"] = False

        # Set additional taxonomic information
        for rank in ["superkingdom", "phylum", "genus", "species"]:
            if rank in full_lineage.index:
                blast.at[index, rank] = full_lineage.loc[rank, "Lineage"]
            else:
                blast.at[index, rank] = ""

    blast = blast.drop(rows_to_drop)
    blast = blast.sort_values(by=["% identity"], ascending=False)
    blast = blast.reset_index(drop=True)

    return blast


def read_blast(blast_file: Union[str, os.PathLike, BinaryIO, TextIO]) -> pd.DataFrame:
    """
    Read in BLAST/DIAMOND files and pre-format the data frame with essential info
    """
    blast = pd.read_csv(blast_file, sep="\t", comment="#", header=None)
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
    blast["q. coverage"] = abs(blast["q. end"] - blast["q. start"]) / blast["query length"].max()
    blast["s. coverage"] = abs(blast["s. end"] - blast["s. start"]) / blast["subject length"]

    blast = blast[blast["subject tax ids"].notna()]
    blast = blast.reset_index(drop=True)

    return blast


def _trim_overlapping(blast: pd.DataFrame):
    """
    Remove any hits that are completely overlapped by another, higher-quality hit.
    """
    # set start to the lowest coordinate and end to the highest to avoid confusion
    blast = shift_hits_pos_strand(blast)

    # if any multispecies hits contain regulated pathogens, put the regulated up top
    if "regulated" in blast:
        blast = blast.sort_values(by=["regulated"], ascending=False)

    # rank hits by percent identity, then bit score
    blast = blast.sort_values(by=["% identity", "bit score"], ascending=False)
    blast = blast.reset_index(drop=True)

    blast2 = blast
    # only keep top-ranked hits that don't overlap
    for query in blast["query acc."].unique():
        df = blast[blast["query acc."] == query]
        for i in df.index:  # run through each hit from the top
            for j in df.index[(i + 1) :]:  # compare to each below
                if j in blast2.index:
                    # if beginning and end of the higher-rank hit both overlap or extend further
                    # than the beginning and end of lower-ranked hit, discard the lower-ranked hit
                    if (
                        df.loc[i, "q. start"] <= df.loc[j, "q. start"]
                        and df.loc[i, "q. end"] >= df.loc[j, "q. end"]
                    ):
                        # Unless the hits have the same coordinates and % identity
                        if (
                            df.loc[i, "q. start"] < df.loc[j, "q. start"]
                            or df.loc[i, "q. end"] > df.loc[j, "q. end"]
                            or df.loc[i, "% identity"] > df.loc[j, "% identity"]
                        ):
                            blast2 = blast2.drop([j])
    blast2 = blast2.reset_index(drop=True)

    return blast2

def shift_hits_pos_strand(blast):
    for j in blast.index:
        if blast.loc[j, "q. start"] > blast.loc[j, "q. end"]:
            start = blast.loc[j, "q. end"]
            end = blast.loc[j, "q. start"]
            blast.loc[j, "q. start"] = start
            blast.loc[j, "q. end"] = end
    return blast

def _trim_edges(df):

    # Use this enumeration, instead of below, for pylint errors - needs testing first.
    #for top, i in enumerate(df.index):  # run through each hit from the top
    #    for next, j in enumerate(df.index[top + 1:], start=top + 1):  # compare to each below

    for top in range(len(df.index)):  # run through each hit from the top
        i = df.index[top]
        for next in range(top + 1, len(df.index)):  # compare to each below
            j = df.index[next]
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
                elif i_end == j_end and df.loc[j, "% identity"] == df.loc[i, "% identity"]:
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
                elif i_start == j_start and df.loc[j, "% identity"] == df.loc[i, "% identity"]:
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


def get_top_hits(blast: pd.DataFrame):
    """
    Trim BLAST results down to the top hit for each base.
    """
    top_hits = _trim_overlapping(blast)
    top_hits = top_hits.sort_values("% identity", ascending=False)

    # only keep coordinates of each hit that are not already covered by a better hit
    for query in top_hits["query acc."].unique():
        df = top_hits[top_hits["query acc."] == query]

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


def get_high_identity_matches(blast_output_file, threshold=90):
    """Read all hits with high sequence identity from a BLAST results file."""
    hits = read_blast(blast_output_file)
    hits = _trim_overlapping(hits)
    return hits[hits["% identity"] >= threshold]
