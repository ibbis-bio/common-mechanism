#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Module for a hidden markov model handler, specifically for calling hmmscan command line interface.
Additional methods for reading hmmscan output, readhmmer, which returns a pandas database.
Instantiate a HmmerHandler, with input local database, input fasta, and output file.
Throws if inputs are invalid. Creates a temporary log file, which is deleted on completion.
"""
import re
import subprocess
import pandas as pd
import itertools
from commec.tools.search_handler import SearchHandler, SearchToolVersion
from commec.utils.coordinates import convert_protein_to_nucleotide_coords


class HmmerHandler(SearchHandler):
    """A Database handler specifically for use with Hmmer files for commec screening."""

    def _search(self):
        command = [
            "hmmscan",
            "--cpu",
            str(self.threads),
            "--domtblout",
            self.out_file,
            self.db_file,
            self.input_file,
        ]
        self.run_as_subprocess(command, self.temp_log_file)

    def get_version_information(self) -> SearchToolVersion:
        """
        At the moment this is just grabbing some basic header info of the
        first entrant of the hmm database. Not really a true version control.
        But better than nothing at the moment. There may be some way to return
        some version information from hmmcan itself, look into that.
        """
        database_info: str = None
        try:
            with open(self.db_file, "r", encoding="utf-8") as file:
                for line in file:
                    if line.startswith("HMMER3/f"):
                        database_info = line.split("[", maxsplit=1)
                        continue
                    # Early exit if data has been found
                    if database_info:
                        break

            tool_version_result = subprocess.run(
                ["hmmscan", "-h"], capture_output=True, text=True, check=True
            )
            tool_info: str = tool_version_result.stdout.splitlines()[1]
            return SearchToolVersion(tool_info, database_info)

        except subprocess.CalledProcessError:
            return None


def readhmmer(fileh):
    """
    Read in HMMER output files
    """
    columns = [
        "target name",
        "accession",
        "tlen",
        "query name",
        " accession",
        "qlen",
        "E-value",
        "score",
        "bias",
        "hit #",
        "of",
        "c-Evalue",
        "i-Evalue",
        "score2",
        "bias",
        "hmm from",
        "hmm to",
        "ali from",
        "ali to",
        "env from",
        "env to",
        "acc",
        "description of target",
    ]

    hmmer = []

    with open(fileh, "r", encoding="utf-8") as f:
        for line in f:
            if "# Program:         hmmscan" in line:
                break
            if "#" in line:
                continue
            bits = re.split(r"\s+", line)
            description = " ".join(bits[22:])
            bits = bits[:22]
            bits.append(description)
            hmmer.append(bits)
    hmmer = pd.DataFrame(hmmer, columns=columns)
    hmmer["E-value"] = pd.to_numeric(hmmer["E-value"])
    hmmer["score"] = pd.to_numeric(hmmer["score"])
    hmmer["ali from"] = pd.to_numeric(hmmer["ali from"])
    hmmer["ali to"] = pd.to_numeric(hmmer["ali to"])
    hmmer["qlen"] = pd.to_numeric(hmmer["qlen"])
    # Extract the frame information.
    hmmer["frame"] = hmmer["query name"].str.split('_').str[-1].astype(int)
    return hmmer

def remove_overlaps(hmmer : pd.DataFrame) -> pd.DataFrame:
    """
    Trims verbosity of a HMMER output, 
    by removing weaker hits which are 
    encompassed in their extent by higher scoring hits.

    Note, works to trim nucleotide coordinates relative to the query, 
    not ali from and ali to from the HMMER itself.

    This means it can be used on any DataFrame with the q. start and q. end NT headings.
    (Consider moving to a general coordinates tool function?)
    """
    assert "q. start" in hmmer.columns, ("No \"q. start\" heading in HMMER output dataframe being "
                                         "passed to remove overlaps, ensure that the dataframe has "
                                         "been processed for converstion to nucleotide coordinates.")

    assert "q. end" in hmmer.columns, ("No \"q. end\" heading in HMMER output dataframe being "
                                         "passed to remove overlaps, ensure that the dataframe has "
                                         "been processed for converstion to nucleotide coordinates.")

    trimmed_hmmer = hmmer # Direct Assignment, reassigned later with .drop() for deep-copy.

    # Ensure all logic is performed per unique Query name.
    for query in hmmer["query name"].unique():

        hmmer_for_query = hmmer[hmmer["query name"] == query]
        sorted_values = hmmer_for_query.sort_values(by=["score"], ascending = False)

        for i, j in itertools.combinations(sorted_values.index, 2):
            # If J is encapsulated:
            if (sorted_values.loc[i, "q. start"] <= sorted_values.loc[j, "q. start"]
                and sorted_values.loc[i, "q. end"] >= sorted_values.loc[j, "q. end"]
                and sorted_values.loc[i, "score"] >= sorted_values.loc[j, "score"]):
                if j in trimmed_hmmer.index:
                    trimmed_hmmer = trimmed_hmmer.drop([j])
                    continue
            # If I is encapsulated:
            if (sorted_values.loc[i, "q. start"] >= sorted_values.loc[j, "q. start"]
                and sorted_values.loc[i, "q. end"] <= sorted_values.loc[j, "q. end"]
                and sorted_values.loc[i, "score"] <= sorted_values.loc[j, "score"]):
                if i in trimmed_hmmer.index:
                    trimmed_hmmer = trimmed_hmmer.drop([i])

    # Tidy the output indices.
    trimmed_hmmer = trimmed_hmmer.reset_index(drop=True)

    return trimmed_hmmer

def recalculate_hmmer_query_coordinates(hmmer : pd.DataFrame):
    """
    Recalculate the coordinates of the hmmer database , such that each translated frame
    reverts to original nucleotide coordinates.
    """
    query_start, query_end = convert_protein_to_nucleotide_coords(
        hmmer["frame"],
        hmmer["ali from"],
        hmmer["ali to"],
        hmmer["qlen"]*3) #TODO: Update this to the scalar for each query, else risk off-by-1 NT reporting.

    hmmer["q. start"] = pd.Series(query_start, dtype="int64")
    hmmer["q. end"] = pd.Series(query_end, dtype="int64")