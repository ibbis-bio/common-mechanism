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
from commec.tools.search_handler import SearchHandler, SearchToolVersion


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

    def read_output(self):
        output_dataframe = readhmmer(self.out_file)
        # Standardize the output column names to be like blast:
        output_dataframe = output_dataframe.rename(columns={
            "ali from": "q. start",
            "ali to": "q. end",
            "coverage": "q. coverage",
            "target name": "subject title",
            "qlen":"query length",
            "hmm from":"s. start",
            "hmm end":"s. end",
            'E-value': "evalue",
        })
        return output_dataframe

    def get_version_information(self) -> SearchToolVersion:
        """
        The first line of the HMM database typically contains creation date
        information, and some version information.
        """
        database_info: str = None
        try:
            with open(self.db_file, "r", encoding="utf-8") as file:
                for line in file:
                    if line.startswith("HMMER3/f"):
                        database_info = line.split(";", maxsplit=1)[0].strip()
                        continue
                    # Early exit if data has been found
                    if database_info:
                        break

            tool_version_result = subprocess.run(
                ["hmmscan", "-h"], capture_output=True, text=True, check=True
            )
            tool_info: str = tool_version_result.stdout.splitlines()[1].strip()
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
    return hmmer

def trimhmmer(hmmer):
    """
    Trim hmmer files.

    Don't forget this is a report on 6-frame translations so coordinates will be complicated.
    """
    # rank hits by bitscore
    hmmer = hmmer.sort_values(by=["score"], ascending=False)
    #     hmmer = hmmer.drop_duplicates(subset=['query acc.', 'q. start', 'q. end'], keep='first', ignore_index=True)

    hmmer2 = hmmer
    # only keep  top ranked hits that don't overlap
    for query in hmmer["query name"].unique():
        df = hmmer[hmmer["query name"] == query]
        for i in df.index:
            for j in df.index[(i + 1) :]:
                if (
                    df.loc[i, "ali from"] <= df.loc[j, "ali from"]
                    and df.loc[i, "ali to"] >= df.loc[j, "ali to"]
                ) | (
                    df.loc[i, "ali from"] >= df.loc[j, "ali from"]
                    and df.loc[i, "ali to"] <= df.loc[j, "ali to"]
                ):
                    if j in hmmer2.index:
                        hmmer2 = hmmer2.drop([j])
        hmmer2 = hmmer2.reset_index(drop=True)
    return hmmer2
