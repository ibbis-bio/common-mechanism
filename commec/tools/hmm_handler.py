#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Database handler for Hidden Markov Model type databases.
"""
import re
import pandas as pd
from commec.tools.database_handler import DatabaseHandler, DatabaseVersion

class HMMDataBase(DatabaseHandler):
    """ A Database handler specifically for use with Hmmer files for commec screening. """
    def screen(self):
        command = [
            "hmmscan", 
            "--domtblout",
            self.out_file,
            self.db_file,
            self.input_file
            ]
        self.run_as_subprocess(command, self.temp_log_file)

    def get_version_information(self) -> DatabaseVersion:
        """ 
        At the moment this is just grabbing some basic header info of the 
        first entrant of the hmm database. Not really a true version control.
        But better than nothing at the moment. There may be some way to return
        some version information from hmmcan itself, look into that. 
        """
        version : str = None
        date : str = None
        comment : str = None
        try:
            with open(self.db_file, 'r', encoding = "utf-8") as file:
                for line in file:
                    if line.startswith("NAME"):
                        comment = line.split(maxsplit=1)[1].strip()
                        continue
                    if line.startswith("DATE"):
                        date = line.split(maxsplit=1)[1].strip()
                        continue
                    if line.startswith("HMMER3/f"):
                        version = line.split("[",maxsplit=1)[1].split("|")[0].strip()
                        continue
                    # Early exit if data has been found
                    if version and date and comment:
                        break
            return DatabaseVersion(version, date, comment)

        except RuntimeError:
            return super().get_version_information()


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