#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Cmscan search handler, and calling cmscan command line interface.
Additional methods for reading handler output, readcmscan, which returns a pandas database.
Instantiate a CmscanHandler, with input local database, input fasta, and output file.
Throws if inputs are invalid. Creates a temporary log file, which is deleted on completion.
"""
import subprocess
import re
import pandas as pd
from commec.tools.search_handler import SearchHandler, SearchToolVersion


class CmscanHandler(SearchHandler):
    """A Database handler specifically for use with Hmmer files for commec screening."""

    def _search(self):
        command = [
            "cmscan",
            "--cpu",
            str(self.threads),
            "--tblout",
            self.out_file,
            self.db_file,
            self.input_file,
        ]
        self.run_as_subprocess(command, self.temp_log_file)
    
    def read_output(self):
        output_dataframe = readcmscan(self.out_file)
        # Standardize the output column names to be like blast:
        output_dataframe = output_dataframe.rename(columns={
            "seq from": "q. start",
            "seq to": "q. end",
            "coverage": "q. coverage",
            "target name": "subject title",
            "mdl from": "s. start",
            "mdl to" : "s. end",
            'E-value': "evalue",
        })
        return output_dataframe
 

    def get_version_information(self) -> SearchToolVersion:
        try:
            database_info = None
            with open(self.db_file, "r", encoding="utf-8") as file:
                for line in file:
                    if line.startswith("INFERNAL1/a"):
                        database_info = line.strip()
                        continue
                    # Early exit if data has been found
                    if database_info:
                        break

            result = subprocess.run(
                ["cmscan", "-h"], capture_output=True, text=True, check=True
            )
            tool_info = result.stdout.splitlines()[1].strip()[2:]

            return SearchToolVersion(tool_info, database_info)

        except subprocess.CalledProcessError:
            return None


def readcmscan(fileh):
    """
    Read in cmscan output files
    """
    columns = [
        "target name",
        "accession",
        "query name",
        "accession",
        "mdl",
        "mdl from",
        "mdl to",
        "seq from",
        "seq to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "E-value",
        "inc",
        "description of target",
    ]

    cmscan = []

    with open(fileh, "r", encoding="utf-8") as f:
        for line in f:
            if "# Program:         cmscan" in line:
                break
            if "#" in line:
                continue
            bits = re.split(r"\s+", line)
            description = " ".join(bits[17:])
            bits = bits[:17]
            bits.append(description)
            cmscan.append(bits)
    cmscan = pd.DataFrame(cmscan, columns=columns)
    cmscan["E-value"] = pd.to_numeric(cmscan["E-value"])
    cmscan["score"] = pd.to_numeric(cmscan["score"])
    cmscan["seq from"] = pd.to_numeric(cmscan["seq from"])
    cmscan["seq to"] = pd.to_numeric(cmscan["seq to"])

    return cmscan
