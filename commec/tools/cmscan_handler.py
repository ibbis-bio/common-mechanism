#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `Input and Output Screen Parameters` class, and associated dataclasses.
"""
import re
import pandas as pd
from commec.tools.database_handler import DatabaseHandler, DatabaseVersion

class CmscanHandler(DatabaseHandler):
    """ A Database handler specifically for use with Hmmer files for commec screening. """
    def screen(self):
        command = [
            "cmscan", 
            "--tblout",
            self.out_file,
            self.db_file,
            self.input_file
            ]
        self.run_as_subprocess(command, self.temp_log_file)

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

    with open(fileh, "r", encoding = "utf-8") as f:
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

    #    print(cmscan)
    return cmscan
