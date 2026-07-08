#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Container for search handlers used throughout the Commec screen workflow.
Sets and alters defaults based on input parameters.
"""

import logging
import os
from commec.config.screen_io import ScreenIO
from commec.tools.blastn import BlastNHandler
from commec.tools.blastx import BlastXHandler
from commec.tools.diamond import DiamondHandler
from commec.tools.cmscan import CmscanHandler
from commec.tools.hmmer import HmmerHandler

logger = logging.getLogger(__name__)

class ScreenTools:
    """
    Using parameters and filenames in `ScreenIo`, set up the tools needed to search datbases.
    """

    def __init__(self, params: ScreenIO):
        self.biorisk: HmmerHandler = None
        self.regulated_protein : BlastXHandler | DiamondHandler = None
        self.regulated_nt: BlastNHandler = None
        self.low_concern_hmm: HmmerHandler = None
        self.low_concern_blastn: BlastNHandler = None
        self.low_concern_cmscan: CmscanHandler = None

        # Paths for annotations still used in the Biorisk, and Low Concern databases.
        self.biorisk_annotations = params.config["databases"]["biorisk"]["annotations"]
        self.low_concern_annotations = params.config["databases"]["low_concern"]["annotations"]

        # Database tools for Biorisks / Protein and NT screens / Benign screen:
        self.biorisk = HmmerHandler(
            params.config["databases"]["biorisk"]["path"],
            params.aa_path,
            f"{params.output_prefix}.biorisk.hmmscan",
            threads=params.config["threads"],
            force=params.config["force"],
        )

        if params.should_do_protein_screening:
            if params.config["protein_search_tool"] == "blastx":
                self.regulated_protein = BlastXHandler(
                    params.config["databases"]["best_match"]["protein"]["path"],
                    input_file=params.nt_path,
                    out_file=f"{params.output_prefix}.nr.blastx",
                    threads=params.config["threads"],
                    force=params.config["force"],
                )
                self.regulated_protein.arguments_dictionary["-mt_mode"] = params.config["blast_mt_mode"]
            elif params.config["protein_search_tool"] in ("nr.dmnd", "diamond"):
                self.regulated_protein = DiamondHandler(
                    params.config["databases"]["best_match"]["protein"]["path"],
                    input_file=params.nt_path,
                    out_file=f"{params.output_prefix}.nr.dmnd",
                    threads=params.config["threads"],
                    force=params.config["force"],
                )
                self.regulated_protein.jobs = params.config["diamond_jobs"]
                if params.config["protein_search_tool"] == "nr.dmnd":
                    logger.info(
                        "Using old \"nr.dmnd\" keyword for search tool will not be supported"
                        " in future releases,consider using \"diamond\" instead."
                    )
            else:
                raise RuntimeError('Search tool not defined as "blastx" or "diamond"')

        if params.should_do_nucleotide_screening:
            self.regulated_nt = BlastNHandler(
                params.config["databases"]["best_match"]["nucleotide"]["path"],
                input_file=params.nc_path,
                out_file=f"{params.output_prefix}.nt.blastn",
                threads=params.config["threads"],
                force=params.config["force"],
            )
            self.regulated_nt.arguments_dictionary["-mt_mode"] = params.config["blast_mt_mode"]

        if params.should_do_low_concern_screening:
            self.low_concern_hmm = HmmerHandler(
                params.config["databases"]["low_concern"]["protein"]["path"],
                input_file=params.aa_path,
                out_file=f"{params.output_prefix}.low_concern.hmmscan",
                threads=params.config["threads"],
                force=params.config["force"],
            )
            self.low_concern_blastn = BlastNHandler(
                params.config["databases"]["low_concern"]["dna"]["path"],
                input_file=params.nt_path,
                out_file=f"{params.output_prefix}.low_concern.blastn",
                threads=params.config["threads"],
                force=params.config["force"],
            )
            self.low_concern_blastn.arguments_dictionary["-mt_mode"] = params.config["blast_mt_mode"]
            self.low_concern_cmscan = CmscanHandler(
                params.config["databases"]["low_concern"]["rna"]["path"],
                input_file=params.nt_path,
                out_file=f"{params.output_prefix}.low_concern.cmscan",
                threads=params.config["threads"],
                force=params.config["force"],
            )