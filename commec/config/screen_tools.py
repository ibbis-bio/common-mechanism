#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

"""
Container for search handlers used throughout the Commec screen workflow.
Sets and alters defaults based on input parameters.
"""

import logging
from typing import Union
from commec.config.io_parameters import ScreenIOParameters
from commec.tools.blastn import BlastNHandler
from commec.tools.blastx import BlastXHandler
from commec.tools.diamond import DiamondHandler
from commec.tools.cmscan import CmscanHandler
from commec.tools.hmmer import HmmerHandler

class ScreenTools:
    """
    Using a set of `ScreenIoParameters`, set up the tools needed to search datbases.
    """

    def __init__(self, params: ScreenIOParameters):
        self.biorisk_hmm: HmmerHandler = None
        self.regulated_protein : Union[BlastXHandler, DiamondHandler] = None
        self.regulated_nt: BlastNHandler = None
        self.benign_hmm: HmmerHandler = None
        self.benign_blastn: BlastNHandler = None
        self.benign_cmscan: CmscanHandler = None

        config_file = params.yaml_configuration

        self.biorisk_hmm = HmmerHandler(
            config_file["databases"]["biorisk_hmm"]["path"],
            params.query.aa_path,
            f"{params.output_prefix}.biorisk.hmmscan",
            threads=params.config.threads,
            force=params.config.force,
        )

        if params.should_do_protein_screening:
            if params.config.protein_search_tool == "blastx":
                self.regulated_protein = BlastXHandler(
                    config_file["databases"]["regulated_protein"]["blast"]["path"],
                    input_file=params.query.nt_path,
                    out_file=f"{params.output_prefix}.nr.blastx",
                    threads=params.config.threads,
                    force=params.config.force,
                )
            elif params.config.protein_search_tool in ("nr.dmnd", "diamond"):
                self.regulated_protein = DiamondHandler(
                    config_file["databases"]["regulated_protein"]["diamond"]["path"],
                    input_file=params.query.nt_path,
                    out_file=f"{params.output_prefix}.nr.dmnd",
                    threads=params.config.threads,
                    force=params.config.force,
                )
                self.regulated_protein.jobs = params.config.diamond_jobs
                if params.config.protein_search_tool == "nr.dmnd":
                    logging.info(
                        "Using old \"nr.dmnd\" keyword for search tool will not be supported"
                        " in future releases,consider using \"diamond\" instead."
                    )
            else:
                raise RuntimeError('Search tool not defined as "blastx" or "diamond"')

        if params.should_do_nucleotide_screening:
            self.regulated_nt = BlastNHandler(
                config_file["databases"]["regulated_nt"]["path"],
                input_file=f"{params.output_prefix}.noncoding.fasta",
                out_file=f"{params.output_prefix}.nt.blastn",
                threads=params.config.threads,
                force=params.config.force,
            )

        if params.should_do_benign_screening:
            self.benign_hmm = HmmerHandler(
                config_file["databases"]["benign"]["hmm"]["path"],
                input_file=params.query.aa_path,
                out_file=f"{params.output_prefix}.benign.hmmscan",
                threads=params.config.threads,
                force=params.config.force,
            )
            self.benign_blastn = BlastNHandler(
                config_file["databases"]["benign"]["fasta"]["path"],
                input_file=params.query.nt_path,
                out_file=f"{params.output_prefix}.benign.blastn",
                threads=params.config.threads,
                force=params.config.force,
            )
            self.benign_cmscan = CmscanHandler(
                config_file["databases"]["benign"]["cm"]["path"],
                input_file=params.query.nt_path,
                out_file=f"{params.output_prefix}.benign.cmscan",
                threads=params.config.threads,
                force=params.config.force,
            )
