#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Defines the `Input and Output Screen Parameters` class, and associated dataclasses.
"""
import os
from commec.io_parameters import ScreenIOParameters
from commec.databases.blastn_db import BlastNDataBase
from commec.databases.blastx_db import BlastXDataBase
from commec.databases.blastdmnd_db import DiamondDataBase
from commec.databases.cmscan_db import CmscanDataBase
from commec.databases.hmm_db import HMMDataBase

class CommecDatabases():
    """
    Consolidation and initialisation of all databases needed for the commec workflow.
    This could be maybe split into a Biorisk_handler, Taxonomy handler, and Benign handler
    class objects respectively.
    """
    def __init__(self, params : ScreenIOParameters):
        # Biorisk
        self.biorisk_db : HMMDataBase = None
        # Taxonomy
        self.protein_db : BlastXDataBase = None
        self.nucleotide_db : BlastNDataBase = None
        # Benign
        self.benign_hmm : HMMDataBase = None
        self.benign_blastn : BlastNDataBase = None
        self.benign_cmscan : CmscanDataBase = None

        if params.should_do_biorisk_screening:
            self.biorisk_db = HMMDataBase(
                os.path.join(params.db_dir, "biorisk_db"),
                os.path.join(params.db_dir, "biorisk_db/biorisk.hmm"),
                params.query.fasta_aa_filepath,
                f"{params.output_prefix}.biorisk.hmmscan",
            )

        if params.should_do_protein_screening:
            if params.inputs.search_tool == "blastx":
                self.protein_db = BlastXDataBase(
                    os.path.join(params.db_dir, "nr_blast"),
                    os.path.join(params.db_dir, "nr_blast/nr"),
                    input_file = params.query.fasta_aa_filepath,
                    out_file = f"{params.output_prefix}.nr.blastx"
                )
            elif params.inputs.search_tool == "nr.dmnd":
                self.protein_db = DiamondDataBase(
                    os.path.join(params.db_dir, "nr_dmnd"),
                    os.path.join(params.db_dir, "nr_blast/nr.dmnd"),
                    input_file = params.query.fasta_aa_filepath,
                    out_file = f"{params.output_prefix}.nr.dmnd"
                )
            else:
                raise RuntimeError("Search tool not defined as \"blastx\" or \"nr.dmnd\"")

        if params.should_do_nucleotide_screening:
            self.nucleotide_db = BlastNDataBase(
                os.path.join(params.db_dir, "nt_blast"),
                os.path.join(params.db_dir, "nt_blast/nt"),
                input_file = f"{params.output_prefix}.noncoding.fasta",
                out_file = f"{params.output_prefix}.nt.blastn"
            )

        if params.should_do_benign_screening:
            self.benign_hmm = HMMDataBase(
                os.path.join(params.db_dir, "benign_db"),
                os.path.join(params.db_dir, "benign_db/benign.hmm"),
                input_file = params.query.cleaned_fasta_filepath,
                out_file = f"{params.output_prefix}.benign.hmmscan"
            )
            self.benign_blastn = BlastNDataBase(
                os.path.join(params.db_dir, "benign_db"),
                os.path.join(params.db_dir, "benign_db/benign.fasta"),
                input_file = params.query.cleaned_fasta_filepath,
                out_file = f"{params.output_prefix}.benign.blastn"
            )
            self.benign_cmscan = CmscanDataBase(
                os.path.join(params.db_dir, "benign_db"),
                os.path.join(params.db_dir, "benign_db/benign.cm"),
                input_file = params.query.cleaned_fasta_filepath,
                out_file = f"{params.output_prefix}.benign.cmscan"
            )
