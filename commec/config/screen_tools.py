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
from commec.setup import get_latest_commec_database_release_tag

logger = logging.getLogger(__name__)

class ScreenTools:
    """
    Using parameters and filenames in `ScreenIo`, set up the tools needed to search datbases.
    """

    def __init__(self, params: ScreenIO):
        self.biorisk_hmm: HmmerHandler = None
        self.regulated_protein : BlastXHandler | DiamondHandler = None
        self.regulated_nt: BlastNHandler = None
        self.benign_hmm: HmmerHandler = None
        self.benign_blastn: BlastNHandler = None
        self.benign_cmscan: CmscanHandler = None

        self.taxonomy_path: str | os.PathLike = None
        self.benign_taxid_path: str | os.PathLike = None
        self.biorisk_taxid_path: str | os.PathLike = None

        self.biorisk_annotations_csv: str | os.PathLike = None

        # Check the existance of a version.txt file. ... 
        validate_commec_database_versions(
            params.config["base_paths"]["default"] + "commec-db-version.txt")

        # Paths for vaxid, taxids, and taxonomy directory, used for check_regulated_pathogens
        # (Declared this way for backwards compatibility to old database structure at this stage)
        value = params.config.get("databases", {}).get("taxonomy", {}).get("path")
        self.taxonomy_path = value if value is not None else params.db_dir + "/taxonomy/"

        value = params.config.get("databases", {}).get("taxonomy", {}).get("regulated_taxids")
        self.biorisk_taxid_path = value if value is not None else os.path.join(
            params.config["databases"]["biorisk_hmm"]["path"],"reg_taxids.txt")

        value = params.config.get("databases", {}).get("taxonomy", {}).get("benign_taxids")
        self.benign_taxid_path = value if value is not None else os.path.join(
            params.config["databases"]["benign"]["hmm"]["path"],"vax_taxids.txt")
        
        self.biorisk_annotations_csv = params.config["databases"]["biorisk_hmm"]["annotations"]

        # Database tools for Biorisks / Protein and NT screens / Benign screen:
        self.biorisk_hmm = HmmerHandler(
            params.config["databases"]["biorisk_hmm"]["path"],
            params.aa_path,
            f"{params.output_prefix}.biorisk.hmmscan",
            threads=params.config["threads"],
            force=params.config["force"],
        )

        if params.should_do_protein_screening:
            if params.config["protein_search_tool"] == "blastx":
                self.regulated_protein = BlastXHandler(
                    params.config["databases"]["regulated_protein"]["blast"]["path"],
                    input_file=params.nt_path,
                    out_file=f"{params.output_prefix}.nr.blastx",
                    threads=params.config["threads"],
                    force=params.config["force"],
                )
            elif params.config["protein_search_tool"] in ("nr.dmnd", "diamond"):
                self.regulated_protein = DiamondHandler(
                    params.config["databases"]["regulated_protein"]["diamond"]["path"],
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
                params.config["databases"]["regulated_nt"]["path"],
                input_file=params.nc_path,
                out_file=f"{params.output_prefix}.nt.blastn",
                threads=params.config["threads"],
                force=params.config["force"],
            )

        if params.should_do_benign_screening:
            self.benign_hmm = HmmerHandler(
                params.config["databases"]["benign"]["hmm"]["path"],
                input_file=params.aa_path,
                out_file=f"{params.output_prefix}.benign.hmmscan",
                threads=params.config["threads"],
                force=params.config["force"],
            )
            self.benign_blastn = BlastNHandler(
                params.config["databases"]["benign"]["fasta"]["path"],
                input_file=params.nt_path,
                out_file=f"{params.output_prefix}.benign.blastn",
                threads=params.config["threads"],
                force=params.config["force"],
            )
            self.benign_cmscan = CmscanHandler(
                params.config["databases"]["benign"]["cm"]["path"],
                input_file=params.nt_path,
                out_file=f"{params.output_prefix}.benign.cmscan",
                threads=params.config["threads"],
                force=params.config["force"],
            )

def validate_commec_database_versions(version_filepath: str | os.PathLike):
    """
    Validate whether the local COMMEC database version is up-to-date.

    This function reads the local version of the COMMEC database from a `version.txt`
    file located at `version_filepath`. It also retrieves the most recent available
    version from GitHub via `get_latest_commec_database_release_tag()`.

    If either the local or remote version cannot be determined, or if the local version
    is out of date, a warning will be logged.

    Parameters
    ----------
    version_filepath : str or os.PathLike
        Path to the local `commec-db-version.txt` file containing the COMMEC database version.

    Logs
    ----
    - A warning if the local version file does not exist or cannot be read.
    - A warning if the latest version cannot be fetched (e.g., due to lack of internet).
    - A warning if the local version is not the most recent.
    """
    local_version: str | None = None
    if os.path.exists(version_filepath):
        try:
            with open(version_filepath, "r", encoding="utf-8") as f:
                local_version = f.readline().strip()
        except (FileNotFoundError, PermissionError, OSError, UnicodeDecodeError) as e:
            logger.warning("Failed to read local version file '%s': %s",
                           version_filepath, e)

    most_recent_version, _reason = get_latest_commec_database_release_tag()

    if not local_version:
        logger.warning("Local COMMEC database version not found or unreadable at '%s'.", version_filepath)
    if not most_recent_version:
        logger.warning("Unable to determine the latest COMMEC database version, %s ", _reason)

    if most_recent_version != local_version:
        logger.warning(
            "COMMEC database outdated. Local version: '%s', Latest version: '%s'. "
            "\nIt is strongly recommend to download the latest version by running commec setup, "
            "or directly downloading the latest from the commec-databases repository: "
            "[https://github.com/ibbis-bio/commec-databases/releases].",
            local_version, most_recent_version
        )
