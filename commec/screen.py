#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Run Common Mechanism screening on an input FASTA.

Screening involves (up to) four steps:

  1. Biorisk scan:      HMM-based scan for matches to a custom database of biorisk sequences.
  2. Protein search:    protein homology search for best matches to regulated pathogens.
  3. Nucleotide search: nucleotide homology search for best matches to regulated pathogens.
  4. Benign scan:       three different scans (against conserved proteins, housekeeping RNAs, and
                        synbio parts) to see if hits identified in homology search can be cleared.

In "fast" mode, only the biorisk scan is run. By default, all four steps are run, but the nucleotide
search is only run for regions that do not have any protein hits with a high sequence identity. The
benign search is not permitted to clear biorisk scan hits, only protein or nucleotide hits. Whether
or not a homology scan hit is from a regulated pathogen is determined by referencing the taxonomy
ids assoicated with each accession that returns a hit, then looking at their lineages.

positional arguments:
  fasta_file            FASTA file to screen

options:
  -h, --help            show this help message and exit
  -d DATABASE_DIR, --databases DATABASE_DIR
                        Path to directory containing reference databases (e.g. taxonomy, protein, HMM)
  -y CONFIG_YAML, --config CONFIG_YAML
                        Configuration for screen run in YAML format, including custom database paths

Screen run logic:
  -f, --fast            Run in fast mode and skip protein and nucleotide homology search
  -p {blastx,diamond}, --protein-search-tool {blastx,diamond}
                        Tool for protein homology search to identify regulated pathogens
  -n, --skip-nt         Skip nucleotide search (regulated pathogens will only be identified based on
                        protein hits)

Parallelisation:
  -t THREADS, --threads THREADS
                        Number of CPU threads to use. Passed to search tools.
  -j DIAMOND_JOBS, --diamond-jobs DIAMOND_JOBS
                        Diamond-only: number of runs to do in parallel on split Diamond databases

Output file handling:
  -o OUTPUT_PREFIX, --output OUTPUT_PREFIX
                        Prefix for output files. Can be a string (interpreted as output basename) or
                        a directory (files will be output there, names determined from input FASTA)
  -c, --cleanup         Delete intermediate output files for run
  -F, --force           Overwrite any pre-existing output for run (cannot be used with --resume)
  -R, --resume          Re-use any pre-existing output run (cannot be used with --force)

"""
import argparse
import datetime
import time
import logging
import os
import shutil
import sys
import pandas as pd
from Bio import SeqIO

from commec.utils.file_utils import file_arg, directory_arg
from commec.utils.json_html_output import generate_html_from_screen_data
from commec.config.io_parameters import ScreenIOParameters
from commec.config.query import Query
from commec.config.screen_tools import ScreenTools

from commec.screeners.check_biorisk import check_biorisk, update_biorisk_data_from_database
from commec.screeners.check_benign import update_benign_data_from_database
from commec.screeners.check_reg_path import (
    check_for_regulated_pathogens,
    update_taxonomic_data_from_database
)

from commec.tools.fetch_nc_bits import calculate_noncoding_regions_per_query

from commec.config.result import (
    ScreenResult,
    ScreenStep,
    QueryResult,
    ScreenStatus,
)

from commec.config.json_io import encode_screen_data_to_json

DESCRIPTION = "Run Common Mechanism screening on an input FASTA."

class ScreenArgumentParser(argparse.ArgumentParser):
    """
    Argument parser that returns a `user_specified_args` namespace item, 
    which helps selectively override other configuration (e.g. provided via YAML)
    i.e. for only when it has explicitly been used as an argument in CLI.

    Importantly, this iterates over all sub-parsers too, required for the 
    cli entry point of Commec. However to do this we access various private 
    parser attributes - which is naughty - but its better than writing our own argsparse.
    """
    def parse_args(self, args=None, namespace=None):     
        # Get argument strings; in most cases, args and sys.argv[1:] will be the same
        cli_strings = args if args is not None else sys.argv[1:]
        user_specified_args = set()

        def collect_user_actions(parser : ScreenArgumentParser):
            """ 
            Recursively collect all actions, including subparsers. 
            _actions has every argument provide to the parser, and
            has every SubParserActions instances.
            """
            for action in parser._actions:
                if isinstance(action, argparse._SubParsersAction):
                    # Recurse into each subparser
                    for _sub_name, subparser in action.choices.items():
                        collect_user_actions(subparser)
                else:
                    for arg_string in action.option_strings:
                        #print("Testing:", arg_string)
                        if arg_string in cli_strings:
                            #print("added!")
                            user_specified_args.add(action.dest)

        # Collect arguments from main parser and all subparsers
        collect_user_actions(self)
        
        ns = super().parse_args(args, namespace)
        setattr(ns,"user_specified_args", user_specified_args)
        return ns

def add_args(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """
    Add module arguments to an ArgumentParser object.
    """

    parser.add_argument(dest="fasta_file", type=file_arg, help="FASTA file to screen")
    parser.add_argument(
        "-d",
        "--databases",
        dest="database_dir",
        type=directory_arg,
        default=None,
        help="Path to directory containing reference databases (e.g. taxonomy, protein, HMM)",
    )
    parser.add_argument(
        "-y",
        "--config",
        dest="config_yaml",
        help="Configuration for screen run in YAML format, including custom database paths",
        default="",
    )
    screen_logic_group = parser.add_argument_group("Screen run logic")
    screen_logic_group.add_argument(
        "-f",
        "--fast",
        dest="in_fast_mode",
        action="store_true",
        help="Run in fast mode and skip protein and nucleotide homology search",
    )
    screen_logic_group.add_argument(
        "-p",
        "-protein-search-tool",
        dest="protein_search_tool",
        choices=["blastx", "diamond"],
        help="Tool for protein homology search to identify regulated pathogens",
    )
    screen_logic_group.add_argument(
        "-n",
        "--skip-nt",
        dest="skip_nt_search",
        action="store_true",
        help="Skip nucleotide search (regulated pathogens will only be identified based on protein hits)",
    )
    parallel_group = parser.add_argument_group("Parallelisation")
    parallel_group.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        help="Number of CPU threads to use. Passed to search tools.",
    )
    parallel_group.add_argument(
        "-j",
        "--diamond-jobs",
        dest="diamond_jobs",
        type=int,
        help="Diamond-only: number of runs to do in parallel on split Diamond databases",
    )
    output_handling_group = parser.add_argument_group("Output file handling")
    output_exclusive_group = output_handling_group.add_mutually_exclusive_group()
    output_handling_group.add_argument(
        "-o",
        "--output",
        dest="output_prefix",
        help="Prefix for output files. Can be a string (interpreted as output basename) or a"
        + " directory (files will be output there, names will be determined from input FASTA)",
        default = ""
    )
    output_handling_group.add_argument(
        "-c",
        "--cleanup",
        dest="do_cleanup",
        action="store_true",
        help="Delete intermediate output files for this Screen run",
    )
    output_exclusive_group.add_argument(
        "-F",
        "--force",
        dest="force",
        action="store_true",
        help="Overwrite any pre-existing output for this Screen run (cannot be used with --resume)",
    )
    output_exclusive_group.add_argument(
        "-R",
        "--resume",
        dest="resume",
        action="store_true",
        help="Re-use any pre-existing output for this Screen run (cannot be used with --force)",
    )
    return parser

class Screen:
    """
    Handles the parsing of input arguments, the control of databases, and
    the logical flow of the screening process for commec.
    """

    def __init__(self):
        self.params : ScreenIOParameters = None
        self.queries : dict[str, Query] = None
        self.database_tools : ScreenTools = None
        self.screen_data : ScreenResult = ScreenResult()
        self.start_time = time.time()

    def __del__(self):
        """ 
        Before we are finished, we attempt to write a JSON and HTML output.
        Doing this in the destructor means that sometimes this will complete
        successfully, despite exceptions.
        """
        time_taken = (time.time() - self.start_time)
        hours, rem = divmod(time_taken, 3600)
        minutes, seconds = divmod(rem, 60)
        self.screen_data.commec_info.time_taken = f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"
        self.screen_data.update()
        encode_screen_data_to_json(self.screen_data, self.params.output_json)
        generate_html_from_screen_data(self.screen_data, self.params.output_prefix+"_summary")
        
        if self.params.config["do_cleanup"]:
            self.params.clean()

    def setup(self, args: argparse.Namespace):
        """Instantiates and validates parameters, and databases, ready for a run."""
        self.params: ScreenIOParameters = ScreenIOParameters(args)
        self.params.setup()

        # Set up logging
        logging.basicConfig(
            level=logging.INFO,
            format="%(message)s",
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler(self.params.output_screen_file, "a"),
                logging.FileHandler(self.params.tmp_log, "a"),
            ],
        )
        logging.basicConfig(
            level=logging.DEBUG,
            format="%(message)s",
            handlers=[logging.FileHandler(self.params.tmp_log, "a")],
        )

        logging.info(" Validating Inputs...")
        self.params.setup()
        self.database_tools: ScreenTools = ScreenTools(self.params)

        # Add the input contents to the log
        shutil.copyfile(self.params.input_fasta_path, self.params.tmp_log)

        # Initialize the queries
        self.queries = self.params.parse_input_fasta()
        for query in self.queries.values():
            query.translate(self.params.nt_path, self.params.aa_path)
            qr = QueryResult(query.original_name,
                                 len(query.seq_record),
                                 str(query.seq_record.seq))
            self.screen_data.queries[query.name] = qr
            query.result_handle = qr
        
        # Initialize the version info for all the databases
        _tools = self.database_tools
        _info = self.screen_data.commec_info
        _info.biorisk_database_info = _tools.biorisk_hmm.get_version_information()
        if self.params.should_do_protein_screening:
            _info.protein_database_info = _tools.regulated_protein.get_version_information()
        if self.params.should_do_nucleotide_screening:
            _info.nucleotide_database_info = _tools.regulated_nt.get_version_information()
        if self.params.should_do_benign_screening:
            _info.benign_protein_database_info = _tools.benign_hmm.get_version_information()
            _info.benign_rna_database_info = _tools.benign_blastn.get_version_information()
            _info.benign_synbio_database_info = _tools.benign_cmscan.get_version_information()

        # Store start time.
        _info.date_run = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def run(self, args : argparse.Namespace):
        """
        Wrapper so that args be parsed in main() or commec.py interface.
        """
        # Perform setup steps.
        self.setup(args)

        self.params.output_yaml(self.params.input_prefix + "_config.yaml")

        # Biorisk screen
        logging.info(">> STEP 1: Checking for biorisk genes...")
        self.screen_biorisks()
        logging.info(
            " STEP 1 completed at %s",
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )

        # Taxonomy screen (Protein)
        if self.params.should_do_protein_screening:
            logging.info(" >> STEP 2: Checking regulated pathogen proteins...")
            self.screen_proteins()
            logging.info(
                " STEP 2 completed at %s",
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
        else:
            logging.info(" SKIPPING STEP 2: Protein search")
            self.reset_protein_recommendations(ScreenStatus.SKIP)

        # Taxonomy screen (Nucleotide)
        if self.params.should_do_nucleotide_screening:
            logging.info(" >> STEP 3: Checking regulated pathogen nucleotides...")
            self.screen_nucleotides()
            logging.info(
                " STEP 3 completed at %s",
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
        else:
            logging.info(" SKIPPING STEP 3: Nucleotide search")
            self.reset_nucleotide_recommendations(ScreenStatus.SKIP)

        # Benign Screen
        if self.params.should_do_benign_screening:
            logging.info(
                ">> STEP 4: Checking any pathogen regions for benign components..."
            )
            self.screen_benign()
            logging.info(
                ">> STEP 4 completed at %s",
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
        else:
            logging.info(" SKIPPING STEP 4: Benign search")
            self.reset_benign_recommendations(ScreenStatus.SKIP)

        logging.info(
            ">> COMPLETED AT %s", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )


    def screen_biorisks(self):
        """
        Call hmmscan` and `check_biorisk.py` to add biorisk results to `screen_file`.
        """
        logging.debug("\t...running hmmscan")
        self.database_tools.biorisk_hmm.search()
        logging.debug("\t...checking hmmscan results")
        check_biorisk(
            self.database_tools.biorisk_hmm.out_file,
            self.database_tools.biorisk_hmm.db_directory,
            self.queries
        )
        update_biorisk_data_from_database(self.database_tools.biorisk_hmm, self.screen_data, self.queries)

    def screen_proteins(self):
        """
        Call `run_blastx.sh` or `run_diamond.sh` followed by `check_reg_path.py` to add regulated
        pathogen protein screening results to `screen_file`.
        """
        logging.debug("\t...running %s", self.params.config["protein_search_tool"])
        self.database_tools.regulated_protein.search()
        if not self.database_tools.regulated_protein.check_output():
            self.reset_protein_recommendations(ScreenStatus.ERROR)
            raise RuntimeError(
                "ERROR: Expected protein search output not created: "
                + self.database_tools.regulated_protein.out_file
            )

        logging.debug(
            "\t...checking %s results", self.params.config["protein_search_tool"]
        )
        #reg_path_coords = f"{self.params.output_prefix}.reg_path_coords.csv"

        # Delete any previous check_reg_path results for this search
        #if os.path.isfile(reg_path_coords):
        #    os.remove(reg_path_coords)

        #check_for_regulated_pathogens(
        #    self.database_tools.regulated_protein.out_file,
        #    self.params.db_dir,
        #    str(self.params.config.threads),
        #)

        update_taxonomic_data_from_database(self.database_tools.regulated_protein,
                                            self.database_tools.benign_taxid_path,
                                            self.database_tools.biorisk_taxid_path,
                                            self.database_tools.taxonomy_path,
                                            self.screen_data,
                                            self.queries,
                                            ScreenStep.TAXONOMY_AA,
                                            self.params.config["threads"])


    def screen_nucleotides(self):
        """
        Call `fetch_nc_bits.py`, search noncoding regions with `blastn` and
        then `check_reg_path.py` to screen regulated pathogen nucleotides in
        noncoding regions (i.e. that would not be found with protein search).
        """
        # Only screen nucleotides in noncoding regions
        #fetch_noncoding_regions(
        #    self.database_tools.regulated_protein.out_file, self.screen_io.nt_path
        #)

        # Calculate non-coding information for each Query.
        calculate_noncoding_regions_per_query(
            self.database_tools.regulated_protein.out_file,
            self.queries)
        
        # Generate the non-coding fasta.
        nc_fasta_sequences = ""
        for query in self.queries.values():
            nc_fasta_sequences += query.get_non_coding_regions_as_fasta()
        
        # Skip if there is no non-coding information.
        if nc_fasta_sequences == "":
            logging.debug(
                "\t...skipping nucleotide search since no noncoding regions fetched"
            )
            self.reset_nucleotide_recommendations(ScreenStatus.SKIP)
            return

        # Create a non-coding fasta file.
        with open(self.params.nc_path, "w", encoding="utf-8") as output_file:
            output_file.writelines(nc_fasta_sequences)

        # Only run new blastn search if there are no previous results
        self.database_tools.regulated_nt.search()

        if not self.database_tools.regulated_nt.check_output():
            self.reset_nucleotide_recommendations(ScreenStatus.ERROR)
            raise RuntimeError(
                "ERROR: Expected nucleotide search output not created: "
                + self.database_tools.regulated_nt.out_file
            )
        
        # Note: Currently noncoding coo rdinataes are converted within update_taxonomic_data_from_database,
        # It may be prudent to instead explictly convert them in the output file itself, or during import.

        logging.debug("\t...checking blastn results")
        #check_for_regulated_pathogens(
        #    self.database_tools.regulated_nt.out_file,
        #    self.params.db_dir,
        #    str(self.params.config["threads"]),
        #)

        update_taxonomic_data_from_database(self.database_tools.regulated_nt,
                                            self.database_tools.benign_taxid_path,
                                            self.database_tools.biorisk_taxid_path,
                                            self.database_tools.taxonomy_path,
                                            self.screen_data,
                                            self.queries,
                                            ScreenStep.TAXONOMY_NT,
                                            self.params.config["threads"])

    def screen_benign(self):
        """
        Call `hmmscan`, `blastn`, and `cmscan` and then pass results
        to `check_benign.py` to identify regions that can be cleared.
        """
        # Start by checking if there are any hits that require clearing...
        hits_to_clear : bool = False
        for _query, hit in self.screen_data.hits():
            if hit.recommendation.status in {ScreenStatus.WARN, ScreenStatus.FLAG}:
                hits_to_clear = True
                break

        if not hits_to_clear:
            logging.info("\t...no regulated regions to clear\n")
            self.reset_benign_recommendations(ScreenStatus.SKIP)
            return

        # Run the benign tools:
        logging.info("\t...running benign hmmscan")
        self.database_tools.benign_hmm.search()
        logging.info("\t...running benign blastn")
        self.database_tools.benign_blastn.search()
        logging.info("\t...running benign cmscan")
        self.database_tools.benign_cmscan.search()


        # Update Screen Data with benign outputs.
        benign_desc = pd.read_csv(
            self.database_tools.benign_hmm.db_directory + "/benign_annotations.tsv",
            sep="\t",
        )

        update_benign_data_from_database(
            self.database_tools.benign_hmm,
            self.database_tools.benign_cmscan,
            self.database_tools.benign_blastn,
            self.queries,
            benign_desc
        )

    def reset_benign_recommendations(self, new_recommendation : ScreenStatus):
        """ Helper function 
        apply a single recommendation to the whole benign step 
        for every query."""
        for query in self.screen_data.queries.values():
            query.recommendation.benign_status = new_recommendation

    def reset_protein_recommendations(self, new_recommendation : ScreenStatus):
        """ Helper function
        apply a single recommendation to the whole protein taxonomy step 
        for every query."""
        for query in self.screen_data.queries.values():
            query.recommendation.protein_taxonomy_status = new_recommendation

    def reset_nucleotide_recommendations(self, new_recommendation : ScreenStatus):
        """ Helper function:
        apply a single recommendation to the whole nucleotide taxonomy step 
        for every query."""
        for query in self.screen_data.queries.values():
            query.recommendation.nucleotide_taxonomy_status = new_recommendation

def run(args: argparse.Namespace):
    """
    Entry point from commec main. Passes args to Screen object, and runs.
    """
    my_screen: Screen = Screen()
    my_screen.run(args)


def main():
    """
    Main function. Passes args to Screen object, which then runs.
    """
    parser = ScreenArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    try:
        main()
    except RuntimeError as e:
        print(f"Runtime error: {e}", file=sys.stderr)
        sys.exit(1)
