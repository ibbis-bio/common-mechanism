#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Parse all .screen files in a directory and create three CSV files reflecting the results of screening.

positional arguments:
  directory            top-level directory containing screen files

options:
  -r, --recursive      show this help message and exit
  -o, --output OUTPUT_DIR
                       path where output CSVs should be written

The screen_pipline_status covers the status of the 4 screening steps, while flags.csv has just
flag outcomes for various things that commec can flag, and flags_recommended just has an overall
flag recommendation for each screen file.
"""

import argparse
import glob
import os
import re
from enum import StrEnum
import pandas as pd
from commec.utils.file_utils import directory_arg

DESCRIPTION = (
    "Parse all .screen files in a directory and create two CSVs file of flags raised"
)

class Outcome(StrEnum):
    """Possible outcomes for each step of the screening process."""

    FLAG = "flag"  # the sequence was flagged in this step
    PASS = "pass"  # the sequence passed in this step
    SKIP = "skip"  # step was intentionally not run
    ERROR = "error"  # an error occurred during this step
    NOT_RUN = "-"  # step was not run due to an error, interrupt, or other unexpected outcome
    # (protein / nucloetide only) equally-good best match to regulated and non-regulated organisms
    MIX = "mix"
    # (biorisk only) significant hit to a virulence factor, but these are often shared between
    # regulated and non-regulated organisms
    WARN = "warn"
    CLEARED = "cleared"  # (benign only) all earlier flags were cleared
    NOT_CLEARED = "not-cleared"  # (benign only) not all earlier flags were cleared


class Flag(StrEnum):
    """Possible values for the flags.csv summary."""

    FLAG = "F"
    PASS = "P"
    ERROR = "Err"
    NOT_RUN = "-"


def add_args(parser: argparse.ArgumentParser):
    """
    Add module arguments to an ArgumentParser object.
    """
    parser.add_argument(
        "directory",
        action="store",
        type=directory_arg,
        help="Directory containing .screen files to summarize",
    )
    parser.add_argument(
        "-r",
        "--recursive",
        dest="recursive",
        action="store_true",
        help="Search directory recursively for screen files",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=directory_arg,
        default=None,
        help="Output directory name (defaults to directory if not provided)",
    )
    return parser


def process_file(file_path) -> dict[str, str | set[str] | bool]:
    """
    Read input screen file, split into steps, and prepare dict of results for CSV output.
    """
    with open(file_path, "r", encoding="utf-8") as file:
        content = file.read()

    filename_without_extension = os.path.splitext(os.path.basename(file_path))[0]

    # Don't capture any of the file before STEP 1
    first_step = re.search(r">> STEP 1:", content)
    # Split by STEP
    if first_step:
        content = content[first_step.start() :]
        steps = re.split(r"(?=(?:>> STEP \d|SKIPPING STEP \d))", content)
    # Strip whitespace and pad out to 4 steps in total
    steps = [step.strip() for step in steps if step.strip()]
    steps = steps + [""] * (4 - len(steps))

    results = {
        "name": filename_without_extension,
        "filepath": file_path,
        "flag": None,
        "biorisk": get_biorisk_outcome(steps[0]),
        "protein": get_protein_outcome(steps[1]),
        "nucleotide": get_nucleotide_outcome(steps[2]),
        "benign": get_benign_outcome(steps[3]),
    }

    # Add regulated flags - only check protein, since nucleotide may have been skipped due to
    # having no noconding regions
    if step_ran_successfully(results["protein"]):
        results.update(get_regulated_taxa(content))

    # Add info about benign sub-steps
    if step_ran_successfully(results["benign"]):
        results.update(get_benign_substeps(steps[3]))

    results["flag"] = get_overall_flag(results)

    return results


def get_biorisk_outcome(step_content: str) -> set[str]:
    """Process biorisk scan step from .screen file."""
    outcomes = set()
    if "FLAG" in step_content:
        outcomes.add(Outcome.FLAG)
    if "Virulence factor found" in step_content:
        outcomes.add(Outcome.WARN)
    if (
        "Biorisks: no hits detected, PASS" in step_content
        or "Biorisks: no significant hits detected, PASS" in step_content
    ):
        outcomes.add(Outcome.PASS)
    if "ERROR" in step_content:
        outcomes.add(Outcome.ERROR)
    return outcomes if outcomes else {Outcome.NOT_RUN}


def get_protein_outcome(step_content: str) -> set[str]:
    """Process protein scan step from .screen file."""
    outcomes = set()
    if "FAST MODE: Skipping" in step_content or "SKIPPING STEP 2" in step_content:
        outcomes.add(Outcome.SKIP)

    if "Best match to sequence(s)" in step_content and "FLAG" in step_content:
        outcomes.add(Outcome.FLAG)
    if "found in both regulated and non-regulated organisms" in step_content:
        outcomes.add(Outcome.MIX)
    if (
        "no top hit exclusive to a regulated pathogen: PASS" in step_content
        and "mix" not in outcomes
    ):
        outcomes.add(Outcome.PASS)
    if "ERROR" in step_content:
        outcomes.add(Outcome.ERROR)
    return outcomes if outcomes else {Outcome.NOT_RUN}


def get_nucleotide_outcome(step_content: str) -> set[str]:
    """Process nucleotide scan step from .screen file."""
    outcomes = set()
    if (
        "skipping nt scan" in step_content
        or "skipping nucleotide search" in step_content
        or "FAST MODE: Skipping" in step_content
        or "SKIPPING STEP 3" in step_content
    ):
        outcomes.add(Outcome.SKIP)

    if "Best match to sequence(s)" in step_content and "FLAG" in step_content:
        outcomes.add(Outcome.FLAG)
    if "no top hit exclusive to a regulated pathogen: PASS" in step_content:
        outcomes.add(Outcome.PASS)
    if "ERROR" in step_content:
        outcomes.add(Outcome.ERROR)
    return outcomes if outcomes else {Outcome.NOT_RUN}


def get_benign_outcome(step_content) -> set[str]:
    """Process benign scan step from .screen file."""
    outcomes = set()
    if "no regulated regions to clear" in step_content:
        outcomes.add(Outcome.SKIP)

    if "Regulated region at bases" in step_content and "failed to clear: FLAG" in step_content:
        outcomes.add(Outcome.NOT_CLEARED)
    if "all regulated regions cleared: PASS" in step_content:
        outcomes.add(Outcome.CLEARED)
    if "ERROR" in step_content:
        outcomes.add(Outcome.ERROR)

    return outcomes if outcomes else {Outcome.NOT_RUN}


def get_regulated_taxa(content) -> dict[str, bool]:
    """
    Check for regulated virus, bacteria, and eukaryote flags in the content.
    """
    return {
        "virus_flag": True if "FLAG (virus)" in content else False,
        "bacteria_flag": True if "FLAG (bacteria)" in content else False,
        "eukaryote_flag": True if "FLAG (eukaryote)" in content else False,
    }


def get_benign_substeps(step_content) -> dict[str, bool]:
    """
    Check for whether benign protiein, RNA or DNA (synbio sequences) is in the content.
    """
    return {
        "benign_protein": True if re.search(r"-->.*?protein.*?PASS", step_content) else False,
        "benign_rna": True if re.search(r"-->.*?RNA.*?PASS", step_content) else False,
        "benign_dna": True if re.search(r"-->.*?Synbio sequence.*?PASS", step_content) else False,
    }


def step_ran_successfully(results: set) -> bool:
    """Step did not have an ERROR, SKIP or NOT_RUN outcome."""
    return (
        Outcome.ERROR not in results
        and Outcome.SKIP not in results
        and Outcome.NOT_RUN not in results
    )


def get_overall_flag(results: dict[str]) -> str:
    """
    Get a single overall recommendation (pass, flag or error) based on the flags.

    We don't take into account virulence factor flags here; those are just a warning, since
    we've found that some virulence factors are housekeeping genes (and are poorly-supported
    Victors entries).
    """
    # Errors mean we can't determine overall pass or flag
    if any(
        [
            Outcome.ERROR in results.get(key)
            for key in ["biorisk", "protein", "nucleotide", "benign"]
        ]
    ):
        return Outcome.ERROR

    # If a biorisk is flagged, flag the whole thing
    if Outcome.FLAG in results["biorisk"]:
        return Outcome.FLAG
    # If flags were cleared, that's a pass
    if Outcome.CLEARED in results["benign"]:
        return Outcome.PASS
    # Flag regulated pathogens unless already cleared in the benign screen
    if Outcome.FLAG in results["protein"] or Outcome.FLAG in results["nucleotide"]:
        return Outcome.FLAG
    # Other outcomes pass
    return Outcome.PASS


def get_flags_from_status(results: dict[str, str | set[str] | bool]) -> dict[str]:
    """
    Map step results to flags format for flags.csv output.
    """
    # Columns determined by biorisk step
    biorisk = Flag.NOT_RUN
    virulence_factor = Flag.NOT_RUN
    # Columns determined by protein and nucleotide step
    regulated_virus = Flag.NOT_RUN
    regulated_bacteria = Flag.NOT_RUN
    regulated_eukaryote = Flag.NOT_RUN
    mixed_regulated_non_reg = Flag.NOT_RUN
    # Column determined by benign step
    benign = Outcome.NOT_RUN

    # Set biorisk and virtulence factor flags based on biorisk search
    if step_ran_successfully(results["biorisk"]):
        biorisk = Flag.FLAG if Outcome.FLAG in results["biorisk"] else Flag.PASS
        virulence_factor = Flag.FLAG if Outcome.WARN in results["biorisk"] else Flag.PASS

    # Error overrides anything else
    if Outcome.ERROR in results["biorisk"]:
        biorisk = virulence_factor = Flag.ERROR

    # Set regulated pathogen flags based on protein and nucleotide search
    if step_ran_successfully(results["protein"]):
        regulated_virus = Flag.FLAG if results["virus_flag"] else Flag.PASS
        regulated_bacteria = Flag.FLAG if results["bacteria_flag"] else Flag.PASS
        regulated_eukaryote = Flag.FLAG if results["eukaryote_flag"] else Flag.PASS
        mixed_regulated_non_reg = (
            Flag.FLAG
            if (Outcome.MIX in results["protein"] or Outcome.MIX in results["nucleotide"])
            else Flag.PASS
        )

    # Error overrides anything else
    if Outcome.ERROR in results["protein"] or Outcome.ERROR in results["nucleotide"]:
        regulated_virus = regulated_bacteria = regulated_eukaryote = mixed_regulated_non_reg = (
            Flag.ERROR
        )

    # Set benign flag based on benign search
    if Outcome.CLEARED in results["benign"]:
        benign = Flag.PASS
    if Outcome.NOT_CLEARED in results["benign"]:
        benign = Flag.FLAG
    if Outcome.ERROR in results["benign"]:
        benign = Flag.ERROR

    return {
        "filename": results["filepath"],
        "biorisk": biorisk,
        "virulence_factor": virulence_factor,
        "regulated_virus": regulated_virus,
        "regulated_bacteria": regulated_bacteria,
        "regulated_eukaryote": regulated_eukaryote,
        "mixed_regulated_and_non_reg": mixed_regulated_non_reg,
        "benign": benign,
    }


def get_recommendation(flag: str) -> str:
    """
    Map between flag column in screen_pipeline_status and older flags format.
    """
    flag_map = {
        Outcome.FLAG: Flag.FLAG,
        Outcome.PASS: Flag.PASS,
        Outcome.ERROR: Flag.ERROR,
    }
    return flag_map.get(flag, Flag.NOT_RUN)


def write_output_csvs(output_dir, status, flags, recommendations):
    """Write results to 3 output CSV files."""
    status_file = os.path.join(output_dir, "screen_pipeline_status.csv")
    flags_file = os.path.join(output_dir, "flags.csv")
    recommendations_file = os.path.join(output_dir, "flags_recommended.csv")

    status_df = pd.DataFrame(status)
    status_df = status_df.map(lambda x: ";".join(sorted(x)) if isinstance(x, set) else x)
    status_df.to_csv(status_file, index=False)
    print(f"Pipeline step status written to {status_file}")

    pd.DataFrame(flags).to_csv(flags_file, index=False)
    print(f"Screen results written to {flags_file}")

    summary = pd.DataFrame(recommendations, columns=["filename", "recommend_flag_or_pass"])
    summary.to_csv(recommendations_file, index=False, header=None)
    print(f"Flag recommendations written to {recommendations_file}")

    print("Flags: ", (summary["recommend_flag_or_pass"] == "F").sum(), "/", len(summary))
    print("Errors: ", (summary["recommend_flag_or_pass"] == "Err").sum())


def run(args: argparse.Namespace):
    """
    Wrapper so that args be parsed in main() or commec.py interface.

    Read all files that end with .screen in the input directory, then summarize their outcomes in
    three CSVs: flags.csv, flags_recommended.csv, and screen_pipeline_status.csv.
    """
    # Find .screen files
    search_dir = args.directory
    search_recursive = args.recursive
    output_dir = args.output or os.path.dirname(search_dir)

    search_pattern = "**/*.screen" if search_recursive else "*.screen"
    screen_paths = glob.glob(os.path.join(search_dir, search_pattern), recursive=search_recursive)
    if not screen_paths:
        raise FileNotFoundError(f"No .screen files were found in directory: {search_dir}")

    screen_status = [process_file(file_path) for file_path in sorted(screen_paths)]
    screen_flags = [get_flags_from_status(status) for status in screen_status]
    recommendations = [
        (status["filepath"], get_recommendation(status["flag"])) for status in screen_status
    ]

    write_output_csvs(output_dir, screen_status, screen_flags, recommendations)


def main():
    """Main function. Passes args to `run`."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
