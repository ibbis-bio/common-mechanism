#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Parse all .screen files in a directory and create three CSV files reflecting the results of screening.

Each line in each CSV corresponds to a .screen file.

The flags.csv file will have the following columns and values:

- filename:                 .screen file basename
- biorisk:                  "F" if flagged, "P" if no flags
- virulence_factor:         "F" if flagged, "P" if no flags
- regulated_virus:          "F" if flagged, "P" if no flags, "Err" if error logged
- regulated_bacteria:       "F" if flagged, "P" if no flags, "Err" if error logged
- regulated_eukaryote:      "F" if flagged, "P" if no flags, "Err" if error logged
- mixed_regulated_non_reg:  "F" if flagged, "P" if no flags, "Err" if error logged
- benign:                   "F" if not cleared, "P" if all cleared, "-" if not run

The flags_recommended CSV just has two columns, "filename" and "recommend_flag_or_pass".

The status.csv file will have the following columns and values:

* flag          the sequence was flagged in this step
* pass          the sequence passed in this step
* skip          this step was intentionally not run
* error         an error occurred during this step
* -             this step was not run due to an error, interrupt, or other unexpected outcome
* mix           (protein only) the best match is to a mix of regulated- and non-regulated organisms
* warn          (biorisk only) found a significant hit to a virulence not from a regulated pathogen
* cleared       (benign only) earlier flags were cleared
* not cleared   (benign only) earlier flags were not cleared

Additionally, it includes three columns indicating whether the sequence was flagged as a regulated
virus, bacteria, or eukaryote, and the full paths to the files are also provided.

Command-line usage:
    flag.py /path/to/directory/with/.screen 
"""
import argparse
import glob
import os
import re
import pandas as pd
from commec.utils.file_utils import directory_arg

DESCRIPTION = (
    "Parse all .screen files in a directory and create CSVs file of flags raised"
)

def add_args(parser: argparse.ArgumentParser):
    """
    Add module arguments to an ArgumentParser object.
    """
    parser.add_argument(
        "directory",
        action="store",
        type=directory_arg,
        help="Directory containing .screen files to summarize"
    )
    parser.add_argument(
        "-r",
        "--recursive",
        dest="recursive",
        action="store_true",
        help="Search directory recursively for screen files"
    )
    return parser


def process_step(step_content: str, step_number: int):
    """
    Process the .screen file output to determine the outcome of the step.
    """
    step_processors = {
        1: get_biorisk_outcome,
        2: get_protein_outcome,
        3: get_nucleotide_outcome,
        4: get_benign_outcome,
    }
    return step_processors.get(step_number, lambda _: "-")(step_content)


def get_biorisk_outcome(step_content: str) -> set[str]:
    """Process biorisk scan step from .screen file."""
    outcomes = set()
    if "FLAG" in step_content:
        outcomes.add("flag")
    if "Virulence factor found" in step_content:
        outcomes.add("warn")
    if (
        "Biorisks: no hits detected, PASS" in step_content
        or "Biorisks: no significant hits detected, PASS" in step_content
    ):
        outcomes.add("pass")
    if "ERROR" in step_content:
        outcomes.add("error")
    return outcomes if outcomes else {"-"}

def get_protein_outcome(step_content: str) -> set[str]:
    """Process protein scan step from .screen file."""
    outcomes = set()
    if "FAST MODE: Skipping" in step_content or "SKIPPING STEP 2" in step_content:
        outcomes.add("skip")
    if "Best match to sequence(s)" in step_content and "FLAG" in step_content:
        outcomes.add("flag")
    if "found in both regulated and non-regulated organisms" in step_content:
        outcomes.add("mix")
    if "no top hit exclusive to a regulated pathogen: PASS" in step_content:
        outcomes.add("pass")
    if "ERROR" in step_content:
        outcomes.add("error")
    return outcomes if outcomes else {"-"}

def get_nucleotide_outcome(step_content: str) -> set[str]:
    """Process nucleotide scan step from .screen file."""
    outcomes = set()
    if (
        "no noncoding regions >= 50 bases found, skipping nt scan" in step_content
        or "FAST MODE: Skipping" in step_content
        or "SKIPPING STEP 3" in step_content
    ):
        outcomes.add("skip")
    if "Best match to sequence(s)" in step_content and "FLAG" in step_content:
        outcomes.add("flag")
    if "no top hit exclusive to a regulated pathogen: PASS" in step_content:
        outcomes.add("pass")
    if "ERROR" in step_content:
        outcomes.add("error")
    return outcomes if outcomes else {"-"}

def get_benign_outcome(step_content) -> set[str]:
    """Process benign scan step from .screen file."""
    outcomes = set()
    if "no regulated regions to clear" in step_content:
        outcomes.add("skip")
    if (
        "Regulated region at bases" in step_content
        and "failed to clear: FLAG" in step_content
    ):
        outcomes.add("not cleared")
    if "all regulated regions cleared: PASS" in step_content:
        outcomes.add("cleared")
    if "ERROR" in step_content:
        outcomes.add("error")
    return outcomes if outcomes else {"-"}


def check_regulated_flags(content):
    """
    Check for regulated virus, bacteria, and eukaryote flags in the content.
    """
    return {
        "virus_flag": "true" if "FLAG (virus)" in content else "false",
        "bacteria_flag": "true" if "FLAG (bacteria)" in content else "false",
        "eukaryote_flag": "true" if "FLAG (eukaryote)" in content else "false",
    }


def process_file(file_path) -> list[str]:
    """
    Read input screen file, split into steps, and prepare dict of results for CSV output.
    """
    with open(file_path, "r", encoding="utf-8") as file:
        content = file.read()

    first_step = re.search(r'>> STEP 1:', content)
    if first_step:
        content = content[first_step.start():]
        steps = re.split(r'(?=(?:>> STEP \d|SKIPPING STEP \d))', content)   
    steps = [step.strip() for step in steps if step.strip()]

    filename = os.path.basename(file_path)
    filename_without_extension = os.path.splitext(filename)[0]

    results = {
        "filename": filename_without_extension,
        "location": file_path,
        "biorisk": process_step(steps[0] if len(steps) > 0 else "-", 1),
        "protein": process_step(steps[1] if len(steps) > 1 else "-", 2),
        "nucleotide": process_step(steps[2] if len(steps) > 2 else "-", 3),
        "benign": process_step(steps[3] if len(steps) > 3 else "-", 4),
    }
    # Add regulated flags
    results.update(check_regulated_flags(content))

    return results


def get_flags_from_status(results: dict[str, str | set[str] | bool]) -> list[str]:
    """
    Map step results to flags format, which has the following columns:
        - filename
        - biorisk:                      set by biorisk step
        - virulence_factor:             set by biorisk step
        - regulated_virus:              set by protein or nucleotide step
        - regulated_bacteria:           set by protein or nucleotide step
        - regulated_eukaryote:          set by protein or nucleotide step
        - mixed_regulated_and_non_reg:  set by protein or nucleotide step
        - benign:                       set by benign step
    
    And can take the following values:
        - "F" if flagged (or, for benign, not cleared)
        - "P" if no flags (or, for benign, all cleared)
        - "Err" if error logged
        - "-" if not run (skipped OR not reached due to an error)
    """
    filename = results["location"]

    # Columns determined by biorisk step
    biorisk = "-"
    virulence_factor = "-"

    if "error" in results["biorisk"]:
        biorisk = virulence_factor = "Err"
    elif "flag" in results["biorisk"]:
        biorisk = "F"
        virulence_factor = "P" if "warn" not in results["biorisk"] else "F"
    elif "warn" in results["biorisk"]:
        biorisk = "P"
        virulence_factor = "F"
    elif "pass" in results["biorisk"]:
        biorisk = "P" 
        virulence_factor = "P"

    # Columns determined by protein and nucleotide step
    regulated_virus = "-"
    regulated_bacteria = "-"
    regulated_eukaryote = "-"
    mixed_regulated_non_reg = "-"

    if "error" in results["protein"]:
        regulated_virus = "Err"
        regulated_bacteria = "Err"
        regulated_eukaryote = "Err"
        mixed_regulated_non_reg = "Err"
    elif "skip" in results["protein"] or "-" in results["protein"]:
        pass # leave values as defaults
    elif "pass" in results["protein"] and "pass" in results["nucleotide"]:
        regulated_virus = "P"
        regulated_bacteria = "P"
        regulated_eukaryote = "P"
        mixed_regulated_non_reg = "F" if ("mix" in results["protein"] or "mix" in results["nucleotide"]) else "P"
    elif "flag" in results["protein"] or "flag" in results["nucleotide"]:
        regulated_virus = "F" if results["virus_flag"] == "true" else "P"
        regulated_bacteria = "F" if results["bacteria_flag"] == "true" else "P"
        regulated_eukaryote = "F" if results["eukaryote_flag"] == "true" else "P"
        mixed_regulated_non_reg = "F" if ("mix" in results["protein"] or "mix" in results["nucleotide"]) else "P"

    # Column determined by benign step
    benign = "-"
    if "error" in results["benign"]:
        benign = "Err"
    if "cleared" in results["benign"]:
        benign = "P"
    if "not cleared" in results["benign"]:
        benign = "F"

    return [
        filename,
        biorisk,
        virulence_factor,
        regulated_virus,
        regulated_bacteria,
        regulated_eukaryote,
        mixed_regulated_non_reg,
        benign,
    ]


def get_recommendation(flags: list[str]) -> str:
    """
    Get a single overall recommendation (P, F, or Err) based on the flags.

    We don't take into account virulence factor flags here; those are just a warning, since
    we've found that some virulence factors are housekeeping genes (and are poorly-supported
    Victors entries).
    """
    _, risk, _, reg_vir, reg_bac, reg_euk, _, benign = flags

    # Errors mean we can't determine overall pass or flag
    if risk == "Err" or reg_bac == "Err" or benign == "Err":
        return "Err"
    # if a biorisk is flagged, flag the whole thing
    if risk == "F":
        return "F"
    # Flag regulated pathogens unless cleared in the benign screen
    elif (reg_vir == "F" or reg_bac == "F" or reg_euk == "F") and benign == "F":
        return  "F"
    else:
        return  "P"


def write_output_csvs(output_dir, status, flags, recommendations):
    status_file = os.path.join(output_dir, "screen_pipeline_status.csv")
    flags_file = os.path.join(output_dir, "flags.csv")
    recommendations_file = os.path.join(output_dir, "flags_recommended.csv")

    status_df = pd.DataFrame(status)
    status_df.columns = [
        "name",
        "filepath",
        "biorisk",
        "protein",
        "nucleotide",
        "benign",
        "virus_flag",
        "bacteria_flag",
        "eukaryote_flag",
    ]
    status_df = status_df.apply(lambda x: x.apply(lambda y: ';'.join(y) if isinstance(y, set) else y))
    columns_to_print = [
        "name",
        "filepath",
        "biorisk",
        "protein",
        "nucleotide",
        "benign"
    ]
    status_df[columns_to_print].to_csv(status_file, index=False)
    print(f"Pipeline step status written to {status_file}")

    flags_df = pd.DataFrame(flags)
    flags_df.columns = [
        "filename",
        "biorisk",
        "virulence_factor",
        "regulated_virus",
        "regulated_bacteria",
        "regulated_eukaryote", 
        "mixed_regulated_and_non_reg",
        "benign"
    ]
    flags_df.to_csv(flags_file, index=False)
    print(f"Screen results written to {flags_file}")

    summary = pd.DataFrame(recommendations)
    summary.columns = ("filename", "recommend_flag_or_pass")
    summary.to_csv(recommendations_file, index=False, header=None)
    print(f"Flag recommendations written to {recommendations_file}")

    print(
        "Flags: ", (summary["recommend_flag_or_pass"] == "F").sum(), "/", len(summary)
    )
    print("Errors: ", (summary["recommend_flag_or_pass"] == "Err").sum())


def run(args : argparse.Namespace):
    """
    Wrapper so that args be parsed in main() or commec.py interface.
    
    Read all files that end with .screen in the input directory, then summarize their outcomes in
    three CSVs: flags.csv, flags_recommended.csv, and screen_pipeline_status.csv.
    """
    # Check that screen files can be found
    search_dir = args.directory
    search_recursive = args.recursive
    search_pattern = "**/*.screen" if search_recursive else "*.screen"
    screen_paths = glob.glob(os.path.join(search_dir, search_pattern), recursive=search_recursive)
    if not screen_paths:
        raise FileNotFoundError(f"No .screen files were found in directory: {search_dir}")

    screen_status = []
    for file_path in screen_paths:
        screen_status.append(process_file(file_path))

    screen_flags = [get_flags_from_status(status) for status in screen_status]
    recommendations = [(flags[0], get_recommendation(flags)) for flags in screen_flags]

    write_output_csvs(os.path.dirname(search_dir), screen_status, screen_flags, recommendations)


def main():
    """
    Main function. Passes args to `run`.
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
