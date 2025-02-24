#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Parse all .screen files in a directory and create two CSV files of flags raised during screening.

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
    "Parse all .screen files in a directory and create two CSVs file of flags raised"
)

def add_args(parser):
    """
    Add module arguments to an ArgumentParser object.
    """
    parser.add_argument(
        action="store",
        type=directory_arg,
        dest="screen_dir",
        help="Directory containing .screen files to summarize",
    )
    return parser


def add_flags(lines, flag_list, pattern, regex=False):
    """
    Parse a subset of input lines from a screen file and add to a list of flags
    """
    if regex:
        screen_lines = [s for s in lines if re.search(pattern, s)]
    else:
        screen_lines = [s for s in lines if pattern in s]

    if len(screen_lines) > 0:
        flag_list.append("F")
    else:
        flag_list.append("P")


def get_flag_list(screen_dir):
    """
    Search the input directory for .screen files, then search those files for lines indicating
    whether each kind of flag (biorisk, virulence factor, regulated pathogens, benign regions)
    was found.
    """
    # columns in the summary file, to be filled with the checks below
    filenames = []
    biorisk_flags = []
    vf_flags = []
    virus_flags = []
    bacteria_flags = []
    eukaryote_flags = []
    reg_nonreg_flags = []
    benign_flags = []

    for screen_path in glob.glob(os.path.join(screen_dir, "*.screen")):
        filenames.append(os.path.basename(screen_path))
        with open(screen_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
            add_flags(
                lines,
                biorisk_flags,
                r"Biorisks: Regulated gene in bases \d+ to \d+: FLAG",
                regex=True,
            )
            add_flags(
                lines,
                vf_flags,
                r"Virulence factor found in bases \d+ to \d+",
                regex=True,
            )
            add_flags(
                lines, virus_flags, "found in only regulated organisms: FLAG (virus)"
            )
            add_flags(
                lines,
                bacteria_flags,
                "found in only regulated organisms: FLAG (bacteria)",
            )
            add_flags(
                lines,
                eukaryote_flags,
                "found in only regulated organisms: FLAG (eukaryote)",
            )
            add_flags(
                lines,
                reg_nonreg_flags,
                r"found in both regulated and non-regulated organisms",
            )

            # All homology-related flags should be replaced with "Err" if search failed
            homol_fail = [s for s in lines if "ERROR: Homology search has failed" in s]
            if len(homol_fail) > 0:
                virus_flags[-1:] = ["Err"]
                bacteria_flags[-1:] = ["Err"]
                eukaryote_flags[-1:] = ["Err"]
                reg_nonreg_flags[-1:] = ["Err"]

            # All homology-related flags should be replaced with "-" if run in fast mode
            fast_mode = [s for s in lines if "FAST MODE" in s]
            if len(fast_mode) > 0:
                virus_flags[-1:] = ["-"]
                bacteria_flags[-1:] = ["-"]
                eukaryote_flags[-1:] = ["-"]
                reg_nonreg_flags[-1:] = ["-"]

            # benign screen
            # 1 means a regulated region failed to clear, 0 means benign coverage and clear
            allpass = [s for s in lines if "all regulated regions cleared: PASS" in s]
            anyfail = [s for s in lines if "failed to clear: FLAG" in s]
            clear = [s for s in lines if "no regulated regions to clear" in s]
            # if any region failed to clear, keep flag
            if len(allpass) > 0:
                benign_flags.append("P")
            # if none failed and passes are observed, drop flag
            elif len(anyfail) > 0:
                benign_flags.append("F")
            elif len(clear) > 0:
                benign_flags.append("-")
            else:
                benign_flags.append("Err")

    return list(
        zip(
            filenames,
            biorisk_flags,
            vf_flags,
            virus_flags,
            bacteria_flags,
            eukaryote_flags,
            reg_nonreg_flags,
            benign_flags,
        )
    )


def write_flag_files(screen_dir):
    """
    Writes flags.csv (all flags) and flags_recommended.csv (single recommended flag), with one row
    for each .screen file found in the input directory.
    """
    output_dir = os.path.dirname(screen_dir)
    detail_file = os.path.join(output_dir, "flags.csv")
    summary_file = os.path.join(output_dir, "flags_recommended.csv")

    flags = get_flag_list(screen_dir)
    summary = []
    for name, risk, vf, reg_vir, reg_bac, eukaryote_flags, _, benign in flags:
        # if a biorisk is flagged, flag the whole thing
        if risk == "Err" or reg_bac == "Err" or benign == "Err":
            summary.append((name, "Err"))
        else:
            if risk == "F":
                summary.append((name, "F"))
            # Flag
            elif reg_vir == "F" and benign == "F":
                summary.append((name, "F"))
            elif reg_vir == "F" and benign == "P":
                summary.append((name, "P"))
            # if it's a regulated bacterial pathogen but a known benign gene, clear it
            elif (reg_bac == "F" and benign == "P" and vf == "P") == 1:
                summary.append((name, "P"))
            # though VFs are housekeeping genes, these seem to be poorly supported Victors genes
            elif (reg_bac == "F" and benign == "P" and vf == "F") == 1:
                summary.append((name, "P"))
            elif (eukaryote_flags == "F" and benign == "P") == 1:
                summary.append((name, "P"))
            # if it's a regulated bacterial hit, flag it
            elif reg_bac == "F":
                summary.append((name, "F"))
            elif eukaryote_flags == "F":
                summary.append((name, "F"))
            else:
                summary.append((name, "P"))
    summary = pd.DataFrame(summary)
    summary.columns = ("filename", "recommend_flag_or_pass")
    summary.to_csv(summary_file, index=False, header=None)

    flags = pd.DataFrame(flags)
    flags.columns = (
        "filename",
        "biorisk",
        "virulence_factor",
        "regulated_virus",
        "regulated_bacteria",
        "regulated_eukaryote",
        "mixed_regulated_and_non_reg",
        "benign",
    )
    flags.to_csv(detail_file, index=False)

    print(
        "Flags: ", (summary["recommend_flag_or_pass"] == "F").sum(), "/", len(summary)
    )
    print("Errors: ", (summary["recommend_flag_or_pass"] == "Err").sum())


def run(parsed_args):
    """
    Wrapper so that args be parsed in main() or commec.py interface.
    """
    write_flag_files(parsed_args.screen_dir)


def main():
    """
    Main function. Passes args to `run`.
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
