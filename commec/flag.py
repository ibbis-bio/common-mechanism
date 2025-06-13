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
import pandas as pd
from commec.utils.file_utils import directory_arg
from commec.config.result import ScreenStatus, ScreenResult, ScreenStep
from commec.config.json_io import get_screen_data_from_json, IoVersionError

DESCRIPTION = "Parse all .screen, or .json files in a directory and create CSVs of flags raised"

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

def read_flags_from_json(file_path) -> list[dict[str, str | set[str] | bool]]:
    """
    Read JSON screen output and prepare screen_pipline_status CSV.
    """
    results = []
    try:
        screen_data : ScreenResult = get_screen_data_from_json(file_path)
    except (KeyError, AttributeError):
        return []
    except IoVersionError as e:
        print(f"The following json was not a compatible version ({file_path}): {e}")
        return []

    for name, query in screen_data.queries.items():
        # Defaults are false, and are overwritten if a single hit occurs.
        virus_flag = False
        bacteria_flag = False
        eukaryote_flag = False
        benign_protein = False
        benign_rna = False
        benign_synbio = False

        qr = query.recommendation

        if (qr.protein_taxonomy_status
            not in [ScreenStatus.SKIP, ScreenStatus.ERROR, ScreenStatus.NULL]):
            for hit in query.hits.values():
                if hit.recommendation.from_step == ScreenStep.TAXONOMY_AA:
                    for info in hit.annotations["regulated_taxonomy"]:
                        bacteria_flag  |= (int(info["regulated_bacteria"]) > 0)
                        virus_flag     |= (int(info["regulated_viruses"]) > 0)
                        eukaryote_flag |= (int(info["regulated_eukaryotes"]) > 0)

        # Which forms of benign hits are present?
        if (qr.benign_status
            not in [ScreenStatus.SKIP, ScreenStatus.ERROR, ScreenStatus.NULL]):
            for hit in query.hits.values():
                match hit.recommendation.from_step:
                    case ScreenStep.BENIGN_PROTEIN:
                        benign_protein = True
                    case ScreenStep.BENIGN_RNA:
                        benign_rna = True
                    case ScreenStep.BENIGN_DNA:
                        benign_synbio = True
                    case _:
                        continue

        mixed_aa_taxonomy = False
        mixed_nt_taxonomy = False
        for hit in query.hits.values():
            match hit.recommendation.from_step:
                case ScreenStep.TAXONOMY_AA:
                    if hit.recommendation.status in [ScreenStatus.FLAG, ScreenStatus.WARN]:
                        reg_dicts = hit.annotations["regulated_taxonomy"]
                        for r in reg_dicts:
                            if len(r["non_regulated_taxids"]) > 0:
                                mixed_aa_taxonomy = True

                case ScreenStep.TAXONOMY_NT:
                    if hit.recommendation.status in [ScreenStatus.FLAG, ScreenStatus.WARN]:
                        reg_dicts = hit.annotations["regulated_taxonomy"]
                        for r in reg_dicts:
                            if len(r["non_regulated_taxids"]) > 0:
                                mixed_nt_taxonomy = True
                case _:
                    continue

        protein_status = "Mixed" if mixed_aa_taxonomy else qr.protein_taxonomy_status
        nucleotide_status = "Mixed" if mixed_nt_taxonomy else qr.nucleotide_taxonomy_status

        results.append({
        "name": name,
        "filepath": file_path,
        "flag": query.recommendation.screen_status,
        "biorisk": query.recommendation.biorisk_status,
        "protein": protein_status,
        "nucleotide": nucleotide_status,
        "benign": query.recommendation.benign_status,
        "virus_flag": virus_flag,
        "bacteria_flag": bacteria_flag,
        "eukaryote_flag": eukaryote_flag,
        "benign_protein": benign_protein,
        "benign_rna": benign_rna,
        "benign_dna": benign_synbio
        })

    return results

def write_output_csv(output_dir : str | os.PathLike, status : dict):
    """
    Write flag data results to an output CSV file.
    -----
    `output_dir` [str | os.PathLike] : Desired output directory to write the output csv into.
    `status`     [dict-like object]  : converted to pandas dataframe, sorted, and output to csv.
    """
    status_file = os.path.join(output_dir, "screen_pipeline_status.csv")
    status_df = pd.DataFrame(status)
    status_df = status_df.map(lambda x: ";".join(sorted(x)) if isinstance(x, set) else x)
    status_df.to_csv(status_file, index=False)
    print(f"Pipeline step status written to {status_file}")


def run(args: argparse.Namespace):
    """
    Wrapper so that args be parsed in main() or commec.py interface.

    Read all files that end with .screen in the input directory, then summarize their outcomes in
    three CSVs: flags.csv, flags_recommended.csv, and screen_pipeline_status.csv.
    """
    search_dir = args.directory
    search_recursive = args.recursive
    output_dir = args.output or os.path.dirname(search_dir)

    search_pattern = "**/*.json" if search_recursive else "*.json"
    screen_paths = glob.glob(os.path.join(search_dir, search_pattern), recursive=search_recursive)

    if not screen_paths:
        raise FileNotFoundError(f"No .json files were found in directory: {search_dir}")

    screen_status = []
    for file_path in sorted(screen_paths):
        result = read_flags_from_json(file_path)
        screen_status.extend(result)

    write_output_csv(output_dir, screen_status)


def main():
    """Main function. Passes args to `run`."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    add_args(parser)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
