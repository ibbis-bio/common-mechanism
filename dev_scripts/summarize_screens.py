#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Summarize the screening for all .screen files in a directory. This is intended to be useful for
debugging the pipeline, rather than for interpreting the outputs.

Produces a CSV (name set by -o, defaults to 'output.csv') which contains the outcome for each step
of the pipeline:

* flag      the sequence was flagged in this step
* pass      the sequence passed in this step
* skip      this step was intentionally not run
* error     an error occurred during this step
* -         this step was not run due to an error, interrupt, or other unexpected outcome
* mix       (protein only) the best match is to a mix of regulated- and non-regulated organisms
* warn      (biorisk only) found a significant hit to a virulence not from a regulated pathogen

Each line in the CSV corresponds to a .screen file. The full paths to the files are also provided.
"""
import os
import csv
import argparse
import re

def process_step(step_content, step_number):
    """
    Process the .screen file output to determine the outcome of the step.
    """
    if step_number == 1:  # Biorisk
        if 'FLAG' in step_content:
            return 'flag'
        elif 'Virulence factor found' in step_content:
            return 'warn'
        elif 'Biorisks: no hits detected, PASS' in step_content or 'Biorisks: no significant hits detected, PASS' in step_content:
            return 'pass'
        elif 'ERROR:' in step_content:
            return 'error'
    elif step_number == 2:  # Protein
        if 'Best match to sequence(s)' in step_content and 'FLAG' in step_content:
            return 'flag'
        elif 'found in both regulated and non-regulated organisms' in step_content:
            return 'mix'
        elif 'no top hit exclusive to a regulated pathogen: PASS' in step_content:
            return 'pass'
        elif 'ERROR:' in step_content:
            return 'error'
    elif step_number == 3:  # Nucleotide
        if 'no noncoding regions >= 50 bases found, skipping nt scan' in step_content:
            return 'skip'
        elif 'Best match to sequence(s)' in step_content and 'FLAG' in step_content:
            return 'flag'
        elif 'no top hit exclusive to a regulated pathogen: PASS' in step_content:
            return 'pass'
        elif 'ERROR:' in step_content:
            return 'error'
    elif step_number == 4:  # Benign
        if 'no regulated regions to clear' in step_content:
            return 'skip'
        elif (
            'Regulated region at bases' in step_content and
            'failed to clear: FLAG' in step_content
        ):
            return 'flag'
        elif 'all regulated regions cleared: PASS' in step_content:
            return 'pass'
        elif 'ERROR:' in step_content:
            return 'error'
    return '-'

def process_file(file_path):
    """
    Read input screen file, split into steps, and prepare dict of results for CSV output.
    """
    with open(file_path, 'r', encoding='utf-8') as file:
        content = file.read()

    filename = os.path.basename(file_path)
    filename_without_extension = os.path.splitext(filename)[0]

    steps = re.split(r'>> STEP \d:', content)
    steps = [step.strip() for step in steps if step.strip()]

    results = {
        'filename': filename_without_extension,
        'location': file_path,
        'biorisk': process_step(steps[0] if len(steps) > 0 else '-', 1),
        'protein': process_step(steps[1] if len(steps) > 1 else '-', 2),
        'nucleotide': process_step(steps[2] if len(steps) > 2 else '-', 3),
        'benign': process_step(steps[3] if len(steps) > 3 else '-', 4)
    }

    return results

def main(directory, output_file):
    results = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.screen'):
                file_path = os.path.join(root, file)
                results.append(process_file(file_path))

    # Write results to CSV
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['filename', 'location', 'biorisk', 'protein', 'nucleotide', 'benign']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for result in results:
            writer.writerow(result)

    print(f"Results written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process .screen files and output results to CSV.")
    parser.add_argument("directory", help="Directory to search for .screen files")
    parser.add_argument(
        "-o", "--output", default="output.csv", help="Output CSV file name (default: output.csv)"
    )

    args = parser.parse_args()

    main(args.directory, args.output)
