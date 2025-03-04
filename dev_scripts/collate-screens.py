#!/usr/bin/env python3
"""
Collate screen files located in any subdirectories of an input directory, then match description
fields between FASTAs in adjacent _input directories and FASTAs in another directory, renaming
collated screen files with the names of matched FASTA files. Created as a workaround for breaking
changes introduced in the commec output format in v0.3.

Required inputs:
    -i, --input-dir     input directory to recursively search for .screen files
    -o, --output-dir    output directory where screen files should be collated
    -f, --fasta-dir     directory to search recursively for FASTAs input to commec screen

Example:
$ python collate-screens.py -i . -o ./test-collate-pls -f ../functional-json-test/
"""
import argparse
import csv
import os
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

def clean_header(header: str) -> str:
    """
    Clean a FASTA header by replacing whitespace and special characters with underscores.
    
    Args:
        header: Original FASTA header
        
    Returns:
        Cleaned header string
    """
    return "".join(
        "_" if c.isspace() or c == "\xc2\xa0" or c == "#" else c
        for c in header
    )

def parse_fasta_header(fasta_path: Path) -> Tuple[str, str]:
    """
    Extract filename and cleaned FASTA header from a FASTA file.
    
    Args:
        fasta_path: Path to the FASTA file
        
    Returns:
        Tuple of (filename, cleaned_header)
    """
    filename = fasta_path.name
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()
                cleaned_header = clean_header(header)
                return filename, cleaned_header
    raise ValueError(f"No valid FASTA header found in {fasta_path}")

def build_fasta_mapping(fasta_dir: Path) -> Dict[str, str]:
    """
    Build mapping of FASTA headers to filenames from all FASTA files in directory.
    
    Args:
        fasta_dir: Directory containing FASTA files
        
    Returns:
        Dictionary mapping FASTA headers to original filenames
    """
    mapping = {}
    for root, _, files in os.walk(fasta_dir):
        for file in files:
            if (file.endswith('.fasta') and 
                not file.endswith('.noncoding.fasta') and 
                not file.endswith('.cleaned.fasta')):
                fasta_path = Path(root) / file
                try:
                    filename, cleaned_header = parse_fasta_header(fasta_path)
                    mapping[cleaned_header] = filename
                except (ValueError, IOError) as e:
                    print(f"Warning: Could not process {fasta_path}: {e}")
    return mapping

def find_screen_files(input_dir: Path) -> List[Path]:
    """
    Find all .screen files in the input directory.
    """
    screen_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.screen'):
                screen_files.append(Path(root) / file)
    return screen_files

def get_matching_fasta_header(screen_path: Path) -> str:
    """
    Get the FASTA header from the input screen file.
    """
    # For a file named "something.screen", look for "input_something/something.cleaned.fasta"
    screen_name = screen_path.stem  # removes .screen extension
    fasta_path = screen_path.parent / f"input_{screen_name}" / f"{screen_name}.cleaned.fasta"
    if not fasta_path.exists():
        raise FileNotFoundError(f"Expected FASTA not found at {fasta_path}")
    
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()
                return clean_header(header)
    raise ValueError(f"No valid FASTA header found in {fasta_path}")

def main():
    parser = argparse.ArgumentParser(description='Collate and rename screen files based on FASTA headers')
    parser.add_argument('-i', '--input-dir', required=True, help='Input directory containing screen files')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory for renamed screen files')
    parser.add_argument('-f', '--fasta-dir', required=True, help='Directory containing original FASTA files')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    fasta_dir = Path(args.fasta_dir).resolve()
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Build mapping of FASTA headers to original filenames
    print(f"Building FASTA header mapping based on files found in {fasta_dir}...")
    header_to_filename = build_fasta_mapping(fasta_dir)
    num_fastas_mapped = len(header_to_filename)
    if num_fastas_mapped == 0:
        print("Could not find any FASTAs to map! Note that .cleaned and .noncoding are filtered out.\nExiting...")
        exit(0)
    print(f"Found {num_fastas_mapped} FASTAS for mapping...")

    # Find all screen files
    print(f"Finding screen files in {input_dir}...")
    screen_files = find_screen_files(input_dir)
    print(f"Processing {len(screen_files)} screen files...")

    mappings = []
    for screen_path in screen_files:
        try:
            # Get the FASTA header for this screen file
            fasta_header = get_matching_fasta_header(screen_path)
            
            # Look up the matching filename
            if fasta_header not in header_to_filename:
                print(f"Warning: No matching FASTA file found for {screen_path}")
                continue
                
            matching_filename = header_to_filename[fasta_header]
            new_filename = f"{Path(matching_filename).stem}.screen"
            output_path = output_dir / new_filename
            
            # Copy the screen file with the new name
            shutil.copy2(screen_path, output_path)
            
            # Record the mapping
            mappings.append({
                'screen': str(screen_path),
                'matched_fasta': matching_filename,
                'renamed_screen': new_filename
            })
            
            print(f"Processed: {screen_path} -> {output_path}")
            
        except (FileNotFoundError, ValueError, IOError) as e:
            print(f"Error processing {screen_path}: {e}")
    
    # Write mapping CSV
    csv_path = output_dir / 'screen_mappings.csv'
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['screen', 'matched_fasta', 'renamed_screen'])
        writer.writeheader()
        writer.writerows(mappings)
    
    print(f"\nProcessed {len(mappings)} screen files")
    print(f"Mapping saved to {csv_path}")

if __name__ == '__main__':
    main()
