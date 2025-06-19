#!/usr/bin/env bash
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
#####################################################################
# run_diamond.sh runs DIAMOND against a specified NCBI nr database.
# DIAMOND citation: 
#   Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at 
#   tree-of-life scale using DIAMOND", 
#   Nature Methods 18, 366–368 (2021). 
#   doi:10.1038/s41592-021-01101-x
#####################################################################

set -eu

# Default values
JOBS=""
THREADS=1
DB_PATH=""
INPUT=""
OUTPUT="out"

usage() {
    echo "Usage: run_diamond.sh -d MY_DB -i INPUT_FILE [-o OUTPUT_FILE] [-j JOBS] [-t THREADS]"
    echo "  MY_DB           location of NCBI nr database (required)"
    echo "  INPUT_FILE      input file to align to each database (required)" 
    echo "  OUTPUT_FILE     output prefix for alignments (default: out)" 
    echo "  JOBS            number of diamond runs to do in parallel (default: # CPUs / THREADS)"
    echo "  THREADS         number of threads for each diamond run (default: 1)"
    exit 1
}

# Parse command line arguments
while getopts "j:t:d:i:o:" opt; do
    case $opt in
        j) JOBS=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        d) DB_PATH=$OPTARG ;;
        i) INPUT=$OPTARG ;;
        o) OUTPUT=$OPTARG ;;
        \?) usage ;;
    esac
done

# Check for required arguments
if [[ -z "$DB_PATH" || -z "$INPUT" ]]; then
    usage
fi

# Validate input
if [[ ! -d "$DB_PATH" ]]; then
    echo "ERROR: nr diamond database folder $DB_PATH does not exist" >&2
    exit 1
fi

shopt -s failglob
if ! files=("${DB_PATH}"/nr*.dmnd); then
    echo "ERROR: No nr diamond database files (nr*.dmnd) found in $DB_PATH" >&2
    exit 1
fi
shopt -u failglob

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file $INPUT does not exist" >&2
    exit 1
fi

# Set JOBS if not specified by user
if [[ -z "$JOBS" ]]; then
    CPU_COUNT=$(parallel --number-of-cpus)
    JOBS=$((CPU_COUNT / THREADS))
fi

echo "Running diamond protein search..."
echo "Using $JOBS job(s) in parallel with $THREADS thread(s) each"

# Run diamond
ls ${DB_PATH}/nr*.dmnd | parallel --will-cite --use-cpus-instead-of-cores --jobs ${JOBS} \
    diamond blastx --quiet \
    -d {} \
    --threads ${THREADS} \
    -q ${INPUT} \
    -o ${OUTPUT}.{/.}.tsv \
    --outfmt 6 qseqid stitle sseqid staxids evalue bitscore pident qlen qstart qend slen sstart send \
    --frameshift 15 --range-culling

# Combine results and clean up
cat ${OUTPUT}.*.tsv > ${OUTPUT}.dmnd
rm ${OUTPUT}.*.tsv

echo "Diamond protein search completed. Results are in ${OUTPUT}.dmnd"
