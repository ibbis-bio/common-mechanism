#!/usr/bin/env bash
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
#####################################################################
#run_blastx.sh runs blastx against a specified database.
#####################################################################
#Usage: run_blastx.sh -d DB -q QUERY -o OUTPUT [-t THREADS -f FURTHEROPT]

#set -eux #debug mode
set -eu
DB=""
QUERY=""
OUTPUT="out" 
THREADS=1
FURTHEROPT=""

#Get options from user
while getopts "t:d:q:o:f:" OPTION
    do 
        case $OPTION in
            t)
                THREADS=$OPTARG
                ;;
            d) 
                DB=$OPTARG
                ;;
            q)
                QUERY=$OPTARG
                ;;
            o) 
                OUTPUT=$OPTARG
                ;;
            f)
                FURTHEROPT=$OPTARG
                ;;
            \?)
                echo "Usage: run_blastx.sh -d DB -q QUERY -s OUTPUT [-t THREADS -f FURTHEROPT]"
                echo "  DB              full path to database (required)"
                echo "  QUERY           query file to align to each database (required)"
                echo "  OUTPUT          output prefix for alignments (default: out)"
                echo "  THREADS         number of threads for each database run (default: 1)"
                echo "  FURTHEROPT      any further options to specify"
                exit
                ;;  
        esac 
    done

#Check for values
if [ "$DB" == "" ] && [ "$INPUT" == "" ]
then
    echo "Usage: run_blastx.sh -d DB -q QUERY -s OUTPUT [-t THREADS -f FURTHEROPT]"
        echo "  DB              full path to database (required)"
        echo "  QUERY           query file to align to each database (required)"
        echo "  OUTPUT          output prefix for alignments (default: out)"
        echo "  THREADS         number of threads for each database run (default: 1)"
        echo "  FURTHEROPT      any further options to specify"
    exit
fi

#Check for database
echo -e "\t...checking for valid options..."
if [ ! -f "$DB".pal ] # this is v. blast specific
then
    echo " ERROR: blastx database $DB does not exist"
    exit
fi
  
#Check for input file 
if [ ! -f  "$QUERY" ]
then
    echo " ERROR: input file $QUERY does not exist"
    exit
fi      

echo -e "\t...running protein search..."
blastx -db ${DB} -query ${QUERY} -out ${OUTPUT}.blastx -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-10 -word_size 6 -threshold 21 -max_target_seqs 5000 -culling_limit 50 -window_size 40 -matrix BLOSUM62 -gapopen 11 -gapextend 1 -seg yes -num_threads ${THREADS} ${FURTHEROPT}
