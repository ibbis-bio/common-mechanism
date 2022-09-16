#! /usr/bin/env bash

##############################################################################
#run_pipeline.sh runs the Common Mechanism against a specified QUERY file.
##############################################################################
# Usage: src/run_pipeline_diamond.sh test_folder/[$name].fasta

# parameters for the run
set -eu
PROCESSES=6         #number of database chunks
THREADS=1           #threads per process
QUERY=""
OUTPUT=""
CLEANUP=0

#Get options from user
while getopts "p:t:q:o:c:" OPTION
    do
        case $OPTION in
            p)
                PROCESSES=$OPTARG
                ;;
            t)
                THREADS=$OPTARG
                ;;
            d)
                DB_PATH=$OPTARG
                ;;
            q)
                QUERY=$OPTARG
                ;;
            o)
                OUTPUT=$OPTARG
                ;;
            c)
                CLEANUP=$OPTARG
                ;;
            \?)
                echo "Usage: src/run_pipeline.sh -q QUERY -s OUTPUT [-p PROCESSES -t THREADS]"
                echo "  QUERY           query file to align to each database (required)"
                echo "  OUTPUT          output prefix for alignments (default: query prefix)"
                echo "  PROCESSES       number of databases to evaluate (default: 6)"
                echo "  THREADS         number of threads for each database run (default: 1)"
                echo "  CLEANUP         tidy up intermediate screening files afterward?"
                exit
                ;;
        esac
    done

#Check for values
if [ "$QUERY" == "" ]
then
    echo "Usage: src/run_pipeline.sh -q QUERY -s OUTPUT [-p PROCESSES -t THREADS]"
        echo "  QUERY           query file to align to each database (required)"
        echo "  OUTPUT          output prefix for alignments (default: query prefix)"
        echo "  PROCESSES       number of databases to evaluate (default: 6)"
        echo "  THREADS         number of threads for each database run (default: 1)"
        echo "  CLEANUP         tidy up intermediate screening files afterward? 0 = no, 1 = y"
    exit
fi


if [ "$OUTPUT" == "" ]
then
    OUTPUT=$(basename "$QUERY" | cut -d. -f1)
fi

#Check for database
echo " >> Checking for valid options..."

#Check for input file
if [ ! -f  $QUERY ]
then
    echo " ERROR: input file $QUERY does not exist"
    exit
fi

#################################################################

CM_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#SCRIPT_DIR="$( dirname "$0" )"

export CM_DIR=$CM_DIR
export PYTHONPATH=$CM_DIR

source ${CM_DIR}/../config

# Step 1: biorisk DB scan
date
echo " >> Running biorisk HMM scan..."
transeq $QUERY ${OUTPUT}.faa -frame 6 -clean &>${QUERY}_tmp
hmmscan --domtblout ${OUTPUT}.biorisk.hmmsearch -E 1e-10 biorisk/biorisk.hmm ${OUTPUT}.faa &>${OUTPUT}_tmp
python ${CM_DIR}/check_biorisk.py ${OUTPUT}

# Step 2: taxon ID
# protein screening
date
echo " >> Running taxid screen for regulated pathogens..."
${CM_DIR}/run_diamond.sh -d $NR_DB_DMND -i $QUERY -o ${OUTPUT}.nr -t $THREADS -p $PROCESSES

python ${CM_DIR}/check_reg_path_dmnd_prot.py ${OUTPUT}

# nucleotide screening
date
echo " >> Running taxid screen for regulated pathogen nucleotides..."
python ${CM_DIR}/fetch_nc_bits_dmnd.py ${OUTPUT} ${QUERY}
if [ -f "${OUTPUT}"_nc.fasta ]
then blastn -query ${OUTPUT}_nc.fasta -db $NT_DB -out ${OUTPUT}.nt.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 500 -num_threads 8 -culling_limit 5 -evalue 30
python ${CM_DIR}/check_reg_path_nt.py ${OUTPUT}
fi

### IF A HIT TO A REGULATED PATHOGEN, PROCEED, OTHERWISE CAN FINISH HERE ONCE TESTING IS COMPLETE ####

# Step 3: benign DB scan
echo " >> Checking any pathogen regions for benign components..."
hmmscan --domtblout ${OUTPUT}.benign.hmmsearch -E 1e-10 benign/benign.hmm ${OUTPUT}.faa &>/dev/null
blastn -db ${PFAMDB}/benign/benign.fasta -query $QUERY -out ${OUTPUT}.benign.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-5

python ${CM_DIR}/check_benign.py ${OUTPUT} ${QUERY} # added the original file path here to fetch sequence length, can tidy this

echo " >> Done with screening!"

# Visualising outputs; functional characterization

#python ${CM_DIR}/viz_outputs.py ${OUTPUT} # turning off until file write permissions are resolved

if [ "$CLEANUP" == 1 ]
then
    if [ -f "${OUTPUT}".reg_path_coords.csv ]
    then
        rm ${OUTPUT}.reg_path_coords.csv
    fi
    if [ -f "${OUTPUT}".reg_path_coords_nt.csv ]
    then
        rm ${OUTPUT}.reg_path_coords_nt.csv
    fi
    rm ${OUTPUT}.*hmmsearch ${OUTPUT}.*dmnd ${OUTPUT}.*blastn
fi

rm ${OUTPUT}*_tmp
