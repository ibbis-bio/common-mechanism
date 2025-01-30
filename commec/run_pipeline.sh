#! /usr/bin/env bash
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
##############################################################################
#run_pipeline.sh runs the Common Mechanism against a specified QUERY file.
#
#Copyright (C) 2022-2023 NTI|Bio 
#This file is part of the CommonMechanism 
##############################################################################
# Usage:
#   src/run_pipeline.sh
#       -q INPUT_FILE
#       -d DATBASE_FOLDER 
#       -t THREADS (default: 1) 
#       -o OUTPUT (output prefix, default: query name) 
#   Optional parameters:
#       -b = use blast instead of diamond (default: diamond) 
#       -f = use fast mode (default: false)
#       -n = run a nucleotide search if no protein hits are found in a region (default: true)
#       -m = run multiple queries (default: false)
#       -c = clean up intermediate files (default: no cleanup)
# Example Usage: src/run_pipeline.sh -q test_folder/[$name].fasta -d databases/ -p 5 -t 1 -o out_prefix 
##############################################################################
# set parameters for the run
#set -eu 
QUERY=""            #query input file 
DB_PATH="" 
OUTPUT=""           #output prefix
THREADS=1           #threads available
JOBS=""             #number of DIAMOND jobs to run at once
BLAST=0             #blast or diamond (default: diamond) 
FAST_MODE=0         #fast mode (default: off)
NT_SEARCH="on"      #nucleotide search (default: on)
MULTIQUERY=0        #multiple query files (default: off)
CLEANUP=0           #cleanup files (default: off)" 

function usage() {
    echo ""
    echo " Usage: src/run_pipeline.sh -q QUERY -d DB_PATH/ [-o OUTPUT -t THREADS -j JOBS -b -f -n -m -c]"
    echo "    QUERY           query file to align to each database (required)"
    echo "    DB_PATH         path to databases (required)"
    echo "    OUTPUT          output prefix for alignments (default: query prefix)"
    echo "    THREADS         threads available (default: 1)"
    echo "    JOBS            number of diamond runs to do in parallel (default: # CPUs / THREADS)"
    echo " OPTIONAL FLAGS" 
    echo "    -b              run blast for protein screen (default: diamond)"
    echo "    -f              turn on fast mode (no bast match function) (default: off)"
    echo "    -n              run a nucleotide search if no protein hits are found in a region (default: on)"
    echo "    -m              input contains multiple queries (default: off)"
    echo "    -c              tidy up intermediate screening files afterward (default: off)"
    
}

#Get options from user
while getopts "q:d:o:t:j:bfn:mc" OPTION
    do
        case $OPTION in
            q)
                QUERY=$OPTARG
                ;;
            d)
                DB_PATH=$OPTARG
                ;;
            o)
                OUTPUT=$OPTARG
                ;;
            t)
                THREADS=$OPTARG
                ;;
            j)
                JOBS=$OPTARG
                ;;
            b)
                BLAST=1
                ;;
            f)
                FAST_MODE=1
                ;;
            n)
                NT_SEARCH=$OPTARG
                ;;
            m)
                MULTIQUERY=1
                ;;
            c)
                CLEANUP=1
                ;;
            \?) usage ;;
        esac
    done

#INITIALIZING: Check for values
#Check input query 
echo " >> STEP 0: Checking for valid options..." | tee -a ${OUTPUT}.screen
if [ "$QUERY" == "" ]; then
    echo " ERROR: no input query file specified" | tee -a ${OUTPUT}.screen
    print_usage
    exit 1 
elif [ ! -f "$QUERY" ]; then 
    echo " ERROR: specified input query file $QUERY does not exist" | tee -a ${OUTPUT}.screen
    print_usage
    exit 1
fi

#Split files if multiquery is true
if [ "$MULTIQUERY" == 1 ]; then
    echo " >> Splitting multiquery file $QUERY" | tee -a ${OUTPUT}.screen
    python3 ${CM_DIR}/split_query.py -f $QUERY > split_list.txt
    echo "Output split list: split_list.txt"
    echo "Please run the pipeline on each file in split_list.txt with your chosen options" 
    # alternative: have the script run each file according to options chosen - ignores any job management system
    exit 0
fi

#If output not specified, set value 
if [ "$OUTPUT" == "" ]; then
    OUTPUT=${QUERY%.*}
    if [ ${#OUTPUT} -ge 200 ]; then
        OUTPUT=`echo $OUTPUT | cut -c1-200`
    fi
    echo "Output handle = $OUTPUT"
fi

#Check input database folder 
if [ "$DB_PATH" == "" ]; then
    echo " ERROR: screening database path not specified" | tee -a ${OUTPUT}.screen
    print_usage
    exit 1
elif [ ! -d "$DB_PATH" ]; then
    echo " ERROR: specified database folder does not exist" | tee -a ${OUTPUT}.screen
    print_usage
    exit 1
else
    #Test database downloads 
    biorisk_test=${DB_PATH}/biorisk_db/biorisk.hmm
    benign_test=${DB_PATH}/benign_db/benign.hmm
    split_nr_test=${DB_PATH}/nr_dmnd/nr.1.dmnd
    if [ ! -f "$biorisk_test" ]; then
        echo " ERROR: database folder does not contain biorisk database" | tee -a ${OUTPUT}.screen
        exit 1 
    fi 
    if [ ! -f "$benign_test" ]; then
        echo " ERROR: database folder does not contain benign database"  | tee -a ${OUTPUT}.screen
        exit 1 
    fi 
fi
# export DB_PATH=${DB_PATH}

#################################################################
CM_DIR="$( dirname "$0" )"

start_time=$(date)
echo -e " >> STARTED AT $start_time" | tee -a ${OUTPUT}.screen
echo -e " >> Screening $QUERY" | tee -a ${OUTPUT}.screen

# remove any spaces or blank characters
cp ${QUERY} ${OUTPUT}.tmp
cat ${OUTPUT}.tmp | sed -E 's/[[:space:]]|\xc2\xa0/_/g' > ${OUTPUT}.fasta

# Step 1: biorisk DB scan
echo " >> STEP 1: Checking for biorisk genes..."  | tee -a ${OUTPUT}.screen

echo -e "\t...running transeq" 
transeq $QUERY ${OUTPUT}.faa -frame 6 -clean >> ${OUTPUT}.tmp 2>&1
if [ ! -f "${OUTPUT}".faa ]; then
    echo -e "\t ERROR: transeq failed" | tee -a ${OUTPUT}.screen
fi

#python3 ${CM_DIR}/concat_seqs.py ${OUTPUT}.faa

echo -e "\t...running hmmscan" 
hmmscan --domtblout ${OUTPUT}.biorisk.hmmscan ${DB_PATH}/biorisk_db/biorisk.hmm ${OUTPUT}.faa >> ${OUTPUT}.tmp 2>&1
echo -e "\t...checking hmmscan results"
python3 ${CM_DIR}/check_biorisk.py -i ${OUTPUT}.biorisk.hmmscan --database ${DB_PATH}/biorisk_db/ | tee -a ${OUTPUT}.screen

s1_time=$(date)
echo -e "    STEP 1 completed at $s1_time\n" | tee -a ${OUTPUT}.screen

# if [ "$FAST_MODE" = 0 ]; then
    # Step 2: taxon ID/protein screening
    echo " >> STEP 2: Checking regulated pathogen proteins..." | tee -a ${OUTPUT}.screen

    if [ "$BLAST" = 1 ]; then
        if [ ! -f "${OUTPUT}".nr.blastx ]; then
            echo -e "\t...running run_blastx.sh"
            ${CM_DIR}/run_blastx.sh -d $DB_PATH/nr_blast/nr -q ${OUTPUT}.fasta -o ${OUTPUT}.nr -t $THREADS # use the shortened filename rather than the original
        fi
        echo -e "\t...checking blast results"
        if [ -f "${OUTPUT}".reg_path_coords.csv ]; then 
            rm "${OUTPUT}".reg_path_coords.csv
        fi
        python3 ${CM_DIR}/check_reg_path.py -i ${OUTPUT}.nr.blastx -d $DB_PATH -t $THREADS | tee -a ${OUTPUT}.screen
    else 
        if [ ! -f "${OUTPUT}".nr.dmnd ]; then
            echo -e "\t...running run_diamond.sh"
            cmd=(
                "${CM_DIR}/run_diamond.sh"
                -d "${DB_PATH}/nr_dmnd/"
                -i "${OUTPUT}.fasta"
                -o "${OUTPUT}.nr"
                -t "${THREADS}"
            )
            # Add -j option only if JOBS is set and non-empty
            if [[ -n "${JOBS:-}" ]]; then
                cmd+=(-j "${JOBS}")
            fi
            # Execute the command
            "${cmd[@]}"
        fi
        echo -e "\t...checking diamond results"
        if [ -f "${OUTPUT}".reg_path_coords.csv ]; then 
            rm "${OUTPUT}".reg_path_coords.csv
        fi
        python3 ${CM_DIR}/check_reg_path.py -i ${OUTPUT}.nr.dmnd --database $DB_PATH -t $THREADS | tee -a ${OUTPUT}.screen
    fi

    s2_time=$(date)
    echo -e "    STEP 2 completed at $s2_time\n" | tee -a ${OUTPUT}.screen

    # Step 3: nucleotide screening
    if [ "$NT_SEARCH" = "on" ]; then
        echo " >> STEP 3: Checking regulated pathogen nucleotides..." | tee -a ${OUTPUT}.screen

        echo -e "\t...fetching noncoding regions"
        if [ "$BLAST" = 1 ]; then
            python3 ${CM_DIR}/fetch_nc_bits.py ${OUTPUT}.nr.blastx ${OUTPUT}.fasta | tee -a ${OUTPUT}.screen
        else
            python3 ${CM_DIR}/fetch_nc_bits.py ${OUTPUT}.nr.dmnd ${OUTPUT}.fasta | tee -a ${OUTPUT}.screen
        fi

        if [ -f "${OUTPUT}".noncoding.fasta ]; then 
            if [ ! -f "${OUTPUT}".nt.blastn ]; then
                echo -e "\t...running blastn"
                blastn -query ${OUTPUT}.noncoding.fasta -db ${DB_PATH}/nt_blast/nt -out ${OUTPUT}.nt.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -max_target_seqs 50 -num_threads 8 -culling_limit 5 -evalue 10
            fi
            echo -e "\t...checking blastn results"
            python3 ${CM_DIR}/check_reg_path.py -i ${OUTPUT}.nt.blastn -d $DB_PATH -t $THREADS | tee -a ${OUTPUT}.screen
        else 
            echo -e "\t...skipping nucleotide search"
        fi

        s3_time=$(date)
        echo -e "    STEP 3 completed at $s3_time\n" | tee -a ${OUTPUT}.screen
    else
        echo " >> SKIPPING STEP 3: Nucleotide search" | tee -a ${OUTPUT}.screen
    fi

else
    echo " >> FAST MODE: Skipping steps 2-3" | tee -a ${OUTPUT}.screen

fi

# Step 4: benign DB scan
#date
echo -e " >> STEP 4: Checking any pathogen regions for benign components..." | tee -a ${OUTPUT}.screen

hmmscan --domtblout ${OUTPUT}.benign.hmmscan ${DB_PATH}/benign_db/benign.hmm ${OUTPUT}.faa >> ${OUTPUT}.tmp 2>&1
blastn -db ${DB_PATH}/benign_db/benign.fasta -query ${OUTPUT}.fasta -out ${OUTPUT}.benign.blastn -outfmt "7 qacc stitle sacc staxids evalue bitscore pident qlen qstart qend slen sstart send" -evalue 1e-5
cmscan --tblout ${OUTPUT}.benign.cmscan ${DB_PATH}/benign_db/benign.cm ${OUTPUT}.fasta >> ${OUTPUT}.tmp 2>&1

python3 ${CM_DIR}/check_benign.py -i ${OUTPUT} -d ${DB_PATH}/benign_db/ | tee -a ${OUTPUT}.screen

finish_time=$(date)
echo -e " >> COMPLETED AT $finish_time" | tee -a ${OUTPUT}.screen

# Visualising outputs; functional characterization

#python3 ${CM_DIR}/viz_outputs.py ${OUTPUT} # turning off until file write permissions are resolved

if [ "$CLEANUP" == 1 ]
then
    if [ -f "${OUTPUT}".reg_path_coords.csv ]
    then
        rm ${OUTPUT}.reg_path_coords.csv
    fi

    rm ${OUTPUT}.*hmmscan ${OUTPUT}.*blastn ${OUTPUT}.faa

    if ["$BLAST" = 1 ]
    then
        rm ${OUTPUT}.*blastx 
    else 
        rm ${OUTPUT}.*dmnd
    fi
fi

rm ${OUTPUT}*.tmp
