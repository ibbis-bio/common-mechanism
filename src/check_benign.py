#! /usr/bin/env python

#################################################################################
#check_benigh.py checks the output from hmmscan and prints to screen the results
#
#Copyright (C) 2022-2023 NTI|Bio
#This file is part of the CommonMechanism
#################################################################################
#Usage: 
#   python check_benign.py -i INPUT -s SEQUENCE -d DATABASE FOLDER 
#       -i, --input = input sample name (will check for sample.benign.hmmscan file)
#       -s, --sequence = input sequence file
#       -d, --database = database folder location/path (will check for benign_annotations.csv) 
#################################################################################
from utils import *
import os, sys, argparse 
import pandas as pd

def check_for_benign(query, coords):
        
        cleared = [0] * coords.shape[0]
        
        # PROTEIN HITS
        # for each set of hits, need to pull out the coordinates covered by benign entries
        hmmscan = query + ".benign.hmmscan"
        if check_blastfile(hmmscan) == 2:
            sys.stdout.write("\t...no housekeeping protein hits\n")
        else:
            hmmer = readhmmer(hmmscan)
    #        print(hmmer)
            for region in range(0, coords.shape[0]): # for each regulated pathogen region
                # look at only the hmmer hits that overlap with it
                htrim = hmmer[~((hmmer['ali from'] > coords['q. end'][region]) & (hmmer['ali to'] > coords['q. end'][region])) & ~((hmmer['ali from'] < coords['q. start'][region]) & (hmmer['ali to'] < coords['q. start'][region]))]
                if htrim.shape[0] > 0:
                    htrim = htrim.assign(coverage = abs(htrim['ali to'] - htrim['ali from']) / htrim['qlen'])
                    if any(htrim['coverage'] > 0.90):
                        htrim = htrim[htrim['coverage'] > 0.90]
                        htrim = htrim.reset_index(drop=True)
                        descriptions = []
                        for row in range(htrim.shape[0]):
                            hit = htrim['target name'][row]
                            hit = hit.replace(".faa.final_tree.fa", "")
                            hit = hit.replace(".faa.final_tree.used_alg.fa", "")
                            descriptions.append(benign_desc['Annotation'][benign_desc['ID'] == hit])
                        annot_string = "\n".join(str(v) for v in descriptions)
                        sys.stdout.write("\t...Housekeeping proteins - >90% coverage of bases " + str(coords['q. start'][region]) + " to " + str(coords['q. end'][region]) + " achieved = PASS\n")
                        sys.stdout.write(annot_string)
                        cleared[region] = 1
                    else:
                        sys.stdout.write("\t...Housekeeping proteins - <90% coverage achieved = FAIL\n")
                    
        # RNA HITS
        # for each set of hits, need to pull out the coordinates covered by benign entries
        cmscan = query + ".benign.cmscan"
        
        if check_blastfile(cmscan) == 2:
            sys.stdout.write("\t...no benign gene hits\n")
        else:
            cmscan = readcmscan(cmscan)
    #        print(cmscan)
            for region in range(0, coords.shape[0]): # for each regulated pathogen region
                # look at only the cmscan hits that overlap with it
                qlen = abs(coords['q. start'][region] - coords['q. end'][region])
                htrim = cmscan[~((cmscan['seq from'] > coords['q. start'][region]) & (cmscan['seq to'] > coords['q. end'][region])) & ~((cmscan['seq from'] < coords['q. end'][region]) & (cmscan['seq to'] < coords['q. start'][region]))]
                if htrim.shape[0] > 0:
                    htrim = htrim.assign(coverage = abs(htrim['seq to'] - htrim['seq from']) / qlen)
                    if any(htrim['coverage'] > 0.90):
                        htrim = htrim[htrim['coverage'] > 0.90]
                        htrim = htrim.reset_index(drop=True)
                        descriptions = []
                        for row in range(htrim.shape[0]):
                            hit = htrim['target name'][row]
                            descriptions.append(hit)
                        annot_string = "\n\t...".join(str(v) for v in descriptions)
                        sys.stdout.write("\t...Housekeeping RNAs - >90% coverage of bases " + str(coords['q. start'][region]) + " to " + str(coords['q. end'][region]) + " achieved: PASS\n")
                        sys.stdout.write("\t...RNA family: " + annot_string + "\n")
                        cleared[region] = 1
                    else:
                        sys.stdout.write("\t...Housekeeping RNAs - <90% coverage achieved = FAIL\n")

        # SYNBIO HITS
        # annotate and clear benign nucleotide sequences
        blast = query + ".benign.blastn"
        if check_blastfile(blast) == 2:
            sys.stdout.write("\t...no Synbio sequence hits\n")
        else:
            blastn = readblast(blast) # synbio parts
            blastn = trimblast(blastn)
            blastn = tophits(blastn)
            for region in range(0, coords.shape[0]): # for each regulated pathogen region
                htrim = blastn[~((blastn['q. start'] > coords['q. end'][region]) & (blastn['q. end'] > coords['q. end'][region])) & ~((blastn['q. start'] < coords['q. start'][region]) & (blastn['q. end'] < coords['q. start'][region]))]
                # print(htrim)
                if any(htrim['q. coverage'] > 0.90):
                    htrim = htrim[htrim['q. coverage'] > 0.90]
                    htrim = htrim.reset_index(drop=True)
                    descriptions = []
                    for row in range(htrim.shape[0]):
                        hit = htrim['subject title'][row]
                        descriptions.append(hit)
                    annot_string = "\n\t...".join(str(v) for v in descriptions)
                    sys.stdout.write("\t...Synbio sequences - >90% coverage achieved = PASS\n")
                    sys.stdout.write("\t...Synbio parts: " + annot_string + "\n")
                    cleared[region] = 1
                else:
                    sys.stdout.write("\t...Synbio sequences - <90% coverage achieved = FAIL\n")
                
            for region in range(0, coords.shape[0]):
                if cleared[region] == 0:
                    sys.stdout.write("\t...Regulated region at bases " + coords.iloc[region, 0] + " to "  + coords.iloc[region, 1] + "failed to clear: FLAG\n")
            if sum(cleared) == len(cleared):
                sys.stdout.write("\n\t...all regulated regions cleared: PASS\n")

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", dest="sample_name",
        required=True, help="Sample name")
    parser.add_argument("-s","--sequence", dest="seq_file",
        required=True, help="FASTA sequence file")
    parser.add_argument("-d","--database", dest="db",
        required=True, help="Benign HMM database folder (must contain benign_annotations.csv)")
    args=parser.parse_args()

    #check input files
    if (not os.path.exists(args.seq_file)):
        sys.stderr.write("\t...sequence file does not exist\n")
        exit(1)
    if (not os.path.exists(args.db + "/benign_annotations.csv")):
        sys.stderr.write("\t...benign_annotations.csv does not exist\n")
        exit(1) 
    
    #Read in database file
    pd.set_option('max_colwidth',200)
    benign_desc = pd.read_csv(args.db + "/benign_annotations.csv")
    
    #Check for file - if exists, check for benign 
    if os.path.exists(args.sample_name + ".reg_path_coords.csv"):
        coords = pd.read_csv(args.sample_name + ".reg_path_coords.csv")
        check_for_benign(args.sample_name, coords)
    else:
        sys.stdout.write("\t...no regulated regions to clear\n")
    
if __name__ == "__main__":
    main()
