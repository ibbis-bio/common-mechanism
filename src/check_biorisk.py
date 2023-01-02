#! /usr/bin/env python

##############################################################################
#check_biorisk.py checks the output from hmmscan and prints to screen the results
#
#Copyright (C) 2022-2023 NTI|Bio 
#This file is part of the CommonMechanism 
##############################################################################
# Usage:
#  python check_biorisk.py -i INPUT.biorisk.hmmsearch -d databases/biorisk_db/ 
##############################################################################
from utils import *
import os, sys, argparse 
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", dest="in_file",
        required=True, help="Input file - hmmscan output file")
    parser.add_argument("-d","--database", dest="db",
        required=True, help="HMM folder (must contain biorisk_lookup.csv)")
    args = parser.parse_args()
    
    #check input files
    if (not os.path.exists(args.in_file)):
        sys.stderr.write("\t...input file does not exist\n") 
        exit(1) 
    if (not os.path.exists(args.db + "/biorisk_lookup.csv")):
        sys.stderr.write("\t...biorisk_lookup.csv does not exist\n")
        exit(1)
    
    #Specify input file and read in database file 
    in_file = args.in_file
    sys.stdout.write("\t...checking %s\n" % in_file) 

    lookup = pd.read_csv(args.db + "/biorisk_lookup.csv")

# read in HMMER output and check for valid hits
res = checkfile(file)
if res == 1:
    hmmer = readhmmer(file)
    hmmer = trimhmmer(hmmer)
    hmmer['description'] = ''
    hmmer = hmmer.reset_index(drop=True)
    # hmmer['target name'] = hmmer['target name'].str.replace("\.", "")
    new_names = []
    for model in range(hmmer.shape[0]):
        name_index = [i for i, x in enumerate([lookup['ID'] == hmmer['target name'][model]][0]) if x]
        # print(name_index)
        # hmmer['description'][model] = lookup['Description'][name_index[0]]
        try:
            new_names.append(lookup['Description'][name_index[0]])
        except:
            new_names.append("")
        # print(lookup['Description'][name_index[0]])
    hmmer['description'] = new_names
    keep1 = [i for i, x in enumerate(hmmer['E-value']) if x < 1e-25]
    hmmer = hmmer.iloc[keep1,:]
    if hmmer.shape[0] > 0:
        print("Biorisks: FLAG\n" + "\n".join(set(hmmer['description'])))
    else:
        sys.stdout.write("\t...Biorisks: unexpected outcome\n")

if __name__ == "__main__":
    main()
