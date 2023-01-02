from utils import *
import sys, os
import pandas as pd

#Input parameter error checking 
if len(sys.argv) < 1:
    sys.stdout.write("\tERROR: Please provide a query file\n")
    exit(1)

file = sys.argv[1] + ".biorisk.hmmsearch"
print("File: ", file)

lookup = pd.read_csv(os.environ['DB_PATH'] + '/biorisk/biorisk_lookup.csv')
# print(lookup.head())

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
        print("Biorisks: PASS")
elif res == 2:
	print("Biorisks: PASS")
else:
	print("Biorisks: unexpected outcome")
