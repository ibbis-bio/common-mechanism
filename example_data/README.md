# Example data

This directory contains two example files, `commec-examples-1.fasta` and `commec-examples-2.fasta`, which includes queries illustrating different possible screening outcomes, as well as the results of running `commec screen` on these files.

A guide to interpreting these results is provided in the [Tutorial](https://github.com/ibbis-bio/common-mechanism/wiki/tutorial) on the `commec` wiki.

### General Examples-1

* **M14107.1 Plasmid pHly152 (from E.coli) hlyA, B** ('WARN'): Fragment of the plasmid encoding genes for Hemolysin A and B, known virulence factors of E. coli and other bacteria. A viruence factor does not directly cause harm by itself, which results in a warn outcome.
* **xylanase_zero_shot_des31** (`PASS`):  This sequence is one of the xylanase variants used in the zero-shot enzyme activity prediction challenge problem from [Align Bio’s](https://alignbio.org/) 2023 [Protein Engineering Tournament](https://alignbio.org/tournamentpilot-results-2023). Xylanase is an enzyme that degrades the second-most-abundant polysaccharide and should not be flagged.
* **BBa_K5108009_creA_** (`PASS`): This is a [composite DNA part](https://parts.igem.org/Part:BBa_K5108009) developed by 2024 iGEM team Toulouse-INSA-UPS for space exploration applications. It is an artificial operon composed of four basic parts: creatinase and creatinine amidohydrolase ORFs (creA BBa K5108003, crnA BBa K5108004) and two RBS (BBa K5108006, BBaK5108007) enabling their expression in the plant growth-promoting rhizobacteria, _Pseudomonas fluorescens_, enabling the metabolization of creatinine by this organism.
* **BBa_K209429_A_15261** (`PASS`): This sequence is another [composite DNA part](https://parts.igem.org/Part:BBa_K209429) created by the [igem UCSF team in 2009](https://2009.igem.org/Team:UCSF) with the goal of manipulating signaling pathways to mediate chemotaxis.
* **addgene-plasmid-36399** (`FLAG`): This is the [pVSV Venus VSV-G plasmid from addgene](https://www.addgene.org/36399/), which contains Vesiculovirus glycoprotein, which is widely used but can be used to generate VSV, which is a controlled organism from many jurisdictions.
* **RVFV_Rift_valley_fever** (`FLAG`): The Rift Valley Fever virus sample is successfully flagged during the taxonomic steps as containing extensive regions of material from regulated organisms, namely nucleocapsid proteins from _Phlebovirus riftense_ AKA Rift Valley Fever virus.
* **encrypted** (`PASS`): This DNA sequence contains an encrypted message generated using the [CryptoGErM](https://2016.igem.org/Team:Groningen/Tour) algorithm developed by the 2016 iGEM team from Groningen. It is therefore an entirely artificial sequence, with no biological function or related taxonomy across the domains of life.

### Regional Examples-2

These examples show some differences between control lists tracked by `commec`. We recommend screening the examples yourself using `--region EU, US`, or `--region BR` to see how commec captures and reports on the unique requirements across jurisdictions. The default results have been generated with no jurisdiction filter.

* **Cadang-cadang** (`FLAG` (Brazil)): Cadang-cadang is a disease caused by Coconut cadang-cadang viroid (CCCVd, Cocadviroid cadangi), a lethal viroid of several palms including coconut, African oil palm , anahaw , and buri. It is of concern for Brazil.
* **Paenibacillus** (`FLAG` (European Union)): Bacterium causative agent of American Foulbrood, a disease that is fatal for larval honey bees, and a part of the European Union List of Animal Pathogens.
* **pUC19** (`FLAG` (United Kingdom, China)): A pUC19 vector containing the spike protein from SARS-CoV-2. Although SARS-CoV-2 is a pandemic pathogen, its endemicity and frequent use in research means it is not controlled in every jurisdiction.
* **addgene-plasmid-111773** (`FLAG`): This is the This is the [pGEMT-PTE2A plasimd from addgene](https://www.addgene.org/111773/), which contains 2A sequences from three different organisms, including one from Porcine Teschovirus (P2A), which is controlled by a wide variety of lists under various terms.
