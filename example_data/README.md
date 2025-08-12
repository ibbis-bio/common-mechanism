# Example data

This directory contains a file, `commec-examples.fasta`, which inludes queries illustrating different possible screening outcomes, as well as the results of running `commec screen` on that file.

A guide to interpreting these results is provided in the [Tutorial](https://github.com/ibbis-bio/common-mechanism/wiki/tutorial) on the `commec` wiki.

### Examples included

* **BBa_K5108009_creA_** (`WARN`): This is a [composite DNA part](https://parts.igem.org/Part:BBa_K5108009) developed by 2024 iGEM team Toulouse-INSA-UPS for space exploration applications. It is an artificial operon composed of four basic parts: creatinase and creatinine amidohydrolase ORFs (creA BBa K5108003, crnA BBa K5108004) and two RBS (BBa K5108006, BBaK5108007) enabling their expression in the plant growth-promoting rhizobacteria, _Pseudomonas fluorescens_, enabling the metabolization of creatinine by this organism.
* **encrypted** (`WARN`): This DNA sequence contains an encrypted message generated using the [CryptoGErM](https://2016.igem.org/Team:Groningen/Tour) algorithm developed by the 2016 iGEM team from Groningen. It is therefore an entirely artificial sequence, with no biological function or related taxonomy across the domains of life.
* **xylanase_zero_shot_des31** (`PASS`):  This sequence is one of the xylanase variants used in the zero-shot enzyme activity prediction challenge problem from [Align Bio’s](https://alignbio.org/) 2023 [Protein Engineering Tournament](https://alignbio.org/tournamentpilot-results-2023). Xylanase is an enzyme that degrades the second-most-abundant polysaccharide and should not be flagged.
* **BBa_K209429_A_15261** (`PASS`): This sequence is another [composite DNA part](https://parts.igem.org/Part:BBa_K209429) created by the [igem UCSF team in 2009](https://2009.igem.org/Team:UCSF) with the goal of manipulating signaling pathways to mediate chemotaxis.
* **RVFV_Rift_valley_fever** (`FLAG`): The Rift Valley Fever virus sample is successfully flagged during the taxonomic steps as containing extensive regions of material from regulated organisms, namely nucleocapsid proteins from _Phlebovirus riftense_ AKA Rift Valley Fever virus.
