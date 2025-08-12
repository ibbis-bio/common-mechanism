# commec: a free, open-source, globally available tool for DNA sequence screening

<picture style="max-width: 512; display: inline-block;">
	<source media="(prefers-color-scheme: dark)" srcset="https://ibbis.bio/wp-content/uploads/2025/06/COMMEC_Logo_Horiz_Color_IBBIS_onDark.png">
	<img align="left" alt="commec logo" style="max-width: 512;" src="https://ibbis.bio/wp-content/uploads/2025/06/COMMEC_Logo_Horiz_Color_IBBIS_onWhite.png">
</picture>

The `commec` package is a tool for DNA sequence screening that is part of the
[Common Mechanism for DNA Synthesis screening](https://ibbis.bio/common-mechanism/).

Introduction
============
The Common Mechanism offers several sub-commands through the `commec` entrypoint:

    screen  Run Common Mechanism screening on an input FASTA.
    flag    Parse .screen.json files in a directory and create a CSV file of outcomes
    setup   A command-line helper tool to download the required databases
    split   Split a multi-record FASTA file into individual files, one for each record

The `screen` command runs an input FASTA through four steps:

  1. Biorisk Search (uses a hmmer search against custom databases)
  2. Regulated protein Search (uses a BLASTX or DIAMOND search against NCBI nr)
  3. Regulated nucleotide Search (uses BLASTN against NCBI nt)
  4. Low-concern Search (uses hmmer, cmscan and BLASTN against custom databases)

The `.screen.json` file produced by that pipeline can be passed to `flag` to produce a convenient output CSV.
The **screen_pipeline_status.csv** file has the following values:

| Column Name     | Type   | Description  |
| --------------- | ------ | ------------- |
| name            | str    | Name of screen file |
| filepath        | path   | Path to screen file |
| flag            | status | Outcome of the overall screen |
| biorisk         | status | Outcome of the biorisk search |
| protein         | status | Outcome of the protein taxonomy search |
| nucleotide      | status | Outcome of the nucleotide taxonomy search |
| cleared         | status | Outcome of the low-concern search |
| virus_flag      | bool   | Did the taxonomy search flag sequences from a regulated virus? |
| bacteria_flag   | bool   | ... or a regulated bacterial species? |
| eukaryote_flag  | bool   | ... or a regulated eukaryotic (e.g. fungal) species? |
| cleared_protein  | bool   | Were any flags cleared by the low_concern protein search? |
| cleared_rna      | bool   | ... or the low_concern RNA search? |
| cleared_dna      | bool   | ... or the low_concern DNA search? |

The flag column contains the overall recommendation is based on the following decision flow:

![Flowchart showing decision-making by the common mechanism flag module.](https://ibbis.bio/wp-content/uploads/2025/08/commec-screening-flow-v1.jpg "Decision Flow")

Documentation
=============
The online documentation is located at the
[GitHub Wiki](https://github.com/ibbis-screening/common-mechanism/wiki).

Development
===========
Development dependencies are managed through a conda environment. Install conda, and make sure
that [your channels are configured correctly](http://bioconda.github.io/).
Then create the environment with:

```
conda env create -f environment.yml
conda activate commec-dev
```

From here, you should have an interactive version of the package installed via (`pip -e .`) as well
as the necessary shell dependencies.

About
=====
The Common Mechanism is a project of [IBBIS](https://ibbis.bio), the International Biosecurity and
Biosafety Initiative for Science. From 2021-2023, the software and databases were developed by a
team of technical consultants working with the Nuclear Threat Initiative, led by Dr. Nicole Wheeler
of the University of Birmingham, and including contributions from Brittany Rife Magalis of the
University of Louisville and Jennifer Lu of the Center for Computational Biology at Johns Hopkins
University. In 2024, IBBIS became the home of the project
