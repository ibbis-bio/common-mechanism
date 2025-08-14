# commec: a free, open-source, globally available tool for DNA sequence screening

<picture style="max-width: 512; display: inline-block;">
	<source media="(prefers-color-scheme: dark)" srcset="https://ibbis.bio/wp-content/uploads/2025/06/COMMEC_Logo_Horiz_Color_IBBIS_onDark.png">
	<img align="left" alt="commec logo" style="max-width: 512;" src="https://ibbis.bio/wp-content/uploads/2025/06/COMMEC_Logo_Horiz_Color_IBBIS_onWhite.png">
</picture>

The `commec` package is a tool for DNA sequence screening that is part of the
[Common Mechanism for DNA Synthesis screening](https://ibbis.bio/common-mechanism/). The package offers several sub-commands through the `commec` entrypoint:

    screen  Run Common Mechanism screening on an input FASTA.
    flag    Parse .screen.json files in a directory and create a CSV file of outcomes
    setup   A command-line helper tool to download the required databases
    split   Split a multi-record FASTA file into individual files, one for each record

The `commec screen` command runs an input FASTA through the following screening steps:

1. **Biorisk search**: Fast HMM-based search against curated sequence profiles
2. **Taxonomy Search**: look for best matches to regulated pathogens using a two-step process:
   * **Protein search**: BLASTX/DIAMOND search against NCBI nr
   *  **Nucleotide search**: BLASTN search against NCBI core_nt
3. **Low concern search**: Clear earlier flags based on matches to common or conserved sequences

![Flowchart showing decision-making by the common mechanism flag module.](https://ibbis.bio/wp-content/uploads/2025/08/commec-screening-flow-v1.jpg "Decision Flow")

Information about the databases supporting screening can be found in the [commec-databases](https://github.com/ibbis-bio/commec-databases/) repostiory.

## Documentation
The [GitHub Wiki](https://github.com/ibbis-screening/common-mechanism/wiki) has documentation for this package, including information about installing `commec` and interpreting screening results.

More information about the Common Mechanism project is available on the [IBBIS project page](https://ibbis.bio/common-mechanism/) and [Common Mechanism FAQ](https://ibbis.bio/our-work/common-mechanism/faq/).

## Development
The `commec` package is being actively developed by IBBIS staff. We welcome contributions! To get started, install conda, and make sure
that [your channels are configured correctly](http://bioconda.github.io/). Then create the dev environment with:

```
conda env create -f environment.yml
conda activate commec-dev
```

From here, you should have an interactive version of the package installed via `pip -e .` and the necessary shell dependencies.
