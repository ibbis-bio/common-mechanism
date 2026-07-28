# commec: a free, open-source, globally available tool for DNA sequence screening

<picture style="max-width: 512; display: inline-block;">
	<source media="(prefers-color-scheme: dark)" srcset="https://ibbis.bio/wp-content/uploads/2025/06/COMMEC_Logo_Horiz_Color_IBBIS_onDark.png">
	<img align="left" alt="commec logo" style="max-width: 512;" src="https://ibbis.bio/wp-content/uploads/2025/06/COMMEC_Logo_Horiz_Color_IBBIS_onWhite.png">
</picture>

The `commec` package is a tool for DNA sequence screening that is part of the
[Common Mechanism for DNA Synthesis screening](https://ibbis.bio/common-mechanism/). The package offers several sub-commands through the `commec` entrypoint:

    setup   Download or update the reference databases required for screening
    screen  Run Common Mechanism screening on an input FASTA.
    flag    Parse .screen or .json files in a directory and create CSVs of flags raised
    list    Display information on available annotated control lists

The `commec screen` command runs an input FASTA through the following screening steps:

1. **Biorisk search**: Fast HMM-based search against curated sequence profiles
2. **Best Match / Taxonomy Search**: look for best matches to regulated pathogens using a two-step process:
   * **Protein search**: BLASTX search against custom-filtered UniRef clusters
   * **Nucleotide search**: BLASTN search against curated viral reference genomes
3. **Low concern search**: Clear earlier flags based on matches to common or conserved sequences

![Flowchart of the commec screening steps, from input FASTA through the biorisk, best match / taxonomy, and low-concern screening to a Warning, Flag, or Clear outcome.](docs/images/common-mechanism-screening-flow-v2.jpg "Commec screening decision flow")

The [GitHub Wiki](https://github.com/ibbis-screening/common-mechanism/wiki) has documentation for this package, including information about installing `commec` and interpreting screening results.

More information about the Common Mechanism project is available on the [commec.ibbis.bio](https://commec.ibbis.bio/) and [IBBIS project page](https://ibbis.bio/common-mechanism/).

## Quick start
`commec` and its package dependencies are installed via conda. The reference databases that support screening are downloaded and kept up to date with the `commec setup` command, which retrieves them from [databases.commec.io](https://databases.commec.io).

Once installed, download the reference databases and screen a FASTA file with:

```
commec setup -d /path/to/databases
commec screen -d /path/to/databases input.fasta
```

See the [GitHub Wiki](https://github.com/ibbis-screening/common-mechanism/wiki) for full installation instructions.

## Run with a container

A container image with `commec` and all of its dependencies (BLAST+, HMMER, Infernal) is published
to the GitHub Container Registry:

```
docker pull ghcr.io/ibbis-bio/common-mechanism:latest
```

The examples below use `docker`; `podman` accepts the same commands. The image entrypoint is
`commec`, so sub-commands are passed straight through.

**The reference databases are not included in the image.** They are tens of gigabytes, are
versioned independently of the package, and are fetched by `commec setup`. Keep them on the host
and bind-mount them.

### Download the databases

This writes tens of gigabytes to the host directory, which must be writable:

```
mkdir -p commec-databases
docker run --rm \
  -v "$PWD/commec-databases:/databases" \
  ghcr.io/ibbis-bio/common-mechanism:latest \
  setup -d /databases
```

### Screen a FASTA

```
docker run --rm \
  -v "$PWD/commec-databases:/databases:ro" \
  -v "$PWD:/data" \
  ghcr.io/ibbis-bio/common-mechanism:latest \
  screen -d /databases input.fasta
```

There are two mounts:

* `/databases` — the reference databases, passed to `commec` with `-d`. A read-only mount is fine
  for `screen`, but **not** for `setup`, and not if you set `auto_update_databases: true` in a
  config file.
* `/data` — the image's working directory. Input FASTA files are read from here and `.screen`,
  `.json`, and report output is written back here.

Instead of `-d`, you can supply a full config file with `-y /data/config.yaml`; the paths in
`base_paths.default` must be absolute paths as seen from inside the container.

### Notes

* **File ownership.** The image runs as UID 1000. On Linux, add `--user "$(id -u):$(id -g)"` if
  your host UID differs, so output files are owned by you.
* **Threads.** Pass `screen --threads N` to use more cores, and constrain the container with
  `--cpus N` to match.
* **Interactive shell.** `docker run --rm -it --entrypoint bash ghcr.io/ibbis-bio/common-mechanism:latest`
* **Architecture.** The image is `linux/amd64` only, because bioconda does not publish
  `linux-aarch64` builds of BLAST+. It runs on Apple Silicon and other arm64 hosts under emulation,
  which is noticeably slower.

### Build the image locally

```
docker build --target runtime -t commec:local -f Containerfile .
```

The `test` target adds `pytest` and the source tree, and runs the test suite inside the image:

```
docker build --target test -t commec:test -f Containerfile .
docker run --rm commec:test
```

## Development
The `commec` package is being actively developed by IBBIS staff. We welcome contributions! To get started, install conda, and make sure
that [your channels are configured correctly](http://bioconda.github.io/). Then create the dev environment with:

```
conda env create -f environment.yaml
conda activate commec-dev
```

From here, you should have an interactive version of the package installed via `pip -e .` and the necessary shell dependencies.

## License
`commec` is released under the [MIT License](LICENSE).
