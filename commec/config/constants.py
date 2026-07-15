#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

# SCREENING
MINIMUM_QUERY_LENGTH = 41
MAXIMUM_QUERY_LENGTH = 100000
BAD_ACCESSIONS = [
    "YP_009724390",
]

# CONTROL LISTS
N_NON_REGIONAL_HITS_TO_WARN = 2 # How many non-regional hits you need before warnings start showing up.

# NUCLEOTIDE TAXONOMY SCREEN
# At what percent identity a protein hit must attain, to remove that region from non-coding label.
NON_CODING_REGION_PERCENT_IDENTITY_THRESHOLD = 90

# BIORISK E-VALUE FILTERING
# Sequences shorter than this threshold use a length-dependent E-value cutoff.
BIORISK_SHORT_QUERY_NT_THRESHOLD = 200
# Exponent for the length-dependent E-value formula: E < 1 / (1 + L^EXPONENT)
BIORISK_SHORT_QUERY_EVALUE_EXPONENT = 2.598
# E-value cutoff applied to sequences at or above the length threshold.
BIORISK_LONG_QUERY_EVALUE_THRESHOLD = 1e-20

# SEARCH TOOL THREAD LIMITS
HMMSCAN_MAX_THREAD_LIMIT = 4
CMSCAN_MAX_THREAD_LIMIT = 4

# Valid values for the BLAST -mt_mode argument (0: auto, 1: split by database, 2: split by query)
VALID_BLAST_MT_MODES = (0, 1, 2)

# I/O
DEFAULT_CONFIG_YAML_PATH = "screen-default-config.yaml"
MAXIMUM_FILENAME_SIZE = 255
MAXIMUM_QUERY_NAME_LENGTH = 64