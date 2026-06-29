#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science

# SCREENING
MINIMUM_QUERY_LENGTH = 41
MAXIMUM_QUERY_LENGTH = 100000

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

# I/O
DEFAULT_CONFIG_YAML_PATH = "screen-default-config.yaml"
MAXIMUM_FILENAME_SIZE = 255
MAXIMUM_QUERY_NAME_LENGTH = 64
