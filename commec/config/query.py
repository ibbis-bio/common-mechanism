#!/usr/bin/env python3
# Copyright (c) 2021-2025 International Biosecurity and Biosafety Initiative for Science
import os
import logging
import subprocess
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class Query:
    """
    A query to screen. Contains a sequence record and derived information, such
    as translated sequences.

    At present, we only support nucleotide queries, though we may add suport for
    amino acid queries in future.
    """

    def __init__(self, seq_record: SeqRecord):
        Query.validate_sequence_record(seq_record)
        self.original_name = seq_record.id
        self.name = self.create_id(seq_record.id)
        self.seq_record = seq_record
        self.non_coding_regions : list[tuple[int, int]] = []

    def translate(self, input_path, output_path) -> None:
        """Run command transeq, to translate our input sequences."""

        # TODO: Update line 53-55 of Check_Benign, to ensure that the query filter is using
        # The correct name, when filtering benign components.
        command = ["transeq", input_path, output_path, "-frame", "6", "-clean"]
        result = subprocess.run(command)
        if result.returncode != 0:
            raise RuntimeError(
                f"Input FASTA {input_path} could not be translated:\n{result.stderr}"
            )

    @staticmethod
    def validate_sequence_record(seq_record: SeqRecord) -> None:
        """
        Validate record has non-empty sequence and id. Raises QueryError.
        """
        if not seq_record.id:
            raise QueryValueError(
                "Could not initialize query with an empty sequence id. Is input FASTA valid?"
            )

        if not seq_record.seq:
            raise QueryValueError(
                f"Could not initialize query with id {seq_record.id} because sequence was empty."
                " Is input FASTA valid?"
            )

    @staticmethod
    def create_id(name : str) -> str:
        """
        Parse the Fasta SeqRecord string ID into a 25 digit maximum Unique Identification.
        For internal Commec Screen Use only.
        """
        return name[:25] if len(name) > 24 else name
    

    def get_non_coding_regions(self) -> str:
        """ 
        Return the concatenation of all non-coding regions as a string,
        to be appended to a non_coding fasta file.
        """
        output : str = ""
        for start, end in self.non_coding_regions:
            output += self.seq_record.seq[start-1:end]
        return output

    def nc_to_nt_query_coords(self, index : int) -> int:
        """
        Given an index in non-coding coordinates,
        calculate the nucleotide index in query coordinates.
        """
        nc_pos : int = 0
        for start, end in self.non_coding_regions:
            region_length : int = end - start
            if index < (nc_pos + region_length):
                return index - nc_pos + start
            nc_pos += region_length

        # index was out put bounds of non-coding list of tuples:
        raise QueryValueError(
            f"Non-coding index provided  ({index}) "
             "which is out-of-bounds for any known NC start-end tuple."
            )

class QueryValueError(ValueError):
    """Custom exception for errors when validating a `Query`."""
