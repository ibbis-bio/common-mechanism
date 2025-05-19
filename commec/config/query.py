#!/usr/bin/env python3
# Copyright (c) 2021-2025 International Biosecurity and Biosafety Initiative for Science

import subprocess
from Bio.SeqRecord import SeqRecord
from dataclasses import dataclass
from commec.config.result import QueryResult

class Query:
    """
    A query to screen. Contains a sequence record and derived information, such
    as translated sequences.

    At present, we only support nucleotide queries, though we may add suport for
    amino acid queries in future.
    """

    def __init__(self, seq_record: SeqRecord):
        Query.validate_sequence_record(seq_record)
        self._seq_record = seq_record
        self.name = self.create_id(seq_record.id)
        self.non_coding_regions : list[tuple[int, int]] = [] # 1 based coordinates for Non-Coding Regions.
        self.result_handle : QueryResult = None
        self.translations: list[QueryTranslation] = []

    @property
    def original_name(self) -> str:
        return str(self._seq_record.id)

    @property
    def length(self) -> int:
        return len(self._seq_record.seq)

    @property
    def sequence(self) -> str:
        return str(self._seq_record.seq)

    def translate(self, input_path, output_path) -> None:
        """Run command transeq, to translate our input sequences."""

        # TODO: Update line 53-55 of Check_Benign, to ensure that the query filter is using
        # The correct name, when filtering benign components.
        command = ["transeq", input_path, output_path, "-frame", "6", "-clean"]
        result = subprocess.run(command, check = False, capture_output=True)
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
    def create_id(input_name : str) -> str:
        """
        Parse the Fasta SeqRecord string ID into 
        a 25 digit maximum Unique Identification.
        For internal Commec Screen Use only.
        Original Fasta name IDs are used during JSON output.
        """
        # Protections on split name conventions for translation output:
        name = input_name.split('|')[0]

        # Protection on accidentally ending with _
        if name.endswith("_"):
            name = name[:-1]

        if len(name) < 26:
            return name
        
        tokens = name.split("_")
        output = ""
        for i in range(len(tokens)):
            testname = "_".join(tokens[:i])
            if len(testname) > 25:
                break
            output = testname

        return output
    

    def get_non_coding_regions_as_fasta(self) -> str:
        """ 
        Return the concatenation of all non-coding regions as a string,
        to be appended to a non_coding fasta file.
        """
        if len(self.non_coding_regions) == 0:
            return ""
        heading : str = f">{self.name}"
        sequence : str = ""
        for start, stop in self.non_coding_regions:
            heading+=f" ({start}-{stop})"
            sequence+=f"{self._seq_record.seq[int(start)-1: int(stop)]}"
        return f"{heading}\n{sequence}\n"

    def nc_to_nt_query_coords(self, index : int) -> int:
        """
        Given an index in non-coding coordinates,
        calculate the nucleotide index in query coordinates.
        """
        nc_pos : int = 1
        for start, end in self.non_coding_regions:
            region_length : int = end - start
            if (index <= (nc_pos + region_length) and
                index >= nc_pos):
                return index - nc_pos + start

            nc_pos += region_length + 1

        # index was out put bounds of non-coding list of tuples:
        raise QueryValueError(
            f"Non-coding index provided  ({index}) for {self.name}"
            f"which is out-of-bounds for any known NC start-end tuple: {self.non_coding_regions}"
            )

@dataclass
class QueryTranslation:
    """
    Represents a single frame translation of a nucleotide sequence.

    Attributes:
        sequence (str): The translated amino acid sequence
        frame (int): Frame number following transeq convention (1-3: forward, 4-6: reverse)
    """

    sequence: str
    frame: int

class QueryValueError(ValueError):
    """Custom exception for errors when validating a `Query`."""