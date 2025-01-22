"""
Container class to hold information pertaining to query from an input fasta file,
 as well as derived information, such as translated sequences and whether or not
 the query was derived from AA or NT.
"""

import os
import subprocess
from Bio.SeqRecord import SeqRecord


class Query:
    """
    A query to screen. Contains a sequence record and derived information, such
    as translated sequences.
    
    At present, we only support nucleotide queries, though we may add suport for
    amino acid queries in future.
    """
    def __init__(self, seq_record: SeqRecord):
        Query.validate_sequence_record(seq_record)
        self.name = seq_record.id
        self.seq_record = seq_record

    def translate(self, input_path, output_path) -> None:
        """Run command transeq, to translate our input sequences."""
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

class QueryValueError(ValueError):
    """Custom exception for errors when validating a `Query`."""
