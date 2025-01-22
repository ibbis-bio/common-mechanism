#!/usr/bin/env python3
# Copyright (c) 2021-2025 International Biosecurity and Biosafety Initiative for Science
import os
import subprocess
from Bio import Seq
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

    def _write_six_frame_translation(self):
        """
        Write a file with translations of the query records in all 6 reading frames.
       
        Each record is named with the id of the sequence record, followed by "_index", where
        the frame indexes, following the same format as transeq, are:
          * _1, _2, _3: Forward frames starting at positions 0, 1, 2
          * _4, _5, _6: Reverse frames, starting at positions 0, 1, 2

        As in previous `transeq -clean` command, all stop codons (*) are replaced with (X).
        """
        with open(self.aa_path, "w", encoding="utf-8") as fout:
            for record in self.seq_records:
                seq = str(record.seq)
                seq_rev = Seq.reverse_complement(seq)
                seq_len = len(seq)

                for i in range(3):
                    # Use integer division to get frame length 
                    frame_len = 3 * ((seq_len - i) // 3)
                    # Forward frame
                    protein = Seq.translate(seq[i:i + frame_len], stop_symbol="X")
                    frame = i+1
                    fout.write(f">{record.id}_{frame}\n{protein}\n")
                    # Reverse frame
                    protein = Seq.translate(seq_rev[i:i + frame_len], stop_symbol="X")
                    rev_frame = i+4
                    fout.write(f">{record.id}_{rev_frame}\n{protein}\n")

class QueryValueError(ValueError):
    """Custom exception for errors when validating a `Query`."""

