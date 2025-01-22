#!/usr/bin/env python3
# Copyright (c) 2021-2025 International Biosecurity and Biosafety Initiative for Science
import os
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
        self.seq_record = seq_record
        self.translations: list[QueryTranslation] = []

    @property
    def name(self):
        return self.seq_record.id

    @property
    def length(self):
        return len(self.seq_record.seq)

    @property
    def sequence(self):
        return self.seq_record.seq

    def translate(self) -> None:
        """
        Get the six-frame translations of the query sequence in all 6 reading frames.

        We follow the same naming format as transeq:
          * 1, 2, 3: Forward frames starting at positions 0, 1, 2
          * 4, 5, 6: Reverse frames, starting at positions 0, -1, -2

        Offsets are a little complicated. Taking the 11-nt sequence 'atgtgccatgg' as an example:

            Frame   Pos     Codon split         Translation
            1       0       atg tgc cat gg      MCH
            2       1       a tgt gcc atg g     CAM
            3       2       at gtg cca tgg      VPW
            4       -0      cc atg gca cat      MAH
            5       -1      c cat ggc aca t     HGT
            6       -2      cca tgg cac at      PWH

        As in previous `transeq -clean` command, all stop codons (*) are replaced with (X).
        """
        seq_rev = Seq.reverse_complement(self.sequence)

        for i in range(3):
            # Use integer division to get frame length based on starting offset
            frame_len = 3 * ((self.length - i) // 3)

            # Forward frame translates from beginning to frame length
            f_start = i
            f_end = i + frame_len
            protein = Seq.translate(self.sequence[f_start:f_end], stop_symbol="X")
            self.translations.append(
                QueryTranslation(
                    sequence=protein, frame=i, nt_start=f_start, nt_end=f_end
                )
            )

            # Reverse frame is offset from the end of the sequence
            r_end = self.length - i
            r_start = r_end - frame_len
            protein = Seq.translate(seq_rev[r_start:r_end], stop_symbol="X")
            # I think these indices are wrong actually
            self.translations.append(
                QueryTranslation(
                    sequence=protein, frame=4+i, nt_start=r_start, nt_end=r_end
                )
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


class QueryTranslation:
    """
    Represents a single frame translation of a nucleotide sequence.

    Attributes:
        sequence (str): The translated amino acid sequence
        frame (int): Frame number (1-6) following transeq convention:
          * 1, 2, 3: Forward frames starting at positions 0, 1, 2
          * 4, 5, 6: Reverse frames, starting at positions 0, -1, -2
        nt_start (int): Start position in the nucleotide sequence (0-based)
        nt_end (int): End position in the nucleotide sequence (0-based, exclusive)
    """

    sequence: str
    frame: int  # 1-6
    nt_start: int
    nt_end: int


class QueryValueError(ValueError):
    """Custom exception for errors when validating a `Query`."""
