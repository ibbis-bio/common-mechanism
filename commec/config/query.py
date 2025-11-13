#!/usr/bin/env python3
# Copyright (c) 2021-2025 International Biosecurity and Biosafety Initiative for Science

import subprocess
import os
from Bio import Seq
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
        self.result : QueryResult = None
        self.translations: list[QueryTranslation] = []
        self.no_hits_warning : bool = True # Updated to False whenever any hit is found.

    @property
    def original_name(self) -> str:
        return str(self._seq_record.id)

    @property
    def length(self) -> int:
        return len(self._seq_record.seq)

    @property
    def sequence(self) -> str:
        return str(self._seq_record.seq)

    def translate(self, output_path: str | os.PathLike) -> None:
        """
        Append the six-frame translation of the query to the output file.
        """
        self._translate()
        with open(output_path, "a", encoding="utf-8") as outfile:
            for translation in self.translations:
                outfile.write(f">{self.name}_{translation.frame}\n")
                outfile.write(f"{translation.sequence}\n")

    def _translate(self) -> None:
        """
        Get the six-frame translations of the query sequence in all 6 reading frames.

        Frame numbers follow the same naming convention as transeq:
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

        One *difference* from transeq is that we only translate full codons. So the frame 1 in the
        example above is translated as MCH, rather than MCHG, even though gg will translate to
        glycine (G) no matter what the subsequent nucleotide is.
        """
        self.translations = []
        seq_rev = Seq.reverse_complement(self.sequence)

        # Frames use offset 0, 1, 2
        for offset in range(3):
            frame_len = self._get_frame_length(offset)

            # Forward frame is offset from the start of the sequence
            f_start = offset
            f_end = offset + frame_len
            protein = str(Seq.translate(self.sequence[f_start:f_end], stop_symbol="X"))
            self.translations.append(
                QueryTranslation(sequence=protein, frame=offset + 1)
            )
            # Reverse frame is offset from the end of the sequence
            r_start = self.length - offset - frame_len
            r_end = self.length - offset
            protein = str(Seq.translate(seq_rev[r_start:r_end], stop_symbol="X"))
            self.translations.append(
                QueryTranslation(sequence=protein, frame=offset + 4)
            )

        # Sort the list in frame order
        self.translations = sorted(self.translations, key=lambda x: x.frame)

    def _get_frame_length(self, frame_offset: int) -> int:
        """
        Get total length of nt sequence that will be translated based on frame offset.
        """
        # Use integer division to get frame length based on starting offset
        return 3 * ((self.length - frame_offset) // 3)

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
    
    def mark_as_hit(self):
        """
        Confirm that this query has had a valid hit, and therefore, has some
        sort of homology to something. 
        
        TODO: In the future, this could also be passed
        coordinate information to mark some areas of the query identified compared
        to other areas.
        """
        self.no_hits_warning = False

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
