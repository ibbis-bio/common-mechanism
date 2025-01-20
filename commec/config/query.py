"""
Container class to hold information pertaining to query from an input fasta file,
 as well as derived information, such as translated sequences and whether or not 
 the query was derived from AA or NT.
"""
import os
import subprocess
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

class Query:
    """
    Query to screen, based on an input FASTA. Self-calculates AA version.
    TODO: back translate to NT when given AA too.
    """
    def __init__(self, input_fasta_filepath : str):
        self.input_fasta_path = input_fasta_filepath
        self.nt_path : str = ""
        self.aa_path : str = ""
        self.seq_records : list[SeqRecord] = []
        self.non_coding_regions : list[list[int]] = []

    def setup(self, output_prefix : str) -> None:
        """ 
        Parse the input FASTA file and set up the NT and AA files that will be used as queries.
        """
        self.nt_path = f"{output_prefix}.cleaned.fasta"
        self.aa_path = f"{output_prefix}.faa"
        self._write_clean_fasta()
        self.parse_query_data()
        self._write_six_frame_translation()

    def parse_query_data(self) -> None:
        """
        Populate a list of query names, and associated sequences, for json formatting purposes.
        """
        with open(self.nt_path, "r", encoding = "utf-8") as fasta_file:
            self.seq_records : list[SeqRecord] = list(SeqIO.parse(fasta_file, "fasta"))

    def _write_clean_fasta(self) -> str:
        """
        Write a FASTA in which whitespace (including non-breaking spaces) and 
        illegal characters are replaced with underscores.
        """
        with (
            open(self.input_fasta_path, "r", encoding="utf-8") as fin,
            open(self.nt_path, "w", encoding="utf-8") as fout,
        ):
            for line in fin:
                line = line.strip()
                modified_line = "".join(
                    "_" if c.isspace() or c == "\xc2\xa0" or c == "#" else c
                    for c in line
                )
                fout.write(f"{modified_line}{os.linesep}")


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

    def get_non_coding_regions(self) -> str:
        """ 
        Return the concatenation of all non-coding regions as a string,
        to be appended to a non_coding fasta file.
        """
        output : str = ""
        for start, end in self.non_coding_regions:
            output += self.seq_records[start-1:end]
        return output

    def convert_noncoding_index_to_query_index(self, index : int) -> int:
        """
        Given an index in non-coding space, calculate the index in query space.
        """
        nc_pos : int = 0
        for start, end in self.non_coding_regions:
            region_length : int = end - start
            if index < (nc_pos + region_length):
                return index - nc_pos + start
            nc_pos += region_length
