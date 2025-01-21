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
    Query to screen for biorisks and 
    Query to screen, based on an input FASTA. Self-calculates AA version.
    TODO: back translate to NT when given AA too.
    """
    def __init__(self, input_fasta_filepath : str):
        self.input_fasta_path = input_fasta_filepath
        self.seq_records : list[SeqRecord] = []

    def translate_query(self) -> None:
        """ Run command transeq, to translate our input sequences. """
        command = ["transeq", self.nt_path, self.aa_path, "-frame", "6", "-clean"]
        result = subprocess.run(command)
        if result.returncode != 0:
            raise RuntimeError(f"Input FASTA {self.nt_path} could not be translated:\n{result.stderr}")

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
