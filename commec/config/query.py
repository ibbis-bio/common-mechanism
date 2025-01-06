"""
Container class to hold information pertaining to query from an input fasta file,
 as well as derived information, such as translated sequences, whether or not 
 the query was derived from AA or NT.
"""
import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Query:
    """
    Query to screen, based on an input FASTA. Self-calculates AA version.
    Future direction to back translate to NT when given AA too.
    """
    def __init__(self, input_fasta_filepath : str):
        self.input_fasta_path = input_fasta_filepath
        self.nt_path : str = ""
        self.aa_path : str = ""
        self.raw : list[SeqRecord] = []
        self.non_coding_regions : list[list[int]] = [[]]

    def setup(self, output_prefix : str):
        """ 
        Translate or reverse translate query, to be ready in AA or NT format. 
        """
        self.nt_path = self.get_cleaned_fasta(output_prefix)
        self.aa_path = f"{output_prefix}.transeq.faa"
        self.parse_query_data()

    def translate_query(self):
        """ Run command transeq, to translate our input sequences. """
        command = ["transeq", self.nt_path, self.aa_path, "-frame", "6", "-clean"]
        result = subprocess.run(command, check=True)
        if result.returncode != 0:
            raise RuntimeError("Input FASTA {fasta_to_screen} could not be translated:\n{result.stderr}")
        
    def parse_query_data(self):
        """
        Populate a list of query names, and associated sequences, for json formatting purposes.
        """
        with open(self.nt_path, "r", encoding = "utf-8") as fasta_file:
                self.raw : list[SeqRecord] = list(SeqIO.parse(fasta_file, "fasta"))

    def get_cleaned_fasta(self, out_prefix):
        """
        Return a FASTA where whitespace (including non-breaking spaces) and 
        illegal characters are replaced with underscores.
        """
        cleaned_file = f"{out_prefix}.cleaned.fasta"
        with (
            open(self.input_fasta_path, "r", encoding="utf-8") as fin,
            open(cleaned_file, "w", encoding="utf-8") as fout,
        ):
            for line in fin:
                line = line.strip()
                modified_line = "".join(
                    "_" if c.isspace() or c == "\xc2\xa0" or c == "#" else c
                    for c in line
                )
                fout.write(f"{modified_line}{os.linesep}")
        return cleaned_file
    
    def get_non_coding_regions(self) -> str:
        """ 
        Return the concatenation of all non-coding regions as a string,
        to be appended to a non_coding fasta file.
        """
        output : str = ""
        for start, end in self.non_coding_regions:
            output += self.raw[start-1:end]
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
