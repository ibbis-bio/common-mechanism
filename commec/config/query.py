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
    Query to screen.
    """
    def __init__(self, seq_record : SeqRecord):
        self.name = seq_record.id
        self.seq_record = seq_record

    def translate(self, input_path, output_path) -> None:
        """ Run command transeq, to translate our input sequences. """
        command = ["transeq", input_path, output_path, "-frame", "6", "-clean"]
        result = subprocess.run(command)
        if result.returncode != 0:
            raise RuntimeError(f"Input FASTA {input_path} could not be translated:\n{result.stderr}")

