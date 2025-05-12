"""
Module for communicating between Commec and FoldSeek APIs.
Responsible for the generation of 3Di files from a Aminoacid .fasta.
As well as acting as a Commec Tool for calling foldseek
"""
import os
import logging
from transformers import T5ForConditionalGeneration, T5Tokenizer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from commec.tools.search_handler import SearchHandler, SearchToolVersion

logger = logging.getLogger(__name__)

class FoldseekHandler(SearchHandler):
    """
    Database handler for FoldSeek.
    """

    def search(self):
        command = [
            "foldseek",
            "search",
            self.input_file,
            self.db_file,
            self.out_file,
            self.out_file + "_tmp/"
        ]
        self.run_as_subprocess(command, self.temp_log_file)
    
    def read_output(self):
        ...

    def get_version_information(self):
        return SearchToolVersion("FOLDSEEK", "FOLDSEEK")
    
def _readfoldseek(filename : str | os.PathLike):
    ...

def _generate_3di(model : T5ForConditionalGeneration,
                  tokenizer : T5Tokenizer,
                  aa_sequence : str) -> str:
    input_text = f"translate amino to 3di: {aa_sequence}"
    inputs = tokenizer(input_text, return_tensors = "pt")
    outputs = model.generate(**inputs)
    return tokenizer.decode(outputs[0], skip_special_tokens=True)

def convert_aa_to_3di(fasta_aa : str | os.PathLike, output_3di : str | os.PathLike):
    """
    Takes an input fasta file, containing aminoacid sequences, and converts it using
    the ProstT5 model into 3Di sequences, which represent the structure of the query.
    Outputs are written into a new fasta file.
    """
    logger.info("Generating Structural 3Di data ...")
    try:
        model_name = "Rostlab/ProstT5"
        model = T5ForConditionalGeneration.from_pretrained(model_name)
        tokenizer = T5Tokenizer.from_pretrained(model_name)

        # Read in the input fasta.
        records : list[SeqRecord] = []
        with open(fasta_aa, "r", encoding = "utf-8") as fasta_file:
            records = list(SeqIO.parse(fasta_file, "fasta"))

        for record in records:
            aa_seq = record.seq
            threedi_seq = _generate_3di(model, tokenizer,aa_seq)
            record.seq = threedi_seq
        
        with open(output_3di, "w", encoding = "utf-8") as fasta_output:
            SeqIO.write(records, fasta_output, "fasta")

    except OSError as e:
        logger.error("Error generating 3Di data (likely due to a lack of RAM): %s", str(e))