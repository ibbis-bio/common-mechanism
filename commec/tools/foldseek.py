"""
Module for communicating between Commec and FoldSeek APIs.
Responsible for the generation of 3Di files from a Aminoacid .fasta.
As well as sending those 3Di fasta files through the Foldseek pipeline,
and returning, and parsing the appropriate data.
"""
import os
import logging
from enum import StrEnum
from time import sleep
import gzip
from io import StringIO, BytesIO

import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

from commec.tools.search_handler import (
    SearchHandler,
    SearchToolVersion,
    DatabaseValidationError
)

logger = logging.getLogger(__name__)

class FoldseekStatus(StrEnum):
    """ 
    Possible values for returned foldseek job status
    """
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETE = "COMPLETE"
    ERROR = "ERROR"

class FoldseekHandler(SearchHandler):
    """
    Database handler for FoldSeek.
    """
    api = "https://search.foldseek.com/api/"
    api_ticket = api + "ticket"
    api_result = api + "result/download/"

    def __init__(
        self, database_file: str, input_file: str, out_file: str, **kwargs,
    ):
        super().__init__("", input_file, out_file, **kwargs)
        # We fill this with defaults, however they can always be overridden before screening.
        self.arguments_dictionary = {
                'mode': '3diaa',
                'database[]': ['pdb100'],
            }

    def _validate_db(self):
        """
        Validates that the database directory and file exists. Called on init.

        Foldseek is based online, so we simply return whether we can http the site.
        """
        if not self._check_foldseek_api():
            raise DatabaseValidationError(
                "Could not access Foldseek Webserver."
                " Check your internet connection."
            )

    def _check_foldseek_api(self, timeout=5):
        try:
            url = self.api + "databases"
            logger.debug("Sending request to %s", url)
            response = requests.get(url, timeout=timeout)
            print(response.status_code)
            return response.status_code == 200
        except requests.RequestException:
            return False

    def _write_3di_fasta(self, output_fasta : str | os.PathLike , record : SeqRecord):
        """
        Takes a SeqRecord, and creates a fasta file containing it, and its seq_3DI sequence.
        SeqRecord must have the seq_3di attribute, and returns and stores the output fasta filename as 
        a fasta_file attribute.
        """
        assert hasattr(record, "seq_3di"), "SeqRecord is missing the '_3di' attribute"
        record.fasta_file = output_fasta
        # Write 3DI single fasta file
        with open(record.fasta_file, 'w', encoding = "utf-8") as write_file:
            # Note, Foldseek is VERY finicky about fasta format correctness.
            write_file.writelines([
                ">"+record.name+"\n",
                str(record.seq)+"\n",
                ">3DI"+"\n",
                str(record.seq_3di+"\n")
            ])
        return output_fasta

    def _create_ticket(self, record : SeqRecord) -> int:
        """
        Takes a record, with seq_3di attribute, and creates a foldseek ticket.
        The ticket is stored as an additional attribute within the SeqRecord.
        The return code of the request is returned.
        Furthermore, a downloaded boolean attribute is added to the record for future use.
        """
        assert hasattr(record, "fasta_file"), "SeqRecord is missing expected 3di containing \'fasta_file\' attribute"
        
        #'BFVD', 'afdb50', 'afdb-swissprot', 'afdb-proteome',
        #'bfmd', 'cath50', 'mgnify_esm30', 'pdb100', 'gmgcl_id'

        # Request Foldseek job
        response = requests.post(
            self.api_ticket,
            files={'q': open(record.fasta_file, 'rb')},
            data=self.arguments_dictionary,
            timeout=5,
        )

        logger.debug("Foldseek ticket request status: %s", response.status_code)

        if response.status_code == 200:
            record.foldseek_ticket = response.json()['id']
            record.downloaded = False
            logger.debug("Created Foldseek ticket: %s", record.foldseek_ticket)
            return 200
        
        return response.status_code
    
    def _check_ticket(self, record : SeqRecord) -> bool:
        """
        Checks a ticket, returns true if we should be continuing . i.e:
        * We are completed, and downloaded.
        * Or we have Errored.
        """
        assert hasattr(record, "downloaded"), "SeqRecord is missing expected \'downloaded\' bool attribute"

        # No need to check anything if complete.
        if record.downloaded:
            return True

        progress = requests.get(self.api_ticket + "/" + record.foldseek_ticket, timeout = 5)
        if progress.status_code == 200:
            progress_status = progress.json()['status']
            logger.debug("Checking foldseek ticket %s for %s : %s",
                         record.foldseek_ticket, record.name, progress_status)
            match progress_status:
                case FoldseekStatus.PENDING:
                    return False
                case FoldseekStatus.RUNNING:
                    return False
                case FoldseekStatus.ERROR:
                    return True
                case FoldseekStatus.COMPLETE:
                    if record.downloaded:
                        return True
                    else:
                        self._download_ticket(record)
                        return False
                case _:
                    logger.debug("Unknown mapping for FoldseekStatus: %s", progress_status)
                    return False
        logger.debug("Bad response from Foldseek: %s", progress.status_code)
        return False

    def _download_ticket(self, record : SeqRecord):
        """
        Downloads the result of a foldseek job as a tsv, and saves the file for later use.
        """
        assert hasattr(record, "foldseek_ticket"), "SeqRecord is missing expected \'foldseek_ticket\' attribute"

        record.output_file = self.out_file + "_" + record.name + ".tsv"
        logger.debug("Downloading foldseek result for %s to %s", record.name, record.output_file)

        result_url = self.api_result + record.foldseek_ticket
        result = requests.get(result_url, timeout=3600)

        if result.status_code != 200:
            record.downloaded = False
            logger.error("A completed Foldseek result download failed.")
            return

        # Process the downloaded file into the more managable
        with gzip.open(BytesIO(result.content), 'rt', encoding='utf-8') as f:
            df = pd.read_csv(f, sep='\t')
        df.to_csv(record.output_file, sep="\t", index=False, encoding="utf-8")
        
        # Assume `result.content` is gzip-compressed binary data
        #with gzip.open(BytesIO(result.content), 'rt', encoding='utf-8') as gz_file:
        #    text = gz_file.read()

        # Now write the decompressed text to a regular UTF-8 file
        #with open(record.output_file, 'w', encoding='utf-8') as out_file:
        #    out_file.write(text)
            
        record.downloaded = True

    def _merge_output_data(self, records : list[SeqRecord]):
        """
        Each translation of each query represents a single output .tsv from foldseek.
        We need to merge the data for each query appropriately such that there is a single output file.
        """

    def _search(self):
        """
        Foldseek provides an API for dealing with remote request to their servers.
        Here is an example of an API request:
        
        curl -X POST -F q=@PATH_TO_FILE -F 'mode=3diaa' 
            -F 'database[]=BFVD' 
            -F 'database[]=afdb50' 
            -F 'database[]=afdb-swissprot' 
            -F 'database[]=afdb-proteome' 
            -F 'database[]=bfmd' 
            -F 'database[]=cath50' 
            -F 'database[]=mgnify_esm30' 
            -F 'database[]=pdb100' 
            -F 'database[]=gmgcl_id' 
            https://search.foldseek.com/api/ticket

        This provides a ticket number, which is also used for further API status checking, as well as requesting the download.
        However this only works for a fasta with 3di and aa, we need to generate 3di sequences first.
        This can be done by using another api call copied from the frontend of the webserver that mimics the 'predict' button.

        All requests to foldseek are single queries, thus we need to manage many request for each translation of the provided sequence.
        All the returned files will need to be stitched together, and processed per query.
        """

        records : list[SeqRecord] = []

        # Foldseek does not accept multi-fasta, so we need to treat queries individually.
        # Furthermore, we need to treat reading frames individually.
        # Thankfully, all of these can be separate requests.
        logger.debug("Starting FoldSeek Search Process...")

        logger.debug("Loading protein primary sequences...")
        with open(self.input_file, "r", encoding = "utf-8") as fasta_file:
            records = list(SeqIO.parse(fasta_file, "fasta"))

        # SeqRecords allow arbitrary attribute additions,
        # so we take advantage of this to append some useful information:
        # - seq_3di, the 3DI sequence
        # - fasta_file, the temprary single fasta containing original and 3di sequences.
        # - foldseek_ticket, the ticket ID needed to check results of the request.

        logger.debug("Processing %i primary sequences", len(records))
        for record in records:
            new_fasta_filepath = self.out_file + "_" + record.name + ".3da.fasta"
            if os.path.isfile(new_fasta_filepath) and not self.force:
                logger.warning("Existing 3DI fasta for query %s detected, which will be used. Disable this behaviour with --force", record.name)
                record.fasta_file = new_fasta_filepath
                continue
            record.seq_3di = convert_aa_to_3di(record.seq)
            self._write_3di_fasta(new_fasta_filepath, record)
            sleep(1.0) # Don't overload the servers too much.

        # Create the tickets from the fasta files.
        logger.debug("Creating tickets for Foldseek jobs...")
        for record in records:
            attempts = 5 # How many attempts should be make before giving up?
            wait_time = 1 # Increments with failed attempts to avoid webserver overload.
            while attempts > 0:
                status_code = self._create_ticket(record)
                if status_code == 429:
                    wait_time += 1
                    attempts -= 1
                    logger.debug("Creating ticket failed due to request overload, "
                                 "increasing wait-time to %i seconds. %i remaining attempts.",
                                 wait_time, attempts)
                sleep(wait_time)

                if status_code == 200:
                    break

        # Wait for all responses:
        logger.debug("Waiting for all responses...")
        while True:
            if all(self._check_ticket(record) for record in records):
                break
            sleep(1)
    
    def read_output(self):
        ...

    def get_version_information(self):
        return SearchToolVersion("FOLDSEEK", "FOLDSEEK")

def convert_aa_to_3di(sequence : str, timeout : int = 5, attempts_to_make : int = 3) -> str:
    """
    Foldseek uses an [internal 3DI api](https://github.com/soedinglab/MMseqs2-App/blob/b51c0821a68badaf330a4da0d640c9687bc9203a/frontend/PredictStructureButton.vue#L152) 
    which we can take advantage of, allowing the conversion of primary sequence AA to 3DI near instantly 
    so long as we have an internet connection.
    
    Returns the 3DI sequence for a given AA on success, or an empty string on failure.

    * `sequence` : Input primary sequence as string, consisting only of 20 AA or X.
    * `timeout = 5` : How long (in seconds) before a single request will timeout.
    * `attempts_to_make = 3` : How many times to fail before giving up on the request.
    """
    assert len(sequence) > 0, "Cannot convert a length of 0 primary sequence to 3DI."

    attempts = attempts_to_make

    while attempts:
        attempts -= 1
        api_3di = "https://3di.foldseek.com/predict/" + sequence
        logger.debug("Requesting 3DI sequence (size = %i aa) from foldseek server (attempt %i)...",
                     len(sequence), (attempts_to_make - attempts))
        logger.debug(api_3di, extra = {"no_prefix" : True, "cap":True})
        
        response = requests.get(api_3di, timeout=timeout)
        logger.debug("Response: %s", response)

        if response.status_code != 200:
            logger.debug("Bad response status code: %i", response.status_code)
            continue

        _3di = response.json()
        logger.debug("Successful 3DI request! (size = %i)", len(_3di))
        return _3di

    logger.debug("3DI request was unsuccessful.")
    return ""

def _readfoldseek(filename : str | os.PathLike) -> pd.DataFrame:
    ...