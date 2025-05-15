"""
Module for communicating between Commec and FoldSeek APIs.
Responsible for the generation of 3Di files from a Aminoacid .fasta.
As well as acting as a Commec Tool for calling foldseek
"""
import os
import logging
from requests import get, post
import sys
from time import sleep
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from commec.tools.search_handler import SearchHandler, SearchToolVersion

logger = logging.getLogger(__name__)

class FoldseekHandler(SearchHandler):
    """
    Database handler for FoldSeek.
    """

    def _search(self):
        """
        Foldseek provides an API for dealing with remote request to their servers.
        Here is an example of an API request:
        
        curl -X POST -F q=@PATH_TO_FILE -F 'mode=3diaa' -F 'database[]=BFVD' -F 'database[]=afdb50' -F 'database[]=afdb-swissprot' -F 'database[]=afdb-proteome' -F 'database[]=bfmd' -F 'database[]=cath50' -F 'database[]=mgnify_esm30' -F 'database[]=pdb100' -F 'database[]=gmgcl_id' https://search.foldseek.com/api/ticket

        This provides a ticket number, which is also used for further API status checking, as well as requesting the download.
        """
        api = "https://search.foldseek.com/api/"
        #api = "https://search.mmseqs.com/api/"
        api_ticket = api + "ticket"
        api_result = api + "result"

        response = ""
        data = ""

        #'BFVD', 'afdb50', 'afdb-swissprot', 'afdb-proteome',
        #'bfmd', 'cath50', 'mgnify_esm30', 'pdb100', 'gmgcl_id'

        clean_seq = "MPKIIEAIYENGVFKPLQKVDLKEGE"

        print(f"Requesting 3Di prediction for sequence of length {len(clean_seq)}...")
        url = f"https://3di.foldseek.com/predict/{clean_seq}"
        response = get(url)
        print(response)
        print(response.json())

        response = post(
            api_ticket,
            files={'q': open("input.fasta", 'rb')},
            data={
                'mode': '3diaa',
                'database[]': ['pdb100'],
                'email':'tomichaelbarnett@gmail.com'},
            #timeout=1000,
        )
                
        print(data)
        
        print(response)
        ticket_data = response.json()
        print("Ticket Data:")
        print(ticket_data)
        ticket = ticket_data['id']
        #progress_url = f"https://search.foldseek.com/api/status?ticket={ticket}"
        #status = get('https://search.mmseqs.com/api/ticket/' + ticket['id']).json()
        while True:
            #progress_data = get(progress_url)
            #print(progress_data)
            progress_data = get(api_ticket + "/" + ticket)
            progress = progress_data.json()
            print(progress)
            if progress['status'] == 'COMPLETE':
                break
            elif progress['status'] == 'ERROR':
                print("Response :\n", progress_data)
                raise RuntimeError("Foldseek job failed.")
            print("waiting...")
            sleep(5)

        result_url = f"https://search.foldseek.com/api/result/ticket={ticket}"
        result = get(result_url)

        with open(self.out_file, "wb") as f:
            f.write(result.content)
                
    
    def read_output(self):
        ...

    def get_version_information(self):
        return SearchToolVersion("FOLDSEEK", "FOLDSEEK")


def _readfoldseek(filename : str | os.PathLike) -> pd.DataFrame:
    ...