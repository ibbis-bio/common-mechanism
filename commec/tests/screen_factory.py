"""
Helper module to quickly create a test run for a full commec screen call.
Simply create a ScreenFactory object,
call add_query() to add a query to the screen.
call add_hits() to add hits to queries.
Call run() to get the ScreenResult object from the test.

This abstracts the difficulty in dealing with external database files,
as well as dealing with the various regulated annotations etc.
"""

import os
import math
import json
from dataclasses import asdict
from unittest.mock import patch

import pandas as pd

from commec.screen import run, ScreenArgumentParser, add_args
from commec.config.result import ScreenResult, ScreenStep
from commec.config.json_io import get_screen_data_from_json

def skip_taxonomy_info(
    blast: pd.DataFrame,
    _regulated_taxids: list[str],
    _vaccine_taxids: list[str],
    _db_path: str | os.PathLike,
    _threads: int,
):
    """
    Override the taxonomy database retrieval with our own information at a 
    per-hit basis.
    """
    blast = blast.merge(TAXONOMY, on = "subject acc.", how = "left")
    return blast

def skip_biorisk_annotations(_input_file):
    """
    Overrides the biorisk annotations retrieval to be
    replaced with our own generated annotation data.
    """
    return BIORISK_ANNOTATIONS_DATA

DATABASE_DIRECTORY = os.path.join(os.path.dirname(__file__), "test_dbs/")
TAXONOMY = pd.DataFrame(columns=["subject acc.","regulated", "superkingdom", "phylum", "genus", "species"])
BIORISK_ANNOTATIONS_DATA = pd.DataFrame(columns=["ID", "Description", "Must flag"])

class ScreenTesterFactory:
    """
    Factory class, to easily create the data for a commec screen run,
    so we can more easily construct various testing scenarios without bloated
    code strings.
    Simply create an instance, and call 
        add_query, : Add queries to the screen
        add_hit,    : Add hits to the screen
        add_regulated_taxid, : Include informatiokn as to what taxids are regulated.
    
    and then run(), which will return the ScreenResult object.

    """
    def __init__(self, test_name, tmp_path):
        # Containers to hold input hits.
        self.name = test_name
        self.tmp_path = tmp_path
        self.tx_reg_info = []
        self.queries = {}
        self.biorisks = []
        self.protein_tx = []
        self.nucl_tx = []
        self.lowconcern_protein = []
        self.lowconcern_dna = []
        self.lowconcern_rna = []
        self.input_fasta_path = ""

        # We reset the globals when a new test is made - its just safer.
        TAXONOMY = pd.DataFrame(columns=["subject acc.","regulated", "superkingdom", "phylum", "genus", "species"])
        BIORISK_ANNOTATIONS_DATA = pd.DataFrame(columns=["ID", "Description", "Must flag"])

    def run(self):
        """
        Run the Factory commec screen. Return the ScreenResult Object.
        """

        self._create_temporary_files()

        print("Using the following Taxonomy Information:\n", TAXONOMY.to_string())
        # We patch taxonomic labels to avoid having to make a mini-taxonomy database.
        # We also patch in the desired CLI arguments to avoid an input yaml, and control output.
        with (patch("commec.screeners.check_reg_path.get_taxonomic_labels", new=skip_taxonomy_info), patch(
                "commec.screeners.check_biorisk.read_biorisk_annotations", new=skip_biorisk_annotations), patch(
            "sys.argv",
            [
                "test.py", str(self.input_fasta_path), 
                "-d", str(DATABASE_DIRECTORY), 
                "-o", str(self.tmp_path), 
                "--resume",
                "--verbose"
            ],
        )):
            parser = ScreenArgumentParser()
            add_args(parser)
            args = parser.parse_args()
            run(args)
        
        # return the screen data output for checking:
        json_output_path = self.tmp_path / f"{self.name}.output.json"
        assert os.path.isfile(json_output_path)
        actual_screen_result : ScreenResult = get_screen_data_from_json(json_output_path)

        print(f"Raw {self.name} test output: ")
        print(json.dumps(asdict(actual_screen_result), indent=2))

        return actual_screen_result

    def add_query(self, name, size):
        self.queries[name] = "a"*size

    def add_hit(self,
        to_step,                    # Which step this hit belongs to...
        to_query,                   # Which query this hit belongs to...
        start,                      # Start Nucleotide in Query coords.
        stop,                       # End Nucleotide in Query Coords...
        title = "untitled",
        accession = "",
        taxid = 0,
        species = "unclassified",
        genus = "unclassified",
        superkingdom = "unclassified",
        phylum = "unclassified",
        regulated = False,          # Used to update taxonomy, and also biorisk Must flag.
        score = 1000,               # Used for HMMSCAN
        description = "no description",
        ):

        assert start > 0
        assert stop > 0

        # Ensure that the taxonomy LUT entry for this hit exists.
        if to_step in [ScreenStep.TAXONOMY_AA, ScreenStep.TAXONOMY_NT]:
            TAXONOMY.loc[len(TAXONOMY)] = [
                accession,
                regulated,
                superkingdom,
                phylum,
                genus,
                species
            ]

        query_length = len(self.queries[to_query])
        query_length_aa = math.floor(query_length/3)
        # Calculate frame if appropriate:
        frame = (start - 1) % 3 + 1
        if start > stop:
            # Frame must be reversed.
            frame += 3

            # Flip the start and stop to account for reverse frames.
            start = query_length - start
            stop = query_length - stop

        print(f"Added Hit using frame {frame}")

        start_aa = math.floor(start / 3)
        end_aa = math.floor(stop / 3)
        
        # Used if HMMSCAN i.e. protein based.
        query_name = f"{to_query}_{str(frame)}"

        length = abs(stop - start)
        length_aa = abs(end_aa - start_aa)
        
        if to_step == ScreenStep.BIORISK:
            BIORISK_ANNOTATIONS_DATA.loc[len(BIORISK_ANNOTATIONS_DATA)] = [
                title, description, regulated]
            print(f"Added Biorisk Annotation: {title},{description},{regulated}")
            self.biorisks.append(
                f"{title}    {accession}    {length_aa}    {query_name}    999    {query_length_aa}  0.0   {score}    10.0    1   1   0    0    {score}  10.0  1   {length_aa}  {start_aa}  {end_aa}  1  {length_aa}  1.00  {description}"
            )
            return

        if to_step == ScreenStep.TAXONOMY_AA:
            self.protein_tx.append(
                f"{to_query}\t{title}\t{accession}\t{taxid}\t0.0\tBITSCORE\t99.999\t{query_length}\t{start}\t{stop}\t{length}\t1\t{length}"
            )
            return
        
        if to_step == ScreenStep.TAXONOMY_NT:
            self.nucl_tx.append(
                f"{to_query}\t{title}\t{accession}\t{taxid}\t0.0\tBITSCORE\t99.999\t{query_length}\t{start}\t{stop}\t{length}\t1\t{length}"
            )
            return

        if to_step == ScreenStep.LOW_CONCERN_PROTEIN:
            self.lowconcern_protein.append(
                f"{title}    {accession}    {length_aa}    {query_name}    999    {query_length_aa}  0.0   {score}    10.0    1   1   0    0    {score}  10.0  1   {length_aa}  {start_aa}  {end_aa}  1  {length_aa}  1.00  {description}"
            )
            return

        if to_step == ScreenStep.LOW_CONCERN_RNA:
            self.lowconcern_rna.append(
                f"{title}\t{accession}\t{to_query}\tQA999\t{length}\t1\t{length}\t{start}\t{stop}\tSTRAND\tTRUNC\tPASS\tGC\t10\t{score}\t0.0\t100\t{description}"
            )
            #RNA, mdl = target, seq = query
            #"""\
            #target name         accession query name                accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value  inc description of target
            #------------------- --------- ------------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ ---------  --- ---------------------
            #BENIGNRNA            12346     FCTEST1	                 Q1         50	    100      200       50      150 STRAND TRUNC PASS   GC    10   1000       0.0  100    BenignCMTestOutput
            #"""
            return

        if to_step == ScreenStep.LOW_CONCERN_DNA:
            self.lowconcern_dna.append(
                f"{to_query}\t{title}\t{accession}\t{taxid}\t0.0\tBITSCORE\t99.999\t{query_length}\t{start}\t{stop}\t{length}\t1\t{length}"
            )
            return


    def _create_temporary_files(self):
        
        # Make the paths.
        os.mkdir(self.tmp_path / f"output_{self.name}")
        os.mkdir(self.tmp_path / f"input_{self.name}")
        
        # Create input fasta:
        self.input_fasta_path = self.tmp_path / f"{self.name}.fasta"
        print("Using fasta file: " + str(self.input_fasta_path))

        for key, value in self.queries.items():
            self.input_fasta_path.write_text(f">{key}\n{value}\n")

        # --RESUME FILES::
        # BIORISK FILES
        header = "#tname    accession  tlen qname        accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target\n"
        biorisk_db_output_path = self.tmp_path / f"output_{self.name}/{self.name}.biorisk.hmmscan"
        hmmer_biorisk_to_parse = header + "\n".join(self.biorisks)
        biorisk_db_output_path.write_text(hmmer_biorisk_to_parse)
        print("writing biorisk hmm: \n", hmmer_biorisk_to_parse)

        # TAXONOMY NR FILES
        header = "#query acc.	    title   subject acc.    taxid	evalue	bit score	% identity	    q.len	q.start	q.end	s.len	s. start	s. end\n"
        nr_db_output_path = self.tmp_path / f"output_{self.name}/{self.name}.nr.blastx"
        blastnr_to_parse = header + "\n".join(self.protein_tx)
        nr_db_output_path.write_text(blastnr_to_parse)
        print("writing blast nr: \n", blastnr_to_parse)

        # TAXONOMY NT FILES:
        nt_db_output_path = self.tmp_path / f"output_{self.name}/{self.name}.nt.blastn"
        blastnt_to_parse = header + "\n".join(self.nucl_tx)
        nt_db_output_path.write_text(blastnt_to_parse)
        print("writing blast nt: \n", blastnt_to_parse)

        # LOW CONCERN FILES:
        header = " #tname    accession        tlen qname        accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias    from    to  from    to  from    to  acc description of target\n"
        low_concern_hmm_output_path = self.tmp_path / f"output_{self.name}/{self.name}.low_concern.hmmscan"
        low_concern_hmmscan_to_parse = header + "\n".join(self.lowconcern_protein)
        low_concern_hmm_output_path.write_text(low_concern_hmmscan_to_parse)
        print("writing lowconcern hmm: \n", low_concern_hmmscan_to_parse)

        header = "#target name         accession query name                accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value  inc description of target\n"
        low_concern_cmscan_output_path = self.tmp_path / f"output_{self.name}/{self.name}.low_concern.cmscan"
        low_concern_cmscan_to_parse = header + "\n".join(self.lowconcern_rna)
        low_concern_cmscan_output_path.write_text(low_concern_cmscan_to_parse)
        print("writing lowconcern rna: \n", low_concern_cmscan_to_parse)

        header = "#query acc.	title	subject acc.taxid	evalue	bit score	% identity	    q.len	q.start	q.end	    s.len	s. start	s. end\n"
        low_concern_nt_output_path = self.tmp_path / f"output_{self.name}/{self.name}.low_concern.blastn"
        low_concern_blastnt_to_parse = header + "\n".join(self.lowconcern_dna)
        low_concern_nt_output_path.write_text(low_concern_blastnt_to_parse)
        print("writing lowconcern dna: \n", low_concern_blastnt_to_parse)

