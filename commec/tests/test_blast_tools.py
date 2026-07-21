from io import StringIO
import pytest
import textwrap
import numpy as np
import pandas as pd
from commec.tools.blast_tools import read_blast


@pytest.fixture
def blast_df():
    """
    Return a dataframe containing 3 BLAST hits, 2 with multiple taxids, 1 of which is invalid and 1
    of which is a synthetic taxid
    """
    # outfmt 6 tabular output: data rows only, no comment/header lines.
    blast_to_parse = textwrap.dedent(
        """\
        BT_01	SUBJECT	SUBJECT_ACC	2371;644357	0.0	BITSCORE	99.999	300	101	200	500	1	100
        BT_01	SUBJECT	SUBJECT_ACC	10760;110011001100	0.0	BITSCORE	99.999	300	25	80	500	1	100
        BT_01	SUBJECT	SUBJECT_ACC	32630	0.0	BITSCORE	99.999	300	275	300	500	1	100
        """
    )
    return read_blast(StringIO(blast_to_parse))


@pytest.fixture
def lineage_df():
    """
    Dataframe subsetting columns from the results of pytaxonkit.lineage applied to blast_df
    """
    return pd.DataFrame(
        {
            "TaxID": [2371, 644357, 10760, 110011001100, 32630],
            "Code": [2371, 644357, 10760, -1, 32630],
            "FullLineage": [
                "cellular organisms;Bacteria;Pseudomonadota;Gammaproteobacteria;Lysobacterales;Lysobacteraceae;Xylella;Xylella fastidiosa",
                "cellular organisms;Bacteria;Pseudomonadota;Gammaproteobacteria;Lysobacterales;Lysobacteraceae;Xylella;Xylella fastidiosa;Xylella fastidiosa subsp. multiplex",
                "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Autographiviridae;Studiervirinae;Teseptimavirus;Teseptimavirus T7;Escherichia phage T7",
                np.nan,
                "other entries;other sequences;artificial sequences;synthetic construct",
            ],
            "FullLineageTaxIDs": [
                "131567;2;1224;1236;135614;32033;2370;2371",
                "131567;2;1224;1236;135614;32033;2370;2371;644357",
                "10239;2731341;2731360;2731618;2731619;2731643;2731653;110456;1985738;10760",
                np.nan,
                "2787854;28384;81077;32630",
            ],
            "FullLineageRanks": [
                "no rank;superkingdom;phylum;class;order;family;genus;species",
                "no rank;superkingdom;phylum;class;order;family;genus;species;subspecies",
                "superkingdom;clade;kingdom;phylum;class;family;subfamily;genus;species;no rank",
                np.nan,
                "no rank;no rank;no rank;species",
            ],
        }
    )
