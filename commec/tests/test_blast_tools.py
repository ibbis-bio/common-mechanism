from io import StringIO
import pytest
import textwrap
from unittest.mock import patch
import numpy as np
import pandas as pd
from commec.tools.blast_tools import _split_by_tax_id, read_blast, _get_lineages, get_taxonomic_labels


@pytest.fixture
def blast_df():
    """
    Return a dataframe containing 3 BLAST hits, 2 with multiple taxids, 1 of which is invalid and 1
    of which is a synthetic taxid
    """
    blast_to_parse = textwrap.dedent(
        """\
        # BLASTX 2.15.0+
        # Query: NC_TEST
        # Database: /root/commec-dbs/mock
        #query acc.	subject title	subject acc.	subject tax ids	evalue	bit score	% identity	query length	q. start	q. end	subject length	s. start	s. end
        # 3 hits found
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


def test_split_by_tax_id(blast_df: pd.DataFrame):
    assert len(blast_df) == 3
    split_blast = _split_by_tax_id(blast_df)
    assert len(split_blast) == 5
    expected_tax_ids = {2371, 644357, 10760, 110011001100, 32630}
    assert set(split_blast["subject tax ids"]) == expected_tax_ids


@patch("pytaxonkit.lineage")
def test_get_lineages(mock_lineage, blast_df, lineage_df):
    mock_lineage.return_value = lineage_df
    blast_df = _split_by_tax_id(blast_df)
    lin = _get_lineages(
        blast_df["subject tax ids"], "commec-dbs/taxonomy/", 8
    )
    # Expect the invalid taxid to be filtered out
    expected_tax_ids = {2371, 644357, 10760, 32630}
    assert set(lin["TaxID"]) == expected_tax_ids


@patch("pytaxonkit.lineage")
def test_taxdist(mock_lineage, blast_df, lineage_df):
    mock_lineage.return_value = lineage_df
    # Fake values - should find 1 regulated hit after filtering
    reg_taxids = ['644357', '10760']
    vax_taxids = ['10760']
    reg_df = get_taxonomic_labels(
        blast_df, reg_taxids, vax_taxids, "commec-dbs/taxonomy/", 8
    )
    # Expect the synthetic taxid to be filtered out
    expected_tax_ids = {2371, 644357, 10760}
    assert set(reg_df["subject tax ids"]) == expected_tax_ids

    # Expect only taxid 644357 to be marked as "regulated"
    assert reg_df[reg_df["subject tax ids"] == 2371]["regulated"].iloc[0] == False
    assert reg_df[reg_df["subject tax ids"] == 644357]["regulated"].iloc[0] == True
    assert reg_df[reg_df["subject tax ids"] == 10760]["regulated"].iloc[0] == False
