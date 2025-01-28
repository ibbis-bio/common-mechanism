import pandas as pd
import pytest
from pandas.testing import assert_frame_equal
from commec.utils.coordinates import convert_aa_to_nt_coordinates


def test_convert_aa_to_nt_coords_dataframe():
    """
    Consider the 22nt sequence atgtgccatggatgtgccatgg. We get these translations:

        Frame   Pos     nt with codon split                 Translation
        1       0       atg tgc cat gga tgt gcc atg g       MCHGCAM
        2       1       a tgt gcc atg gat gtg cca tgg       CAMDVPW
        3       2       at gtg cca tgg atg tgc cat gg       VPWMCH
        4       -0      c cat ggc aca tcc atg gca cat       HGTSMAH
        5       -1      cca tgg cac atc cat ggc aca t       PWHIHGT       
        6       -2      cc atg gca cat cca tgg cac at       MAHPWH

    Then consider 6 alignments, one to each of the reading frames, including the
    coordinates given for the Alignment (a.) and Query (q.), noting that all
    coordinates are 1-indexed and *inclusive* (be careful with splicing!):

        Frame   Match   nt start    nt end  a. start   a. end   q. start    q.end  
        1       CHGC    tgc         tgt     2          5        4           15
        2       CAMDV   tgt         gtg     1          5        2           16
        3       CH      tgc         cat     5          6        15          20
        4       HGTSM   cat         atg     1          5        7           21
        5       IHGT    atc         ata     4          7        2           13
        6       AHP     gca         cat     2          4        9           17

    Query coordinates are w.r.t. to the original query sequence, not its reverse complement.
    """
    input_coords = pd.DataFrame(
        {
            "query name": ["F1", "F2", "F3", "R1", "R2", "R3"],
            "frame": [1, 2, 3, 4, 5, 6],
            "ali from": [2, 1, 5, 1, 4, 2],
            "ali to": [5, 5, 6, 5, 7, 4],
            "qlen": [22, 22, 22, 22, 22, 22],
        }
    )

    nt_start_expected = pd.Series([4, 2, 15, 7, 2, 9])
    nt_end_expected = pd.Series([15, 16, 20, 21, 13, 17])

    nt_start, nt_end = convert_aa_to_nt_coordinates(
        input_coords["frame"],
        input_coords["ali from"],
        input_coords["ali to"],
        input_coords["qlen"],
    )
    assert_frame_equal(nt_start_expected, nt_start)
    assert_frame_equal(nt_end_expected, nt_end)


def test_convert_aa_to_nt_coords_ints():
    nt_start, nt_end = convert_aa_to_nt_coordinates(
        frame=5, aa_start=4, aa_end=7, seq_length=22
    )
    assert nt_start == 1
    assert nt_end == 12
