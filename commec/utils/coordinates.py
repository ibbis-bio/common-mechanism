"""
Helper functions associated with the handling of
basepair|amino-acid, nucleotide|peptide coordinate systems.
"""

import numpy as np
import numpy.typing as npt
import pandas as pd


def _convert_single_aa_coordinate_to_nt(
    frame: int, aa_start: int, aa_end: int, seq_length: int
) -> tuple[int, int]:
    """
    Convert protein coordinates to nucleotide coordinates, considering the reading frame.

    Parameters:
        frame: Reading frame (1-6), numbered following the same naming convention as transeq:
            * 1, 2, 3: Forward frames starting at positions 0, 1, 2
            * 4, 5, 6: Reverse frames, starting at positions 0, -1, -2
        protein_start: Start position in protein coordinates, counting from 1.
        protein_end: End position in protein coordinates, counting from 1.
        seq_length: Length of the original nucleotide sequence.

    Returns:
        tuple: (nucleotide_start, nucleotide_end)

    A worked example of the conversion between coordinates:

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
    frame_is_forward = frame <= 3

    # We assume that codons are split the same way in the forward and reverse frames, so offsets
    # are calculated from the start in the forward frames, and from the end in the reverse frames.
    #
    # For example: Frame 1 = atg tgc cat gg, offset 0 from the start
    #              Frame 4 = cc atg gca cat, offset 2 from the end
    frame_offset = frame - 1 if frame_is_forward else (seq_length - (frame - 4)) % 3

    if frame_is_forward:
        # Count from 1 to the start of the aa_start codon (hence aa_start - 1)
        nt_start = 1 + frame_offset + (aa_start - 1) * 3
        # Count from starting offset to the end of the aa_end codon
        nt_end = frame_offset + (aa_end * 3)
    else:
        # Count from length to the end of the aa_end codon
        # (the *end* of the reversed codon = the start in nt)
        nt_start = 1 + seq_length - frame_offset - (aa_end * 3)
        # Count from the end of the sequence to the start of the aa_start codon
        nt_end = seq_length - frame_offset - (aa_start - 1) * 3

    return nt_start, nt_end


_vectorized_convert_aa_coordinate_to_nt = np.vectorize(
    _convert_single_aa_coordinate_to_nt
)


def convert_aa_to_nt_coordinates(
    frame: int | pd.Series,
    aa_start: int | pd.Series,
    aa_end: int | pd.Series,
    seq_length: int | pd.Series,
) -> tuple[int | npt.NDArray[np.int64], int | npt.NDArray[np.int64]]:
    """
    Convert protein coordinates to nucleotide coordinates considering the reading frame.
    The nucleotide coordinates are all 1-indexed, inclusive, and w.r.t. to the original
    query sequence, not its reverse complement.

    Parameters:
        frame: Reading frame (1-6), numbered following the same naming convention as transeq:
            * 1, 2, 3: Forward frames starting at positions 0, 1, 2
            * 4, 5, 6: Reverse frames, starting at positions 0, -1, -2
        protein_start: Start position in protein coordinates, counting from 1.
        protein_end: End position in protein coordinates, counting from 1.
        seq_length: Length of the original nucleotide sequence.

    Returns:
        tuple: (nucleotide_start, nucleotide_end)
    """
    if isinstance(frame, (int, np.integer)):
        return _convert_single_aa_coordinate_to_nt(frame, aa_start, aa_end, seq_length)

    frame = np.asarray(frame)
    aa_start = np.asarray(aa_start)
    aa_end = np.asarray(aa_end)
    seq_length = np.asarray(seq_length)

    return _vectorized_convert_aa_coordinate_to_nt(frame, aa_start, aa_end, seq_length)
