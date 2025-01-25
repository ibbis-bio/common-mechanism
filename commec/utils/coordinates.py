"""
Helper functions associated with the handling of 
basepair|amino-acid, nucleotide|peptide coordinate systems.
"""

import numpy as np
import pandas as pd

def convert_protein_to_nucleotide_coords(frame: int | pd.Series,
                                         aa_start: int | pd.Series,
                                         aa_end: int | pd.Series,
                                         seq_length: int | pd.Series):
    """
    Convert protein coordinates to nucleotide coordinates considering the reading frame.

    Parameters:
        frame: Reading frame (1-6), numbered following the same naming convention as transeq:
            * 1, 2, 3: Forward frames starting at positions 0, 1, 2
            * 4, 5, 6: Reverse frames, starting at positions 0, -1, -2
        protein_start (int, or [int]): Start position in protein coordinates, counting from 1.
        protein_end (int, or [int]): End position in protein coordinates, counting from 1.
        seq_length (int): Length of the original nucleotide sequence.

    Returns:
        tuple: (nucleotide_start, nucleotide_end)

    A worked example of the conversion between coordinates:

    Consider the 22nt sequence atgtgccatggatgtgccatgg. We get these translations:

        Frame   Pos     Codon split                         Translation
        1       0       atg tgc cat gga tgt gcc atg g       MCHGCAM
        2       1       a tgt gcc atg gat gtg cca tgg       CAMDVPW
        3       2       at gtg cca tgg atg tgc cat gg       VPWMCH
        4       -0      c cat ggc aca tcc atg gca cat       HGTSMAH
        5       -1      cca tgg cac atc cat ggc aca t       PWHIHGT       
        6       -2      cc atg gca cat cca tgg cac at       MAHPWH

    Then consider 6 alignments, one to each of the reading frames, including the
    coordinates given for the Alignment (a.) and Query (q.), noting that query
    coordinates are 0-indexed, and alignment coordinates are 1-indexed, and all
    are *inclusive* (which means be careful with splicing!):

        Frame   Match   a. start   a. end   q. start    q.end  
        1       CHGC    2          5        3           15
        2       CAMDV   1          5        1           15
        3       CH      5          6        14          19
        4       HGTSM   1          5        6           20
        5       IHGT    4          7        1           12
        6       AHP     2          4        8           16

    """
    aa_start = np.asarray(aa_start)
    aa_end = np.asarray(aa_end)
    frame = np.asarray(frame)
    seq_length = np.asarray(seq_length)

    frame_is_forward = frame <= 3

    # We assume that codons are split the same way in the forward and reverse frames, so offsets
    # are calculated from the start in the forward frames, and from the end in the reverse frames.
    #
    # For example: Frame 1 = atg tgc cat gg, offset 0 from the start
    #              Frame 4 = cc atg gca cat, offset 2 from the end
    frame_offset = np.where(frame_is_forward, frame - 1, (seq_length - (frame - 4)) % 3)

    # For the starting nt coordinate:
    #   Forward: count from the start to aa_start, accounting for offset
    #            subtract 1 since we want the *start* of the first codon
    #   Reverse: count from end to aa_end, accounting for offset
    #            don't subtract since the *end* of the reversed codon = the start in nt
    nt_start = np.where(frame_is_forward,
                        frame_offset + (aa_start - 1) * 3,
                        seq_length - frame_offset - (aa_end * 3))

    # For the ending nt coordinate:
    #   Forward: start to aa_end, accounting for offset
    #   Reverse: end to aa_start, accounting for offset
    nt_end = np.where(frame_is_forward,
                      frame_offset + (aa_end - 1) * 3,
                      seq_length - frame_offset - (aa_start * 3))
   
    # Convert to back to 1-based coordinates for reporting.
    nt_start += 1
    nt_end += 1

    return nt_start, nt_end
