"""
Helper functions associated with the handling of 
basepair|amino-acid, nucleotide|peptide coordinate systems.
"""

import numpy as np

def convert_protein_to_nucleotide_coords(frame,
                                         protein_start,
                                         protein_end,
                                         seq_length):
    """
    Convert protein coordinates to nucleotide coordinates considering the reading frame.

    Parameters:
    frame (int, or [int]): Reading frame (1-6)
        Frames 1-3: Forward frames starting at positions 0, 1, 2
        Frames 4-6: Reverse frames starting from the end at positions 0, 1, 2
    protein_start (int, or [int]): Start position in protein coordinates, counting from 1.
    protein_end (int, or [int]): End position in protein coordinates, counting from 1.
    seq_length (int): Length of the original sequence, mandatory for reverse frames (4,5,6) only.

    Returns:
    tuple: (nucleotide_start, nucleotide_end)
    """
    # Convert protein coordinates to 0-based, for calculation.
    protein_start = np.asarray(protein_start, dtype=np.int64) - 1
    protein_end = np.asarray(protein_end, dtype=np.int64) - 1
    frame = np.asarray(frame, dtype=np.int64)
    seq_length = np.asarray(seq_length, dtype=np.int64)

    # Reverse frame offsets, for when total length not divisible into codons.
    reverse_offset = seq_length % 3

    # Initialize arrays for nucleotide start and end
    nucleotide_start = np.zeros_like(protein_start, dtype=np.int64)
    nucleotide_end = np.zeros_like(protein_end, dtype=np.int64)

    # Forward frames (1, 2, 3)
    forward_mask = frame <= 3
    nucleotide_start[forward_mask] = (protein_start[forward_mask] * 3) + (frame[forward_mask] - 1)
    nucleotide_end[forward_mask] = (protein_end[forward_mask] * 3) + 2 + (frame[forward_mask] - 1)

    # Reverse frames (4, 5, 6)
    reverse_mask = frame > 3
    reverse_frame = frame[reverse_mask] - 3
    nuc_start_reverse = (protein_start[reverse_mask] * 3) + (reverse_frame - 1)
    nuc_end_reverse = (protein_end[reverse_mask] * 3) + 2 + (reverse_frame - 1)

    nucleotide_start[reverse_mask] = seq_length[reverse_mask] - nuc_end_reverse - 1 + reverse_offset[reverse_mask]
    nucleotide_end[reverse_mask] = seq_length[reverse_mask] - nuc_start_reverse - 1 + reverse_offset[reverse_mask]

    # Convert to back to 1-based coordinates for reporting.
    nucleotide_start += 1
    nucleotide_end += 1

    return nucleotide_start, nucleotide_end
