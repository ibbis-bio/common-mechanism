#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script for clustering a pandas dataframe based on ranges.
"""
from commec.tools.blast_tools import _trim_edges, shift_hits_pos_strand
import pandas as pd

def find_clusters(
    input_data: pd.DataFrame,
    start_heading: str = "q. start",
    end_heading: str = "q. end",
    label_heading: str = "cluster",
) -> tuple[pd.DataFrame, list[tuple]]:
    """
    Groups rows into clusters where any set of rows with collectively overlapping
    ranges (transitively) belong to the same cluster.

    Returns a copy of input_data with a new integer cluster-ID column, and a list
    of (min_start, max_end) tuples — one per cluster.
    """
    assert not input_data.empty, "Empty dataframe given to find_clusters()"
    assert start_heading in input_data.columns
    assert end_heading in input_data.columns

    sorted_df = input_data.sort_values(by=start_heading)

    clusters: list[tuple] = []
    row_cluster_ids: dict = {}
    current_cluster = -1
    current_max_end = None

    for idx, row in sorted_df.iterrows():
        start = row[start_heading] - 1
        end = row[end_heading] + 1

        if current_cluster == -1 or start >= current_max_end:
            current_cluster += 1
            clusters.append((row[start_heading], row[end_heading]))
            current_max_end = end
        else:
            current_max_end = max(current_max_end, end)
            clusters[current_cluster] = (clusters[current_cluster][0], current_max_end)

        row_cluster_ids[idx] = current_cluster

    output_data = input_data.copy()
    output_data[label_heading] = output_data.index.map(row_cluster_ids)

    return output_data, clusters


def trim_overlapping(blast: pd.DataFrame):
    """
    Remove any hits that are completely overlapped by another, higher-quality hit.
    """
    # set start to the lowest coordinate and end to the highest to avoid confusion
    blast = shift_hits_pos_strand(blast)

    # if any multispecies hits contain regulated pathogens, put the regulated up top
    if "regulated" in blast:
        blast = blast.sort_values(by=["regulated"], ascending=False)

    # rank hits by percent identity, then bit score
    blast = blast.sort_values(by=["% identity", "bit score"], ascending=False)
    blast = blast.reset_index(drop=True)

    blast2 = blast
    # only keep top-ranked hits that don't overlap
    #df = blast[blast["query acc."] == query]

    for i in blast.index:  # run through each hit from the top
        for j in blast.index[(i + 1) :]:  # compare to each below
            if j in blast2.index:
                # if beginning and end of the higher-rank hit both overlap or extend further
                # than the beginning and end of lower-ranked hit, discard the lower-ranked hit
                if (
                    blast.loc[i, "q. start"] <= blast.loc[j, "q. start"]
                    and blast.loc[i, "q. end"] >= blast.loc[j, "q. end"]
                ):
                    # Unless the hits have the same coordinates and % identity
                    if (
                        blast.loc[i, "q. start"] < blast.loc[j, "q. start"]
                        or blast.loc[i, "q. end"] > blast.loc[j, "q. end"]
                        or blast.loc[i, "% identity"] > blast.loc[j, "% identity"]
                    ):
                        blast2 = blast2.drop([j])
    blast2 = blast2.reset_index(drop=True)

    return blast2

def get_top_hits(blast: pd.DataFrame):
    """
    Trim BLAST results down to the top hit for each base.
    """

    assert isinstance(blast, pd.DataFrame), "get_top_hits expects a pandas dataframe object."

    if blast.empty:
        logger.debug("Empty dataframe passed to Get Top Hits.")
        return blast

    top_hits = trim_overlapping(blast)
    top_hits = top_hits.sort_values("% identity", ascending=False)

    # only keep coordinates of each hit that are not already covered by a better hit
    df = top_hits

    rerun = 1
    while (
        rerun == 1
    ):  # edges of hits can be moved within a higher scoring hit in the first pass
        df, rerun = _trim_edges(df)

    for j in df.index:
        top_hits.loc[j, "subject length"] = max(
            [df.loc[j, "q. start"], df.loc[j, "q. end"]]
        ) - min([df.loc[j, "q. start"], df.loc[j, "q. end"]])
        top_hits.loc[j, "q. start"] = df.loc[j, "q. start"]
        top_hits.loc[j, "q. end"] = df.loc[j, "q. end"]

    top_hits = top_hits.sort_values("q. start")
    top_hits = top_hits[top_hits["q. start"] != 0]

    # only keep annotations covering 50 bases or more
    top_hits = top_hits[top_hits["subject length"] >= 50]
    top_hits = top_hits.reset_index(drop=True)
    return top_hits
