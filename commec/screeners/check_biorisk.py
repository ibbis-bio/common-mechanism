#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks the output from hmmscan and prints to screen the results

Usage:
 python check_biorisk.py -i INPUT.biorisk.hmmscan -d databases/biorisk_db/ 
"""
import logging
import os
import pandas as pd
from commec.config.query import Query
from commec.tools.hmmer import (
    readhmmer,
    remove_overlaps,
    recalculate_hmmer_query_coordinates,
    append_nt_querylength_info,
    HmmerHandler
)
from commec.config.result import (
    ScreenResult,
    HitResult,
    ScreenStep,
    ScreenStatus,
    HitScreenStatus,
    MatchRange,
    compare
)

logger = logging.getLogger(__name__)

def _guess_domain(search_string : str) -> str:
    """ 
    Given a string description, try to determine 
    which domain of life this has come from. Temporary work around
    until we can retrieve this data directly from biorisk outputs.
    """
    def contains(search_string : str, search_terms):
        for token in search_terms:
            if search_string.find(token) == -1:
                continue
            return True
        return False

    search_token = search_string.lower()
    if contains(search_token, ["vir", "capsid", "RNA Polymerase"]):
        logger.debug("Determined virus from \"%s\"", search_string)
        return "Virus"
    if contains(search_token, ["cillus","bact","coccus","phila","ella","cocci","coli"]):
        logger.debug("Determined bacteria from \"%s\"", search_string)
        return "Bacteria"
    if contains(search_token, ["eukary","nucleus","sona","odium","myces"]):
        logger.debug("Determined Eukaryote from \"%s\"", search_string)
        return "Eukaryote"
    logger.debug("Could not guess domain from \"%s\"", search_string)
    return "not assigned"

def update_biorisk_data_from_database(search_handle : HmmerHandler,
                                      biorisk_annotations_file : str | os.PathLike,
                                      data : ScreenResult,
                                      queries : dict[str, Query]):
    """
    Takes an input database, reads its outputs, and updates the input data to contain
    biorisk hits from the database. Also requires passing of the biorisk annotations CSV file.
    Inputs:
        search : search_handle - The handler which has performed a search on a database.
        biorisk_annotations_csv_file : str - directory/filename of the biorisk annotations provided by Commec.
        data : ScreenResult - The ScreenResult to be updated with information from database, interpeted as Biorisks.
    """
    logger.debug("Starting Update Biorisk data from database...")
    logger.debug("Directory: %s", search_handle.db_directory)
    logger.debug("Directory/file: %s", search_handle.db_file)

    hmm_folder_csv = biorisk_annotations_file #os.path.join(search_handle.db_directory,"biorisk_annotations.csv")
    if not os.path.exists(hmm_folder_csv):
        logger.error("\t...biorisk_annotations.csv does not exist\n %s", hmm_folder_csv)
        return 1
    if not search_handle.validate_output():
        logger.error("\t...database output file does not exist\n %s", search_handle.out_file)
        return 1
    if search_handle.is_empty(search_handle.out_file):
        logger.error("\t...ERROR: biorisk search results empty\n")
        return 1

    # We delay non-debug logging to sort messages via query.
    log_container = {key : [] for key in data.queries.keys()}

    for query in data.queries.values():
        query.status.set_step_status(ScreenStep.BIORISK, ScreenStatus.PASS)

    if not search_handle.has_hits(search_handle.out_file):
        return 0

    # Read in Output, and parse.
    hmmer : pd.DataFrame = readhmmer(search_handle.out_file)
    logger.debug("Biorisk Import: shape: %s preview:\n%s", hmmer.shape, hmmer.head())
    keep1 = [i for i, x in enumerate(hmmer['E-value']) if x < 1e-20]
    hmmer = hmmer.iloc[keep1,:]
    logger.debug("Biorisk Filterd by E-Value: shape: %s preview:\n%s", 
                 hmmer.shape, hmmer.head())
    append_nt_querylength_info(hmmer, queries)
    logger.debug("Appended Query NT length: shape: %s preview:\n%s", 
                 hmmer.shape, hmmer[["query name","nt_qlen"]].head())
    recalculate_hmmer_query_coordinates(hmmer)
    logger.debug("Recalculated AA to NT coords: shape: %s preview:\n%s", 
                 hmmer.shape, hmmer[["ali from","ali to","q. start","q. end"]].head())
    hmmer = remove_overlaps(hmmer)
    logger.debug("Removed overlaps: shape: %s preview:\n%s", hmmer.shape, hmmer.head())

    # Read in annotations.
    lookup : pd.DataFrame = pd.read_csv(hmm_folder_csv)
    lookup.fillna(False, inplace=True)

    # Append description, and must_flag columns from annotations:
    hmmer['description'] = ''
    hmmer['Must flag'] = False
    hmmer = hmmer.reset_index(drop=True)
    for model in range(hmmer.shape[0]):
        name_index = [i for i, x in enumerate([lookup['ID'] == hmmer['target name'][model]][0]) if x]
        hmmer.loc[model, 'description'] = lookup.iloc[name_index[0], 1]
        hmmer.loc[model, 'Must flag'] = lookup.iloc[name_index[0], 2]

    # Update the data state to capture the outputs from biorisk search:
    unique_queries = hmmer['query name'].unique()
    logger.debug("Unique Queries: shape: %s preview:\n%s", unique_queries.shape, unique_queries)
    for affected_query in unique_queries:
        logger.debug("\tProcessing query: %s", affected_query)
        biorisk_overall : ScreenStatus = ScreenStatus.PASS

        query_data = data.get_query(affected_query)
        if not query_data:
            logger.error("Query during hmmscan could not be found! [%s]", affected_query)
            continue

        queries[affected_query[:-2]].confirm_has_hits()

        # Grab a list of unique queries, and targets for iteration.
        unique_query_data : pd.DataFrame = hmmer[hmmer['query name'] == affected_query]
        unique_targets = unique_query_data['target name'].unique()
        
        logger.debug("\tData for %s: shape: %s preview:\n%s", 
                     affected_query, unique_query_data.shape, unique_query_data.head())
        logger.debug("\tTargets (%i) hit by %s: %s", 
                     len(unique_targets), affected_query , unique_targets)

        for affected_target in unique_targets:
            logger.debug("\t\tProcessing target %s", affected_target)
            unique_target_data = unique_query_data[unique_query_data['target name'] == affected_target]
            target_description = ", ".join(set(unique_target_data['description'])) # First should be unique.
            must_flag = unique_target_data['Must flag'].iloc[0] # First should be unique.
            match_ranges = []
            match_string = ""
            for _, region in unique_target_data.iterrows():
                match_range = MatchRange(
                    float(region['E-value']),
                    int(region['hmm from']), int(region['hmm to']),
                    int(region['q. start']), int(region['q. end'])
                )
                match_ranges.append(match_range)
                match_string += f"{match_range.query_start}-{match_range.query_end}, "

            # Remove final ", "
            match_string = match_string[:-2]

            target_recommendation = ScreenStatus.FLAG if must_flag > 0 else ScreenStatus.WARN
            logger.debug("\t\tTarget recommendation from this hit: %s", target_recommendation)

            biorisk_overall = compare(target_recommendation, biorisk_overall)

            # Precalculate some logging information...
            regulation_str : str = "Regulated Gene" if must_flag else "Virulence Factor"
            log_message = f"\t --> {regulation_str} found at coordinates: {match_string}."
            logger.debug(f"{affected_query[:-2]:<25}" + log_message)
            log_container[affected_query[:-2]].append(log_message)

            # Deal with whether this biorisk should collapse into an existing hit or not.
            hit_data : HitResult = query_data.get_hit(affected_target)
            if hit_data:
                logger.debug("\t\tHit already existed! Extending hit data ranges only...")
                hit_data.ranges.extend(match_ranges)
                logger.debug("Updated hit: %s", hit_data)
                continue

            domain : str = _guess_domain(""+str(affected_target)+target_description)

            new_hit : HitResult = HitResult(
                HitScreenStatus(
                    target_recommendation,
                    ScreenStep.BIORISK
                ),
                affected_target,
                target_description,
                match_ranges,
                {"domain" : [domain],"regulated":[regulation_str]},
            )
            logger.debug("\t\tHit was unique, creating new hit result information...\n%s", new_hit)
            query_data.hits[affected_target] = new_hit


        # Update the recommendation for this query for biorisk.
        query_data.status.set_step_status(ScreenStep.BIORISK, biorisk_overall)

    # Do all non-verbose logging in order of query:
    for query_name, log_list in log_container.items():
        if len(log_list) == 0:
            continue
        s = "" if len(log_list) == 1 else "s"
        logger.info(f" Biorisk{s} in {query_name}:")
        for log_text in log_list:
            logger.info(log_text)

    return 0