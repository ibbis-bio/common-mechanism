#!/usr/bin/env python3
# Copyright (c) 2021-2024 International Biosecurity and Biosafety Initiative for Science
"""
Script that checks results for regulated pathogen and prints any matched coordinates. Ignores any
synthetic constructs.

Usage:
  python check_reg_path.py -i INPUT -d database_folder -t threads

"""
import argparse
import logging
import os
import re
import sys
import textwrap
import pandas as pd
from commec.tools.blast_tools import read_blast, get_taxonomic_labels, get_top_hits
from commec.tools.blastn import BlastNHandler
from commec.tools.search_handler import SearchHandler
from commec.config.query import Query
from commec.config.result import (
    ScreenResult,
    HitResult,
    ScreenStep,
    ScreenStatus,
    HitScreenStatus,
    MatchRange,
    compare
)

pd.set_option("display.max_colwidth", 10000)


def _check_inputs(
        search_handle : SearchHandler,
        benign_taxid_path : str | os.PathLike,
        biorisk_taxid_path : str | os.PathLike,
        taxonomy_directory : str | os.PathLike
        ):
    """ 
    Simply check for the existance of files, 
    returns True if it is safe to continue. 
    """
    # check input files
    if not search_handle.check_output():
        logging.info("\t...ERROR: Taxonomic search results empty\n %s", search_handle.out_file)
        return False

    if not os.path.exists(benign_taxid_path):
        logging.error("\t...benign db file %s does not exist\n", benign_taxid_path)
        return False

    if not os.path.exists(biorisk_taxid_path):
        logging.error("\t...biorisk db file %s does not exist\n", biorisk_taxid_path)
        return False
    
    if not os.path.exists(taxonomy_directory):
        logging.error("\t...taxonomy directory %s does not exist\n", taxonomy_directory)
        return False

    if search_handle.is_empty(search_handle.out_file):
        logging.info("\tERROR: Homology search has failed\n")
        return False

    if not search_handle.has_hits(search_handle.out_file):
        logging.info("\t...no hits\n")
        return False
    
    return True

def update_taxonomic_data_from_database(
        search_handle : SearchHandler,
        benign_taxid_path : str | os.PathLike,
        biorisk_taxid_path : str | os.PathLike,
        taxonomy_directory : str | os.PathLike,
        data : ScreenResult,
        queries : dict[str, Query],
        step : ScreenStep,
        n_threads : int
        ):
    """
    Given a Taxonomic database screen output, update the screen data appropriately.
        search_handle : The handle of the search tool used to screen taxonomic data.
        benign/biorisk_handlers : Only used to determine the location of taxid related information
        taxonomy_directory : The location of taxonomy.
        data : the Screen data object, to be updated.
        step : Which taxonomic step this is (Nucleotide, Protein, etc)
        n_threads : maximum number of available threads for allocation.
    """
    logging.debug("Acquiring Taxonomic Data for JSON output:")

    if not _check_inputs(search_handle, benign_taxid_path, 
                         biorisk_taxid_path, taxonomy_directory):
        return 1

    # Read in lists of regulated and benign tax ids
    vax_taxids = pd.read_csv(benign_taxid_path, header=None).squeeze().astype(str).tolist()
    reg_taxids = pd.read_csv(biorisk_taxid_path, header=None).squeeze().astype(str).tolist()

    # The default is to pass, its up to the data to over-write this.
    if step == ScreenStep.TAXONOMY_AA:
        for query in data.queries.values():
            query.recommendation.protein_taxonomy_status = ScreenStatus.PASS
    if step == ScreenStep.TAXONOMY_NT:
        for query in data.queries.values():
            query.recommendation.nucleotide_taxonomy_status = ScreenStatus.PASS

    blast = read_blast(search_handle.out_file)
    blast = get_taxonomic_labels(blast, reg_taxids, vax_taxids, taxonomy_directory, n_threads)
    blast = blast[blast["species"] != ""]  # ignore submissions made above the species level

    # label each base with the top matching hit, but include different taxids attributed to same hit
    top_hits = get_top_hits(blast)

    if top_hits["regulated"].sum() == 0:
        logging.info("\t...no regulated hits\n")
        return 0
    
    # The non-coding fasta appends (start, stop) info to the filenames. This counters that.
    if step == ScreenStep.TAXONOMY_NT:
        top_hits["query acc."] = top_hits["query acc."].str.split(" ").str[0]
    
    # if ANY of the trimmed hits are regulated
    with pd.option_context('display.max_rows', None,
                    'display.max_columns', None,
                    'display.precision', 3,
                    ):

        unique_queries = top_hits['query acc.'].unique()
        for query in unique_queries:
            query_write = data.get_query(query)
            if not query_write:
                logging.debug("Query during %s could not be found! [%s]", str(step), query)
                continue

            unique_query_data : pd.DataFrame = top_hits[top_hits['query acc.'] == query]
            unique_query_data.dropna(subset = ['species'])
            regulated_only_data = unique_query_data[unique_query_data["regulated"] == True]
            regulated_hits = regulated_only_data['subject acc.'].unique()

            for hit in regulated_hits:
                regulated_hit_data : pd.DataFrame = regulated_only_data[regulated_only_data["subject acc."] == hit]
                hit_description = regulated_hit_data['subject title'].values[0]

                n_regulated_bacteria = 0
                n_regulated_virus = 0
                n_regulated_eukaryote = 0
                
                reg_taxids = [] # Regulated Taxonomy IDS
                non_reg_taxids = [] # Non-regulated Taxonomy IDS.
                reg_species = [] # List of species
                domains = [] # List of domains.
                match_ranges = [] # Ranges where hit matches query.

                for _, region in regulated_hit_data.iterrows():
                    match_range = MatchRange(
                        float(region['evalue']),
                        int(region['s. start']), int(region['s. end']),
                        int(region['q. start']), int(region['q. end'])
                    )

                    # Convert from non-coding to nt query coordinates if we're doing a NT taxonomy step.
                    if step == ScreenStep.TAXONOMY_NT:
                        match_range.query_start = queries[query].nc_to_nt_query_coords(match_range.query_start)
                        match_range.query_end = queries[query].nc_to_nt_query_coords(match_range.query_end)

                    match_ranges.append(match_range)

                    # Filter shared_site based on 'q. start' or 'q. end' (Previously only shared starts were used)
                    shared_site = top_hits[(top_hits['q. start'] == region['q. start']) | (top_hits['q. end'] == region['q. end'])]

                    # Filter for regulated and non-regulated entries
                    regulated = shared_site[shared_site["regulated"] == True]
                    non_regulated = shared_site[shared_site["regulated"] == False]

                    # Count domain information.
                    domain = region['superkingdom']
                    if domain == "Viruses":
                        n_regulated_virus += 1
                    if domain == "Bacteria":
                        n_regulated_bacteria +=1
                    if domain == "Eukaryota":
                        n_regulated_eukaryote+=1
                    domains.append(domain)

                    # Collect unique species from both regulated and non-regulated
                    reg_species.extend(regulated["species"])
                    # JSON serialization requires int, not np.int64, hence the map()
                    reg_taxids.extend(map(str, regulated["subject tax ids"]))
                    non_reg_taxids.extend(map(str, non_regulated["subject tax ids"]))

                    # These are the old values, now we simply count the size of the regulated, and non_regulated taxid arrays.
                    #n_reg += (top_hits["regulated"][top_hits['q. start'] == region['q. start']] != False).sum()
                    #n_total += len(top_hits["regulated"][top_hits['q. start'] == region['q. start']])

                # Uniquefy.
                reg_species = list(set(reg_species))
                reg_taxids = list(set(reg_taxids))
                non_reg_taxids = list(set(non_reg_taxids))
                match_ranges = list(set(match_ranges))

                reg_species.sort()
                reg_taxids.sort()
                non_reg_taxids.sort()

                recommendation : ScreenStatus = ScreenStatus.FLAG

                # TODO: Currently, we recapitulate old behaviour,
                # # " no top hit exclusive to a regulated pathogen: PASS" 
                #  however in the future:
                # if all hits are in the same genus n_reg > 0, and n_total > n_reg, WARN, or other logic.
                # the point is, this is where you do it.

                if len(non_reg_taxids) > 0:
                    recommendation = ScreenStatus.PASS

                # Update the query level recommendation of this step.
                if step == ScreenStep.TAXONOMY_AA:
                    query_write.recommendation.protein_taxonomy_status = compare(
                        query_write.recommendation.protein_taxonomy_status,
                        recommendation)
                if step == ScreenStep.TAXONOMY_NT:
                    query_write.recommendation.nucleotide_taxonomy_status = compare(
                        query_write.recommendation.nucleotide_taxonomy_status,
                        recommendation)

                regulation_dict = {"number_of_regulated_taxids" : str(len(reg_taxids)),
                                   "number_of_unregulated_taxids" : str(len(non_reg_taxids)),
                                   "regulated_eukaryotes": str(n_regulated_eukaryote),
                                   "regulated_bacteria": str(n_regulated_bacteria),
                                   "regulated_viruses": str(n_regulated_virus),
                                   "regulated_taxids": reg_taxids,
                                   "non_regulated_taxids" : non_reg_taxids,
                                   "regulated_species" : reg_species}

                # Append our hit information to Screen data.
                new_hit = HitResult(
                    HitScreenStatus(
                        recommendation,
                        step
                    ),
                    hit,
                    hit_description,
                    match_ranges,
                    {"domain" : [domain],"regulation":[regulation_dict]},
                )

                if query_write.add_new_hit_information(new_hit):
                    write_hit = query_write.get_hit(hit)
                    if write_hit:
                        write_hit.annotations["domain"] = domains # Always overwrite, better than our guess from biorisk.
                        write_hit.annotations["regulation"].append(regulation_dict)
                        write_hit.recommendation.status = compare(write_hit.recommendation.status, recommendation)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="in_file",
        required=True,
        help="Input query file (e.g. QUERY.nr.dmnd)",
    )
    parser.add_argument(
        "-d",
        "--database",
        dest="db",
        required=True,
        help="database folder (must contain vax_taxids and reg_taxids file)",
    )
    parser.add_argument("-t", "--threads", dest="threads", required=True, help="number of threads")
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(message)s",
        handlers=[logging.StreamHandler(sys.stderr)],
    )

    exit_code = check_for_regulated_pathogens(args.in_file, args.db, args.threads)
    sys.exit(exit_code)

def check_for_regulated_pathogens(input_file: str, input_database_dir: str, n_threads: int):
    """
    Check an input file (output from a database search) for regulated pathogens, from the benign and
    biorisk database taxids.
    """
    # Check input file
    if not os.path.exists(input_file):
        logging.error("\t...input query file %s does not exist\n", input_file)
        return 1
    sample_name = re.sub(r"\.nr.*|\.nt\.blastn", "", input_file)

    # Read in lists of regulated and benign tax ids
    benign_taxid_path = f"{input_database_dir}/benign_db/vax_taxids.txt"
    if not os.path.exists(benign_taxid_path):
        logging.error("\t...benign db file %s does not exist\n", benign_taxid_path)
        return 1
    vax_taxids = pd.read_csv(benign_taxid_path, header=None).squeeze().astype(str).tolist()

    biorisk_taxid_path = f"{input_database_dir}/biorisk_db/reg_taxids.txt"
    if not os.path.exists(biorisk_taxid_path):
        logging.error("\t...biorisk db file %s does not exist\n", biorisk_taxid_path)
        return 1
    reg_taxids = pd.read_csv(biorisk_taxid_path, header=None).squeeze().astype(str).tolist()

    # if there are already regulated regions written to file for this query, add to them
    hits1 = None
    if os.path.exists(sample_name + ".reg_path_coords.csv"):
        hits1 = pd.read_csv(sample_name + ".reg_path_coords.csv")

    if BlastNHandler.is_empty(input_file):
        logging.info(
            "\tERROR: Cannot check for regulated pathogens in empty or non-existent file: %s\n",
            input_file,
        )
        return 1

    if not BlastNHandler.has_hits(input_file):
        logging.info("\t... Skipping regulated pathogens check, no hits in: %s\n", input_file)
        return 0

    blast = read_blast(input_file)
    blast = get_taxonomic_labels(blast, reg_taxids, vax_taxids, input_database_dir + "/taxonomy/", n_threads)
    blast = blast[blast["species"] != ""]  # ignore submissions made above the species level

    # label each base with the top matching hit, but include different taxids attributed to same hit
    blast2 = get_top_hits(blast)

    reg_bac = 0
    reg_vir = 0
    reg_fung = 0

    # if this is the nucleotide screen, check if any weak protein flags can be negated with strong non-regulated nt ones
    if re.findall(".nt.blastn", input_file):
        if hits1 is not None:
            for region in range(0, hits1.shape[0]):  # for each regulated pathogen region
                # look at only the hits that overlap it
                htrim = blast2[
                    ~(
                        (blast2["q. start"] > hits1["q. end"][region])
                        & (blast2["q. end"] > hits1["q. end"][region])
                    )
                    & ~(
                        (blast2["q. start"] < hits1["q. start"][region])
                        & (blast2["q. end"] < hits1["q. start"][region])
                    )
                ]
                species_list = textwrap.fill(", ".join(set(htrim["species"])), 100).replace(
                    "\n", "\n\t\t     "
                )
                taxid_list = textwrap.fill(
                    ", ".join(map(str, set(htrim["subject tax ids"]))), 100
                ).replace("\n", "\n\t\t     ")
                percent_ids = " ".join(map(str, set(htrim["% identity"])))
                if htrim.shape[0] > 0:
                    if any(htrim["q. coverage"] > 0.90):
                        htrim = htrim[htrim["q. coverage"] > 0.90]
                        htrim = htrim.reset_index(drop=True)
                        descriptions = []
                        for row in range(htrim.shape[0]):
                            hit = htrim["subject title"][row]
                            descriptions.append(hit)
                        annot_string = "\n\t...".join(str(v) for v in descriptions)
                        logging.info(
                            "\t...Regulated protein region at bases "
                            + str(int(hits1["q. start"][region]))
                            + " to "
                            + str(int(hits1["q. end"][region]))
                            + " overlapped with a nucleotide hit\n"
                        )
                        logging.info(
                            "\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n",
                            species_list,
                            taxid_list,
                            percent_ids,
                        )

    if blast2["regulated"].sum():  # if ANY of the trimmed hits are regulated
        hits = pd.DataFrame(columns=["q. start", "q. end"])
        with pd.option_context(
            "display.max_rows",
            None,
            "display.max_columns",
            None,
            "display.precision",
            3,
        ):
            # for each hit (subject acc) linked with at least one regulated taxid
            for site in set(blast2["q. start"][blast2["regulated"] != False]):
                subset = blast2[(blast2["q. start"] == site)]
                subset = subset.sort_values(by=["regulated"], ascending=False)
                subset = subset.reset_index(drop=True)
                org = ""

                blast2 = blast2.dropna(subset=["species"])
                n_reg = (blast2["regulated"][blast2["q. start"] == site] != False).sum()
                n_total = len(blast2["regulated"][blast2["q. start"] == site])
                gene_names = ", ".join(set(subset["subject acc."]))
                end = blast2["q. end"][blast2["q. start"] == site].max()
                coordinates = str(int(site)) + " - " + str(int(end))

                species_list = textwrap.fill(
                    ", ".join(set(blast2["species"][blast2["q. start"] == site])), 100
                ).replace("\n", "\n\t\t     ")
                desc = blast2["subject title"][blast2["q. start"] == site].values[0]
                taxid_list = textwrap.fill(
                    ", ".join(
                        map(
                            str,
                            set(blast2["subject tax ids"][blast2["q. start"] == site]),
                        )
                    ),
                    100,
                ).replace("\n", "\n\t\t     ")
                percent_ids = " ".join(
                    map(str, set(blast2["% identity"][blast2["q. start"] == site]))
                )
                reg_ids = " ".join(
                    map(
                        str,
                        set(
                            blast2["regulated"][
                                (blast2["q. start"] == site) & (blast2["regulated"] != False)
                            ]
                        ),
                    )
                )

                # if some of the organisms with this sequence aren't regulated, say so
                if n_reg < n_total:
                    logging.info(
                        "\t\t --> Best match to sequence(s) %s at bases %s found in both regulated and non-regulated organisms\n"
                        % (gene_names, coordinates)
                    )
                    logging.info(
                        "\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n"
                        % (species_list, taxid_list, percent_ids)
                    )
                    logging.info("\t\t     Description: %s\n" % (desc))
                    # could explicitly list which are and aren't regulated?
                # otherwise, raise a flag and say which superkingdom the flag belongs to
                elif n_reg == n_total:
                    if subset["superkingdom"][0] == "Viruses":
                        reg_vir = 1
                        org = "virus"
                    elif subset["superkingdom"][0] == "Bacteria":
                        reg_bac = 1
                        org = "bacteria"
                    elif "superkingdom" in subset:
                        if subset["superkingdom"][0] == "Eukaryota":
                            org = "eukaryote"
                            reg_fung = 1

                    new_hits = subset[["q. start", "q. end"]].dropna()
                    if not new_hits.empty and not hits.empty:
                        hits = pd.concat([hits, new_hits], ignore_index=True)
                    elif not new_hits.empty:
                        hits = new_hits.copy()
                    logging.info(
                        "\t\t --> Best match to sequence(s) %s at bases %s found in only regulated organisms: FLAG (%s)\n"
                        % (gene_names, coordinates, org)
                    )
                    logging.info(
                        "\t\t     Species: %s (taxid(s): %s) (%s percent identity to query)\n"
                        % (species_list, taxid_list, percent_ids)
                    )
                    logging.info("\t\t     Description: %s\n" % (desc))
                else:  # something is wrong, n_reg > n_total
                    logging.info("\t...gene: %s\n" % gene_names)
                    logging.info("%s\n" % (blast["regulated"][blast["subject acc."] == gene_names]))
        hits = hits.drop_duplicates()
        # Create output file
        if hits1 is not None:
            hits = pd.concat([hits1, hits])
        hits.to_csv(sample_name + ".reg_path_coords.csv", index=False)

    if reg_vir == 0 and reg_bac == 0 and reg_fung == 0 and reg_fung == 0:
        logging.info("\t\t --> no top hit exclusive to a regulated pathogen: PASS\n")

    return 0


if __name__ == "__main__":
    main()
