# The Common Mechanism - How screening works
This document outlines the logical flow, and intended behaviour of commec screen. It is to go in-depth in the current state of the code as of writing (2.0.0) and updated as changes occur. This document focuses on `commec screen`, as the other cli options are more simple, and of less concern than the screening function of commec.

## Step 0 - Initialisation
The commec entry point starts at cli.py. Where the primary mode is determined by the initial argument, commec arg0, where arg0 is:
* `screen` Use of commec as a DNA synthesis screening tool.
* `list` Reporting on Control List annotations.
* `flag` Summarise multiple screen runs in a csv.
* `setup` Tool for downloading required commec databases.
* `gui` Tool for running a local server that can be accessed to run commec screen via graphical user interface.

`screen` is the major functionality of commec here, and requires additional input and parsing of that input.

All configurable input parameters for commec are derived from configuration yaml file. The default configuration is encoded in screen-default-config.yaml and represents the contract for how the input data is laid out. A second .yaml can be ingested during run time, which can override any or all of these values. Finally some config options can be further overriden via CLI inputs. (commec/config/screen_io.py)

a `.fasta` file is a required input, and must be in a valid .fasta format. Internally, this .fasta file is split in its individual queries, any Query that is too small or too large in size is identified to be skipped. A copy of the input fasta is made with all valid queries and stored in the `/inputs_` folder. Each valid query is initiated as a dictionary of Query objects. (commec/config/query.py)
Each Query then handles its own 6-frame translation, each frame is stored in a second fasta, and labelled query_name_[1-6] with _X component representing the translation frame. The two new .fasta files represent the nucleotide, and protein sequences of all input queries for the screen and are passed as inputs into all subsequent tool calls. A final non-coding region fasta is also generated during step, which is also handled by each Query object, the specifics of which will be covered in Step 3.

Some of the configuration points to the location of expected database files, used by all subsequent steps. Search Handler objects are initialised using these settings. (commec/tools/search_handler.py). During initialisation, the expected directory location, and entry file for each tool are checked during object initialisation. Inheritance is utilised for customised behaviour for the differences in search handlers for BLASTX, or HMSCAN for example. Search handlers handle the sub process calls to external tools, as well as all of their parameters, their output files, and the ingestion of those output files into a pandas DataFrame object. Any agnostic to screening logic data ingestion sanitisation is the role of the Search Handler.

The Control Lists database is not imported with a Search Handler, and instead is imported using the control_list module (commec/control_list/). Import for all control lists follows a recursive search in the given directory. region_definitions.json is a requirement at the top-most level. Control list data is merged into a pandas DataFrame object, where accessions (currently represent taxid, but also could be any string i.e. uniprot id etc) are used as indexes, with allowed duplicate entries.

Finally, the primary output of commec screen is the output json. The structure of this is contracted in commec/config/results.py, and consists of a hard coded hierarchy of dataclass objects. json_io.py handles the dynamic import and export of these dataclasses even under changes. Breaking changes should increment the version of the json output to avoid undefined behaviour. The structure of dataclasses is initialised based of all queries input, and represents the current State of the commec program. Theorhetically, at any point in time, the json can be output (and is on unexpected failure). The update of the fields within this structure occurs over the commec screen runtime. The first updates to the results is the version information of commec, its databases, and the start time of the run.

All of the above logic discussed above occurs in the setup function of the Screen class in commec/screen.py. They are executed as of writing in the current order:
* Init console logging
* Parse input parameters (cli and yaml)
* Init file logging.
* Check for database updates (if applicable)
* Validate all input databases that use Searh Handlers
* Validate and import Control List Data
* Import all query information, generate Query objects, add query to results, translate query.
* Populate results with database revision info, search tool version info, global query info, start time.

After the setup function is successful, Screen.run continues, which simply calls the functions for steps 1,2,3, and 4. The logic for each step is handled in the commec/screeners/ scripts.

## Step 1 - Biorisk Search

The Biorisk step starts with asking the HMM search handler to perform its search by externally calling hmmer on the input translated query .fasta file. This calls cli arguments are handled by the search handler, and include number of threads, which is also capped at 4 for hmmer. This is done due to extremely diminishing returns on additional thread use on hmmer past 4 cores.

The screen handler, as well as the biorisk annotations .csv file location, and the output data json structure, and the current query information, is all passed to `commec/screeners/check_biorisk.py`, via the parse_biorisk_hits entry point.

We check the outputs of the Biorisk step for validity, and if everything is okay, we proceed.

We set all queries to PASS in the biorisk step. They will be set to False if hits are found.

We check for the existance of _any_ hits first, to allow for early return.

The hmmer output file is imported via the search handlers read methods. The data is imported as a pandas DataFrame object. First, we append the information of the total query length measured in nucleotides to each row. Note that all Biorisk data is in AA coordinates, or translated coordinates. The nucleotide length is important information to be used to calculate the nucleotide coordinates, as well as the frame info. Furthermore, appended nucleotide info is used during biorisk e-value filtering.

Next we filter the biorisk hits based on e-value, the exact filter is dependant on whether the query is considered short ot long, the constant of which is currently at 200 nucleotides (constants.py). 

If the query is long, the cut-off for e-value is 1e-20.
If the query is short, the cut-off for e-value is calculated dynamically as the inverse  of the query length to the power of 2.598 (BIORISK_SHORT_QUERY_EVALUE_EXPONENT). These constants were determined experimentally.

After filtering, the nucleotide coordinates are back calculated, to ensure that the hmmer output is correctly mapped onto the users input query nucleotide coordinates.

Then we remove overlapping hmmer hits, this has the effect of checking every nucleotide, and keeping the top hit (best on hmmer score output) per nucleotide. This allows overlapping hits, so long as they have at least 1 nucleotide outside of the overlapping region. It also allows hits that are completely encompassed within another hit, if they have a better score.

The annotations are imported, a simple mapping between the biorisk database hmmer profile names, and hand curated annotation information, as well as the "must_flag" column, which simple separates virulence factors from pathogens of concern. The must_flag and description information are added to the hmmer DataFrame.

For each query, the list of hits pertaining to that query is taken from the DataFrame. A HitResult object is constructed for each hit, and the status is set to WARN if the must_flag for that hit is False. Otherwise, the status for that his is FLAG. The hit object is then added to the Hits list for that query in the json output data structure.

After appending all hits for a query, its overall biorisk status is set based on all of its biorisk hits.

Finally, the logs for all saved hits are displayed.

The output of the biorisk step is now stored as a series of HitResult objects per query as a state in the results.py output dataclasses structure object.

## Step 2 - Best Match Protein Search

Screening the best match protein search is done so long as skip-taxonomy is not true, as passed by cli arguments. Like the previous biorisk step, a Blastx search handler is invoked to run the external tool for best match protein search, which is in this case, a Blastx call. It takes the nucleotide query fasta (not the translated version), as translation is handled within the blastx tool itself. The exact cli arguments used to call Blastx are present in the Search Handler itself, some are hard coded, and others like thread count is derived from input cli/yaml configurations. Once the blastx call is complete, we proceed with the the parsing of the output using the parse_taxonomy_hits() entry point from commec/sceeners/check_reg_path.py. In a parallel to biorisk, it also takes the search handler, the output data ScreenResult object, the query info, the exact screen step.

The output from the Blastx call is imported as a pandas DataFrame, the details of which are handled by the search handler, as a glorified csv document.

We check that the call to blast resulted in a valid output, and after which set all queries best match protein step to PASS, for subsequent update.

We early return on no hits from this step.

We start by immediately removing "Bad accessions" which is a small list of hand curated accessions which are highly problematic. They are identified by matching to the subject acc. column.

We then label the imported data with control information. This process starts with splitting the multiple taxid entries in a logical way:
Any blast hit labelled "Multispecies" in its subject title, or has more than 100 taxids in its taxid field has its taxid set to the keyword "MULTISPECIES" for special processing. Multispecies are treated as non-regulated/mixed entries.

Any remaining multi-taxid rows are now split into multiple rows, each representing a single taxid.

We then remove any taxids that are in the control list modules "Ignored taxids" list. These are calculated during control list construction and typically involved hypothetical proteins etc.

Next, we create a new column for the DataFrame "regulated" and set it to True/False depending on whether each taxid is detected to be on any list in the control list module. This is a rapid lookup, and only outputs True/False to save time, full control list information is imported later.

Then we get the unique list of taxids, and calculate the "parent of all parents" for that entry. This effectively maps all taxids into the same cluster if they are related (i.e. one is the child of another) in the context of what is present in the control list hierarchy. The parental cluster, which is itself the valid regulated taxid of the parent, is stored in the "control_hash" column.

We now check if there were any regulated hits, and exist early if there are none.

We now go through the remaining DataFrame hits at the per-Query level. (Which makes the following operations on the query specific part of the dataframe cheaper)

For each query, we get the top hits. First the dataframe coordinates are updated such that q. start is always smaller than q. end. (i.e. reverse hits are simply restored with the start and end points swapped). Then we sort the dataframe such that regulated hits are first, and then by percent identity, and bit score. After sorting, the top ranked hits based on "% identity" are kept so long as they have no overlaps with another hit.
We then sort the remaining "top hits" by "% identity" before iteratively trimming the edges of the hits, which alters that start and end points of each hit that is weaker than an overlapping stronger hit to match the start and end points. (This step is important later for when we check for mixed results)
Finally, all hits are filtered based on being > 50 bps.

> Note: _trim_overlapping() form blast_tools has some old code, and should be revisted.

Then the remaining hits are iterated over by unique control list cluster, each unique cluster can have 1 or more "hits" associated with it. Each set of clustered data is always used to generate a single HitResult object.

We use the parent taxid (the control cluster hash) as the name for the HitResult object.

Any duplicate "Subject title" entries are removed, these are rare, but are always useless duplicates. More of an issue when we matched against the full nr database.

We check what overlapping non-regulated data from the blast search exists for this cluster, only the top 10 are kept. They are detected with shared start and end sites, which is conveniently the case after get_top_hits() has updated the start and stop sites with trim edges.

Now for each regulated taxid in this cluster, we grab all the control list information (A list of which lists it has hit, as well as the information of the list item itself) We store the regulated information, and also append the annotation information into separate buckets depending on whether each list it hits is regionally controlled, fully controlled, just a warning etc. This is derived from the ListMode, a flag that was set on the control list info on import, which depended on the regional context at import. What that regional context was is currently not important at this stage.

After capturing every regulated taxids control list info per list mode outcome, we act differently depending on whether:
* There are any fully controlled lists. Here we treat other list modes as mixed / passed.
* There are only warning compliance controlled items, again treat other modes as mixed.
* There are only conditionally controlled list items, this will be a pass, but still recorded if there are 2 or more lists involved.

Depending on these outcomes, the data is used to populate the rest of a HitResult object, its final screen status for the hit updated, and most importantly, the list annotation information is copied into the free-form annotations field of the HitResult object.

This processed is done for each query, for each unique set of cluster taxids, and for each regulated taxid. The result for each is a HitResult is is appended to the Hits list for each relevant query. The Query is then updated in its best match protein status based on the most important outcome here.

Finally, cached logs are output at the per query level, separated from the logic above.

## Step 3 - Best Match Nucleotide Search

The best match search uses exactly the same logic as above, however there are some important differences in the input data, which will be the focus on this section.

Like the previous steps, all external calls is performed by a search handler, which in this instance calls Blastn. The cli arguments of which are of course stored in the search handlers run function, and some of which will be updated via cli/config options.

Importantly, before Blastn can be called, we need to detect non-coding regions for each query. The nucleotide best match step is **only** run on regions of the sequence that have no protein best match hit. Importantly, the threshold is based on "% identity" and all of the raw blastx data is used to find **any** hit that meets the threshold (NON_CODING_REGION_PERCENT_IDENTITY_THRESHOLD in constants, currently 80 %). This is unrelated to any of the hit matching or filtering of the protein best match step.

It is therefore possible for a weaker hit to be present in the best match protein output, but not be strong enough to prevent nucleotide best match hits from overlapping with it. This is by design to be redundant rather than over-sure.

The non-coding regions are calculated as any region of the query sequence not covered by those protein hits which pass the % identity threshold, which also covers a total size of greater than 50 nucleotides. The non-coding region information is stored in the Query objects themselves, and are used to create a third .fasta file, the non-coding.fasta. Each region is labelled with an integer "_X" suffix. In much the same way that the 6 frame translations are encoded on the name of the query in the amino acid .fasta.

Because of this, all blastn data at this step has to be prepared:
* The query name is interrogated, and suffix _X gleamed to understand which non-coding region the hit is for.
* The offset for that non-coding region is applied to all q. start and q. end coordinates, to represent the data in query sequence space, rather than non-coding sequence space.

Apart from this initial setup, the following per query logic for handling hits is identical to the best match protein step as outlined above.

## Step 4 - Low Concern Search

We now have a ScreenResult data structure, fill of Querys with 0 or many hits from the previous steps, however some of the hits may be "cleared" as low concern, if they overlap with known common parts. The low concern step consists of 3 tool calls: blastn, cmscan, and hmmer, Each checking for low concern dna, rna, and protein respectively.

We start with an early exit if there are no hits that require clearing.
Then each of the 3 screen searchers are run sequentially, hmmer, blastn, and then cmscan.
We then read the low concern annotation csv, which contains a simple ID,Description column mapping. The 3 search handlers for each low concern step, as well as these descriptions, and the query info is used at the parse_low_concern_hits() entrypoint of check_low_concern.py.

The protein data is read from the hmmer output via search handler. We then filter all results based on the low concern protein e-value cutoff, a constant (1e-20). The protein data is further modified to also contain the total query length. Which is also used to recalculate the protein AA coordinates into nucleotide coordinates.

The RNA data is read from the cmscan handler, this is the only current cmscan use in commec. The data is simply read in via the search handler.

The DNA data is read from the blastn handler.

> Note: Currently no additional filtering is occuring for the RNA and DNA outputs, The DNA blastn is already somewhat filtered in outputs due to the search handlers blastn call cli arguments. But the CMSCAN output doesn't have this, and may require investigating whether it can output poorly matching hits.

Each query is then interrogated, if there are not hits, the low concern step is skipped, and we move onto the next query,

We create filtered copies of the input data to be processed at a per query level. We then iterate over every Hit in the query. If the hit is not FLAG or WARN, we ignore it.

If the Hit is a Flag or Warn, we then interrogate each set of low concern data specific for that query:

The protein data is trimmed to the region of the Hit. If there are no remaining data, we early exit. We calculate the coverage of the low concern data over the hit region in question. The coverage calculations contain both a coverage_nt (a raw count of how many nucleotides overlap) and a coverage_ratio (The percentage of overlap with the Hit region.)
We further filter the protein hits based on having at least 50 (constant MINIMUM_PEPTIDE_COVERAGE) nucleotides in coverage.
We then filter the protein hits based on having at least 80% (constant MINIMUM_QUERY_COVERAGE_FRACTION) coverage_ratio.

We then use the top remaining hit to create a HitResult for the low concern hit. We also update the FLAG/WARN hit to be FLAG (Cleared) or WARN (Cleared).

The RNA low concern data for the query is also trimmed to the hit region. The same coverage calculations also are done for coverage_nt and coverage ratio.  Only the coverage_nt is checked to be less than 50 (constant MINIMUM_RNA_BASEPAIR_COVERAGE). A HitResult for the low concern is created by the top entry in the remaining data, and the original hit status updated to be cleared.

The DNA low concern data for the query is also trimmed to the region, and coverage calculated. Similar to RNA, only the ratio of coverage is checked against the constant MINIMUM_SYNBIO_COVERAGE_FRACTION. (Currently at 80%) The top remaining hit is then used to create a HitResult, and the original hit cleared.

> Note, potentialy issues in that the output low concern data is never sorted() before grabbing the literal "top" hit. i.e. that hit at the 0th location in the remaining filtered dataframe. It is likely there is only 1 hit here anyway, but it may be prudent to perform a sort on e-value at these steps.

At each of the previous steps, the resulting low concern hit results are returned. They are then added to a list to be added to the query data, and the total number of clears is compared against the total number of flagged / warn hit results.

> Note: Much of this counting region / cleared_regions is redundant, and a hangover from when HitResults had multiple regions associated with them, as well as that this _was_ the place that you would clear the whole query, now the whole query will be automatically cleared because we overwrite the status with Flag (Cleared) rather than keeping the Flag status. Allowing the ScreenResult data to be updated at any time.

If all regions have been cleared, we log this info, and update the Query screens status. We take the full list of new low concern hits, ensure it is fully of unique entries only, and then log the outcome and add the information to the screen result json output data.

## Step 5 - Conclude

Now commec is completed, there is some minor tidy up. The complete time is logged, the output ScreenResult data object is updated(), which ensures that all the hits and their screen statuses propagate to the Querys status, and that all rationales are updated, so that the data is ready to be output to JSON. Some summary data is printed to the log, and the lifetime of the Screen class object expires.

As the Screen object is deleted, the __del__ method is called. This calculates the time taken, and updates that field to the result data. It is then converted into a json, and output. If the screen was a success, the HTML report is also generated and output. This all occurs in the delete function as this can be called even under some exception circumstances during python shutdown, and allows us to get a JSON in partial amounts depending on what error has occured.