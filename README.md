# GWAS summary statistics sharing project

This repository is intended to store the data and results from the *Sharing GWAS summary statistics results in more citations: evidence from the GWAS catalog* paper, currently submitted and available as a preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.09.27.509657v1)!

After the first round of revision (December 2022), we have updated our manuscript and code

**Code** contains the following (updated) scripts:

* **R1_preparing-data_v2.R**: contains the steps we followed to clean our data prior to analysis. These included matching journal names between different datasets to get citation and impact factor together, as well as fetching individual-level paper citations per year, adding per-year SJR, and updating GWAS sharing classification according to our new search.
* **R1_analysis_v1.R** contains all analyses and figures generated in this work.
* **PG_NG_nonsharers_selection.R** contains the code used to select a random 50% of articles classified as non-sharers at GWAS catalog published in two of the journals with the most published GWAS: Nature Genetics and PLoS Genetics.
* **R1_response_to_reviewers_v2.R** contains the scripts used in our first revision.

**Data** contains the downloaded and generated datasets for analysis. 

See `code/` for details on file origins and generation.

* **R1_Main_data_20221207.tsv** contains the latest dataset used in our analysis.


For the code and data used in the extended search for shared files, see [here](https://github.com/chr1swallace/data-sharing-search).


Guillermo Reales
2022/12/12
