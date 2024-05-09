# GWAS summary statistics sharing project

This repository is intended to store the data and results from the *Sharing GWAS summary statistics results in more citations: evidence from the GWAS catalog* paper, currently submitted and available as a preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.09.27.509657v1)!

After the first round of revision (December 2022), we have updated our manuscript and code.

Find below the description of the main scripts and datasets used in this work. Note that this repository contains the main analysis only; see [here](https://github.com/chr1swallace/data-sharing-search) for the code and data used in the extended search for shared files.


**Code** contains the following (updated) scripts:

* **01_Preparing_data.R**: contains the steps we followed to clean our data prior to analysis. These included matching journal names between different datasets to get citation and impact factor together, as well as fetching individual-level paper citations per year, adding per-year SJR, and updating GWAS sharing classification according to our new search. Check this file for links to the original data sources and how files used in analyses were generated.
* **02_Analysis.R** contains all analyses and figures generated in this work. 
* **03_PG_NG_nonsharers_selection.R** contains the code used to select a random 50% of articles classified as non-sharers at GWAS catalog published in two of the journals with the most published GWAS: Nature Genetics and PLoS Genetics.
* **04_response_to_reviewers_v2.R** contains the scripts used in our response to referees during the first peer-review round.


**Data** contains the downloaded and generated datasets for analysis. 

Downloaded datasets:

* **gwas-catalog-v1.0.3-studies-r2022-05-17.tsv**: List of published studies in GWAS catalog. Downloaded on late May 2022 from https://www.ebi.ac.uk/gwas/api/search/downloads/studies_new.
* **list_gwas_summary_statistics_24_May_2022.csv**: List of studies with summary statistics in GWAS catalog. Downloaded on late May 2022 from https://www.ebi.ac.uk/gwas/downloads/summary-statistics.
* **Citation_info_20220525.tsv**: Citation information for articles in GWAS catalog, by PMID, extracted using iCiteR on 2022-05-25.
* **J_Entrez.txt**: Full NLM catalog in PubMed, used for journal name matching between SCimago and GWAS catalog. Downloaded from https://ftp.ncbi.nih.gov/pubmed/J_Entrez.txt
* **scimagojr_[2005-2021].tsv**: SCimagoJR journal impact factors for individual years. Manually downloaded from https://www.scimagojr.com/journalrank.php.

Generated datasets (See `code/01_Preparing_data.R` for code to generate them):

* **R1_Main_data_20221207.tsv**: Latest dataset used in our analysis. Other files named similary but with different dates correspond to previous data versions. See scripts in `code/` to learn more about their generation/usage.
* **Full_scimago_info.tsv**: Merged SCimagoJR dataset with journal impact factors for years 2005-2021, after merging scimagojr_*.tsv files.
* **NLM_journals.tsv**: NLM catalog with journal names, ISSN, and abbreviations.
* **Citations_per_year_20220611.tsv**: Dataset on the number of citations obtained per article, per year since publication.
* **PG_NG_nonsharers.tsv**: A list of 353 randomly selected papers published in PLoS Genetics or Nature Genetics, and classified by GWAS catalog as non-sharers, used for our first sharing status chek-up. The result from this analysis is shown in Supplementary Table 2.
* **Non_sharers_all.tsv**: List of articles classified as non-sharers in GWAS catalog, input for our extended sharing search.
* **data_found.txt**: List of articles reclassified as sharers by our extended sharing search.

We also provide the datasets used to generate the figures in their respective tab-separated files (ie. Figure1.tsv, Figure2.tsv, Figure3a.tsv, and Figure3b.tsv).



Guillermo Reales
2023/01/17
