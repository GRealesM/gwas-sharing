# GWAS summary statistics sharing project

This repository is intended to store the data and results from the GWAS data sharing project started on summer 2022.

**Code** contains the following scripts:

* **preparing-data.R**: contains the steps we followed to clean our data prior to analysis. These included matching journal names between different datasets to get citation and impact factor together, as well as fetching individual-level paper citations per year, among other steps.
* **Process_journal_names_dict.R** contains a wrapper to extract a journal name/journal abbreviation dictionary to match journal names across datasets.
* **analysis-and-figures-guille.R** contains all analyses and figures generated in this work.
* **PG_NG_nonsharers_selection.R** contains the code used to select a random 50% of articles classified as non-sharers at GWAS catalog published in two of the journals with the most published GWAS: Nature Genetics and PLoS Genetics.

**Data** contains the downloaded and generated datasets for analysis. 
Datasets were downloaded on the last week of May, 2022. **GWAS catalog** list of studies was dated from 2022/05/17, whereas the list of published summary statistics dated from 2022/05/24.
**SCImago** dataset on journal impact factor and other metrics correspond to 2021.
See `code/` for details on file origins and generation.


Guillermo Reales
2022/08/03
