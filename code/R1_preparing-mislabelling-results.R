################################################################################
###  Preparing results after new sharing classification status analysis      ###   
###  for GWAS summary statistics sharing project. Revision 1                 ###
################################################################################

# Background: We webscraped text from GWAS articles classified as non-sharers by 
# GWAS catalog from 2017, to double-check if they actually shared data.
# To that end, we used a Ruby robot that used keywords to classify articles.
# We'll now extract the results for further analysis steps.


library(data.table)
library(magrittr)

d="../data" #"~/rds/rds-cew54-basis/Projects/gwas-sharing/data/"

# This file is the result of the analysis to find articles that published data outside GWAS catalog. Can be found at https://raw.githubusercontent.com/chr1swallace/data-sharing-search/main/data_found.txt
f <- fread(file.path(d,"data_found.txt"), header = F)$V1 %>% gsub("data/([0-9]+).xml", "\\1", ., perl = TRUE) %>% unique %>% as.numeric
f <- data.table(PMID = f, Sharer_updated = "Y")


### Load main data source, which includes studies in GWAS catalog, whether they share summary statistics or not, journal information, etc.
md <- fread(file.path(d,"R1_Main_data_20221020.tsv"))

md <- merge(md, f, all.x = TRUE, by = "PMID")
table(md[ !is.na(Sharer_updated), year]) # Distribution of newly found sharers

# Use GWAS catalog classification for the rest
md[ is.na(Sharer_updated), Sharer_updated:=Public_ss]

  
# Save it
fwrite(md ,file.path(d, "R1_Main_data_20221207.tsv"), sep = "\t")


# Update non-sharer table to generate table S5
ns <- fread(file.path(d, "Nonsharers_all.tsv"))
ns <- merge(ns, f, by="PMID", all.x = TRUE)
ns[ is.na(Sharer_updated), Sharer_updated:="N"]
table(ns$Sharer_updated) # 217 new sharers

fwrite(ns , "../tables/R1_S5_mislabelled_updated.tsv", sep = "\t")
