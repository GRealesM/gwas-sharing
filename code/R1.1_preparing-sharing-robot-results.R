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

f <- fread(file.path(d,"data_found.txt"), header = F)$V1 %>% gsub("data/([0-9]+).xml", "\\1", ., perl = TRUE) %>% unique %>% as.numeric
nf <- fread(file.path(d, "data_notfound.txt"), header = F)$V1 %>% gsub("data/([0-9]+).xml", "\\1", ., perl = TRUE) %>% unique %>% as.numeric

intersect(f,nf) # Should eval zero

f <- data.table(PMID = f, Sharer_17 = "Y")
nf <- data.table(PMID = nf, Sharer_17 ="N")
ar <- rbind(f, nf)


### Load main data source, which includes studies in GWAS catalog, whether they share summary statistics or not, journal information, etc.
md <- fread(file.path(d,"R1_Main_data_20221020.tsv"))

md <- merge(md, ar, all.x = TRUE, by = "PMID")
table(md[ !is.na(Sharer_17), year]) # 2017 onwards only

# Use GWAS catalog classification for the rest
md[ is.na(Sharer_17), Sharer_17:=Public_ss]

  
# Save it
fwrite(md ,file.path(d, "R1_Main_data_20221128.tsv"), sep = "\t")
  