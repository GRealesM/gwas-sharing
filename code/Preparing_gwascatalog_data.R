########################################
####   Preparing GWAS catalog data  ####
########################################

# Background: This script is intended to pre-process GWAS catalog data on published studies with and without summary statistics.

# Date: 2022-05-24
# Author: Guillermo Reales

# Load packages
library(data.table)
library(ggplot2)
library(iCiteR)

# List of publications and available summary statistics downloaded on 2022-05-24
pubs <- fread("../data/gwas-catalog-v1.0.3-studies-r2022-05-17.tsv", sep="\t") # Published studies
ss <- fread("../data/list_gwas_summary_statistics_24_May_2022.csv")  # Studies with summary statistics

# Replace key column names by something more useful
setnames(pubs, c("FIRST AUTHOR", "PUBMED ID", "STUDY ACCESSION", "ASSOCIATION COUNT"), c("First_Author", "PMID", "accession", "association_count")) 
setnames(ss, c("PubMed ID", "Study accession"), c("PMID", "accession"))

# Label studies with summary statistics and check how many of them in each study
pubs[, Public_ss := ifelse(PMID %in% ss$PMID, "Y","N")]
ds_num <- ss[, .(ds_count = .N), by=PMID] # Number of published summary statistics datasets
pubs <- merge(pubs, ds_num, by="PMID", all.x=TRUE)
pubs[is.na(ds_count), ds_count:=0]

# Update changed PMID. This PMID was removed because it was duplicated in PubMed, so we update it to the current ID. It didn't have summary statistics.  
pubs[PMID == 24513584, PMID:=24516880]

# Subset columns to keep key information per paper
m <- unique(pubs[, .(PMID, Public_ss, ds_count)])


# Get citation info
citinfo <- lapply(0:5, function(x){
	start_idx = (x * 1000) + 1
	if(x == 0) start_idx = 1
	end_idx = start_idx + 999
	if(end_idx > nrow(m)) end_idx = nrow(m)
	
	
	pmids <- m$PMID[start_idx:end_idx]
	message("Getting metrics for papers ", start_idx, " to ", end_idx) 
	papers <- get_metrics(pmids)
})
citinfo <- rbindlist(citinfo)

# Extract key data for analysis
cit <- citinfo[, .(pmid, journal, year, is_research_article, relative_citation_ratio, nih_percentile, citation_count, citations_per_year, expected_citations_per_year, field_citation_rate)] 

m <- merge(m, cit, by.x="PMID", by.y="pmid")

jn <- unique(m$journal)
get_impactfactor(jn)

# Store data
fwrite(pubs, "../data/Publications_20220525.tsv", sep="\t")
fwrite(citinfo, "../data/Citation_info_20220525.tsv", sep="\t")
fwrite(m, "../data/Summary_info_20220525.tsv", sep="\t")

# Questions about the data:

# 1 - What's the proportion of studies with shared or not shared data?
# Most studies don't have shared summary statistics in GWAS catalog
ggplot(m, aes(x=Public_ss))+ geom_bar(stat="count") + xlab("Summary statistics available?")  

# 2 - How did sharing change over time?
ggplot(m, aes(x=year, fill=Public_ss))+ geom_bar(position = "dodge", stat="count") + xlab("Summary statistics available?") 

# 3 - Do sharing result in higher citations?
ggplot(m, aes(x = Public_ss, y=citations_per_year)) + geom_boxplot() + xlab("Summary statistics available?") + ylab("Citations per year")

# 4 - Does this change over time?
ggplot(m, aes(x=as.factor(year), y=citations_per_year, fill=Public_ss))+ geom_boxplot() + ylab("Citations per year")

# 5 - Are there differences between journals?
j_50 <- unique(m[, .N, by=journal][N > 50, journal])
#hist(m[, .N, by=journal]$N, breaks=200)
mj <- m[journal %in% j_50 & journal != "Nature"]
ggplot(mj, aes(x=journal, y=citations_per_year, fill=Public_ss))+ geom_boxplot() + xlab("Journal (20+ papers)") + ylab("Citations per year") + theme(axis.text.x = element_text(angle = -45, hjust=0))

# 6 - 


