#####################################################################
#### Preparing data for GWAS summary statistics sharing project  ####
#####################################################################

# Background: We collected data on studies and datasets published and publicly
# available at the GWAS catalog, as well as data on SCImago for journal impact
# factor, and NLM catalog data on journals to match journal names.
# Here we'll process and prepare data for further analyses.

# Date: 2022/06/10

## Load packages
library(data.table)
library(magrittr)
library(ggplot2)
library(iCiteR)


## List of publications and available summary statistics downloaded on 2022-05-24
pubs <- fread("../data/gwas-catalog-v1.0.3-studies-r2022-05-17.tsv", sep="\t") # Published studies. Downloaded from: https://www.ebi.ac.uk/gwas/api/search/downloads/studies_new
ss <- fread("../data/list_gwas_summary_statistics_24_May_2022.csv")  # Studies with summary statistics. Downloaded from https://www.ebi.ac.uk/gwas/downloads/summary-statistics 
SCR <- fread("../data/scimagojr_2021.csv") # SCimagoJR dataset with journal impact factors for year 2021. Downloaded from https://www.scimagojr.com/journalrank.php for 2021, without additional filters.
NLM <- fread("../data/NLM_journals.tsv") # NLM catalog with journal names, ISSN, and abbreviations. See Process_journal_names_dict.R for code to generate it.

## Replace key column names by something more useful
setnames(pubs, c("FIRST AUTHOR", "PUBMED ID", "STUDY ACCESSION", "ASSOCIATION COUNT"), c("First_Author", "PMID", "accession", "association_count")) 
setnames(ss, c("PubMed ID", "Study accession"), c("PMID", "accession"))

## Label studies with summary statistics and check how many of them in each study
pubs[, Public_ss := ifelse(PMID %in% ss$PMID, "Y","N")]
ds_num <- ss[, .(ds_count = .N), by=PMID] # Number of published summary statistics datasets
pubs <- merge(pubs, ds_num, by="PMID", all.x=TRUE)
pubs[is.na(ds_count), ds_count:=0]

## Update changed PMID. This PMID was removed because it was duplicated in PubMed, so we update it to the current ID. It didn't have summary statistics.  
pubs[PMID == 24513584, PMID:=24516880]

## Subset columns to keep key information per paper
m <- unique(pubs[, .(PMID, Public_ss, ds_count, DATE, First_Author)])

# Get citation info. Retrieving is capped at 1000 articles, so we'll run it 6 times
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

# Extract key data for analysis and merge both datasets
m <- merge(m, citinfo, by.x="PMID", by.y="pmid")

# Match Journal Abbreviations from GWAS catalog with journal names 

jab <- unique(m$journal) # unique abbreviations
NLM.g <- NLM[ MedAbbr %in% jab] 
# There are some duplications, let's get rid of them
dup.nlmg <- NLM.g[duplicated(MedAbbr), MedAbbr]
NLM.g[MedAbbr %in% dup.nlmg]
NLM.g <- NLM.g[!NlmId %in% c("7803060", "0074126")] # Checked duplicated in SCR and the latter two didn't have an entry by ISSN in SCR, so we'll remove them

# Now we need to match SCR info with NLM.g. We'll do it by matching ISSN. 
# SCR file has two types of ISSN (Online and Print) in one comma-separated string, and the order is not constant (????) so we'll try to separate and match
NLM.g[, Issn.print:=gsub("-", "", `ISSN (Print)`)][, Issn.online:=gsub("-", "", `ISSN (Online)`)] # SCR Issn are missing the hyphen, so let's remove it
mSCR <- copy(SCR)
mSCR[, c("V1", "V2", "V3", "V4"):=tstrsplit(Issn, ", ")]
guide.SCR <- melt(mSCR, id.vars = "Title", measure.vars =  c("V1", "V2", "V3", "V4"))
guide.SCR <- na.omit(guide.SCR)
guide.NLM <- melt(NLM.g, id.vars = "MedAbbr", measure.vars = c("Issn.print", "Issn.online"))
m.guides <- merge(guide.SCR, guide.NLM, by="value")
m.guides <- unique(m.guides[, .(Title, MedAbbr)])
setnames(m.guides, "Title", "SCR.title")

m <- merge(m, m.guides, by.x = "journal", by.y = "MedAbbr", all.x = TRUE) # Merge with guides to get the journal title
sSCR <- SCR[ Title %in% m$SCR.title, .(Title, SJR, `SJR Best Quartile`, `H index`, Publisher)][, SJR:=as.numeric(gsub(",", ".", SJR))] # Select columns to use. Also replace comma for numeric SJR
m <- merge(m, sSCR, by.x = "SCR.title", by.y="Title", all.x = TRUE) 
setnames(m, c("H index", "SJR Best Quartile", "journal", "DATE", "year", "title", "authors"), c("H_index", "SJR_Best_Quartile","Journal", "Date", "Year", "Title", "Authors"))

# Order m
setcolorder(m, c("PMID", "Public_ss", "ds_count", "Year", "First_Author", "Title", "Authors", "Date", "citation_count","citations_per_year", "relative_citation_ratio", "expected_citations_per_year", "cited_by", "references", "doi", "Journal", "SJR", "H_index", "SJR_Best_Quartile", "Publisher"))

# We want to get the citation history of the articles, we'll extract this info in the same fashion as we did previously. We'll extract citing articles from the `cited_by` column and retrieve their year, then we'll be able to count number of citations each year.
cith <- m[,.(PMID, cited_by)]
cith[, paste0("V", 1:6770):= tstrsplit(cited_by, " ")] # Max citations is 6770, so we'll need that number of columns
cith <- melt(cith, id.vars = "PMID", measure.vars = paste0("V", 1:6770))
cith <- na.omit(cith)
cith[, variable:=NULL]
lup <- unique(cith$value) # Lookup items


citing.items <- lapply(0:188, function(x){
	start_idx = (x * 1000) + 1
	if(x == 0) start_idx = 1
	end_idx = start_idx + 999
	if(end_idx > length(lup)) end_idx = length(lup)
	
	pmids <- lup[start_idx:end_idx]
	message("Getting metrics for papers ", start_idx, " to ", end_idx) 
	papers <- get_metrics(pmids)
})
citing.items <- rbindlist(citing.items)
citing.items <- citing.items[,.(pmid, year)]
citing.items$pmid <- as.character(citing.items$pmid)
cith <- merge(cith, citing.items, by.x = "value", by.y = "pmid")
# Now create a summary
sum.cith <- cith[, .(count = .N), by=c("PMID", "year")]


# Store data
fwrite(pubs, "../data/Publications_20220525.tsv", sep="\t")
fwrite(citinfo, "../data/Citation_info_20220525.tsv", sep="\t")
fwrite(m, "../data/Main_data_20220610.tsv", sep="\t", na = NA)
fwrite(sum.cith, "../data/Citations_per_year_20220611.tsv", sep="\t")
