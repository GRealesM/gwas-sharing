################################################################################
### Preparing data for GWAS summary statistics sharing project. Revision 1  ####
################################################################################

# Background: We collected data on studies and datasets published and publicly
# available at the GWAS catalog, as well as data on SCImago for journal impact
# factor, and NLM catalog data on journals to match journal names.
# Here we'll process and prepare data for further analyses.

# Date: 2022/10/20


## Load packages
library(data.table)
library(magrittr)
library(ggplot2)
library(iCiteR)


#---- Data preprocessing

# We need to pre-process some of the helper files prior to moving on. If those files are already available, you can skip this part.

# Journal info from both GWAS catalog and iCiteR use journal abbreviations while SCR data use full journal names and ISSN. We'll use data from NLM catalog to match them.

system("wget https://ftp.ncbi.nih.gov/pubmed/J_Entrez.txt -O ../data/J_Entrez.txt")
system('grep -v "\\-\\-\\-\\-" ../data/J_Entrez.txt > ../data/J_Entrez.tsv')
jm <- fread("../data/J_Entrez.tsv", sep="", header = FALSE)
jm[, V1:=sub(":", "\t", V1)]
jm[, c("label", "value"):=tstrsplit(V1, "\t")][,V1:=NULL]
jm[, num:=rep(1:table(jm$label)[1], each =7)] 
jm <- dcast(jm, num ~ label,value.var = "value")
jm[, num:=NULL]
jm <- jm[!is.na(IsoAbbr) | !is.na(MedAbbr)] # Remove NA from abbreviations
all(jm$IsoAbbr == jm$MedAbbr) # TRUE, so we can remove one of them
jm[, IsoAbbr:=NULL]
fwrite(jm, "../data/NLM_journals.tsv", sep="\t")

# Get SCImago journal impact factor data together

# Load journal info
jinf <- lapply(2005:2021, function(x){
  j <- fread(paste0("../data/scimagojr_",x,".csv"))
  j[, Year:=x]
  setnames(j, paste0("Total Docs. (",x,")"), "Total Docs (Year)")
}) %>% rbindlist

# Keep only relevant columns
jinf <- jinf[, .(Title, Issn, SJR, Year, `SJR Best Quartile`, `H index`, Publisher)][, SJR:=as.numeric(gsub(",", ".", SJR))] 
setnames(jinf, c("SJR Best Quartile", "H index", "Title"), c( "SJR_Best_Quartile", "H_index", "SCR.title") )

fwrite(jinf, "../data/Full_scimago_info.tsv", sep="\t")



#---- Load datasets

## List of publications and available summary statistics downloaded on 2022-05-24
pubs <- fread("../data/gwas-catalog-v1.0.3-studies-r2022-05-17.tsv", sep="\t") # Published studies. Downloaded from: https://www.ebi.ac.uk/gwas/api/search/downloads/studies_new
ss <- fread("../data/list_gwas_summary_statistics_24_May_2022.csv")  # Studies with summary statistics. Downloaded from https://www.ebi.ac.uk/gwas/downloads/summary-statistics 
citinfo <- fread("../data/Citation_info_20220525.tsv") # Citation info for PMID, extracted using iCiteR on 2022-05-25. See below for the code used to generate it.
SCR <- fread("../data/Full_scimago_info.tsv") # SCimagoJR dataset with journal impact factors for year 2021. Downloaded from https://www.scimagojr.com/journalrank.php for years 2005-2021. 
NLM <- fread("../data/NLM_journals.tsv") # NLM catalog with journal names, ISSN, and abbreviations. See Process_journal_names_dict.R for code to generate it.


#---- Tidy up the dataset

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
# fwrite(pubs, "../data/Publications_20220525.tsv", sep="\t")

#----  Get citation info. Retrieving is capped at 1000 articles, so we'll run it 6 times

# Note: The result of the function below is time-dependent. We ran it initially on 2022-05-25 and saved it, so now the code is commented and we'll use the saved file

# citinfo <- lapply(0:5, function(x){
#   start_idx = (x * 1000) + 1
#   if(x == 0) start_idx = 1
#   end_idx = start_idx + 999
#   if(end_idx > nrow(m)) end_idx = nrow(m)
#   
#   pmids <- m$PMID[start_idx:end_idx]
#   message("Getting metrics for papers ", start_idx, " to ", end_idx) 
#   papers <- get_metrics(pmids)
# })
# citinfo <- rbindlist(citinfo)
# fwrite(citinfo, "../data/Citation_info_20220525.tsv", sep="\t")

# Extract key data for analysis and merge both datasets
m <- merge(m, citinfo, by.x="PMID", by.y="pmid")

# We'll now make the distinction between the year of online publication and year of publication in print
# For citation purposes, we'll use the year they became available online, rather than when they got in print
setnames(m,"year","year_inprint")
m[,year:=as.numeric(substr(DATE,1,4))]

# Check that year > year_inprint
table(m$year_inprint - m$year)
# -1    0    1    2    3 
#  1 4698 1014   42    1
# For one paper, its "online" year is posterior to its "inprint" year. We'll coerce both years to match.
m[ year_inprint < year, year_inprint:=year]


#----  Match Journal Abbreviations from GWAS catalog with journal names and add SCImago SJR impact factor info

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
guide.SCR <- melt(mSCR, id.vars = c("SCR.title", "SJR", "Year", "SJR_Best_Quartile", "H_index", "Publisher"), measure.vars =  c("V1", "V2", "V3", "V4")) %>% unique
guide.SCR <- na.omit(guide.SCR)
guide.NLM <- melt(NLM.g, id.vars = "MedAbbr", measure.vars = c("Issn.print", "Issn.online"))
m.guides <- merge(guide.SCR, guide.NLM, by="value") # Match by ISSN, either online and print
m.guides <- unique(m.guides[, .(SCR.title, MedAbbr, SJR, Year, SJR_Best_Quartile, H_index, Publisher)]) # Select variables and remove duplicates arisen from a journal having both online and printed ISSN.

m <- merge(m, m.guides, by.x = c("journal","year"), by.y = c("MedAbbr", "Year"), all.x = TRUE) # Merge to get SJR info. Note: We're using year of publishing, not year in print


#---- Clean up the dataset

setnames(m, c("DATE", "First_Author", "Publisher"), c("date", "first_author", "publisher"))
# Order columns
setcolorder(m, c("PMID", "journal", "Public_ss", "ds_count",  "relative_citation_ratio", "SJR", "year", "first_author", "title", "authors", "date", "year_inprint", "citation_count","citations_per_year", "expected_citations_per_year", "cited_by", "references", "doi", "H_index", "SJR_Best_Quartile", "publisher"))

#---- Save dataset

# Note that this is a revised dataset, so the date and name changed
fwrite(m, "../data/R1_Main_data_20221020.tsv", sep="\t", na = NA)



#---- Get citation history, year per year of the articles

# We'll extract this info in the same fashion as we did previously. We'll extract citing articles from the `cited_by` column and retrieve their year, then we'll be able to count number of citations each year.

cith <- m[,.(PMID, cited_by)]
max(sapply(strsplit(m$cited_by, " "), length)) # Max citations is 6770, so we'll need that number of columns
cith[, paste0("V", 1:6770):= tstrsplit(cited_by, " ")]
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
fwrite(sum.cith, "../data/Citations_per_year_20220611.tsv", sep="\t")





