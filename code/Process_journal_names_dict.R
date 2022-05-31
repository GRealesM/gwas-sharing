#################################
### Processing Journal Names ####
#################################

# Journal info from both GWAS catalog and iCiteR use journal abbreviations while SCR data use full journal names and ISSN. We'll use data from NLM catalog to match them.

library(data.table)

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
