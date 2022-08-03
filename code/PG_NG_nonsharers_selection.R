#####################################################################
####     Preparing data for analysis of non-sharer articles      ####
#####################################################################

# Background: We'll select a random 50% of Plos Genet and Nat Genet papers classified as non-sharers in GWAS catalog for further analysis.
# Date: 2022/06/10

library(data.table)
set.seed(1)

x <- fread("../data/Main_data_20220610.tsv")
xp <- x[Journal == "PLoS Genet" & Public_ss == "N", .(Journal, Date, PMID, First_Author, Title)][, .SD[sample(.N, length(PMID)/2)]]
xn <- x[Journal == "Nat Genet" & Public_ss == "N", .(Journal, Date, PMID, First_Author, Title)][, .SD[sample(.N, length(PMID)/2)]]
xj <- rbind(xn, xp)
fwrite(xj, "../data/PG_NG_nonsharers.tsv", sep="\t")
