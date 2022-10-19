#####################################################################
#### Preparing data for GWAS summary statistics sharing project  ####
#####################################################################

# Background: We collected data on studies and datasets published and publicly
# available at the GWAS catalog, as well as data on SCImago for journal impact
# factor, and NLM catalog data on journals to match journal names.
# Here we'll process and prepare data for further analyses.

# Upon revision, it was suggested if including SJR data from 2021 only was enough or appropriate.
# We downloaded SJR info from 2005-2021 and we'll prepare a new main dataset with the appropriate SJR data by journal and year the paper was published. 

# Date: 2022/10/11

## Load packages
library(data.table)
library(magrittr)

# Load main data to get the journal abbreviations
m <- fread("../data/Main_data_20220610.tsv")

# Remove old SJR columns
m[, c("SJR", "H_index", "SJR_Best_Quartile", "Publisher"):=NULL]

# Load journal info
jinf <- lapply(2005:2021, function(x){
          j <- fread(paste0("../data/scimagojr_",x,".csv"))
          j[, Year:=x]
          setnames(j, paste0("Total Docs. (",x,")"), "Total Docs (Year)")
}) %>% rbindlist

# Keep only relevant columns
jinf <- jinf[, .(Title, SJR, Year, `SJR Best Quartile`, `H index`, Publisher)][, SJR:=as.numeric(gsub(",", ".", SJR))] 
setnames(jinf, c("SJR Best Quartile", "H index", "Title"), c( "SJR_Best_Quartile", "H_index", "SCR.title") )

fwrite(jinf, "../data/Full_scimago_info.tsv", sep="\t")

m[jinf , SJR := i.SJR, on=.(SCR.title, Year)][jinf , SJR_Best_Quartile := i.SJR_Best_Quartile, on=.(SCR.title, Year)][jinf , H_index := i.H_index, on=.(SCR.title, Year)] # Update values per year

fwrite(m, "../data/R1_Main_data_20221011.tsv", sep="\t", na = NA)
