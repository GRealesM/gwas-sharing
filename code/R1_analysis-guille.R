## Revision 1 - Extra analyses

# Here we'll address two questions
# 1) Is molcel important? How much does the result change if we remove it?
# 2) Is 2021 SJR data enough? 
# 3) How does mislabelling change over time?

# For the first question, we'll apply the same dataset we used for the main analyses, for the second, we'll use an updated dataset with yearly SJR scores.

# Load packages

library(data.table)
library(xlsx)
library(magrittr)
library(ggplot2)
# devtools::install_github("chr1swallace/seaborn_palettes")
library(seaborn)
library(cowplot)
theme_set(theme_cowplot())


#------------ Preparing data


## Load and process data
d="../data" #"~/rds/rds-cew54-basis/Projects/gwas-sharing/data/"
list.files(d)

### Load main data source, which includes studies in GWAS catalog, whether they share summary statistics or not, journal information, etc.
x=file.path(d,"Main_data_20220610.tsv") %>% fread()
dim(x)

# Check journals and SJR info
x[ , .(Journal, SJR)] %>% unique %>% nrow
# 723
x[ !is.na(SJR), .(Journal, SJR)] %>% unique %>% nrow
# 696

setnames(x,c("Journal","Year"),c("journal","year")) # to match previous data columns
## Date is first published date (perhaps online), Year is published in print
## year. Some papers have citations after their Date, but before their Year. Use
## truncated Date as year of publication.
setnames(x,"year","year_inprint")
x[,year:=as.numeric(substr(Date,1,4))]
## exclude 2022 as incomplete and 2006 and before because no shared studies
with(x,table(year, Public_ss))
x=x[year < 2022 & year > 2006]

summary(x$relative_citation_ratio)
hist((x$relative_citation_ratio))
hist(log(x$relative_citation_ratio)) # but missing values for 0
## offset by minimum observed non-zero value
mincit=min(x[!is.na(relative_citation_ratio) & relative_citation_ratio > 0 , relative_citation_ratio])
x[relative_citation_ratio==0, relative_citation_ratio:=relative_citation_ratio+mincit]
x[,y:=log(relative_citation_ratio)] # Our Response Variable
hist(x$y,breaks=100)
x[,y_ss:=Public_ss=="Y"]

#------------ Question 1 - Is molcel important? How much does the result change if we remove it?

### Model fitting


# Prepare dataset for modelling, we'll remove missing data for relevant variables. We'll also consider years 2007-2020
xmodel <- x[!is.na(SJR) & !is.na(y) & year <=2020]
## Exploration by journal
xjournal <- xmodel[, .(articles = .N, shared_ss = sum(y_ss)), by= "journal"][order(articles, decreasing = TRUE)]
#fwrite(xjournal, "../tables/Journals_table.tsv", sep="\t")

## look at most common journals - do we see the same story within journal?
tt=table(xmodel$journal)
sort(tt,decreasing = TRUE) %>% head(., 20)
jkeep=sort(tt,decreasing = TRUE) %>% names() %>% intersect(., x[Public_ss=="Y"]$journal) %>% head(.,20)
# ggplot(x[journal %in% jkeep], aes(x=year,fill=Public_ss)) +
#   geom_histogram(binwidth=1,position="dodge",col="grey") +
#   scale_fill_seaborn() +
#   facet_grid(journal ~ .)

xmodel[,sjournal:=ifelse(journal %in% jkeep,journal,"Other")]             # Use top 20 most common journals with at least one shared dataset, pool the rest as "other"
xmodel[,sjournal:=factor(sjournal)][,sjournal:=relevel(sjournal,"Other")]


# Model without ss, see analysis-and-figures-guille.R for information on how we reached these models
m9=glm(y~year + log(SJR) + sjournal + molecular_cellular, data=xmodel)   # m9: molecular_cellular. Original model used 

m5=glm(y~year + log(SJR) + sjournal, data=xmodel)   # m5: year, log journal impact factor, and journal, using top20 and "other". Same as m9 but without
BIC(m5, m9)
#    df      BIC
# m5 24 12884.21
# m9 25 12866.83
# BIC gives lower score to m9, but we'll use m5 to see how it changes the result.

m5a=glm(y~year + log(SJR) + sjournal + Public_ss, data=xmodel)
BIC(m5,m5a)
#     df      BIC
# m5  24 12884.21
# m5a 25 12743.70
# Adding Public_ss has a huge impact on model fit

summary(m5a)
# Call:
# glm(formula = y ~ year + log(SJR) + sjournal + Public_ss, data = xmodel)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -3.9180  -0.5204   0.0327   0.5433   3.3051  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                   147.796092   7.630790  19.368  < 2e-16 ***
#   year                                           -0.073537   0.003786 -19.423  < 2e-16 ***
#   log(SJR)                                        0.730144   0.020258  36.043  < 2e-16 ***
#   sjournalAm J Hum Genet                          0.027753   0.086997   0.319  0.74973    
#   sjournalAm J Med Genet B Neuropsychiatr Genet   0.341875   0.119247   2.867  0.00416 ** 
#   sjournalAnn Rheum Dis                          -0.081598   0.144320  -0.565  0.57183    
#   sjournalDiabetes                                0.164485   0.127597   1.289  0.19743    
#   sjournalEur J Hum Genet                         0.040081   0.132656   0.302  0.76255    
#   sjournalFront Genet                             0.010231   0.141415   0.072  0.94233    
#   sjournalGastroenterology                        0.087036   0.166260   0.523  0.60066    
#   sjournalHum Mol Genet                           0.216801   0.055670   3.894 9.98e-05 ***
#   sjournalJ Allergy Clin Immunol                 -0.247862   0.129120  -1.920  0.05496 .  
#   sjournalJ Hum Genet                            -0.200439   0.130214  -1.539  0.12379    
#   sjournalJ Med Genet                            -0.147504   0.161029  -0.916  0.35971    
#   sjournalMol Psychiatry                          0.369783   0.079535   4.649 3.42e-06 ***
#   sjournalNat Commun                              0.240348   0.061243   3.924 8.81e-05 ***
#   sjournalNat Genet                               0.071593   0.059714   1.199  0.23061    
#   sjournalNature                                  0.734996   0.123824   5.936 3.12e-09 ***
#   sjournalPLoS Genet                              0.418624   0.062050   6.747 1.69e-11 ***
#   sjournalPLoS One                                0.162764   0.058223   2.796  0.00520 ** 
#   sjournalProc Natl Acad Sci U S A                0.436411   0.156928   2.781  0.00544 ** 
#   sjournalSci Rep                                 0.158102   0.077910   2.029  0.04248 *  
#   sjournalTransl Psychiatry                       0.147252   0.104713   1.406  0.15972    
#   Public_ssY                                      0.566643   0.046181  12.270  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.7427586)
# 
# Null deviance: 7793.4  on 4941  degrees of freedom
# Residual deviance: 3652.9  on 4918  degrees of freedom
# AIC: 12581
# 
# Number of Fisher Scoring iterations: 2

## conclude sharing has highly significant effect (t=12.270, p<2e-16, coef=0.566643), this is slightly better than our earlier model!
summary(m5a)$coeff["Public_ssY",,drop=FALSE]
cat("Ratio of RCR in publications that shared vs did not = ",
    exp(coef(m5a)["Public_ssY"]),
    "\n95% CI: ",
    exp(coef(m5a)["Public_ssY"]+ c(-1,1) * 1.96 * sqrt(vcov(m5a)["Public_ssY","Public_ssY"])),
    "\n")
# Ratio of RCR in publications that shared vs did not =  1.762342 
# 95% CI:  1.609829 1.929303

# The model now predicts slightly more effect on sharing on citations!

# Table of results
predictors = rownames(summary(m5a)$coeff)
sigtable1 <- data.table(preds=predictors,
                        exp_estimate = exp(coef(m5a)),
                        low.ci = sapply(predictors, function(x){exp(coef(m5a)[x]+ -1 * 1.96 * sqrt(vcov(m5a)[x,x]))}),
                        hi.ci =  sapply(predictors, function(x){exp(coef(m5a)[x]+ 1 * 1.96 * sqrt(vcov(m5a)[x,x]))}))

#  fwrite(sigtable1, "../tables/Sigtable_m5a.tsv", sep = "\t")

######################################
### Some additional visualisations ###
######################################

## How is molcel related to log(RCR) and sharing status?

mcp <- ggplot(xmodel, aes(molecular_cellular, y, col=y_ss)) +
  geom_point(aes(group=y_ss))+
  geom_smooth()+
  scale_colour_seaborn("Shared Sumstats    ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes"))+
  labs(x="Molcel score",y="log(RCR)")
mcp

hp <- ggplot(xmodel, aes(human, y, col=y_ss)) +
  geom_point(aes(group=y_ss))+
  geom_smooth()+
  scale_colour_seaborn("Shared Sumstats    ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes"))+
  labs(x="Human score",y="log(RCR)")
hp

hmcp <- ggplot(xmodel, aes(human, molecular_cellular, col=y_ss)) +
  geom_point(aes(group=y_ss))+
  geom_smooth()+
  scale_colour_seaborn("Shared Sumstats    ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes"))+
  labs(x="Human score",y="Molecular cellular score")
hmcp



#------------ Question 2 - Is 2021 SJR data enough?

# We'll load an updated dataset containing yearly SJR scores for each journal, and apply the same data cleaning procedure to it.

### Load main data source, which includes studies in GWAS catalog, whether they share summary statistics or not, journal information, etc.
xb=file.path(d,"R1_Main_data_20221011.tsv") %>% fread()
dim(xb)

# Check journals and SJR info
xb[ , .(Journal)] %>% unique %>% nrow
# 723
xb[ !is.na(SJR), .(Journal)] %>% unique %>% nrow
# 674
# Some journals are lost. This is because they had SJR data for 2021, but not for the years the papers in GWAS catalog were published, so effectively we can't count them.

setnames(xb,c("Journal","Year"),c("journal","year")) # to match previous data columns
## Date is first published date (perhaps online), Year is published in print
## year. Some papers have citations after their Date, but before their Year. Use
## truncated Date as year of publication.
setnames(xb,"year","year_inprint")
xb[,year:=as.numeric(substr(Date,1,4))]
## exclude 2022 as incomplete and 2006 and before because no shared studies
with(xb,table(year, Public_ss))
xb=xb[year < 2022 & year > 2006]
# NOTE: Here SJR scores make reference to the old "year", which is now year_inprint.


summary(xb$relative_citation_ratio)
hist((xb$relative_citation_ratio))
hist(log(xb$relative_citation_ratio)) # but missing values for 0
## offset by minimum observed non-zero value
mincit=min(xb[!is.na(relative_citation_ratio) & relative_citation_ratio > 0 , relative_citation_ratio])
xb[relative_citation_ratio==0, relative_citation_ratio:=relative_citation_ratio+mincit]
xb[,y:=log(relative_citation_ratio)] # Our Response Variable
hist(xb$y,breaks=100)
xb[,y_ss:=Public_ss=="Y"]

## As we saw in Figure 3, it takes 2-3 years for sharing effect to stabilise. Exclude 2021 because not enough time has passed
xb=xb[year<=2020]
# Prepare dataset for modelling, we'll remove missing data for relevant variables. We'll also consider years 2007-2020
xbmodel <- xb[!is.na(SJR) & !is.na(y)]
## Exploration by journal
xbjournal <- xbmodel[, .(articles = .N, shared_ss = sum(y_ss)), by= "journal"][order(articles, decreasing = TRUE)]
#fwrite(xjournal, "../tables/Journals_table.tsv", sep="\t")

## look at most common journals - do we see the same story within journal?
ttb=table(xbmodel$journal)
sort(ttb,decreasing = TRUE) %>% head(., 20)
jbkeep=sort(ttb,decreasing = TRUE) %>% names() %>% intersect(., xb[Public_ss=="Y"]$journal) %>% head(.,20)

xbmodel[,sjournal:=ifelse(journal %in% jbkeep,journal,"Other")]             # Use top 20 most common journals with at least one shared dataset, pool the rest as "other"
xbmodel[,sjournal:=factor(sjournal)][,sjournal:=relevel(sjournal,"Other")]

# Evolution of SJR over time for the top 20 journals

xjevo <- xbmodel[ sjournal != "Other", .(year_inprint, journal, SJR)] %>% unique


jevop <- ggplot(xjevo, aes(x=year_inprint, y=SJR, group=journal, colour=journal))+
        geom_point()+
        geom_line()
jevop


# Model fitting

# We'll now implement the same model building process as we did for the previous model

## start building a model without Public_ss
mb0=glm(y~year, data=xbmodel)                 # m0: year only, as numeric
mb1=glm(y~factor(year), data=xbmodel)          # m1: year only, as factor
BIC(mb1,mb0) # relationship with year better modelled as linear

mb2=glm(y~year + sjournal, data=xbmodel)      # m2: year and journal, using top20 and "other"
mb2a=glm(y~year + journal, data=xbmodel)      # m2a: year and journals, without pooling
BIC(mb0,mb2,mb2a) # sjournal is useful (m2), but journal has too many levels to be useful

mb3=glm(y~year + log(SJR), data=xbmodel)              # m3: year and log journal impact factor
mb4=glm(y~year + log(SJR) + journal, data=xbmodel)    # m4: year, log journal impact factor, and journals, without pooling
mb5=glm(y~year + log(SJR) + sjournal, data=xbmodel)   # m5: year, log journal impact factor, and journal, using top20 and "other"
BIC(mb2,mb3,mb4,mb5) # SJR has bigger impact than sjournal. Interestingly, this time mb3 has lower BIC than mb5, meaning that we might be better off without sjournal in our model.
#      df      BIC
# mb2  23 13921.04
# mb3   4 12819.85
# mb4 618 16646.78
# mb5  24 12838.07

# The next "human", "molecular_cellular", and "animal" variables correspond to estimations of study relevance for translation human, molecular, and animal
# studies, from MeSH terms contained in the article.
mb6=glm(y~year + log(SJR) + human , data=xbmodel)               # m6: include the previous variables + human score
mb7=glm(y~year + log(SJR) + animal, data=xbmodel)               # m7: animal
mb8=glm(y~year + log(SJR) + human + animal, data=xbmodel)       # m8: human + animal
mb9=glm(y~year + log(SJR) + molecular_cellular, data=xbmodel)   # m9: molecular_cellular
mb10=glm(y~year + log(SJR) + human + molecular_cellular, data=xbmodel)          # m10 : human + molecular_cellular
mb11=glm(y~year + log(SJR) + human + animal + molecular_cellular, data=xbmodel) # m11: human + animal + molecular_cellular
mb12=glm(y~year + log(SJR) + is_clinical, data=xbmodel)                         # m12: is the study categorised clinical?
BIC(mb3, mb6, mb7, mb8, mb9, mb10, mb11, mb12) # adding molecular_cellular (mb9) improves fit somewhat
#      df      BIC
# mb3   4 12819.85
# mb6   5 12816.74
# mb7   5 12827.73
# mb8   6 12820.27
# mb9   5 12812.74
# mb10  6 12821.21
# mb11  7 12828.61
# mb12  5 12825.59

## check all vars still useful - yes
drop1(mb9)
# Single term deletions
# 
# Model:
#   y ~ year + log(SJR)
# Df Deviance   AIC
# <none>        3905.4 12794
# year      1   4017.7 12931
# log(SJR)  1   6890.1 15571

## go with m9
mb9a=glm(y~year + log(SJR) + molecular_cellular + Public_ss, data=xbmodel)# %>% summary()
BIC(mb9,mb9a)
# Adding Public_ss has a substantial effect
#      df      BIC
# mb9   5 12812.74
# mb9a  6 12632.08

summary(mb9a)
# 
# Call:
#   glm(formula = y ~ year + log(SJR) + molecular_cellular + Public_ss, 
#       data = xbmodel)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -4.0917  -0.5341   0.0208   0.5497   3.5935  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        119.115635   7.534203  15.810  < 2e-16 ***
#   year                -0.059398   0.003737 -15.895  < 2e-16 ***
#   log(SJR)             0.784345   0.013768  56.968  < 2e-16 ***
#   molecular_cellular  -0.229779   0.060344  -3.808 0.000142 ***
#   Public_ssY           0.641786   0.046238  13.880  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.7659339)
# 
# Null deviance: 7674.7  on 4894  degrees of freedom
# Residual deviance: 3745.4  on 4890  degrees of freedom
# AIC: 12593
# 
# Number of Fisher Scoring iterations: 2

## conclude sharing has highly significant effect (t=12.270, p<2e-16, coef=0.566643), this is slightly better than our earlier model!
summary(mb9a)$coeff["Public_ssY",,drop=FALSE]
cat("Ratio of RCR in publications that shared vs did not = ",
    exp(coef(mb9a)["Public_ssY"]),
    "\n95% CI: ",
    exp(coef(mb9a)["Public_ssY"]+ c(-1,1) * 1.96 * sqrt(vcov(mb9a)["Public_ssY","Public_ssY"])),
    "\n")
# Ratio of RCR in publications that shared vs did not =  1.899871 
# 95% CI:  1.735263 2.080093 

# The model now predicts larger effect of sharing on citations!

# Table of results
predictors = rownames(summary(mb9a)$coeff)
sigtableb1 <- data.table(preds=predictors,
                        exp_estimate = exp(coef(mb9a)),
                        low.ci = sapply(predictors, function(x){exp(coef(mb9a)[x]+ -1 * 1.96 * sqrt(vcov(mb9a)[x,x]))}),
                        hi.ci =  sapply(predictors, function(x){exp(coef(mb9a)[x]+ 1 * 1.96 * sqrt(vcov(mb9a)[x,x]))}))


#######################################################################################################################

#---- Sensitivity Analysis

# We'll now remove each variable from the original model and see who coefficients for sharing change
# For this, we'll use the original dataset (xmodel), with SJR year == 2021


ms1=glm(y~year + log(SJR) + sjournal + molecular_cellular + Public_ss, data=xmodel) # Original model
ms2=glm(y~year + log(SJR) + sjournal + Public_ss, data=xmodel)                      # No molcel
ms3=glm(y~year + log(SJR) + molecular_cellular + Public_ss, data=xmodel)            # No sjournal
ms4=glm(y~year + sjournal + molecular_cellular + Public_ss, data=xmodel)            # No SJR
ms5=glm(y~ log(SJR) + sjournal + molecular_cellular + Public_ss, data=xmodel)       # No year
BIC.ms <- BIC(ms1,ms2,ms3,ms4,ms5)
#    df      BIC
#ms1 26 12729.46
#ms2 25 12743.70
#ms3  6 12705.49
#ms4 25 13889.52
#ms5 25 13095.02

sens.table <- data.table(pred_removed=c("None", "molecular_cellular", "sjournal", "log(SJR)", "Year"),
                         exp_estimate = sapply(1:5, function(i) {exp(coef(get(paste0("ms", i)))["Public_ssY"])}),
                         low.ci = sapply(1:5, function(i){exp(coef(get(paste0("ms", i)))["Public_ssY"]+ -1 * 1.96 * sqrt(vcov(get(paste0("ms", i)))["Public_ssY","Public_ssY"]))}),
                         hi.ci  = sapply(1:5, function(i){exp(coef(get(paste0("ms", i)))["Public_ssY"]+ 1 * 1.96 * sqrt(vcov(get(paste0("ms", i)))["Public_ssY","Public_ssY"]))}),
                         df     = BIC.ms$df,
                         BIC    = BIC.ms$BIC)
#           pred_removed exp_estimate   low.ci    hi.ci df      BIC
# 1:               None     1.749873 1.598681 1.915363 26 12729.46
# 2: molecular_cellular     1.762342 1.609829 1.929303 25 12743.70
# 3:           sjournal     1.883189 1.721459 2.060114  6 12705.49
# 4:           log(SJR)     2.108419 1.905758 2.332631 25 13889.52
# 5:               Year     1.359551 1.242375 1.487778 25 13095.02


#------------ Question 3 - How does mislabelling evolve over time?

mlb <- read.xlsx("../tables/Tables_submission_bioRxiv.xlsx", sheetIndex = 1, startRow = 2) %>% data.table
mlb[, year:=as.numeric(sapply(strsplit(Date, "-"), `[`, 1))]
setnames(mlb, "Full.summary.statistics.available.", "Public_ss")
mlb[Public_ss %in% c("Yes, authors' website", "Yes, author's website", "Yes, authors' website, but not advertised in manuscript"), Public_ss:="Yes"]
mlb[Public_ss %in% c("Yes, authors' website, but link does not contain data", "Authors' website, but link broken"), Public_ss:="No"] # Simplify, coerce to no

ggplot(mlb, aes(x=year,fill=Public_ss)) +
  geom_histogram(binwidth=1,colour="grey",position="dodge") +
  scale_fill_seaborn("Shared outside GWAS catalog  ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes")) +
  labs(x="Publication year",y="Observations") +
  #scale_fill_discrete(guide = guide_legend(reverse=TRUE), breaks = rev(levels(x$Public_ss)), labels=c("Yes", "No")) +
  background_grid(major="y")+
  theme(legend.position=c(0.5,0.96),legend.justification=c(0,1),legend.direction = "horizontal")


xms <- x[!PMID %in% mlb$PMID, .(year, PMID, Public_ss)] # Remove sample from full dataset, x minus sample
cmp <- xms[Public_ss == "Y", .(ss_count_large = .N), by = year] # Counts of shared ss per year
mlbsum <- mlb[Public_ss == "Yes", .(ss_count_sample = .N), by = year]

cmp <- merge(cmp, mlbsum, by="year", all = TRUE)
cmp[is.na(ss_count_sample), ss_count_sample:=0]

cor(cmp$ss_count_large, cmp$ss_count_sample)
# 0.8456454



## Extract nonsharers from 2017 onwards

xj17 <- x[Public_ss =="N" & year >= 2017, .(PMID, doi, First_Author, Title, journal)]
xj17[, c("PubMed_url","doi_url"):=list(paste0("https://pubmed.ncbi.nlm.nih.gov/", PMID), paste0("https://doi.org/", doi))]

fwrite(xj17, "../data/Nonsharers_2017.tsv", sep="\t")


