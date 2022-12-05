########################################################
######            Response to reviewers           ######
########################################################

# Date: 2022-10-21

# Here we'll address reviewer questions about our analyses
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



#------------ Question 1 - Is molcel important? How much does the result change if we remove it?


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

# Compare both models with Public_ss
m9a=glm(y~year + log(SJR) + sjournal + molecular_cellular + Public_ss, data=xmodel) # Model used in the first version of the paper.

BIC.59a <- BIC(m9a, m5a)
BIC.59a
#     df      BIC
# m5a 25 12743.70
# m9a 26 12729.46

comparison.table <- data.table(pred_removed=c("Original", "without_molecular_cellular"),
                         exp_estimate = sapply(c(9,5), function(i) {exp(coef(get(paste0("m", i, "a")))["Public_ssY"])}),
                         low.ci = sapply(c(9,5), function(i){exp(coef(get(paste0("m", i, "a")))["Public_ssY"]+ -1 * 1.96 * sqrt(vcov(get(paste0("m", i, "a")))["Public_ssY","Public_ssY"]))}),
                         hi.ci  = sapply(c(9,5), function(i){exp(coef(get(paste0("m", i, "a")))["Public_ssY"]+ 1 * 1.96 * sqrt(vcov(get(paste0("m", i, "a")))["Public_ssY","Public_ssY"]))}),
                         df     = BIC.59a$df,
                         BIC    = BIC.59a$BIC)
#                  pred_removed exp_estimate   low.ci    hi.ci df      BIC
# 1:                   Original     1.749873 1.598681 1.915363 26 12729.46
# 2: without_molecular_cellular     1.762342 1.609829 1.929303 25 12743.70


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

# The model now predicts slightly larger effect on sharing on citations!

# Table of results
predictors = rownames(summary(m5a)$coeff)
data.table(preds=predictors,
                        exp_estimate = exp(coef(m5a)),
                        low.ci = sapply(predictors, function(x){exp(coef(m5a)[x]+ -1 * 1.96 * sqrt(vcov(m5a)[x,x]))}),
                        hi.ci =  sapply(predictors, function(x){exp(coef(m5a)[x]+ 1 * 1.96 * sqrt(vcov(m5a)[x,x]))}))
#                                            preds exp_estimate       low.ci        hi.ci
# 1:                                    (Intercept) 1.538251e+64 4.915496e+57 4.813790e+70
# 2:                                           year 9.291015e-01 9.222325e-01 9.360217e-01
# 3:                                       log(SJR) 2.075379e+00 1.994590e+00 2.159441e+00
# 4:                         sjournalAm J Hum Genet 1.028141e+00 8.669615e-01 1.219287e+00
# 5:  sjournalAm J Med Genet B Neuropsychiatr Genet 1.407585e+00 1.114216e+00 1.778196e+00
# 6:                          sjournalAnn Rheum Dis 9.216427e-01 6.945686e-01 1.222954e+00
# 7:                               sjournalDiabetes 1.178786e+00 9.179559e-01 1.513729e+00
# 8:                        sjournalEur J Hum Genet 1.040895e+00 8.025794e-01 1.349976e+00
# 9:                            sjournalFront Genet 1.010283e+00 7.657168e-01 1.332963e+00
# 10:                      sjournalGastroenterology 1.090936e+00 7.875457e-01 1.511202e+00
# 11:                         sjournalHum Mol Genet 1.242097e+00 1.113699e+00 1.385297e+00
# 12:                sjournalJ Allergy Clin Immunol 7.804677e-01 6.059621e-01 1.005228e+00
# 13:                           sjournalJ Hum Genet 8.183712e-01 6.340297e-01 1.056309e+00
# 14:                           sjournalJ Med Genet 8.628590e-01 6.293170e-01 1.183069e+00
# 15:                        sjournalMol Psychiatry 1.447421e+00 1.238492e+00 1.691596e+00
# 16:                            sjournalNat Commun 1.271692e+00 1.127848e+00 1.433882e+00
# 17:                             sjournalNat Genet 1.074218e+00 9.555716e-01 1.207596e+00
# 18:                                sjournalNature 2.085474e+00 1.636077e+00 2.658312e+00
# 19:                            sjournalPLoS Genet 1.519869e+00 1.345823e+00 1.716423e+00
# 20:                              sjournalPLoS One 1.176759e+00 1.049850e+00 1.319010e+00
# 21:              sjournalProc Natl Acad Sci U S A 1.547145e+00 1.137499e+00 2.104316e+00
# 22:                               sjournalSci Rep 1.171286e+00 1.005414e+00 1.364524e+00
# 23:                     sjournalTransl Psychiatry 1.158646e+00 9.436634e-01 1.422606e+00
# 24:                                    Public_ssY 1.762342e+00 1.609829e+00 1.929303e+00


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

### Load new data source, which includes updated per-year SJR
x=file.path(d,"R1_Main_data_20221020.tsv") %>% fread()
dim(x)
# [1] 5756   34

# Check journals and SJR info
j.total <- x[ , journal] %>% unique
# 723 journals in total
j.with.sjr <- x[ !is.na(SJR), journal] %>% unique
# 691 journals with SJR data for at least one year
setdiff(j.total, j.with.sjr)

## exclude 2022 as incomplete and 2005-2006 and before because no shared studies
with(x,table(year, Public_ss))
x=x[year < 2022 & year > 2006]

## Visualise RCR and offset it by minimum non-zero value
summary(x$relative_citation_ratio)
hist((x$relative_citation_ratio))
hist(log(x$relative_citation_ratio)) # but missing values for 0
## offset by minimum observed non-zero value
mincit=min(x[!is.na(relative_citation_ratio) & relative_citation_ratio > 0 , relative_citation_ratio])
x[relative_citation_ratio==0, relative_citation_ratio:=relative_citation_ratio+mincit]

# Create our response variable: log(RCR)
x[,y:=log(relative_citation_ratio)] # Our Response Variable
hist(x$y,breaks=100)
x[,y_ss:=Public_ss=="Y"]

# Prepare dataset for modelling, we'll remove missing data for relevant variables. 
# We'll also consider years 2007-2020 since it takes 2-3 years for sharing effect to stabilise and not enough time has passed since 2021

xmodel <- x[!is.na(SJR) & !is.na(y) & year<=2020]

## Add sjournal variable, which contains individual categories for each of the 
# 20 journals that share the most GWAS (with at least one sharing study), plus "Other" for the rest.
xjournal <- xmodel[, .(articles = .N, shared_ss = sum(y_ss)), by= "journal"][order(articles, decreasing = TRUE)]

# Save Table S4
#fwrite(xjournal, "../tables/Journals_table.tsv", sep="\t")

jkeep <- xjournal[shared_ss > 0][1:20, journal]
xmodel[,sjournal:=ifelse(journal %in% jkeep,journal,"Other")]             # Use top 20 most common journals with at least one shared dataset, pool the rest as "other"
xmodel[,sjournal:=factor(sjournal)][,sjournal:=relevel(sjournal,"Other")]

# Evolution of SJR over time for the top 20 journals

xjevo <- xmodel[ sjournal != "Other", .(year, journal, SJR)] %>% unique

jevop <- ggplot(xjevo, aes(x=year, y=log(SJR), group=journal, colour=journal))+
  geom_point()+
  geom_line()
jevop



## Start building a model without Public_ss
m0=glm(y~ year, data=xmodel)                       # m0: year only, as numeric
m1=glm(y~ factor(year), data=xmodel)               # m1: year only, as factor
BIC(m1,m0) # relationship with year better modelled as **factor**
#    df      BIC
# m1 15 15756.77
# m0  3 15766.09

m2=glm(y~factor(year) + sjournal, data=xmodel)    # m2: year and journal, using top20 and "other"
m2a=glm(y~factor(year) + journal, data=xmodel)            # m2a: year and journals, without pooling
BIC(m1,m2,m2a) # sjournal is useful (m2), but journal has too many levels to be useful
#      df      BIC
# m1   15 15756.77
# m2   35 14158.46
# m2a 638 16927.22

m3=glm(y~factor(year) + log(SJR), data=xmodel)              # m3: year and log journal impact factor
m4=glm(y~factor(year) + log(SJR) + journal, data=xmodel)    # m4: year, log journal impact factor, and journals, without pooling
m5=glm(y~factor(year) + log(SJR) + sjournal, data=xmodel)   # m5: year, log journal impact factor, and journal, using top20 and "other"
BIC(m2,m3,m4,m5) # SJR has bigger impact than sjournal, and after adding log(SJR), we no longer need sjournal in the model. We'll go with m3
#     df      BIC
# m2  35 14158.46
# m3  16 12952.25
# m4 639 16913.02
# m5  36 12984.40

# The next "human", "molecular_cellular", and "animal" variables correspond to estimations of study relevance for translation human, molecular, and animal
# studies, from MeSH terms contained in the article.
m6=glm(y~factor(year) + log(SJR)  + human , data=xmodel)                              # m6: include the previous variables + human score
m7=glm(y~factor(year) + log(SJR)  + animal, data=xmodel)                              # m7: animal
m8=glm(y~factor(year) + log(SJR)  + human + animal, data=xmodel)                      # m8: human + animal
m9=glm(y~factor(year) + log(SJR)  + molecular_cellular, data=xmodel)                  # m9: molecular_cellular
m10=glm(y~factor(year) + log(SJR) + human + molecular_cellular, data=xmodel)          # m10 : human + molecular_cellular
m11=glm(y~factor(year) + log(SJR) + human + animal + molecular_cellular, data=xmodel) # m11: human + animal + molecular_cellular
m12=glm(y~factor(year) + log(SJR) + is_clinical, data=xmodel)                         # m12: is the study categorised clinical?
BIC(m3, m6, m7, m8, m9, m10, m11, m12) # molecular_cellular seems to improve the fit a bit. We'll go with m9
#     df      BIC
# m3  16 12952.25
# m6  17 12947.20
# m7  17 12960.31
# m8  18 12950.75
# m9  17 12943.53
# m10 18 12952.04
# m11 19 12959.21
# m12 17 12959.43


## check all vars still useful - yes
drop1(m9)
# Single term deletions
# 
# Model:
#   y ~ factor(year) + log(SJR) + molecular_cellular
#                    Df Deviance   AIC
# <none>                  3834.1 12833
# factor(year)       13   4056.2 13086
# log(SJR)            1   6774.5 15654
# molecular_cellular  1   3847.5 12848


## go with m9
m9a=glm(y~factor(year) + log(SJR) + molecular_cellular + Public_ss, data=xmodel) # %>% summary()
BIC(m9,m9a)
#     df      BIC
# m9  17 12943.53
# m9a 18 12775.37
# Sharing has a powerful effect on BIC


summary(m9a)
# Call:
#   glm(formula = y ~ factor(year) + log(SJR) + molecular_cellular + 
#         Public_ss, data = xmodel)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -4.0642  -0.5207   0.0214   0.5456   3.7389  
# 
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)         0.59044    0.09498   6.216 5.51e-10 ***
#   factor(year)2008   -0.54834    0.11844  -4.630 3.76e-06 ***
#   factor(year)2009   -0.67876    0.10843  -6.260 4.17e-10 ***
#   factor(year)2010   -0.86225    0.10409  -8.284  < 2e-16 ***
#   factor(year)2011   -1.03540    0.10216 -10.135  < 2e-16 ***
#   factor(year)2012   -1.06082    0.10222 -10.378  < 2e-16 ***
#   factor(year)2013   -1.12082    0.10232 -10.955  < 2e-16 ***
#   factor(year)2014   -1.23359    0.10339 -11.932  < 2e-16 ***
#   factor(year)2015   -1.21981    0.10312 -11.829  < 2e-16 ***
#   factor(year)2016   -1.25569    0.10344 -12.140  < 2e-16 ***
#   factor(year)2017   -1.27408    0.10153 -12.549  < 2e-16 ***
#   factor(year)2018   -1.27358    0.10019 -12.711  < 2e-16 ***
#   factor(year)2019   -1.28576    0.10014 -12.839  < 2e-16 ***
#   factor(year)2020   -1.45259    0.10193 -14.251  < 2e-16 ***
#   log(SJR)            0.78147    0.01359  57.485  < 2e-16 ***
#   molecular_cellular -0.23729    0.05935  -3.998 6.48e-05 ***
#   Public_ssY          0.61184    0.04570  13.388  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.7485294)
# 
# Null deviance: 7721  on 4959  degrees of freedom
# Residual deviance: 3700  on 4943  degrees of freedom
# AIC: 12658
# 
# Number of Fisher Scoring iterations: 2

## conclude sharing has highly significant effect (t=12.131, p<2e-16, coef=0.559133)
summary(m9a)$coeff["Public_ssY",,drop=FALSE]
cat("Ratio of RCR in publications that shared vs did not = ",
    exp(coef(m9a)["Public_ssY"]),
    "\n95% CI: ",
    exp(coef(m9a)["Public_ssY"]+ c(-1,1) * 1.96 * sqrt(vcov(m9a)["Public_ssY","Public_ssY"])),
    "\n")
# Ratio of RCR in publications that shared vs did not =  1.843827 
# 95% CI:  1.685849 2.016609

# Table of results
predictors = rownames(summary(m9a)$coeff)
data.table(preds=predictors,
                        exp_estimate = exp(coef(m9a)),
                        low.ci = sapply(predictors, function(x){exp(coef(m9a)[x]+ -1 * 1.96 * sqrt(vcov(m9a)[x,x]))}),
                        hi.ci =  sapply(predictors, function(x){exp(coef(m9a)[x]+ 1 * 1.96 * sqrt(vcov(m9a)[x,x]))}))
#                  preds exp_estimate    low.ci     hi.ci
# 1:         (Intercept)    1.8047733 1.4982037 2.1740747
# 2:    factor(year)2008    0.5779067 0.4581833 0.7289139
# 3:    factor(year)2009    0.5072474 0.4101330 0.6273573
# 4:    factor(year)2010    0.4222111 0.3442932 0.5177630
# 5:    factor(year)2011    0.3550832 0.2906479 0.4338035
# 6:    factor(year)2012    0.3461725 0.2833225 0.4229645
# 7:    factor(year)2013    0.3260107 0.2667715 0.3984047
# 8:    factor(year)2014    0.2912443 0.2378230 0.3566653
# 9:    factor(year)2015    0.2952854 0.2412463 0.3614291
# 10:   factor(year)2016    0.2848781 0.2326016 0.3489035
# 11:   factor(year)2017    0.2796868 0.2292180 0.3412678
# 12:   factor(year)2018    0.2798280 0.2299346 0.3405476
# 13:   factor(year)2019    0.2764397 0.2271727 0.3363912
# 14:   factor(year)2020    0.2339637 0.1915950 0.2857016
# 15:           log(SJR)    2.1846844 2.1272420 2.2436780
# 16: molecular_cellular    0.7887621 0.7021386 0.8860724
# 17:         Public_ssY    1.8438267 1.6858487 2.0166085




#######################################################################################################################

#---- Sensitivity Analysis

# We'll now remove each variable from the original model and see who coefficients for sharing change

ms1=glm(y~factor(year) + log(SJR) + molecular_cellular + Public_ss, data=xmodel) # Original model
ms2=glm(y~factor(year) + log(SJR) + Public_ss, data=xmodel)                      # No molcel
ms3=glm(y~factor(year) + molecular_cellular + Public_ss, data=xmodel)            # No SJR
ms4=glm(y~ log(SJR) + molecular_cellular + Public_ss, data=xmodel)               # No year

BIC.ms <- BIC(ms1,ms2,ms3,ms4)
#     df      BIC
# ms1 18 12775.37
# ms2 17 12782.87
# ms3 17 15306.06
# ms4  5 13035.97
# Dropping any variable will have a negative effect on BIC

sens.table <- data.table(pred_removed=c("None", "molecular_cellular", "log(SJR)", "Year"),
                         exp_estimate = sapply(1:4, function(i) {exp(coef(get(paste0("ms", i)))["Public_ssY"])}),
                         low.ci = sapply(1:4, function(i){exp(coef(get(paste0("ms", i)))["Public_ssY"]+ -1 * 1.96 * sqrt(vcov(get(paste0("ms", i)))["Public_ssY","Public_ssY"]))}),
                         hi.ci  = sapply(1:4, function(i){exp(coef(get(paste0("ms", i)))["Public_ssY"]+ 1 * 1.96 * sqrt(vcov(get(paste0("ms", i)))["Public_ssY","Public_ssY"]))}),
                         df     = BIC.ms$df,
                         BIC    = BIC.ms$BIC)
#          pred_removed exp_estimate   low.ci    hi.ci df      BIC
# 1:               None     1.843827 1.685849 2.016609 18 12775.37
# 2: molecular_cellular     1.849471 1.690801 2.023031 17 12782.87
# 3:           log(SJR)     3.509796 3.137394 3.926401 17 15306.06
# 4:               Year     1.517310 1.388877 1.657620  5 13035.97

# Removing molecular_cellular would mean an increase in effect, at the cost of a poorer fit.



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

xj17 <- x[Public_ss =="N" & year >= 2017, .(PMID, doi, year, First_Author, Title, journal)]
xj17[, c("PubMed_url","doi_url"):=list(paste0("https://pubmed.ncbi.nlm.nih.gov/", PMID), paste0("https://doi.org/", doi))]

fwrite(xj17, "../data/Nonsharers_2017.tsv", sep="\t")

# Extract all nonsharers

xj.all <- x[Public_ss =="N", .(PMID, doi, year, first_author, title, journal)]
xj.all[, c("PubMed_url","doi_url"):=list(paste0("https://pubmed.ncbi.nlm.nih.gov/", PMID), paste0("https://doi.org/", doi))]

fwrite(xj.all, "../data/Nonsharers_all.tsv", sep="\t")
