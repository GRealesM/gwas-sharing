################################################################################
###         Analyses, figures, and tables for GWAS sharing project          ####
################################################################################

## This script contains the most recent analysis and figures for the GWAS sharing 
## project. In this new version, we build onto the previous version with 
## incorporated SJR scores per year, but with added classification of 2017 onwards
## papers. Then remade figures, tables, and analyses.


#---- Load Packages

library(data.table)
library(magrittr)
library(ggplot2)
# devtools::install_github("chr1swallace/seaborn_palettes")
library(seaborn)
library(cowplot)
theme_set(theme_cowplot())

#---- Load data, extract journal statistics and process it for further analyses

d="../data" #"~/rds/rds-cew54-basis/Projects/gwas-sharing/data/"
list.files(d)

### Load main data source, which includes studies in GWAS catalog, whether they share summary statistics or not, journal information, etc.
x=fread(file.path(d,"R1_Main_data_20221128.tsv"))
dim(x)
# [1] 5756   35 # Just an additional column for the new classification (Sharer_17)

# Check journals and SJR info
j.total <- x[ , journal] %>% unique
# 723 journals in total
j.with.sjr <- x[ !is.na(SJR), journal] %>% unique
# 691 journals with SJR data for at least one year
setdiff(j.total, j.with.sjr)

# Proportion of sharing/non-sharing studies following new classification
x[ , .N, by = "Sharer_17"][, N/sum(N)][2] 
# 0.1370744, or 13.7%
# Up from 10.5% GWAS catalog

# Proportion of sharing/non-sharing studies in 2021
x[ year == 2021, .N, by = "Sharer_17"][, N/sum(N)][2]
# 0.283737, or 28.73%
# Up from 20.93% GWAS catalog

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
x[,y_ss:=Sharer_17=="Y"]


#---- Load and process number of citations per year since publication, for each study.

y=file.path(d,"Citations_per_year_20220611.tsv") %>% fread()
dim(y)
y=merge(y, x[,.(PMID,Sharer_17,pubyear=year)],by="PMID")
y[,years_since:=year - pubyear] # We used online publication year, not in-print as "pubyear", so no problem of discrepancies for being cited "before" being published in print.
summary(y$years_since)          # Some negative values due to errors in citation metadata.  We'll remove them.
table(y$years_since)            
c(sum(y$years_since), sum(y$years_since < 0)) # 87 errors out of 169210 retrieved citations. Some of them are true citations (eg. online early, earlier that online publication date?)
y=y[years_since >= -1] # allow -1 because so many - citations based on online early?!?
y[years_since<0, years_since:=0] # but truncate the 38 remaining with years_since=-1 to 0
y=y[order(PMID,years_since)]



#---- Exploratory Data Analysis



## Figure 1 : observations by publication year

ggplot(x, aes(x=year,fill=Sharer_17)) +
  geom_histogram(binwidth=1,colour="grey",position="dodge") +
  scale_fill_seaborn("Shared Sumstats   ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes")) +
  labs(x="Publication year",y="Observations") +
  #scale_fill_discrete(guide = guide_legend(reverse=TRUE), breaks = rev(levels(x$Public_ss)), labels=c("Yes", "No")) +
  background_grid(major="y")+
  theme(legend.position=c(0.05,0.96),legend.justification=c(0,1),legend.direction = "horizontal")
# ggsave("../figures/R1.1_figure1.tiff",height=6,width=7,units="in", bg="white")
# ggsave("../figures/R1.1_figure1.png",height=6,width=7,units="in", bg="white")

# total observations increases with year, and sharing increases.
# Need to consider year as a potential confounder.


## Figure 2 : Relative citations have decreased for GWAS papers over the years, with a recent trend up for 2021

x[,year_ss:=paste(year,Sharer_17)]

top=ggplot(x, aes(x=year,y=y)) + geom_boxplot(aes(group=year)) + geom_smooth() +
  background_grid(major="y") +
  labs(x="Publication year",y="log relative citation ratio")

bot=ggplot(x, aes(x=year,y=y,col=Sharer_17)) +
  geom_boxplot(aes(group=year_ss)) +
  # geom_jitter(width=0.2,height=0) +
  geom_smooth() +
  scale_colour_seaborn("Shared Sumstats    ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes")) +
  background_grid(major="y") +
  theme(legend.position=c(0.08,1),legend.justification=c(0,1),legend.direction = "horizontal") +
  labs(x="Publication year",y="log relative citation ratio")
plot_grid(top,bot,nrow=1, labels = "AUTO")

# ggsave("../figures/R1.1_figure2.tiff",height=6,width=12,units="in", bg="white")
# ggsave("../figures/R1.1_figure2.png",height=6,width=12,units="in", bg="white")


## Figure 3 : Shared/unshared citation ratio over the years, by year of publication, and mean citation count of shared and unshared datasets over the years, by year of publication.

y_avg=y[,.(mean_count=mean(count),n=.N),by=c("year","years_since","Sharer_17","pubyear")][order(years_since)]
y_avg=y_avg[pubyear >= 2010 & pubyear<2019 & year >=pubyear & year < 2022]
y_avg[,max_years:=max(years_since),by="pubyear"]
y_rel=dcast(y_avg, year + pubyear +  years_since ~ Sharer_17, value.var="mean_count") %>% as.data.table
y_rel[,rel_mean_count:=Y/N]

top=ggplot(y_rel, aes(x=years_since, y=rel_mean_count)) +
  geom_path() +
  facet_wrap(~pubyear,scales="free_y") +
  scale_colour_seaborn("Shared") +
  scale_x_continuous(breaks=c(0,5,10)) +
  background_grid(major = "y") +
  labs(x="Years since publication",y="Mean citation count ratio, shared / unshared")
bot=ggplot(y_avg[], aes(x=years_since, y=mean_count, col=Sharer_17)) +
  geom_path() +
  geom_label(aes(label=n),data=y_avg[years_since==max_years & pubyear<2019 & pubyear>=2010]) +
  facet_wrap(~pubyear,scales="free_y") +
  scale_colour_seaborn("Shared  ",guide = guide_legend(reverse=TRUE), labels=c("No", "Yes")) +
  scale_x_continuous(breaks=c(0,5,10), limits=c(0,12)) +
  background_grid(major = "y") +
  theme(legend.position=c(0.84,0.24),legend.justification=c(0,1), legend.background = element_rect(fill="white")) +
  ## scale_y_log10() +
  labs(x="Years since publication",y="Mean citation count")
plot_grid(top,bot,nrow=1, labels = "AUTO")

# ggsave("../figures/R1.1_figure3.tiff",height=8,width=14,units="in", bg ="white")
# ggsave("../figures/R1.1_figure3.png",height=8,width=14,units="in", bg ="white")


#---- Model fitting

# Prepare dataset for modelling, we'll remove missing data for relevant variables. 
# We'll also consider years 2007-2020 since it takes 2-3 years for sharing effect to stabilise and not enough time has passed since 2021

xmodel <- x[!is.na(SJR) & !is.na(y) & year<=2020]

## Add sjournal variable, which contains individual categories for each of the 
# 20 journals that share the most GWAS (with at least one sharing study), plus "Other" for the rest.
xjournal <- xmodel[, .(articles = .N, shared_ss = sum(y_ss)), by= "journal"][order(articles, decreasing = TRUE)]

# Save Table S4
# fwrite(xjournal, "../tables/R1.1_Journals_table.tsv", sep="\t")

jkeep <- xjournal[shared_ss > 0][1:20, journal]
xmodel[,sjournal:=ifelse(journal %in% jkeep,journal,"Other")]             # Use top 20 most common journals with at least one shared dataset, pool the rest as "other"
xmodel[,sjournal:=factor(sjournal)][,sjournal:=relevel(sjournal,"Other")]


# ---- Building a logistic model for effects on sharing

## Build logistic linear models to see which factors affect sharing, if any

# We'll use y_ss to dichotomise the variable
ml0=glm(y_ss ~ factor(year), data=xmodel, family = "binomial") # ml0: is year better as linear or as a factor?
ml1=glm(y_ss ~ year, data=xmodel, family = "binomial")         # ml1: does time affect sharing?
BIC(ml0, ml1) # In this case, it's better modelled as linear (again)
#     df      BIC
# ml0 14 3171.931
# ml1  2 3102.029

ml2=glm(y_ss ~ log(SJR), data=xmodel, family = "binomial")                   # ml2: what about IF?
ml3=glm(y_ss ~ year + log(SJR), data=xmodel, family = "binomial")            # ml3: year and IF
ml4=glm(y_ss ~ year + log(SJR) + sjournal, data=xmodel, family = "binomial") # ml4: add sjournal
BIC(ml1, ml2, ml3, ml4) # go with year + IF  (ml3), since sjournal does not contribute
#     df      BIC
# ml1  2 3102.029
# ml2  2 3489.449
# ml3  3 2709.153
# ml4 23 2738.885

# Let's add the MeSH terms
ml5=glm(y_ss ~ year + log(SJR) + human, data=xmodel, family = "binomial")              # ml5: human MeSH terms
ml6=glm(y_ss ~ year + log(SJR) + animal, data=xmodel, family = "binomial")             # ml6: animal MeSH terms
ml7=glm(y_ss ~ year + log(SJR) + molecular_cellular, data=xmodel, family = "binomial") # ml7: molecular_cellular MeSH terms
ml8=glm(y_ss ~ year + log(SJR) + human + animal + molecular_cellular, data=xmodel, family = "binomial") # ml8: MeSH terms all together
ml9=glm(y_ss ~ year + log(SJR) + human + molecular_cellular, data=xmodel, family = "binomial") # ml9: human and molecular_cellular MeSH terms
BIC(ml3, ml5, ml6, ml7, ml8, ml9) # In this case, doesn't help. We'll go with ml3
#     df      BIC
# ml3  3 2709.153
# ml5  4 2717.025
# ml6  4 2717.339
# ml7  4 2716.318
# ml8  6 2730.434
# ml9  5 2724.677

drop1(ml3) # All variables useful
# Single term deletions
# 
# Model:
#   y_ss ~ year + log(SJR)
# Df Deviance    AIC
# <none>        2683.6 2689.6
# year      1   3472.4 3476.4
# log(SJR)  1   3085.0 3089.0


summary(ml3)
#  Call:
#   glm(formula = y_ss ~ year + log(SJR), family = "binomial", data = xmodel)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.7565  -0.4809  -0.2506  -0.1197   3.3238  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -945.56855   42.92717  -22.03   <2e-16 ***
#   year           0.46715    0.02126   21.98   <2e-16 ***
#   log(SJR)       1.09684    0.05949   18.44   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 3603.1  on 4959  degrees of freedom
# Residual deviance: 2683.6  on 4957  degrees of freedom
# AIC: 2689.6
# 
# Number of Fisher Scoring iterations: 6

# Table of results for significant
predictors2 = rownames(summary(ml3)$coeff)
sigtable2 <- data.table(preds=predictors2,
                        OR = exp(coef(ml3)),
                        low.ci = sapply(predictors2, function(x){exp(coef(ml3)[x]+ -1 * 1.96 * sqrt(vcov(ml3)[x,x]))}),
                        hi.ci =  sapply(predictors2, function(x){exp(coef(ml3)[x]+ 1 * 1.96 * sqrt(vcov(ml3)[x,x]))}))

# Save Table S2 complement
fwrite(sigtable2, "../tables/R1.1_S2_complement.tsv", sep = "\t")


  
# ---- Building a linear model for effects on log(RCR)

## Start building a model without Sharer_17. We expect it to be exactly the same as before, since we simple added Sharer_17
m0=glm(y~ year, data=xmodel)                       # m0: year only, as numeric
m1=glm(y~ factor(year), data=xmodel)               # m1: year only, as factor
BIC(m0,m1) # relationship with year better modelled as **factor**
# df      BIC
# m0  3 15766.09
# m1 15 15756.77

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
# m2  35 14136.93
# m3  16 12952.25
# m4 639 16913.02
# m5  36 12971.49

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
m9a=glm(y~factor(year) + log(SJR) + molecular_cellular + Sharer_17, data=xmodel) # %>% summary()
BIC(m9,m9a)
#     df      BIC
# m9  17 12943.53
# m9a 18 12739.34
# Sharing has a powerful effect on BIC


summary(m9a)
# Call:
#   glm(formula = y ~ factor(year) + log(SJR) + molecular_cellular + 
#         Sharer_17, data = xmodel)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -4.0192  -0.5138   0.0196   0.5485   3.7472  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.60380    0.09467   6.378 1.95e-10 ***
#   factor(year)2008   -0.54324    0.11801  -4.603 4.26e-06 ***
#   factor(year)2009   -0.67603    0.10804  -6.257 4.24e-10 ***
#   factor(year)2010   -0.85979    0.10371  -8.290  < 2e-16 ***
#   factor(year)2011   -1.03603    0.10179 -10.178  < 2e-16 ***
#   factor(year)2012   -1.06262    0.10185 -10.433  < 2e-16 ***
#   factor(year)2013   -1.12248    0.10194 -11.011  < 2e-16 ***
#   factor(year)2014   -1.23608    0.10301 -11.999  < 2e-16 ***
#   factor(year)2015   -1.22330    0.10275 -11.906  < 2e-16 ***
#   factor(year)2016   -1.25806    0.10304 -12.209  < 2e-16 ***
#   factor(year)2017   -1.30933    0.10129 -12.927  < 2e-16 ***
#   factor(year)2018   -1.31948    0.10007 -13.186  < 2e-16 ***
#   factor(year)2019   -1.34550    0.10015 -13.435  < 2e-16 ***
#   factor(year)2020   -1.50226    0.10185 -14.750  < 2e-16 ***
#   log(SJR)            0.77170    0.01364  56.568  < 2e-16 ***
#   molecular_cellular -0.22306    0.05916  -3.770 0.000165 ***
#   Sharer_17Y          0.61324    0.04167  14.716  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.7431127)
# 
# Null deviance: 7721.0  on 4959  degrees of freedom
# Residual deviance: 3673.2  on 4943  degrees of freedom
# AIC: 12622
# 
# Number of Fisher Scoring iterations: 2


## conclude sharing has highly significant effect (t=14.716, p<2e-16, coef=0.61324)
summary(m9a)$coeff["Sharer_17Y",,drop=FALSE]
cat("Ratio of RCR in publications that shared vs did not = ",
    exp(coef(m9a)["Sharer_17Y"]),
    "\n95% CI: ",
    exp(coef(m9a)["Sharer_17Y"]+ c(-1,1) * 1.96 * sqrt(vcov(m9a)["Sharer_17Y","Sharer_17Y"])),
    "\n")
# Ratio of RCR in publications that shared vs did not =  1.846397 
# 95% CI:  1.70159 2.003528 

# For R1 (ie. using Public_SS, or GWAS catalog classification) it was
# Ratio of RCR in publications that shared vs did not =  1.843827 
# 95% CI:  1.685849 2.016609

# Table of results
predictors = rownames(summary(m9a)$coeff)
sigtable1 <- data.table(preds=predictors,
                        exp_estimate = exp(coef(m9a)),
                        low.ci = sapply(predictors, function(x){exp(coef(m9a)[x]+ -1 * 1.96 * sqrt(vcov(m9a)[x,x]))}),
                        hi.ci =  sapply(predictors, function(x){exp(coef(m9a)[x]+ 1 * 1.96 * sqrt(vcov(m9a)[x,x]))}))
# Save table S3 complement
fwrite(sigtable1, "../tables/R1.1_S3_complement.tsv", sep = "\t")


