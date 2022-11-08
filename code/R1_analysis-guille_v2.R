################################################################################
###         Analyses, figures, and tables for GWAS sharing project          ####
################################################################################

## This script contains the most recent analysis and figures for the GWAS sharing 
## project. In this new version, we incorporated SJR scores per year, when available, 
## then remade figures, tables, and performed new analyses.


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
x=file.path(d,"R1_Main_data_20221020.tsv") %>% fread()
dim(x)
# [1] 5756   34

# Check journals and SJR info
j.total <- x[ , journal] %>% unique
# 723 journals in total
j.with.sjr <- x[ !is.na(SJR), journal] %>% unique
# 691 journals with SJR data for at least one year
setdiff(j.total, j.with.sjr)


# Proportion of sharing/non-sharing studies
x[ , .N, by = "Public_ss"][, N/sum(N)][2] 
# 0.104934, or 10.5%

# Proportion of sharing/non-sharing studies in 2021
x[ year == 2021, .N, by = "Public_ss"][, N/sum(N)][2]
# 0.2093426, or 20.93%

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


#---- Load and process number of citations per year since publication, for each study.

y=file.path(d,"Citations_per_year_20220611.tsv") %>% fread()
dim(y)
y=merge(y, x[,.(PMID,Public_ss,pubyear=year)],by="PMID")
y[,years_since:=year - pubyear] # We used online publication year, not in-print as "pubyear", so no problem of discrepancies for being cited "before" being published in print.
summary(y$years_since)          # Some negative values due to errors in citation metadata.  We'll remove them.
table(y$years_since)            
c(sum(y$years_since), sum(y$years_since < 0)) # 87 errors out of 169210 retrieved citations. Some of them are true citations (eg. online early, earlier that online publication date?)
y=y[years_since >= -1] # allow -1 because so many - citations based on online early?!?
y[years_since<0, years_since:=0] # but truncate the 38 remaining with years_since=-1 to 0
y=y[order(PMID,years_since)]



#---- Exploratory Data Analysis



## Figure 1 : observations by publication year

ggplot(x, aes(x=year,fill=Public_ss)) +
  geom_histogram(binwidth=1,colour="grey",position="dodge") +
  scale_fill_seaborn("Shared Sumstats   ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes")) +
  labs(x="Publication year",y="Observations") +
  #scale_fill_discrete(guide = guide_legend(reverse=TRUE), breaks = rev(levels(x$Public_ss)), labels=c("Yes", "No")) +
  background_grid(major="y")+
  theme(legend.position=c(0.05,0.96),legend.justification=c(0,1),legend.direction = "horizontal")
#ggsave("../figures/R1_figure1.tiff",height=6,width=7,units="in", bg="white")
#ggsave("../figures/R1_figure1.png",height=6,width=7,units="in", bg="white")

# total observations increases with year, and sharing increases.
# Need to consider year as a potential confounder.


## Figure 2 : Relative citations have decreased for GWAS papers over the years, with a recent trend up for 2021

x[,year_ss:=paste(year,Public_ss)]

top=ggplot(x, aes(x=year,y=y)) + geom_boxplot(aes(group=year)) + geom_smooth() +
  background_grid(major="y") +
  labs(x="Publication year",y="log relative citation ratio")

bot=ggplot(x, aes(x=year,y=y,col=Public_ss)) +
  geom_boxplot(aes(group=year_ss)) +
  # geom_jitter(width=0.2,height=0) +
  geom_smooth() +
  scale_colour_seaborn("Shared Sumstats    ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes")) +
  background_grid(major="y") +
  theme(legend.position=c(0.08,1),legend.justification=c(0,1),legend.direction = "horizontal") +
  labs(x="Publication year",y="log relative citation ratio")
plot_grid(top,bot,nrow=1, labels = "AUTO")

#ggsave("../figures/R1_figure2.tiff",height=6,width=12,units="in", bg="white")
#ggsave("../figures/R1_figure2.png",height=6,width=12,units="in", bg="white")


## Figure 3 : Shared/unshared citation ratio over the years, by year of publication, and mean citation count of shared and unshared datasets over the years, by year of publication.

y_avg=y[,.(mean_count=mean(count),n=.N),by=c("year","years_since","Public_ss","pubyear")][order(years_since)]
y_avg=y_avg[pubyear >= 2010 & pubyear<2019 & year >=pubyear & year < 2022]
y_avg[,max_years:=max(years_since),by="pubyear"]
y_rel=dcast(y_avg, year + pubyear +  years_since ~ Public_ss, value.var="mean_count") %>% as.data.table
y_rel[,rel_mean_count:=Y/N]

top=ggplot(y_rel, aes(x=years_since, y=rel_mean_count)) +
  geom_path() +
  facet_wrap(~pubyear,scales="free_y") +
  scale_colour_seaborn("Shared") +
  scale_x_continuous(breaks=c(0,5,10)) +
  background_grid(major = "y") +
  labs(x="Years since publication",y="Mean citation count ratio, shared / unshared")
bot=ggplot(y_avg[], aes(x=years_since, y=mean_count, col=Public_ss)) +
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

#ggsave("../figures/R1_figure3.tiff",height=8,width=14,units="in", bg ="white")
#ggsave("../figures/R1_figure3.png",height=8,width=14,units="in", bg ="white")


#---- Model fitting

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


# ---- Building a logistic model for effects on sharing

## Build logistic linear models to see which factors affect sharing, if any

# We'll use y_ss to dichotomise the variable
ml0=glm(y_ss ~ factor(year), data=xmodel, family = "binomial") # ml0: is year better as linear or as a factor?
ml1=glm(y_ss ~ year, data=xmodel, family = "binomial")         # ml1: does time affect sharing?
BIC(ml0, ml1) # In this case, it's better modelled as linear
#     df      BIC
# ml0 14 2797.897
# ml1  2 2717.718

ml2=glm(y_ss ~ log(SJR), data=xmodel, family = "binomial")                   # ml2: what about IF?
ml3=glm(y_ss ~ year + log(SJR), data=xmodel, family = "binomial")            # ml3: year and IF
ml4=glm(y_ss ~ year + log(SJR) + sjournal, data=xmodel, family = "binomial") # ml4: add sjournal
BIC(ml1, ml2, ml3, ml4) # go with year + IF  (ml3), since sjournal does not contribute
#     df      BIC
# ml1  2 2717.718
# ml2  2 2889.515
# ml3  3 2400.554
# ml4 23 2448.544

# Let's add the MeSH terms
ml5=glm(y_ss ~ year + log(SJR) + human, data=xmodel, family = "binomial")              # ml5: human MeSH terms
ml6=glm(y_ss ~ year + log(SJR) + animal, data=xmodel, family = "binomial")             # ml6: animal MeSH terms
ml7=glm(y_ss ~ year + log(SJR) + molecular_cellular, data=xmodel, family = "binomial") # ml7: molecular_cellular MeSH terms
ml8=glm(y_ss ~ year + log(SJR) + human + animal + molecular_cellular, data=xmodel, family = "binomial") # ml8: MeSH terms all together
ml9=glm(y_ss ~ year + log(SJR) + human + molecular_cellular, data=xmodel, family = "binomial") # ml9: human and molecular_cellular MeSH terms
BIC(ml3, ml5, ml6, ml7, ml8, ml9) # In this case, doesn't help. We'll go with ml3
# df      BIC
# ml3  3 2400.554
# ml5  4 2409.011
# ml6  4 2409.045
# ml7  4 2408.926
# ml8  6 2425.558
# ml9  5 2417.402

drop1(ml3) # All variables useful
# Single term deletions
# 
# Model:
#   y_ss ~ year + log(SJR)
# Df Deviance    AIC
# <none>        2375.0 2381.0
# year      1   2872.5 2876.5
# log(SJR)  1   2700.7 2704.7


summary(ml3)
# Call:
#   glm(formula = y_ss ~ year + log(SJR), family = "binomial", data = xmodel)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.4721  -0.4309  -0.2581  -0.1400   3.2971  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -790.72856   42.69691  -18.52   <2e-16 ***
#   year           0.39024    0.02115   18.45   <2e-16 ***
#   log(SJR)       1.05953    0.06281   16.87   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 3004  on 4959  degrees of freedom
# Residual deviance: 2375  on 4957  degrees of freedom
# AIC: 2381
# 
# Number of Fisher Scoring iterations: 6

# Table of results for significant
predictors2 = rownames(summary(ml3)$coeff)
sigtable2 <- data.table(preds=predictors2,
                        OR = exp(coef(ml3)),
                        low.ci = sapply(predictors2, function(x){exp(coef(ml3)[x]+ -1 * 1.96 * sqrt(vcov(ml3)[x,x]))}),
                        hi.ci =  sapply(predictors2, function(x){exp(coef(ml3)[x]+ 1 * 1.96 * sqrt(vcov(ml3)[x,x]))}))

# Save Table S2 complement
fwrite(sigtable2, "../tables/R1_S2_complement.tsv", sep = "\t")



# ---- Building a linear model for effects on log(RCR)

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
sigtable1 <- data.table(preds=predictors,
                        exp_estimate = exp(coef(m9a)),
                        low.ci = sapply(predictors, function(x){exp(coef(m9a)[x]+ -1 * 1.96 * sqrt(vcov(m9a)[x,x]))}),
                        hi.ci =  sapply(predictors, function(x){exp(coef(m9a)[x]+ 1 * 1.96 * sqrt(vcov(m9a)[x,x]))}))
# Save table S3 complement
fwrite(sigtable1, "../tables/R1_S3_complement.tsv", sep = "\t")


