library(data.table)
library(magrittr)
library(ggplot2)
# devtools::install_github("chr1swallace/seaborn_palettes")
library(seaborn)
library(cowplot)
theme_set(theme_cowplot())

## Load and process data
d="../data" #"~/rds/rds-cew54-basis/Projects/gwas-sharing/data/"
list.files(d)

### Load main data source, which includes studies in GWAS catalog, whether they share summary statistics or not, journal information, etc.
x=file.path(d,"Main_data_20220610.tsv") %>% fread()
dim(x)

# Check journals and SJR info
x[ , .(Journal, SJR)] %>% unique %>% nrow
x[ !is.na(SJR), .(Journal, SJR)] %>% unique %>% nrow

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


# Load and process number of citations per year since publication, for each study.
y=file.path(d,"Citations_per_year_20220611.tsv") %>% fread()
dim(y)
y=merge(y, x[,.(PMID,Public_ss,pubyear=year)],by="PMID")
y[,years_since:=year - pubyear]
y=y[years_since >= -1] # allow -1 because so many - citations based on online early?!?
y[years_since<0, years_since:=0] # but truncate the 38 remaining with years_since=-1 to 0
y=y[order(PMID,years_since)]


### Exploratory Data Analysis

# Proportion of sharing/non-sharing studies in 2021
x[ year == 2021, .N, by = "Public_ss"][, N/sum(N)]

##' Figure: observations by publication year
ggplot(x, aes(x=year,fill=Public_ss)) +
  geom_histogram(binwidth=1,colour="grey",position="dodge") +
  scale_fill_seaborn("Shared Sumstats   ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes")) +
  labs(x="Publication year",y="Observations") +
  #scale_fill_discrete(guide = guide_legend(reverse=TRUE), breaks = rev(levels(x$Public_ss)), labels=c("Yes", "No")) +
  background_grid(major="y")+
  theme(legend.position=c(0.05,0.96),legend.justification=c(0,1),legend.direction = "horizontal")
#ggsave("../figures/fig-obs-by-pub-year.tiff",height=6,width=7,units="in", bg="white")
#ggsave("../figures/fig-obs-by-pub-year.png",height=6,width=7,units="in", bg="white")
## total observations increases with year, and sharing increases.
## Need to consider year as a potential confounder.

##' Figure: Relative citations have decreased for GWAS papers over the years, with a recent trend up for 2021
top=ggplot(x, aes(x=year,y=y)) + geom_boxplot(aes(group=year)) + geom_smooth() +
  background_grid(major="y") +
  labs(x="Publication year",y="log relative citation ratio")
x[,year_ss:=paste(year,Public_ss)]
bot=ggplot(x[year<=2022], aes(x=year,y=y,col=Public_ss)) +
  geom_boxplot(aes(group=year_ss)) +
 # geom_jitter(width=0.2,height=0) +
  geom_smooth() +
  scale_colour_seaborn("Shared Sumstats    ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes")) +
  background_grid(major="y") +
  theme(legend.position=c(0.08,1),legend.justification=c(0,1),legend.direction = "horizontal") +
  labs(x="Publication year",y="log relative citation ratio")
plot_grid(top,bot,nrow=1, labels = "AUTO")
#ggsave("../figures/fig-rel-cit-ratio-per-year.tiff",height=6,width=12,units="in", bg="white")
#ggsave("../figures/fig-rel-cit-ratio-per-year.png",height=6,width=12,units="in", bg="white")


## odd that the gap closes in the final year.
## look at citations per year since publication
## individual profiles too many to discern patterns
ggplot(y, aes(x=years_since, y=count, col=Public_ss,group=PMID)) + geom_path() + facet_wrap(~pubyear,scales="free_x")


## Figure: Shared/unshared citation ratio over the years, by year of publication, and mean citation count of shared and unshared datasets over the years, by year of publication.
y_avg=y[,.(mean_count=mean(count),n=.N),by=c("year","years_since","Public_ss","pubyear")][order(years_since)]
y_avg=y_avg[pubyear >= 2010 & pubyear<2019 & year >=pubyear & year < 2022]
y_avg[,max_years:=max(years_since),by="pubyear"]
y_rel=dcast(y_avg, year + pubyear +  years_since ~ Public_ss, value.var="mean_count")
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
#ggsave("../figures/fig-mean-citation-count-per-year.tiff",height=8,width=14,units="in", bg ="white")
#ggsave("../figures/fig-mean-citation-count-per-year.png",height=8,width=14,units="in", bg ="white")



### Model fitting

## As we saw in Figure 3, it takes 2-3 years for sharing effect to stabilise. Exclude 2021 because not enough time has passed
x=x[year<=2020]
glm(y ~ factor(year) + Public_ss, data=x) %>% summary()
## headline: exp(1.25)=3.49 - presumably some of this can be explained by other variables


# Prepare dataset for modelling, we'll remove missing data for relevant variables. We'll also consider years 2007-2020
xmodel <- x[!is.na(SJR) & !is.na(y)]

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


## start building a model without Public_ss
m0=glm(y~year, data=x)                       # m0: year only, as numeric
m1=glm(y~factor(year), data=xmodel)          # m1: year only, as factor
BIC(m1,m0) # relationship with year better modelled as linear

m2=glm(y~year + sjournal, data=xmodel)      # m2: year and journal, using top20 and "other"
m2a=glm(y~year + journal, data=xmodel)      # m2a: year and journals, without pooling
BIC(m0,m2,m2a) # sjournal is useful (m2), but journal has too many levels to be useful

m3=glm(y~year + log(SJR), data=xmodel)              # m3: year and log journal impact factor
m4=glm(y~year + log(SJR) + journal, data=xmodel)    # m4: year, log journal impact factor, and journals, without pooling
m5=glm(y~year + log(SJR) + sjournal, data=xmodel)   # m5: year, log journal impact factor, and journal, using top20 and "other"
BIC(m2,m3,m4,m5) # SJR has bigger impact than sjournal, but sjournal still useful in addition (m5)

# The next "human", "molecular_cellular", and "animal" variables correspond to estimations of study relevance for translation human, molecular, and animal
# studies, from MeSH terms contained in the article.
m6=glm(y~year + log(SJR) + sjournal + human , data=xmodel)               # m6: include the previous variables + human score
m7=glm(y~year + log(SJR) + sjournal + animal, data=xmodel)               # m7: animal
m8=glm(y~year + log(SJR) + sjournal + human + animal, data=xmodel)       # m8: human + animal
m9=glm(y~year + log(SJR) + sjournal + molecular_cellular, data=xmodel)   # m9: molecular_cellular
m10=glm(y~year + log(SJR) + sjournal + human + molecular_cellular, data=xmodel)          # m10 : human + molecular_cellular
m11=glm(y~year + log(SJR) + sjournal + human + animal + molecular_cellular, data=xmodel) # m11: human + animal + molecular_cellular
m12=glm(y~year + log(SJR) + sjournal + is_clinical, data=xmodel)                         # m12: is the study categorised clinical?
BIC(m5, m6, m7, m8, m9, m10, m11, m12) # adding human (m6) also helpful, but molecular_cellular (m9) has a larger effect

## check all vars still useful - yes
drop1(m9)
# Single term deletions
# 
# Model:
#   y ~ year + log(SJR) + sjournal + molecular_cellular
# Df Deviance   AIC
# <none>                  3745.0 12704
# year                1   3957.6 12975
# log(SJR)            1   4802.1 13931
# sjournal           20   3890.0 12852
# molecular_cellular  1   3764.7 12728


## go with m9
m9a=glm(y~year + log(SJR) + sjournal + molecular_cellular + Public_ss, data=xmodel)# %>% summary()
BIC(m9,m9a)

summary(m9a)
# Call:
#   glm(formula = y ~ year + log(SJR) + sjournal + molecular_cellular + 
#         Public_ss, data = xmodel)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -3.9377  -0.5169   0.0298   0.5434   3.4694  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)                                   149.477550   7.622213  19.611  < 2e-16 ***
#   year                                           -0.074357   0.003782 -19.663  < 2e-16 ***
#   log(SJR)                                        0.732220   0.020218  36.216  < 2e-16 ***
#   sjournalAm J Hum Genet                          0.038881   0.086837   0.448  0.65436    
#   sjournalAm J Med Genet B Neuropsychiatr Genet   0.325891   0.119032   2.738  0.00621 ** 
#   sjournalAnn Rheum Dis                          -0.074202   0.144011  -0.515  0.60640    
#   sjournalDiabetes                                0.161516   0.127319   1.269  0.20464    
#   sjournalEur J Hum Genet                         0.046050   0.132371   0.348  0.72794    
#   sjournalFront Genet                             0.076391   0.141787   0.539  0.59007    
#   sjournalGastroenterology                        0.081094   0.165899   0.489  0.62499    
#   sjournalHum Mol Genet                           0.227154   0.055591   4.086 4.46e-05 ***
#   sjournalJ Allergy Clin Immunol                 -0.252773   0.128841  -1.962  0.04983 *  
#   sjournalJ Hum Genet                            -0.181562   0.129989  -1.397  0.16255    
#   sjournalJ Med Genet                            -0.141193   0.160680  -0.879  0.37960    
#   sjournalMol Psychiatry                          0.354671   0.079424   4.466 8.17e-06 ***
#   sjournalNat Commun                              0.247920   0.061129   4.056 5.08e-05 ***
#   sjournalNat Genet                               0.093893   0.059766   1.571  0.11624    
#   sjournalNature                                  0.800134   0.124307   6.437 1.34e-10 ***
#   sjournalPLoS Genet                              0.432826   0.061986   6.983 3.28e-12 ***
#   sjournalPLoS One                                0.162814   0.058095   2.803  0.00509 ** 
#   sjournalProc Natl Acad Sci U S A                0.481150   0.156865   3.067  0.00217 ** 
#   sjournalSci Rep                                 0.160260   0.077740   2.061  0.03931 *  
#   sjournalTransl Psychiatry                       0.134156   0.104519   1.284  0.19936    
#   molecular_cellular                             -0.282608   0.059347  -4.762 1.97e-06 ***
#   Public_ssY                                      0.559543   0.046104  12.137  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.7394993)
# 
# Null deviance: 7793.4  on 4941  degrees of freedom
# Residual deviance: 3636.1  on 4917  degrees of freedom
# AIC: 12560
# 
# Number of Fisher Scoring iterations: 2

## conclude sharing has highly significant effect (t=12.131, p<2e-16, coef=0.559133)
summary(m9a)$coeff["Public_ssY",,drop=FALSE]
cat("Ratio of RCR in publications that shared vs did not = ",
exp(coef(m9a)["Public_ssY"]),
"\n95% CI: ",
exp(coef(m9a)["Public_ssY"]+ c(-1,1) * 1.96 * sqrt(vcov(m9a)["Public_ssY","Public_ssY"])),
"\n")
# Ratio of RCR in publications that shared vs did not =  1.749873 
# 95% CI:  1.598681 1.915363 

# Table of results
predictors = rownames(summary(m9a)$coeff)
sigtable1 <- data.table(preds=predictors,
                       OR = exp(coef(m9a)),
                       low.ci = sapply(predictors, function(x){exp(coef(m9a)[x]+ -1 * 1.96 * sqrt(vcov(m9a)[x,x]))}),
                       hi.ci =  sapply(predictors, function(x){exp(coef(m9a)[x]+ 1 * 1.96 * sqrt(vcov(m9a)[x,x]))}))

#  fwrite(sigtable1, "../tables/Sigtable_m9a.tsv", sep = "\t")

 
## Build logistic linear models to see which factors affect sharing, if any

# We'll use y_ss to dichotomise the variable?
ml0=glm(y_ss ~ factor(year), data=xmodel, family = "binomial") # ml0: is year better as linear or as a factor?
ml1=glm(y_ss ~ year, data=xmodel, family = "binomial") # ml1: does time affect sharing?
BIC(ml0, ml1) # Better modelled as linear
ml2=glm(y_ss ~ log(SJR), data=xmodel, family = "binomial") # ml2: what about IF?
ml3=glm(y_ss ~ year + log(SJR), data=xmodel, family = "binomial") # ml3: year and IF
ml4=glm(y_ss ~ year + log(SJR) + sjournal, data=xmodel, family = "binomial") # ml4: add sjournal
BIC(ml1, ml2, ml3, ml4) # go with year + IF + sjournal (ml4)

# Let's add the MeSH terms
ml5=glm(y_ss ~ year + log(SJR) + sjournal + human, data=xmodel, family = "binomial") # ml5: human MeSH terms
ml6=glm(y_ss ~ year + log(SJR) + sjournal + animal, data=xmodel, family = "binomial") # ml6: animal MeSH terms
ml7=glm(y_ss ~ year + log(SJR) + sjournal + molecular_cellular, data=xmodel, family = "binomial") # ml7: molecular_cellular MeSH terms
ml8=glm(y_ss ~ year + log(SJR) + sjournal + human + animal + molecular_cellular, data=xmodel, family = "binomial") # ml8: MeSH terms all together
ml9=glm(y_ss ~ year + log(SJR) + sjournal + human + molecular_cellular, data=xmodel, family = "binomial") # ml9: human and molecular_cellular MeSH terms
BIC(ml4, ml5, ml6, ml7, ml8, ml9) # Adding mesh terms doesn't really help. let's keep ml4

drop1(ml4) # All variables still important
# Single term deletions
# 
# Model:
#   y_ss ~ year + log(SJR) + sjournal
# Df Deviance    AIC
# <none>        2238.3 2284.3
# year      1   2692.0 2736.0
# log(SJR)  1   2322.2 2366.2
# sjournal 20   2387.7 2393.7


summary(ml4)
# Call:
#   glm(formula = y_ss ~ year + log(SJR) + sjournal, family = "binomial", 
#       data = xmodel)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.9620  -0.3926  -0.2413  -0.1194   3.1986  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)                                   -798.31555   45.00877 -17.737  < 2e-16 ***
#   year                                             0.39385    0.02230  17.664  < 2e-16 ***
#   log(SJR)                                         1.02319    0.11417   8.962  < 2e-16 ***
#   sjournalAm J Hum Genet                           0.62290    0.36634   1.700 0.089071 .  
#   sjournalAm J Med Genet B Neuropsychiatr Genet    1.16511    0.75340   1.546 0.121995    
#   sjournalAnn Rheum Dis                           -0.04556    0.56535  -0.081 0.935765    
#   sjournalDiabetes                                 2.12796    0.43732   4.866 1.14e-06 ***
#   sjournalEur J Hum Genet                          1.49580    0.51644   2.896 0.003775 ** 
#   sjournalFront Genet                              0.55730    0.75206   0.741 0.458676    
#   sjournalGastroenterology                        -1.40663    1.05101  -1.338 0.180778    
#   sjournalHum Mol Genet                            1.56355    0.25998   6.014 1.81e-09 ***
#   sjournalJ Allergy Clin Immunol                  -0.74808    1.03230  -0.725 0.468653    
#   sjournalJ Hum Genet                              0.98841    0.76671   1.289 0.197341    
#   sjournalJ Med Genet                              0.39431    1.04186   0.378 0.705084    
#   sjournalMol Psychiatry                           0.38990    0.34798   1.120 0.262515    
#   sjournalNat Commun                               1.25073    0.17994   6.951 3.63e-12 ***
#   sjournalNat Genet                                0.67673    0.25180   2.688 0.007197 ** 
#   sjournalNature                                   1.54730    0.42143   3.672 0.000241 ***
#   sjournalPLoS Genet                               1.82746    0.25590   7.141 9.25e-13 ***
#   sjournalPLoS One                                 1.82731    0.33221   5.501 3.79e-08 ***
#   sjournalProc Natl Acad Sci U S A                -0.18753    1.05763  -0.177 0.859266    
#   sjournalSci Rep                                  0.94854    0.37143   2.554 0.010658 *  
#   sjournalTransl Psychiatry                       -0.77840    0.73048  -1.066 0.286600    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 3000.5  on 4941  degrees of freedom
# Residual deviance: 2238.3  on 4919  degrees of freedom
# AIC: 2284.3
# 
# Number of Fisher Scoring iterations: 6

# Table of results for significant
predictors2 = rownames(summary(ml4)$coeff)
sigtable2 <- data.table(preds=predictors2,
                       OR = exp(coef(ml4)),
                       low.ci = sapply(predictors2, function(x){exp(coef(ml4)[x]+ -1 * 1.96 * sqrt(vcov(ml4)[x,x]))}),
                       hi.ci =  sapply(predictors2, function(x){exp(coef(ml4)[x]+ 1 * 1.96 * sqrt(vcov(ml4)[x,x]))}))

# fwrite(sigtable2, "../tables/Sigtable_ml4.tsv", sep = "\t")


################################################################################################################################
### Some additional visualisations

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



