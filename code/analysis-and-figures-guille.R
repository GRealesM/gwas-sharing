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


##' Figure: observations by publication year
ggplot(x, aes(x=year,fill=Public_ss)) +
  geom_histogram(binwidth=1,colour="grey",position="dodge") +
  scale_fill_seaborn("Shared Sumstats   ", guide = guide_legend(reverse=TRUE), labels=c("No", "Yes")) +
  labs(x="Publication year",y="Observations") +
  #scale_fill_discrete(guide = guide_legend(reverse=TRUE), breaks = rev(levels(x$Public_ss)), labels=c("Yes", "No")) +
  background_grid(major="y")+
  theme(legend.position=c(0.05,0.96),legend.justification=c(0,1),legend.direction = "horizontal")
ggsave("../figures/fig-obs-by-pub-year.tiff",height=6,width=7,units="in", bg="white")
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
plot_grid(top,bot,nrow=1)
ggsave("../figures/fig-rel-cit-ratio-per-year.tiff",height=6,width=12,units="in", bg="white")



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
ggsave("../figures/fig-mean-citation-count-per-year.tiff",height=8,width=14,units="in", bg ="white")

## it takes 2-3 years for sharing effect to stabilise. Exclude 2021 because not enough time has passed
x=x[year<=2020]

glm(y ~ factor(year) + Public_ss, data=x) %>% summary()
## headline: exp(1.25)=3.49 - presumably some of this can be explained by other variables

## Exploration by journal
## look at most common journals - do we see the same story within journal?
tt=table(x$journal)
sort(tt,decreasing = TRUE) %>% head(., 20)
jkeep=sort(tt,decreasing = TRUE) %>% names() %>% intersect(., x[Public_ss=="Y"]$journal) %>% head(.,20)
ggplot(x[journal %in% jkeep], aes(x=year,fill=Public_ss)) +
  geom_histogram(binwidth=1,position="dodge",col="grey") +
  scale_fill_seaborn() +
  facet_grid(journal ~ .)
# 
# 
# 
# results=lapply(seq_along(jkeep), function(j) {
#   m=glm(y ~ factor(year) + Public_ss, data=x[year<2021 & journal==jkeep[j]])
#   c(coefficients(m)["Public_ssY"],vcov(m)["Public_ssY","Public_ssY"])
#   })  %>% do.call("rbind",.) %>% as.data.frame() %>% as.data.table()
# setnames(results,c("beta","vbeta"))
# results$journal=jkeep
#   m=glm(y ~ factor(year) + Public_ss, data=x[year < 2021 & !(journal %in% jkeep)])
#   results=rbind(results, data.table(beta=coefficients(m)["Public_ssY"],vbeta=vcov(m)["Public_ssY","Public_ssY"],journal="Other"))
# results
# results[,journal:=factor(journal)][,journal:=relevel(journal,"Other")]
# results=merge(results,unique(x[,.(journal,SJR)]),by="journal",all.x=TRUE)
# ggplot(results, aes(x=exp(beta),y=journal,xmin=exp(beta-1.96*sqrt(vbeta)),xmax=exp(beta+1.96*sqrt(vbeta)))) +
#   geom_pointrange() +
#   labs(x="Ratio of RCR in papers with shared vs unshared data",y="Journal") +
#   background_grid(major="x")
# ## plenty of heterogeneity
# 
# 
# # Top 50 journals with the most published GWAS, with at least one sharer. Compare SJR and RCR.
# 
# j50=sort(tt,decreasing = TRUE) %>% names() %>% intersect(., x[Public_ss=="Y"]$journal) %>% head(.,50)
# res50=lapply(seq_along(j50), function(j) {
#   m=glm(y ~ factor(year) + Public_ss, data=x[year<2021 & journal==j50[j]])
#   c(coefficients(m)["Public_ssY"],vcov(m)["Public_ssY","Public_ssY"])
# })  %>% do.call("rbind",.) %>% as.data.frame() %>% as.data.table()
# setnames(res50,c("beta","vbeta"))
# res50$journal=j50
# 
# mt50=glm(y ~ factor(year) + Public_ss, data=x[year < 2021 & !(journal %in% j50)])
# res50=rbind(res50, data.table(beta=coefficients(mt50)["Public_ssY"],vbeta=vcov(mt50)["Public_ssY","Public_ssY"],journal="Other"))
# res50
# res50 <- na.omit(res50)
# res50[,journal:=factor(journal)][,journal:=relevel(journal,"Other")]
# res50=merge(res50,unique(x[,.(journal,SJR)]),by="journal",all.x=TRUE)
# 
# 
# ## trend that higher SJR = lower RCR ratio, but some considerable outliers
# ggplot(res50, aes(x=exp(beta),y=log(SJR),xmin=exp(beta-1.96*sqrt(vbeta)),xmax=exp(beta+1.96*sqrt(vbeta)))) +
#   geom_pointrange() +
#   ## geom_label(aes(label=journal)) +
#   geom_smooth(method="lm") +
#   labs(x="Ratio of RCR in papers with shared vs unshared data",y="SJR") +
#   background_grid(major="x")

### Building linear models to test whether sharing makes a difference in citation, when considering other variables

## start building a model without Public_ss
# Prepare dataset for modelling, we'll remove missing data for relevant variables. We'll also consider years 2007-2020
xmodel <- x[!is.na(SJR) & !is.na(y)]
m0=glm(y~year, data=x)                       # m0: year only, as numeric
m1=glm(y~factor(year), data=xmodel)          # m1: year only, as factor
BIC(m1,m0) # relationship with year better modelled as linear

xmodel[,sjournal:=ifelse(journal %in% jkeep,journal,"Other")]             # Use top 20 most common journals with at least one shared dataset, pool the rest as "other"
xmodel[,sjournal:=factor(sjournal)][,sjournal:=relevel(sjournal,"Other")]
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
# <none>                  3745.1 12702
# year                1   3957.8 12973
# log(SJR)            1   4836.0 13964
# sjournal           19   3890.0 12852
# molecular_cellular  1   3764.8 12726


## go with m9
m9a=glm(y~year + log(SJR) + sjournal + molecular_cellular + Public_ss, data=xmodel)# %>% summary()
BIC(m9,m9a)

summary(m9a)
## conclude sharing has highly significant effect (t=12.131, p<2e-16, coef=0.559961)
summary(m9a)$coeff["Public_ssY",,drop=FALSE]
cat("Ratio of RCR in publications that shared vs did not = ",
exp(coef(m9a)["Public_ssY"]),
"\n95% CI: ",
exp(coef(m9a)["Public_ssY"]+ c(-1,1) * 1.96 * sqrt(vcov(m9a)["Public_ssY","Public_ssY"])),
"\n")
# Ratio of RCR in publications that shared vs did not =  1.749155 
# 95% CI:  1.598061 1.914535 





