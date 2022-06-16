library(data.table)
library(magrittr)
library(ggplot2)
library(seaborn)
library(cowplot)
theme_set(theme_cowplot())
d="data" #"~/rds/rds-cew54-basis/Projects/gwas-sharing/data/"
list.files(d)
x=file.path(d,"Main_data_20220610.tsv") %>% fread()
dim(x)
setnames(x,c("Journal","Year"),c("journal","year")) # to match previous data columns
## Date is first published date (perhaps online), Year is published in print
## year. Some papers have citations after their Date, but before their Year. Use
## truncated Date as year of publication.
setnames(x,"year","year_inprint")
x[,year:=as.numeric(substr(Date,1,4))]

y=file.path(d,"Citations_per_year_20220611.tsv") %>% fread()
dim(y)
y=merge(y, x[,.(PMID,Public_ss,pubyear=year)],by="PMID")
y[,years_since:=year - pubyear]
y=y[years_since >= -1] # allow -1 because so many - citations based on online early?!?
y[years_since<0, years_since:=0] # but truncate the 38 remaining with years_since=-1 to 0
y=y[order(PMID,years_since)]

with(x,table(year, Public_ss))
## exclude 2022 as incomplete and 2006 and before because no shared studies
x=x[year < 2022 & year > 2006]

##' Figure: observations by publication year
ggplot(x, aes(x=year,fill=Public_ss)) +
  geom_histogram(binwidth=1,colour="grey",position="dodge") +
  scale_fill_seaborn("Shared") +
  labs(x="Publication year",y="Observations") +
  background_grid(major="y")
ggsave("figures/fig-obs-by-pub-year.tiff",height=6,width=6,units="in")
## total observations increases with year, and sharing increases.
## Need to consider year as a potential confounder.



summary(x$relative_citation_ratio)
hist((x$relative_citation_ratio))
hist(log(x$relative_citation_ratio)) # but missing values for 0
## offset by minimum observed non-zero value
mincit=min(x$relative_citation_ratio[x$relative_citation_ratio>0])
x[relative_citation_ratio==0, relative_citation_ratio:=relative_citation_ratio+mincit]
x[,y:=log(relative_citation_ratio)]
hist(x$y,breaks=100)

## our dependent variable
x[,y_ss:=Public_ss=="Y"]

## ggplot(x, aes(x=year,y=y)) + geom_point() + geom_smooth()
##' Figure: Relative citations have decreased for GWAS papers over the years, with a recent trend up for 2021
top=ggplot(x, aes(x=year,y=y)) + geom_boxplot(aes(group=year)) + geom_smooth() +
  background_grid(major="y") +
  labs(x="Publication year",y="log relative citation ratio")
x[,year_ss:=paste(year,Public_ss)]
bot=ggplot(x[year<=2022], aes(x=year,y=y,col=Public_ss)) +
  geom_boxplot(aes(group=year_ss)) +
 # geom_jitter(width=0.2,height=0) +
  geom_smooth() +
  scale_colour_seaborn("Shared") +
  background_grid(major="y") +
  theme(legend.position=c(0,1),legend.justification=c(0,1),legend.direction = "horizontal") +
  labs(x="Publication year",y="log relative citation ratio")
plot_grid(top,bot,nrow=1)
ggsave("figures/fig-rel-cit-ratio-per-year.tiff",height=6,width=8,units="in")

## odd that the gap closes in the final year.
## look at citations per year since publication

## individual profiles too many to discern patterns
ggplot(y, aes(x=years_since, y=count, col=Public_ss,group=PMID)) + geom_path() + facet_wrap(~pubyear,scales="free_x")

y_avg=y[,.(mean_count=mean(count),n=.N),by=c("year","years_since","Public_ss","pubyear")][order(years_since)]
y_avg=y_avg[pubyear >= 2010 & pubyear<2019 & year >=pubyear & year < 2022]
y_avg[,max_years:=max(years_since),by="pubyear"]
y_rel=dcast(y_avg, year + pubyear +  years_since ~ Public_ss, value.var="mean_count")
y_rel[,rel_mean_count:=Y/N]
top=ggplot(y_avg[], aes(x=years_since, y=mean_count, col=Public_ss)) +
  geom_path() +
  geom_label(aes(label=n),data=y_avg[years_since==max_years & pubyear<2019 & pubyear>=2010]) +
  facet_wrap(~pubyear,scales="free_y") +
  scale_colour_seaborn("Shared") +
  scale_x_continuous(breaks=c(0,5,10)) +
  background_grid(major = "y") +
  ## scale_y_log10() +
  labs(x="Years since publication",y="Mean citation count")
bot=ggplot(y_rel, aes(x=years_since, y=rel_mean_count)) +
  geom_path() +
  facet_wrap(~pubyear,scales="free_y") +
  scale_colour_seaborn("Shared") +
  scale_x_continuous(breaks=c(0,5,10)) +
  background_grid(major = "y") +
  labs(x="Years since publication",y="Mean citation count ratio, shared / unshared")
plot_grid(top,bot,nrow=1)
ggsave("figures/fig-mean-citation-count-per-year.tiff",height=6,width=8,units="in")

## it takes 2-3 years for sharing effect to stabilise. Exclude 2021 because not enough time has passed
x=x[year<=2020]

glm(y ~ factor(year) + Public_ss, data=x) %>% summary()
## headline: exp(1.17)=3.22 - presumably some of this can be explained by other variables

## look at most common journals - do we see the same story within journal?
tt=table(x$journal)
sort(tt,decreasing = TRUE) %>% head(., 20)
jkeep=sort(tt,decreasing = TRUE) %>% names() %>% intersect(., x[Public_ss=="Y"]$journal) %>% head(.,20)
ggplot(x[journal %in% jkeep], aes(x=year,fill=Public_ss)) +
  geom_histogram(binwidth=1,position="dodge",col="grey") +
  scale_fill_seaborn() +
  facet_grid(journal ~ .)

results=lapply(seq_along(jkeep), function(j) {
  m=glm(y ~ factor(year) + Public_ss, data=x[year<2021 & journal==jkeep[j]])
  c(coefficients(m)["Public_ssY"],vcov(m)["Public_ssY","Public_ssY"])
  })  %>% do.call("rbind",.) %>% as.data.frame() %>% as.data.table()
setnames(results,c("beta","vbeta"))
results$journal=jkeep
  m=glm(y ~ factor(year) + Public_ss, data=x[year < 2021 & !(journal %in% jkeep)])
  results=rbind(results, data.table(beta=coefficients(m)["Public_ssY"],vbeta=vcov(m)["Public_ssY","Public_ssY"],journal="Other"))
results
results[,journal:=factor(journal)][,journal:=relevel(journal,"Other")]
results=merge(results,unique(x[,.(journal,SJR)]),by="journal",all.x=TRUE)
ggplot(results, aes(x=exp(beta),y=journal,xmin=exp(beta-1.96*sqrt(vbeta)),xmax=exp(beta+1.96*sqrt(vbeta)))) +
  geom_pointrange() +
  labs(x="Ratio of RCR in papers with shared vs unshared data",y="Journal") +
  background_grid(major="x")
## plenty of heterogeneity

## trend that higher SJR = lower RCR ratio, but some considerable outliers
ggplot(results, aes(x=exp(beta),y=log(SJR),xmin=exp(beta-1.96*sqrt(vbeta)),xmax=exp(beta+1.96*sqrt(vbeta)))) +
  geom_pointrange() +
  ## geom_label(aes(label=journal)) +
  geom_smooth(method="lm") +
  labs(x="Ratio of RCR in papers with shared vs unshared data",y="SJR") +
  background_grid(major="x")

## start building a model without Public_ss
m0=glm(y~year, data=x[!is.na(SJR) & year<2021])
m1=glm(y~factor(year), data=x[!is.na(SJR) & year<2021])
BIC(m1,m0) # relationship with year better modelled as linear
x[,sjournal:=ifelse(journal %in% jkeep,journal,"Other")]
x[,sjournal:=factor(sjournal)][,sjournal:=relevel(sjournal,"Other")]
m2=glm(y~year + sjournal, data=x[!is.na(SJR) & year<2021])# %>% summary()
m2a=glm(y~year + journal, data=x[!is.na(SJR) & year<2021])# %>% summary()
BIC(m0,m2,m2a) # sjournal is useful, but journal has too many levels to be useful
m3=glm(y~year + log(SJR), data=x[!is.na(SJR) & year<2021])# %>% summary()
m4=glm(y~year + log(SJR) + journal, data=x[!is.na(SJR) & year<2021])# %>% summary()
m5=glm(y~year + log(SJR) + sjournal, data=x[!is.na(SJR) & year<2021])# %>% summary()
BIC(m2,m3,m4,m5) # SJR has bigger impact than sjournal, but sjournal still useful in addition

m6=glm(y~year + log(SJR) + sjournal + human , data=x[!is.na(SJR) & year<2021])# %>% summary()
m7=glm(y~year + log(SJR) + sjournal + animal, data=x[!is.na(SJR) & year<2021])# %>% summary()
m8=glm(y~year + log(SJR) + sjournal + human + animal, data=x[!is.na(SJR) & year<2021])# %>% summary()
m9=glm(y~year + log(SJR) + sjournal + is_clinical, data=x[!is.na(SJR) & year<2021])# %>% summary()
BIC(m5,m6,m7,m8,m9) # adding human also helpful

## check all vars still useful - yes
drop1(m6)
m6a=glm(y~year + log(SJR) + human , data=x[!is.na(SJR) & year<2021])# %>% summary()
BIC(m6,m6a)

## go with m6, lowest BIC
m6a=glm(y~year + log(SJR) + sjournal + human + Public_ss, data=x[!is.na(SJR) & year<2021])# %>% summary()
BIC(m6,m6a)

summary(m6a)
## conclude sharing has highly significant effect (t=12.146, p<2e-16, coef=0.559961)
summary(m6a)$coeff["Public_ssY",,drop=FALSE]
cat("Ratio of RCR in publications that shared vs did not = ",
exp(coef(m6a)["Public_ssY"]),
"\n95% CI: ",
exp(coef(m6a)["Public_ssY"]+ c(-1,1) * 1.96 * sqrt(vcov(m6a)["Public_ssY","Public_ssY"])),
"\n")
