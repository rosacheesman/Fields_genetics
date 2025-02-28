library("devtools")
install.packages("/cluster/projects/p805/software/GenomicSEM-master", repos = NULL, type = "source")
library(GenomicSEM)
library(data.table)

setwd("/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/")

broad_files<-c("edu_EA4sub_sumstats.txt.gz",
"arts_EA4sub_sumstats.txt.gz",
"social_EA4sub_sumstats.txt.gz",
"business_EA4sub_sumstats.txt.gz",
"natural_sci_EA4sub_sumstats.txt.gz",
"ict_EA4sub_sumstats.txt.gz",
"engineering_EA4sub_sumstats.txt.gz",
"agri_EA4sub_sumstats.txt.gz",
"health_EA4sub_sumstats.txt.gz",
"services_EA4sub_sumstats.txt.gz")

broad_trait_names<-c("edu_EA4sub","arts_EA4sub","social_EA4sub","business_EA4sub","natural_sci_EA4sub","ict_EA4sub","engineering_EA4sub","agri_EA4sub","health_EA4sub","services_EA4sub")

# neff for ea-adjusted data
broad_N<-c(102970	,97262	,69123	,261182	,40072	,50819	,317209	,63834	,292929	,168157)


munge(files=broad_files, hm3 = "w_hm3.snplist",trait.names=broad_trait_names,N=broad_N)

# ...........................................................................................................

# LDSC 
# SNP h2s
# rGs with each other 
# rGs with previous EA-adjusted results.

ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"

traits <- c("edu_EA4sub.sumstats.gz","arts_EA4sub.sumstats.gz","social_EA4sub.sumstats.gz","business_EA4sub.sumstats.gz","natural_sci_EA4sub.sumstats.gz","ict_EA4sub.sumstats.gz","engineering_EA4sub.sumstats.gz","agri_EA4sub.sumstats.gz","health_EA4sub.sumstats.gz","services_EA4sub.sumstats.gz")
population.prev <- c(0.0515	,0.0509	,0.0337	,0.146	,0.0186	,0.0268	,0.226	,0.0327	,0.141	,0.0949)
# .5 when using sumneff
sample.prev <- rep(0.5,10)

names<-c("edu_EA4sub","arts_EA4sub","social_EA4sub","business_EA4sub","natural_sci_EA4sub","ict_EA4sub","engineering_EA4sub","agri_EA4sub","health_EA4sub","services_EA4sub")
Fields_wo_EA_LDSCStand <- ldsc(traits, sample.prev, population.prev, ld, wld,stand=T,trait.names=names)
save(Fields_wo_EA_LDSCStand, file="Fields_wo_EA_LDSCStand.RData")

