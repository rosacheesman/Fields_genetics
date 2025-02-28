
# _____________________________________________________________________________________
library("devtools")
library(GenomicSEM)
# _____________________________________________________________________________________
# get LDSC object for ALL fields and EA4. can use this for mod w/o snp effects.
traits <- c("edu.sumstats.gz","arts.sumstats.gz","social.sumstats.gz","business.sumstats.gz",
"natural_sci.sumstats.gz","ict.sumstats.gz","engineering.sumstats.gz","agri.sumstats.gz"
,"health.sumstats.gz","services.sumstats.gz",
"EA4.sumstats.gz")
population.prev <- c(0.0515	,0.0509	,0.0337	,0.146	,0.0186	,0.0268	,0.226	,0.0327	,0.141	,0.0949,
NA)
sample.prev   <- c(rep(0.5,10),NA)
trait_names<-c("edu","arts","social","business","natural_sci","ict","engineering","agri","health","services",
"EA4")
ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld,trait.names=trait_names)
save(LDSCoutput, file="Fields_raw_EA4_LDSCoutput.RData")
# _____________________________________________________________________________________

