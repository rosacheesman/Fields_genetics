
# steps for running genomic SEM GWAS
#.......................................................................................................

# add rsid column to sumstats from METAL
# this is just including the EA-adjusted Fields sumstats used for the factor modellng

for (i in sprintf("%02d", 1:10)) {
  # Read the input file
  filename <- paste0("META_EA_", i, "1.txt.gz")
  
  # Read the input file
  gp <- fread(filename)
  
  # Read in rsids
  rsids <- fread('chrpos_b37_rsids.txt')
  
  # Merge gzfile with rsids file based on MarkerName column
  merged_data <- merge(gp, rsids, by.x = "MarkerName", by.y = "MarkerName", all.x = TRUE)
  # rename markername
  colnames(merged_data)[1]<-'SNP2'
  
  # Make new file
  pop_filename <- gsub("1.txt.gz", "_2.txt", filename)
  write.table(merged_data, pop_filename, col.names = TRUE, row.names = FALSE, quote = FALSE)
}


# gzip

cd /Users/rosacguio.no/Dropbox/PROMENTA/Choice of educational fields/summary_stats/meta
gzip *_2.txt
# .......................................................................................................


# -- munge fields sumstats and ran ldsc in Rstudio
library(devtools)
library(GenomicSEM)
library(data.table)

setwd("/cluster/p/p805/cluster/rosac/fields/output/meta/rsids")
# MarkerName      Allele1 Allele2 Freq1   FreqSE  MinFreq MaxFreq Weight  Zscore  P-value Direction
# 11:79226604     t       c       0.1275  0.0151  0.1183  0.1523  463092.00       0.852   0.394   +-
# 5:29439275      t       c       0.3282  0.0180  0.3173  0.3579  462831.00       -0.151  0.8802  -+

files = c("META_EA_01_2.txt.gz","META_EA_02_2.txt.gz","META_EA_03_2.txt.gz","META_EA_04_2.txt.gz","META_EA_05_2.txt.gz","META_EA_06_2.txt.gz","META_EA_07_2.txt.gz","META_EA_08_2.txt.gz","META_EA_09_2.txt.gz","META_EA_10_2.txt.gz")
ref = "/tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/reference.1000G.maf.0.005.txt"
trait.names = c("edu","arts","soc","bus","natsci","ict","eng","agri","health","serv")
se.logit = c(F,F,F,F,F,F,F,F,F,F)
info.filter = 0.6
maf.filter = 0.01
# When backing out a logistic beta using the linprob argument this requires the sum of effective sample sizes across the cohorts contributing to the GWAS
Ns=c(102970	,97262	,69123	,261182	,40072	,50819	,317209	,63834	,292929	,168157)
linprobs=c(T,T,T,T,T,T,T,T,T,T)

# Whether the phenotype was a dichotomous outcome for which there are only Z-statistics in the summary statistics file -or- it was a dichotomous outcome analyzed using an OLS estimator, as is the case UKB phenotypes such as the 12 used here that were analyzed using the Hail software. In this case, a form of standardization is performed where logistic betas are backed out from a signed effect column. If no information is provided the default is for the function to assume NULL.
# all traits are binary and were analyzed using a linear function so we set linprob to TRUE and supply the sum of neff
# the object needs to be called "sumstats"
sumstats<-sumstats(files, ref, trait.names, se.logit, info.filter, maf.filter, OLS=NULL,linprob=linprobs,N=Ns,betas=NULL)
save(sumstats, file="/tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/Sumstats3.RData")

# run ldsc
traits <- c("edu_EA_adj.sumstats.gz","arts_EA_adj.sumstats.gz","social_EA_adj.sumstats.gz","business_EA_adj.sumstats.gz","natural_sci_EA_adj.sumstats.gz","ict_EA_adj.sumstats.gz","engineering_EA_adj.sumstats.gz","agri_EA_adj.sumstats.gz","health_EA_adj.sumstats.gz","services_EA_adj.sumstats.gz")
population.prev <- c(0.0515	,0.0509	,0.0337	,0.146	,0.0186	,0.0268	,0.226	,0.0327	,0.141	,0.0949)
# When inputting the sum of effective sample sizes, the sample prevalence should then be entered as 0.5 when running ldsc to reflect the fact that effective sample size already corrects for sample ascertainment
sample.prev   <- rep(0.5,20)
ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"
LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld,trait.names = c("edu","arts","soc","bus","natsci","ict","eng","agri","health","serv"))
save(LDSCoutput, file="/tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/Fields_META_EA_adj_rG.RData")


# using ththis pipeline for splitting up sumstats etc.:
# https://github.com/sarahcolbert/gsemGWAS/blob/master/README.md

# -- moved pipeline to folder 


# -- edited config file then activated it
source ./config

# -- split sumstats up into chunks of  SNPs
# array nr is max 4001, and i have 7m SNPs, so use 2000 SNPs per 4000 array sets
# or 5000 SNPs 1300 ish sets 
cd /tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/gsemGWAS
sed -i 's/1000/5000/g' ./scripts/split_sumstats.R
module load R/4.0.0-foss-2020a
Rscript ./scripts/split_sumstats.R
# two outputs from this script: (1) the new summary statistics files saved as "/split_sumstats/sumstats*.txt" 
# and (2) the number of SNP subsets created which is saved as "num_SNP_sets.txt"
cd /tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/gsemGWAS/split_sumstats
cat num_SNP_sets.txt
# 1313 (including 0 )


# -- make GWAS script
cd /tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/gsemGWAS/scripts
echo '
library(GenomicSEM)
library(data.table)

### load the summary statistics RData file in the split form
print(paste("loading summary statistics from set ",Sys.getenv("cc"),"...", sep = ""))
split_sumstats <- read.table(paste(Sys.getenv("sumstats_dir"),"sumstats",Sys.getenv("cc"),".txt", sep = ""), header = TRUE)
print(paste("finished loading summary statistics from set ",Sys.getenv("cc"), sep = ""))

### load the LDSC covariance matrix
print("loading LDSC covariance matrix...")
load(paste(Sys.getenv("ldsc_file")))
print("finished loading LDSC covariance matrix")


# userGWAS
model<-"F1 =~ NA*eng + natsci + ict + edu + health+ serv +agri+ soc+arts
F2 =~ NA*arts + soc + bus + edu +health+ serv+eng+agri+ ict
F1~~0*F2
soc ~~ a*soc
a > .001
F1~SNP
F2~SNP
SNP~~SNP"

#Run the Genomic SEM GWAS
TwoFactors<-userGWAS(covstruc = LDSCoutput, SNPs = split_sumstats, model = model, sub=c("F1~SNP", "F2~SNP"),smooth_check=TRUE,fix_measurement=TRUE,std.lv = TRUE,parallel=FALSE)

print("GWAS completed")

print("writing results to file...")
write.csv(TwoFactors, file=paste(Sys.getenv("results_dir"),Sys.getenv("cc"),".csv", sep = ""), row.names=FALSE)
print(paste("analysis for set ",Sys.getenv("cc")," complete", sep = ""))
' > multi_GWAS.R




# -- Using multi_GWAS.sbatch, you can run a separate job for each set of SNPs
# just need to edit nr subsets in script wrapper /when run job
#SBATCH --array=0-3000%400

cd /tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/gsemGWAS/scripts

echo '#!/bin/bash
#SBATCH --mem-per-cpu=6GB
#SBATCH --partition=bigmem
#SBATCH --account=p805
#SBATCH --time=2:00:00
#SBATCH --output=./outerr/set.%a.out
#SBATCH --error=./outerr/set.%a.err

### display time at beginning and end of script to track how long it takes
date
### display node that job was run on
hostname

### set variable that will be the number of the SNP set read by R
export cc="${SLURM_ARRAY_TASK_ID}"

### run the Rscript
module load R/4.0.0-foss-2020a
Rscript ${project_dir}scripts/multi_GWAS.R

date' > multi_GWAS.sh
sbatch --array=0-1312%100 multi_GWAS.sh



# check output 
cat /tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/gsemGWAS/scripts/outerr/set.2.out
cd /tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/gsemGWAS/results


# -- compile results into multiple sets of summary statistics for multiple factors


cd /tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/gsemGWAS/scripts/
echo '#!/bin/bash
#SBATCH --job-name=combine_SNPsets
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=5GB
#SBATCH --account=p805
#SBATCH --partition=bigmem
#SBATCH --output=./outerr/combine_Factor%a.out
#SBATCH --error=./outerr/combine_Factor%a.err

date
hostname

### create variable that is equal to the number of factors (2)
cc="${SLURM_ARRAY_TASK_ID}"

### concatenate all results files and only keep header from first file
awk "FNR>1 || NR==1" ${results_dir}*.csv > ${results_dir}F_stats.csv

date' > cat_results.sh 
sbatch cat_results.sh

# first 22 cols are F1 results, second 22 are F2 
cd /tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/gsemGWAS/results
awk -F ',' '{print NF; exit}' F_stats.csv

# -- Split into F1 and F2 sumstats 

# Write first 22 columns to output1.csv
awk -F ',' '{
    for(i=1; i<=22; i++) {
        printf "%s%s", $i, (i==22 ? "\n" : ",")
    }
}' F_stats.csv > output1.csv

# Write remaining columns to output2.csv
awk -F ',' '{ for(i=23; i<=NF; i++) printf "%s%s", $i, (i==NF ? "\n" : ",") }' F_stats.csv > output2.csv


# Replace headers of output2.csv with headers of output1.csv
header_line=$(head -n 1 output1.csv)
sed -i "1s/.*/$header_line/" output2.csv


# rename  
mv output1.csv F1_sumstats.csv 
mv output2.csv F2_sumstats.csv


# -- Calculate Neff for latent factors


### Read in summary statistics for latent factorS which you wish to calculate Neff for
setwd("/tsd/p805/data/durable/projects/Rosa_gwas/fields_gwas/cfa_gwas/gsemGWAS/results/")

awk -F ',' '{print NF; exit}' F1_sumstats.csv



library(data.table)
factor <- fread("F1_sumstats.csv")

# factor <- fread("F2_sumstats.csv")
factor$P<-factor$Pval_Estimate

### Restrict to MAF of 40% and 10%
# factor<-subset(factor, factor$MAF <= .4 & factor$MAF >= .1)

### Calculate Neff using equation from Mallard et al.
# Effective_N<-(mean(((factor$Z_Estimate/factor$est)^2)/(2*factor$MAF*(1-factor$MAF))))

# change to txt format

write.table(factor,"F1_sumstats.txt", quote=F,col.names = T, row.names=F)
write.table(factor,"F2_sumstats.txt", quote=F,col.names = T, row.names=F)

# f1
# 10413.49
# f2
# 7352.692






