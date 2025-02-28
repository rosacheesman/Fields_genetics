
# using the genomic sem gwas pipeline 
# https://github.com/sarahcolbert/gsemGWAS/blob/master/README.md

# -- moved pipeline to NEW folder 
# copied CFA GWAS but deleted the big files.

# make directories
cd /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/
# for i in edu arts social business natural_sci ict engineering agri health services
for i in arts social business natural_sci ict engineering agri health services
do 
mkdir ${i}
done

# move the Sumstats and rG matrices in to correct folder


# _____________________________________________________________________________________
# make config files
for i in arts social business natural_sci ict engineering agri health services
# for i in edu
do 
cd /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}
echo '#!/usr/bin/env bash
# Please put the full path to this directory
# e.g. if the config file path is "/home/Projects/my_gsem/config"
#then project_dir="/home/Projects/my_gsem/"
# it is important that you include a backslash at the end of this path
export project_dir="/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/'$i'/"
## put location for sumstats RData file
export sumstats_file="/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/'$i'/EA4_'$i'_Sumstats.RData"
## put location for LDSC matrix file
export ldsc_file="/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/'$i'/EA4_'$i'_LDSCoutput.RData"
## DO NOT MODIFY
export code_dir="${project_dir}code/"
mkdir -p ${code_dir}
export outerr_dir="${project_dir}outerr/"
mkdir -p ${outerr_dir}
export sumstats_dir="${project_dir}split_sumstats/"
mkdir -p ${sumstats_dir}
export results_dir="${project_dir}results/"
mkdir -p ${results_dir}' > config_${i}
chmod +rwx config_${i} 
done


# _____________________________________________________________________________________

# -- split sumstats up into chunks of  SNPs
# Chose 5000 SNPs 1300 ish sets 
for i in arts social business natural_sci ict engineering agri health services
do 
cd /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}
source ./config_${i}
# sed -i 's/1000/5000/g' ./scripts/split_sumstats.R
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0
Rscript /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/scripts/split_sumstats.R
# two outputs from this script: (1) the new summary statistics files saved as "/split_sumstats/sumstats*.txt" 
# and (2) the number of SNP subsets created which is saved as "num_SNP_sets.txt"
cd /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}/split_sumstats/
cat num_SNP_sets.txt
done


# _____________________________________________________________________________________
# create 10 r and bash scripts for gwas by sub
# NB must change the starting values for each field manually

for i in arts social business natural_sci ict engineering agri health services
do

cd /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}/code/

echo '
cat("\n### Preparing workspace \n")
remove(list = ls()); start.time <- Sys.time()
source("/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/edu/code/setup_R_gsem.R")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("\n ... loading packages \n")
library(devtools)
library(R.utils)
library(GenomicSEM,lib.loc=LIB)
library(data.table)

### load the summary statistics RData file in the split form
split_sumstats <- read.table(paste(Sys.getenv("sumstats_dir"),"sumstats",Sys.getenv("cc"),".txt", sep = ""), header = TRUE)
print(paste("finished loading summary statistics from set ",Sys.getenv("cc"), sep = ""))

### load the LDSC covariance matrix
load(paste(Sys.getenv("ldsc_file")))
print("finished loading LDSC covariance matrix")

# userGWAS
# using starting values from model w/o snp effects

# v1 is field, v2 is ea4
model<-"EA=~NA*V2 + start(xxxx)*V2+ start(xxxx)*V1
'$i'=~NA*V1 + start(xxxx)*V1
'$i'~SNP
EA~SNP
'$i'~~1*'$i'
EA~~1*EA
EA~~0*'$i'
V1 ~~ 0*V1
V2~~0*V2
V2~~0*V1
SNP~~SNP"

#Run the Genomic SEM GWAS
outputGWAS<-userGWAS(covstruc = LDSCoutput, SNPs = split_sumstats, model = model, sub="'$i'~SNP",smooth_check=TRUE,fix_measurement=TRUE,std.lv = TRUE,parallel=FALSE)
print("GWAS completed")

write.csv(outputGWAS, file=paste(Sys.getenv("results_dir"),Sys.getenv("cc"),".csv", sep = ""), row.names=FALSE)
print(paste("analysis for set ",Sys.getenv("cc")," complete", sep = ""))
' > /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}/code/Gw_by_sub_${i}.R
done 
# NEED to manually enter starting values afterwards-- automate later!!


# run config then run script

for i in arts social business natural_sci ict engineering agri health services
do
cd /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}/
source ./config_${i}
echo '#!/bin/bash
#SBATCH --mem-per-cpu=6GB
#SBATCH --partition=bigmem
#SBATCH --account=p805
#SBATCH --time=2:00:00
#SBATCH --output=./outerr/set.%a.out
#SBATCH --error=./outerr/set.%a.err
date
hostname
export cc="${SLURM_ARRAY_TASK_ID}"
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
Rscript /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/'$i'/code/Gw_by_sub_'$i'.R
date' > /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}/code/Gw_by_sub_${i}.sh
sbatch --array=0-1317%100 /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}/code/Gw_by_sub_${i}.sh
done



# _____________________________________________________________________________________

# -- compile results into multiple sets of summary statistics for multiple factors


for i in edu arts social business natural_sci ict engineering agri health services
do
cd /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}/
source ./config_${i}
cd /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/${i}/code/
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

### concatenate all results files and only keep header from first file
awk "FNR>1 || NR==1" ${results_dir}*.csv > ${results_dir}'$i'_EA4sub_sumstats.csv

date' > cat_results_${i}.sh 
sbatch cat_results_${i}.sh
done

# _____________________________________________________________________________________

# check length
zcat Cheesman_meta_edu.txt.gz |wc -l
