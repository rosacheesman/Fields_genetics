

echo 'source("/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/scripts/N_weighted_GWAMA.function.1_2_6.R")
# need the supporting functions for GWAMA
# use modified GWAMA code replacing h2 with loadings ie. Changed w = sqrt(N*h2)to w = sqrt(N)*standardised_loadings
source("/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/scripts/GWAMA_function_AF.R")

# Load dependencies
library(data.table)
library(stringr)

# Make an empty list with length equal to the number of input files
dat <- vector("list")

# Set working directory to the GWAMA format data folder
setwd("/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/GWAMA_format/")
# List all GWAMA files available in this folder
gwama_files <- list.files(pattern = "gwama")

# Create a dataframe pairing GWAMA files with their respective file names
# this is what r calls the object read in.
all_files <- as.data.frame(cbind(gwama_files, file_names = c("agri_EA4sub","arts_EA4sub","business_EA4sub",
"edu_EA4sub","engineering_EA4sub","health_EA4sub","ict_EA4sub","natural_sci_EA4sub",
"services_EA4sub","social_EA4sub")))

# Read in the data from each file, storing them in the list 
n <- 0
for (file in 1:nrow(all_files)) {
  n <- n + 1
  print(n)
  dat[[all_files$file_names[n]]] <- fread(all_files$gwama_files[n], data.table = F)
}
# ...........................................................................................................

# Load covariance matrix from LDSC

# got rG object for only the gwsub x10

load("/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/Fields_wo_EA_LDSCStand.RData")


# Select the matrix of LDSC intercepts and ensure it matches data order
order <- names(dat)

dimnames(Fields_wo_EA_LDSCStand$I)[[1]] <- dimnames(Fields_wo_EA_LDSCStand$S)[[2]]
dimnames(Fields_wo_EA_LDSCStand$I)[[2]] <- dimnames(Fields_wo_EA_LDSCStand$S)[[2]]
CTI <- as.matrix(Fields_wo_EA_LDSCStand$I[order, order])

# Load correlation matrix
dimnames(Fields_wo_EA_LDSCStand$S_Stand)[[1]] <- dimnames(Fields_wo_EA_LDSCStand$S)[[2]]
dimnames(Fields_wo_EA_LDSCStand$S_Stand)[[2]] <- dimnames(Fields_wo_EA_LDSCStand$S)[[2]]
cormatrix <- Fields_wo_EA_LDSCStand$S_Stand[order, order]

# Perform eigen decomposition to get eigenvalues and eigenvectors
eigenvectors <- eigen(cormatrix)$vectors
eigenvalues <- eigen(cormatrix)$values

# Calculate standardised loadings for PC1 and PC2
loadings_PC2 <- as.vector(eigenvectors %*% sqrt(diag(eigenvalues))[, 2])

# Run my modified GWAMA function for PC2
my_GWAMA(x = dat,
         cov_Z = CTI,
         h2 = loadings_PC2,
         out = "/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/GWAMA_format/",
         name = "PC2",
         output_gz = T,
         check_columns = F)' > /ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/scripts/gwama_fields_pc2.R

cd /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/scripts/
echo '#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=15GB
#SBATCH --account=p805_tsd
#SBATCH --time=2:00:00
date
hostname
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
Rscript /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/scripts/gwama_fields_pc2.R
date' >  /cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/scripts/gwama_fields_pc2.sh
sbatch gwama_fields_pc2.sh
