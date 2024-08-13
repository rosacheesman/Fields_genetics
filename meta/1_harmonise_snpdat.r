# .......................................................................................................
# ...............RSCRIPT for QC and harmonising of finngen sumstats for educ fields with moba GWAS.......................


#run this on the cluster without EA adjustemnt, and then also wi EA adjustemnt:
# for i in {0..10}; do
#     Rscript finngen_SNPs_harmonise.R "/cluster/projects/p805/rosac/fields/output/finngen/finngen_raw/no_ea_C$(printf "%02d" $i).gz" > finngen_noEA_harmonise_log.txt
# 	done

# for i in {0..10}; do
#     Rscript finngen_SNPs_harmonise.R "/cluster/projects/p805/rosac/fields/output/finngen/finngen_EA/C$(printf "%02d" $i).gz" > finngen_harmonise_log.txt
# done
# harmonise moba sumstats so they match finngen
#  change SNP column to chrbp

# .......................................................................................................
library(dplyr)
library(data.table)
# library(psych)
library(tidyr)
# .......................................................................................................

# import FinnGen sumstats
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if an argument is provided
if (length(args) == 0) {
  stop("No filename provided")
}

# Use the provided filename
filename <- args[1]

# Use the dynamically generated filename
finn <- fread(filename)

# import list of 7.2m FinnGen SNPs with high info, maf 1%, no multiallelic or strand ambiguous (from qc_sumstats_030124.R)
# /Users/rosacguio.no/Dropbox/PROMENTA/Choice of educational fields/summary_stats/cohorts_original/FinnGen_annotated_maf01_info8_b37.txt
# also contains the b37 reference allele data, which is what I will meta-analyse on.
ref5<-fread("/cluster/projects/p805/rosac/fields/output/finngen/FinnGen_annotated_maf01_info8_b37.txt")
colnames(ref5)[7]<-"chrom"
colnames(finn)[1]<-"chrom"

# subset FinnGen to these SNPs
finn_2<-merge(finn,ref5,by = c("chrom","pos","ref","alt"))
finn<-NULL

# .......................................................................................................
# ...............now performing additional QC of FinnGen and harmonising with MoBa.......................
# .......................................................................................................

# read in MoBa file that is a subset of Sumstats just including snps and EAF to use for allele flipping stuff
# MoBa sumstats  will be restricted upon meta-analysis
# write.table(mob,"/Users/rosacguio.no/Dropbox/PROMENTA/Choice of educational fields/summary_stats/cohorts_original/moba_snps_af1.txt",quote=F,col.names = T,row.names = F )
moba<-fread("/cluster/projects/p805/rosac/fields/output/finngen/moba_snps_af1.txt")
data_merge<-merge(moba,finn_2,by = "chrpos_b37")
# should be # 5,526,445

# .......................................................................................................
# PROGRAM to identify strand and sign flips
allele.qc = function(merged_data, a1_col1, a2_col1, a1_col2, a2_col2) {
  # Strand flip, to change the allele representation in the 2nd data-set
  strand_flip = function(ref) {
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip
  }
  
  flip1 = strand_flip(merged_data[[a1_col1]])
  flip2 = strand_flip(merged_data[[a1_col2]])
  
  # Create a dataframe to store the results
  snp = data.frame(merged_data)
  
  # Remove strand ambiguous SNPs (scenario 3)
  snp$keep = !((merged_data[[a2_col1]] == "A" & merged_data[[a2_col2]] == "T") |(merged_data[[a2_col1]] == "T" & merged_data[[a2_col2]] == "A")|(merged_data[[a2_col1]] == "C" & merged_data[[a2_col2]] == "G")|(merged_data[[a2_col1]] == "G" & merged_data[[a2_col2]] == "C"))
  
  # Remove non-ATCG coding
  snp$keep[ merged_data[[a2_col1]] != "A" & merged_data[[a2_col1]] != "T" & merged_data[[a2_col1]] != "G" & merged_data[[a2_col1]] != "C" ] = FALSE
  snp$keep[ merged_data[[a2_col2]] != "A" & merged_data[[a2_col2]] != "T" & merged_data[[a2_col2]] != "G" & merged_data[[a2_col2]] != "C" ] = FALSE
  
  # As long as scenario 1 is involved, sign_flip will return TRUE
  snp$sign_flip = (merged_data[[a2_col1]] == merged_data[[a1_col2]] & merged_data[[a2_col2]] == merged_data[[a1_col1]]) | (merged_data[[a2_col1]] == flip2 & merged_data[[a2_col2]] == flip1)
  
  # As long as scenario 2 is involved, strand_flip will return TRUE
  snp$strand_flip = (merged_data[[a2_col1]] == flip1 & merged_data[[a2_col2]] == flip2) | 
    (merged_data[[a2_col1]] == flip2 & merged_data[[a2_col2]] == flip1)
  
  # Remove other cases, e.g., tri-allelic, one dataset is A C, the other is A G
  exact_match = (merged_data[[a2_col1]] == merged_data[[a1_col1]] &merged_data[[a2_col2]] == merged_data[[a1_col2]]) 
  snp$keep[!(exact_match | snp$sign_flip | snp$strand_flip)] = FALSE
  
  return(snp)
}


# .......................................................................................................
# run it
qc = allele.qc(data_merge, "A1", "alt_b37", "A2", "ref_b37")

# .......................................................................................................
# Function to adjust beta directions (beta) and effect allele frequencies (af_alt) in finngen IF sign flip=T

adjust_summary_statistics <- function(data) {
  # Change beta directions for sign flip - should be other sign.
  data$adjusted_beta <- ifelse(data$sign_flip, -data$beta, data$beta)
  
  # Adjust effect allele frequencies for sign flip- should be freq of what is finngen reference allele
  data$adjusted_af_alt <- ifelse(data$sign_flip, 1-data$af_alt, data$af_alt)
  
  # Subset adjusted data to exclude rows where keep is FALSE
  adjusted_data <- subset(data, keep)
}
# .......................................................................................................

# Apply the function 
qc_noflip <- adjust_summary_statistics(qc)
# .......................................................................................................

# Deal with EAF discrepancies between Finngen and MoBa
# Specify a threshold for EAF difference
eaf_threshold <- 0.2  # Adjust as needed

# Calculate EAF difference
qc_noflip$eaf_difference <-abs(qc_noflip$AF1-qc_noflip$adjusted_af_alt)

# Identify SNPs with large EAF differences
snps_to_exclude <- qc_noflip$eaf_difference > eaf_threshold

# Exclude SNPs with large EAF differences
qc_noflip_sim <- qc_noflip[!snps_to_exclude, ]


# .......................................................................................................

# May still be duplicate SNPs
# keep only 1 per duplicate-- and choose the one with lowest EAF discrepancy
qc_noflip_sim_nodup <- qc_noflip_sim %>%
  group_by(chrpos_b37) %>%
  slice(which.min(eaf_difference))

# .......................................................................................................

# subset to FinnGen sumstats data

clean_cols<-qc_noflip_sim_nodup[,c("chrpos_b37","A1","A2","adjusted_af_alt","adjusted_beta","pval","sebeta")]

# ......................................................................................................

# Write data to the gzipped file 
# Replace "C07.txt.gz" finngen names with the dynamically generated output filename
# basename removes file path and sans etc removes .gz
output_filename <- paste0("/cluster/projects/p805/rosac/fields/output/finngen/harmonised/harmonized_", tools::file_path_sans_ext(basename(filename)))
write.table(clean_cols, file = gzfile(output_filename), sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


