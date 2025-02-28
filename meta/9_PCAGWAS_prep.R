
# format sumstats for PCA:
#  zcat health_EA4sub_sumstats.txt.gz|head
# SNP CHR BP MAF A1 A2 lhs op rhs free label est SE Z_Estimate Pval_Estimate chisq chisq_df chisq_pval AIC error warning Z_smooth Neff
# rs149168804 1 828539 0.0188867 A T health ~ SNP NA 1 -0.0754506966387086 0.0770823488145085 -0.978832350065949 0.327662822592488 1.39362215891318e-05 3 0.999999986163235 6.000013

# restrict to HM3 SNPs
# restrict to columns  --> # SNPID,CHR,BP,EA,OA,EAF,N,Z,P)
# ensure EAF is correct


# List of input files

library(data.table)
input_files <- c("/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/edu/results/edu_EA4sub_sumstats.csv",
"/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/arts/results/arts_EA4sub_sumstats.csv",
"/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/social/results/social_EA4sub_sumstats.csv",
"/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/business/results/business_EA4sub_sumstats.csv",
"/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/natural_sci/results/natural_sci_EA4sub_sumstats.csv",
"/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/ict/results/ict_EA4sub_sumstats.csv",
"/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/engineering/results/engineering_EA4sub_sumstats.csv",
"/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/agri/results/agri_EA4sub_sumstats.csv",
"/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/health/results/health_EA4sub_sumstats.csv",
"/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/services/results/services_EA4sub_sumstats.csv")

munged_files<- c("/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/edu_EA4sub.sumstats.gz",
"/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/arts_EA4sub.sumstats.gz",
"/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/social_EA4sub.sumstats.gz",
"/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/business_EA4sub.sumstats.gz",
"/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/natural_sci_EA4sub.sumstats.gz",
"/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/ict_EA4sub.sumstats.gz",
"/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/engineering_EA4sub.sumstats.gz",
"/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/agri_EA4sub.sumstats.gz",
"/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/health_EA4sub.sumstats.gz",
"/ess/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/services_EA4sub.sumstats.gz")

output_directory<- "/cluster/p/p805/cluster/rosac/fields/output/meta/gwas_by_sub/results/"

# Process each file
for (i in 1:length(input_files)) {
  # Read the original sumstats file
  ss <- fread(input_files[i])

  ss$EAF_A1_prelim<-ss$MAF
  names(ss)[which(names(ss)=="A1")]<-"EA"
  names(ss)[which(names(ss)=="A2")]<-"OA"
  
  # read the munged file w just hm3 SNPs and A1 A2 in the order that's needed.
  ss_munged<-fread(munged_files[i])
  merged<-merge(ss_munged,ss,by = "SNP",all.x=T)

  # create an effect allele frequency column that matches the allele used in the munged file (because we take the frequency from the reference file that might have flipped effect and other alleles)
  merged$EAF<-ifelse(merged$A1 == merged$EA, merged$EAF_A1_prelim, 1-merged$EAF_A1_prelim)

  # only keep required columns
  names(merged)[which(names(merged) == "SNP")]<-"SNPID"
  names(merged)[which(names(merged) == "Z_Estimate")]<-"Z"
  names(merged)[which(names(merged) == "Pval_Estimate")]<-"P"
  merged<-merged[,c("SNPID","CHR","BP","A1","A2","EAF","N","Z","P")] 

  # name according to GWAMA expectations# In genomicSEM munged file A1 is the effect allele
  names(merged)[which(names(merged) == "A1")]<-"EA"
  names(merged)[which(names(merged) == "A2")]<-"OA"
  
  output_file <- file.path(output_directory, paste0(tools::file_path_sans_ext(basename(input_files[i])), "gwama.txt"))

  # Write the gzipped file
  write.table(merged, file = output_file, quote=F,col.names=T,row.names=F)
  # system(paste("gzip -f", shQuote(output_file)))

}
