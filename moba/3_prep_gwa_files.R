# Rosa Cheesman
# SCRIPT TO CREATE FILES FOR SNP-BASED ANALYSES gcta format

library(dplyr)
library(psych)
library(stringr)
library(ggplot2)
library(data.table)
library(knitr)
library(fastDummies)
library(ggplot2)
library(readxl)
library(readr)

# read in the main file with all variables -- only includes genotyped people already. see script 2. 
fields_gencov <- read.csv("N:/durable/projects/Rosa_gwas/fields/ISCED_codes_covs.csv")

  # all# **********************************************************************************
phenos<-fields_gencov[,c("FID","IID", "gen_prog","edu","arts","social","business",
                         "natural_sci","ict","engineering","agri","health","services")]
names(phenos)
describe(phenos)

setwd("N:/durable/projects/Rosa_gwas/fields")
write.table(phenos, "field_phenos.txt", quote=F,col.names = F, row.names=F)

# just male # **********************************************************************************
male_phenos<-subset(fields_gencov, fields_gencov$kjoenn==1)
male_phenos<-male_phenos[,c("FID","IID", "gen_prog","edu","arts","social","business",
                            "natural_sci","ict","engineering","agri","health","services")]
write.table(male_phenos, "male_field_phenos.txt", quote=F,col.names = F, row.names=F)

# just female # **********************************************************************************
female_phenos<-subset(fields_gencov, fields_gencov$kjoenn==2)
female_phenos<-female_phenos[,c("FID","IID", "gen_prog","edu","arts","social","business",
                                "natural_sci","ict","engineering","agri","health","services")]
write.table(female_phenos, "female_field_phenos.txt", quote=F,col.names = F, row.names=F)



# **********************************************************************************

# continuous data:

# only basic
covs<-fields_gencov[,c("FID","IID",
                "age_2023", 
                "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                "PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20",
                "genotyping_center","Plate_id",
                "genotyping_batch","imputation_batch")]

# Eduyears
covs2a<-fields_gencov[,c("FID","IID",
                "age_2023", 
                "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
                "PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20",
                "genotyping_center","Plate_id",
                "genotyping_batch","imputation_batch",
                "EduYears11_2018")]
# **********************************************************************************

# categorical

# basic, just sex:
covs3<-fields_gencov[,c("FID","IID","kjoenn")]
# **********************************************************************************


# Birthplace Komm: create dummy variables
# set munic to NA if only <10 people have it

mobageo <- fields_gencov %>%
  group_by(foede_kommnune) %>%
  mutate(count = n()) %>%
  mutate(komm = ifelse(count < 10, NA, foede_kommnune)) %>%
  ungroup() %>%
  select(-count)


as.data.frame(table(mobageo$komm))
# 216 municiaplities
describe(mobageo$komm)

komm_dummies<-fastDummies::dummy_cols(mobageo, select_columns = "komm")

names(komm_dummies)

# table(komm_dummies$komm_NA)
# 5123 
komm_dummies$komm_NA<-NULL

covs4<-komm_dummies[,c(72,2,98:313)]
names(covs4)

# **********************************************************************************

# Birthplace komm dummies + parent field dummies
covs5<-komm_dummies[,c(72,2,98:313,48:69)]
names(covs5)


# **********************************************************************************


# continuous
write.table(covs, "Contin_cov.txt", quote=F,col.names =F, row.names=F)
write.table(covs2a, "Contin_cov_EA.txt", quote=F,col.names =F, row.names=F)
# categories
write.table(covs3, "Categ_cov.txt", quote=F,col.names =F, row.names=F)
write.table(covs4, "Categ_cov_komm.txt", quote=F,col.names =F, row.names=F)
write.table(covs5, "Categ_cov_Komm_Parfield_Dummies.txt", quote=F,col.names =F, row.names=F)




