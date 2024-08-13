# script to:
#   convert NUS to isced narrow fields and make dummy variables
#   find individuals who're also genotyped
#   descriptive stats


library(dplyr)
library(psych)
library(stringr)
library(ggplot2)
library(data.table)
library(knitr)
library(fastDummies)
library(ggplot2)
library(readxl)

# **********************************************************************************

# population identified in step 1 
fixed7 <- fread("~/priv_fields/fixed_edu.txt")
head(fixed7$BU_2018)


# **********************************************************************************

# tehse are the NUS codes from SSB website
codes <- read.csv("N:/durable/projects/Rosa_gwas/fields/525.csv", sep=";")
head(codes)
dat<-codes[,c("sourceCode","targetCode","targetName")]
codes_convert <- read.csv("N:/durable/projects/Rosa_gwas/fields/codes_convert.csv", header=T, sep=";")
# to get the broad  codes from the narrow codes.
sample_sizes_130923 <- read.csv("~/priv_fields/sample_sizes_130923.csv", sep=";")
brd<-sample_sizes_130923[,c(1,4)]

codes_convert2<-merge(codes_convert,brd,by = 'Narrow_code',all=T)
codes_convert3<-merge(dat,codes_convert2,by = 'targetName',all=T)
# 0    1    2    3    4    5    6    7    8    9   10 <NA> 
#   77  448 1011  471  453  355  151  828  306  507  379  339
checknas<-subset(codes_convert3,is.na(codes_convert3$Broad_code))
table(checknas$sourceCode,exclude=NULL)

# now fill in the interdisc programmes etc.
# The criterion for determining the leading subject(s) is,
#  the share of learning credits or students' intended learning time
x<-codes_convert3$targetCode
codes_convert3$Broad_field<-ifelse(is.na(x)|x== 9999, NA,
                                   ifelse(!is.na(x) & x == 99|x == 21,0,
                                          ifelse(!is.na(x) & x == 188,1,
                                                 ifelse(!is.na(x) & x == 219|x ==229|x ==288,2,
                                                        ifelse(!is.na(x) & x == 388|x ==319,3,
                                                               ifelse(!is.na(x) & x == 419,4,
                                                                      ifelse(!is.na(x) & x ==588 |x ==599,5,
                                                                             ifelse(!is.na(x) & x == 619|x ==688,6,
                                                                                    ifelse(!is.na(x) & x ==719|x ==799|x ==788 |x ==729,7,
                                                                                           ifelse(!is.na(x) & x == 819|x ==899, 8,
                                                                                                  ifelse(!is.na(x) & x ==919|x ==988 ,9,
                                                                                                         ifelse(!is.na(x) & x == 1088|x ==1019|x ==1022|x ==1039|x==1099,10,  
                                                                                                                NA))))))))))))




# chcekc out remaining missing codes.
codes_convert4<-codes_convert3 %>%
  mutate(Broad_Code = ifelse(!is.na(Broad_code), Broad_code, Broad_field), .keep="unused") 
checknas2<-subset(codes_convert4,is.na(codes_convert4$Broad_Code))
checknas3<-subset(codes_convert4,is.na(codes_convert4$sourceCode))

# remove last three rows: unknown, NA
codes_convert5<-subset(codes_convert4,!is.na(codes_convert4$Broad_Code)&!is.na(codes_convert4$sourceCode))

# now subset to needed variables for making wide datagrame with 
head(codes_convert5)
getwd()
setwd("C:/Users/p805-rosacg/Documents/priv_fields")
write.table(codes_convert5,'NUS_ISCED_Convert.txt',sep = '\t',row.names = F, 
            col.names = T, 
            quote = F)
# **********************************************************************************
# **********************************************************************************
# merge with pop identified b4
fixed7$sourceCode<-fixed7$BU_2018
edu_isced_narrow3<-merge(fixed7,codes_convert5,by ='sourceCode')


# **********************************************************************************

# make broad field dummies 

field_dummies<-fastDummies::dummy_cols(edu_isced_narrow3, select_columns = "Broad_Code")
names(field_dummies)

colnames(field_dummies)[24:34]<-c("gen_prog","edu","arts","social","business",
                                  "natural_sci","ict","engineering","agri","health","services")


# **********************************************************************************
# make field dummies for parents

# get field codes for parents
head(codes_convert5)
codes<-codes_convert5[,c('sourceCode','Broad_Code')]
codes_mor<-codes
codes_far<-codes
colnames(codes_mor)<-c('BU_2018_mor','Broad_Code_mor')
colnames(codes_far)<-c('BU_2018_far','Broad_Code_far')
# merge in
fields_dummies2<-fields_dummies%>%
  left_join(codes_mor, by='BU_2018_mor') %>%
  left_join(codes_far, by='BU_2018_far')

# dummy code

fields_dummies3<-fastDummies::dummy_cols(fields_dummies2, select_columns = "Broad_Code_mor")
fields_dummies4<-fastDummies::dummy_cols(fields_dummies3, select_columns = "Broad_Code_far")
# remove NA dummies
fields_dummies4$Broad_Code_far_NA<-NULL
fields_dummies4$Broad_Code_mor_NA<-NULL
# **********************************************************************************

# Make EA variable to use as covar
fields_dummies4$edu_2018 = as.numeric(str_sub(fields_dummies4$BU_2018,1,1))
#  allow people with missing or unkown education to be counted as 0 rather than NA
# using years of edu coding from GWA consortium instructions
x<-fields_dummies4$edu_2018
fields_dummies4$ISCED11_2018<-ifelse(is.na(x), NA,
                                       ifelse(!is.na(x) & x == 0, 0, # combines ppl with edu up to age 2 & up to age 5
                                              ifelse(!is.na(x) & x == 1, 1, # primary
                                                     ifelse(!is.na(x) & x == 2, 2, # secondary lower
                                                            ifelse(!is.na(x) & x == 4 | x == 3, 3, # secondary upper AND post secondary non tertiary shorter than 2 yrs
                                                                   ifelse(!is.na(x) & x == 5, 5,  # postsecondary vocational 0.5-1.5 or 2 years tertiary--not distinguished so dunno if isced 4 vs 5
                                                                          ifelse(!is.na(x) & x == 6, 6, #undergrad
                                                                                 ifelse(!is.na(x) & x == 7, 7, #MSc/MA
                                                                                        ifelse(!is.na(x) & x == 8, 8, #PhD
                                                                                               ifelse(!is.na(x) & x == 9,0, NA))))))))) )

x<-fields_dummies4$ISCED11_2018
fields_dummies4$EduYears11_2018<-ifelse(is.na(x), NA,
                                          ifelse(!is.na(x) & x == 0, 1, # combines ppl with edu up to age 2 & up to age 5
                                                 ifelse(!is.na(x) & x == 1, 7, # primary
                                                        ifelse(!is.na(x) & x == 2, 10, # secondary lower
                                                               ifelse(!is.na(x) & x == 3, 13,# secondary upper AND post secondary non tertiary shorter than 2 yrs
                                                                      ifelse(!is.na(x) & x == 5, 15,  # postsecondary vocational 0.5-1.5 or 2 years tertiary--not distinguished so dunno if isced 4 vs 5
                                                                             ifelse(!is.na(x) & x == 6, 19, #undergrad
                                                                                    ifelse(!is.na(x) & x == 7, 20, #MSc/MA
                                                                                           ifelse(!is.na(x) & x == 8, 22, #PhD
                                                                                                  NA))))))))) 





# remember TO SCALE EA here




# **********************************************************************************

# merge register data with MoBagenetics IDs

IDs_AUG22_long <- read.csv("N:/durable/projects/Rosa_gwas/IDs_AUG22_long.csv")

ids_gen<-subset(IDs_AUG22_long,select = c("w19_0634_lnr","SENTRIXID","IID","Role"))
ids_gen1<-ids_gen %>% 
  group_by(w19_0634_lnr) %>%
  slice_sample(n = 1)

fields_gen<-merge(fields_dummies4,ids_gen1, by ="w19_0634_lnr",all.x = T )
names(fields_gen)
describe(fields_gen[,c('SENTRIXID',"IID","kjoenn")])

# **********************************************************************************
# subset to geno ppl
cmplt <- fields_gen[complete.cases(fields_gen[ , c('IID')]), ] 

# **********************************************************************************

# add genetic covs

cov <- read.delim("N:/durable/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov-noMoBaIDs.txt")
# use batch dummies!
fields_gencov<-merge(cmplt,cov,by= 'IID',all.x=T)
describe(fields_gencov)
names(fields_gencov)
# **********************************************************************************

# create an overall file with all the info
write.table(fields_gencov, "ISCED_codes_covs.txt", quote=F,col.names = T, row.names=F)
