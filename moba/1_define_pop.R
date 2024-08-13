
library(data.table)
library(dplyr)
library(psych)
library(stringr)
library(ggplot2)
library(data.table)
library(knitr)
library(fastDummies)
library(ggplot2)
library(readxl)

# read in fixed information from pop register
W19_0634_FASTE_OPPLYSNINGER <- read.csv("N:/durable/data/registers/original 2022/csv/W19_0634_FASTE_OPPLYSNINGER.csv",
                                          sep = "," , header = T,na.strings ="", stringsAsFactors= F)
# drop those without parent IDs
pop<-W19_0634_FASTE_OPPLYSNINGER[!(is.na(W19_0634_FASTE_OPPLYSNINGER$lopenr_mor) | W19_0634_FASTE_OPPLYSNINGER$lopenr_mor==""), ]
pop2<-pop[!(is.na(pop$lopenr_far) | pop$lopenr_far==""), ]

# calculate YOB 
fixed<-transform(pop2, YOB = substr(foedsels_aar_mnd, 1, 4), MonthOB = substr(foedsels_aar_mnd, 5, 6))
# calculate YODeath
fixed2<-transform(fixed, YOD = substr(doeds_aar_mnd, 1, 4), MonthOD = substr(doeds_aar_mnd, 5, 6))
fixed2$YOB<-as.numeric(fixed2$YOB)
fixed2$YOD<-as.numeric(fixed2$YOD)

# survive older than 25 years old, meaning we are not including individuals younger then 25 in 2021 (end of the registers)
# keep if survived to age 25
fixed2$early_death<-ifelse(fixed2$YOD-fixed2$YOB<25,1,0)
fixed3<-subset(fixed2,fixed2$early_death==0|is.na(fixed2$early_death))
fixed33<-subset(fixed3,fixed3$age_2018>=25)

# **********************************************************************************

# find adults with available EA data (NUS version, not converting to ISCED at this stage)
edu_dat <- read.csv("N:/durable/data/registers/original/csv/w19_0634_utd_1970_2018_ut.csv")
edu_dat_nodup <- edu_dat[!duplicated(edu_dat$w19_0634_lnr),] 
# **********************************************************************************

# assign Edu to proband mum dad
edu<-edu_dat_nodup[,c("w19_0634_lnr","BU_2018" )]
edu_far<-edu
edu_mor<-edu
colnames(edu_far)<-c("lopenr_far","BU_2018_far")
colnames(edu_mor)<-c("lopenr_mor","BU_2018_mor" )


fixed5<-fixed33%>%
  left_join(edu, by='w19_0634_lnr') %>%
  left_join(edu_far, by='lopenr_far') %>%
  left_join(edu_mor, by='lopenr_mor')

# **********************************************************************************

# remove probands with missing  eduyears
fixed66 <- fixed5[complete.cases(fixed5[ , c('BU_2018')]), ] 
fixed7<-fixed66[,c("w19_0634_lnr","lopenr_mor","lopenr_far","foedsels_aar_mnd"  , 
                  "kjoenn","YOB",
                  "BU_2018","BU_2018_far","BU_2018_mor","foede_kommnune" )]
# **********************************************************************************
# write population file
write.table(fixed7, "fixed_edu.txt", quote=F,col.names = T, row.names=F)


