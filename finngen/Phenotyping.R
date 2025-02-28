# Load libraries
library(data.table)
library(dplyr)
library(stringr)
library(readr)

### Educational fields ###

# 1) Read in the data files
d <- fread("/finngen/pipeline/finngen_R11/socio_register_1.0/data/finngen_R11_socio_register_1.0.txt.gz", header=TRUE)  # Socio-economic data
t <- fread("/finngen/red/vilpin/ea_pheno_cov.txt.gz", header=TRUE)
allID <- fread("/finngen/library-red/finngen_R11/analysis_covariates/R11_COV_V0.FID.txt.gz") # Individuals to be included in GWAS

# 2) Filter only education records
d_edu <- d %>% 
  filter(CATEGORY == "EDUC")

# 3) Keep only the last record for each individual and remove those whose last record was before age 25
d_edu_last <- d_edu %>%
  arrange(FINNGENID, EVENT_AGE) %>%
  group_by(FINNGENID) %>%
  summarise_all(last) %>%
  filter(EVENT_AGE > 25)

# 4) Pad CODE1 to 4 digits, then keep only its first two digits
d_edu_last$C <- str_pad(as.character(d_edu_last$CODE1), 4, pad="0")
d_edu_last$C <- substr(d_edu_last$C, 0, 2)

# 5) Create a file (matrix) where each column is a phenotype for the GWAS
temp <- model.matrix(~ 0 + C, d_edu_last)
temp <- as.data.frame(temp)

# 6) Join the phenotype matrix to the IDs that have AGE_AT_DEATH_OR_END_OF_FOLLOWUP > 25
final <- allID %>%
  filter(AGE_AT_DEATH_OR_END_OF_FOLLOWUP > 25) %>%
  select(IID) %>%
  left_join(temp, by = c("IID" = "FINNGENID"))

# 7) For individuals with missing educational data, set missing values to 0
final[is.na(final)] <- 0

# 8) Read covariates
cov_file <- fread("/finngen/library-red/finngen_R11/analysis_covariates/R11_COV_V0.FID.txt.gz", header=TRUE)

# 9) Merge covariates with the educational fields
merged_df <- merge(cov_file, final, by = "IID")
merged_df <- merged_df %>%
  select(FID, IID, everything())

### Educational attainment ###

# 10) First define a function to check digits and then construct EA from CODE2
check_number <- function(num) {
  # Convert the number to a string
  num_str <- as.character(num)
  
  # Check if the first digit is between 0 and 2
  first_digit_check <- substr(num_str, 1, 1) %in% c('0', '1', '2')
  
  # Check if the entire number is 91
  full_number_check <- num == 91
  
  # Return TRUE if either condition is TRUE
  return(first_digit_check | full_number_check)
}

x <- d_edu_last$CODE2
d_edu_last$EA <- ifelse(
  check_number(x), 10,
  ifelse(substr(x, 1, 1) == '3' | substr(x, 1, 1) == '4', 13,
    ifelse(substr(x, 1, 1) == '5', 15,
      ifelse(substr(x, 1, 1) == '6', 19,
        ifelse(substr(x, 1, 1) == '7', 20,
          ifelse(substr(x, 1, 1) == '8', 22, 10)
        )
      )
    )
  )
)

# 11) Create a data frame for EA
ea_df <- d_edu_last[, c(1, 9)]

### Merge educational fields and attainment ###

# 12) Rename the first column to "IID" for merging
names(ea_df)[1] <- "IID"

final_ea_df <- merge(merged_df, ea_df, by = "IID", all = TRUE)

# Make sure columns are in the correct order
final_ea_df <- final_ea_df %>%
  select("FID", "IID", everything())

# 13) Write the file for GWAS
setwd("/home/ivm/Desktop/first_GWAS/")
write.table(
  final_ea_df,
  file = gzfile("ea_pheno_cov.txt.gz"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = "NA"
)
