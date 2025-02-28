library(GenomicSEM)

#  MODELS WITHOUT SNP EFFECTS FOR ALL FIELDS
# First, define the base model template as a function to allow substitution of variables
create_model <- function(var, field) {
   model <- paste0('EA=~NA*EA4 + ', var, '
                   ', field, '=~NA*', var, '

                   ', field, '~~1*', field, '
                   EA~~1*EA
                   EA~~0*', field, '

                   ', var, ' ~~ 0*', var, '
                   EA4~~0*EA4
                   EA4~~0*', var)
  return(model)
}

# List of variables and their corresponding field components
variables <- c("edu", "arts", "social", "business", "natural_sci", "ict", "engineering", "agri", "health", "services")
field_components <- c("EDU", "ARTS", "SOCIAL", "BUSINESS", "NATURAL_SCI", "ICT", "ENGINEERING", "AGRI", "HEALTH", "SERVICES")

# Initialize a list to hold output results
all_results <- list()
  # load(paste0("Fields_raw_EA4_LDSCoutput", i, ".RData"))  # Assuming each RData file is named as "dat_model1.RData", "dat_model2.RData", etc.
load("Fields_raw_EA4_LDSCoutput.RData")

# Loop through each variable and field_component pair
for (i in 1:length(variables)) {
  
  # Create model based on the variable and field_component
  model <- create_model(variables[i], field_components[i])
  
  # Use the model with the corresponding data
  output <- usermodel(LDSCoutput, estimation="DWLS", model=model)
  all_results[[i]] <- output$results
}

# Combine all results into a single data frame
final_output <- do.call(rbind, all_results)

# Write the output
write.table(final_output, 'gsem_subtract_results.csv', row.names=F,col.names=T,quote=F)
# _____________________________________________________________________________________
# _____________________________________________________________________________________
# _____________________________________________________________________________________
# Prep for gwas by subtraction
# need 10 LDSC objects for EA4 + individual field (w/o EA adjustment)

# Define the list of trait pairs
traits_pairs <- c("edu.sumstats.gz", "arts.sumstats.gz", "social.sumstats.gz", "business.sumstats.gz", "natural_sci.sumstats.gz", "ict.sumstats.gz", "engineering.sumstats.gz", "agri.sumstats.gz", "health.sumstats.gz", "services.sumstats.gz")
ea4_trait <- "EA4.sumstats.gz"

# Define the corresponding population prevalences and sample prevalences for the pairs
population_prevs <- c(0.0515, 0.0509, 0.0337, 0.146, 0.0186, 0.0268, 0.226, 0.0327, 0.141, 0.0949)
ea4_population_prev <- NA
sample_prev <- 0.5
ea4_sample_prev <- NA

# Loop through each trait pair, perform LDSC analysis, and save the results
for (i in 1:length(traits_pairs)) {
  # Define the current traits and prevalences
  current_trait <- traits_pairs[i]
  current_population_prev <- population_prevs[i]
  
  # Construct the traits, population prevalence, and sample prevalence vectors
  traits <- c(current_trait, ea4_trait)
  population.prev <- c(current_population_prev, ea4_population_prev)
  sample.prev <- c(sample_prev, ea4_sample_prev)
  
  # Perform LDSC analysis for the current traits pair
  LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld,trait.names=c(current_trait,"EA4"))
  
  # Generate a filename based on the current trait
  filename <- paste("EA4", sub(".sumstats.gz", "", current_trait), "LDSCoutput.RData", sep="_")
  
  # Save the LDSC output to a file
  save(LDSCoutput, file=filename)
}
# _____________________________________________________________________________________

# Prep sumstats for gwas-by-sub: run sumstats() x10 per field+ea4 combo

# Define the list of trait pairs
traits_pairs <- c(
"Cheesman_meta_edu.txt.gz","Cheesman_meta_arts.txt.gz","Cheesman_meta_social.txt.gz","Cheesman_meta_business.txt.gz",
"Cheesman_meta_natural_sci.txt.gz","Cheesman_meta_ict.txt.gz","Cheesman_meta_engineering.txt.gz","Cheesman_meta_agri.txt.gz",
"Cheesman_meta_health.txt.gz","Cheesman_meta_services.txt.gz")
ea4_trait <- "EA4_excl_23andMe_exclMOBA_2022_04_04.meta.gz"
# Define the corresponding population prevalences and sample prevalences for the pairs
# sum of neff:
traits_N <- c(102970	,97262	,69123	,261182	,40072	,50819	,317209	,63834	,292929	,168157)
ea4_N <- 765283
# define names
field_names <- c("edu", "arts", "social", "business", "natural_sci", "ict", "engineering", "agri", "health", "services")
ea4_name <- "EA4_excl_23andMe_exclMOBA_2022_04_04.meta.gz"


ref = "reference.1000G.maf.0.005.txt"
se.logit = c(F,F)
info.filter = 0.6
maf.filter = 0.01
linprobs=c(T,F)
ols=c(F,T)

# Loop through each trait pair, perform sumstats() 
for (i in 1:length(traits_pairs)) {
  # Define the current traits and prevalences
  current_trait <- traits_pairs[i]
  current_N <- traits_N[i]
  current_name <- field_names[i]
  
  # Construct the traits, population prevalence, and sample prevalence vectors
  traits <- c(current_trait, ea4_trait)
  Ns <- c(current_N, ea4_N)
  names <- c(current_name, ea4_name)
  
  # sumstats() command per field
  sumstats<-sumstats(traits, ref, trait.names=names, se.logit, info.filter, maf.filter, OLS=ols,linprob=linprobs,N=Ns,betas=NULL)
  filename <- paste( "EA4", current_name, "Sumstats.RData", sep="_")
  save(sumstats, file=filename)

}
