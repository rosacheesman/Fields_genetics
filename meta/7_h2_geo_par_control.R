

library(devtools)
# install_github("GenomicSEM/GenomicSEM")
library(GenomicSEM)


broad_files<-c("2_rawfield_komm.fastGWA.gz","3_rawfield_komm.fastGWA.gz","4_rawfield_komm.fastGWA.gz","5_rawfield_komm.fastGWA.gz","6_rawfield_komm.fastGWA.gz","7_rawfield_komm.fastGWA.gz","8_rawfield_komm.fastGWA.gz","9_rawfield_komm.fastGWA.gz","10_rawfield_komm.fastGWA.gz","11_rawfield_komm.fastGWA.gz",
               "2_rawfield_komm_parfield.fastGWA.gz","3_rawfield_komm_parfield.fastGWA.gz","4_rawfield_komm_parfield.fastGWA.gz","5_rawfield_komm_parfield.fastGWA.gz","6_rawfield_komm_parfield.fastGWA.gz","7_rawfield_komm_parfield.fastGWA.gz","8_rawfield_komm_parfield.fastGWA.gz","9_rawfield_komm_parfield.fastGWA.gz","10_rawfield_komm_parfield.fastGWA.gz","11_rawfield_komm_parfield.fastGWA.gz",
               "2_rawfield_komm_par_komm2018.fastGWA.gz","3_rawfield_komm_par_komm2018.fastGWA.gz","4_rawfield_komm_par_komm2018.fastGWA.gz","5_rawfield_komm_par_komm2018.fastGWA.gz","6_rawfield_komm_par_komm2018.fastGWA.gz","7_rawfield_komm_par_komm2018.fastGWA.gz","8_rawfield_komm_par_komm2018.fastGWA.gz","9_rawfield_komm_par_komm2018.fastGWA.gz","10_rawfield_komm_par_komm2018.fastGWA.gz","11_rawfield_komm_par_komm2018.fastGWA.gz")

broad_trait_names<-c("edu_komm","arts_komm","social_komm","business_komm","natural_sci_komm","ict_komm","engineering_komm","agri_komm","health_komm","services_komm",
                     "edu_komm_par","arts_komm_par","social_komm_par","business_komm_par","natural_sci_komm_par","ict_komm_par","engineering_komm_par","agri_komm_par","health_komm_par","services_komm_par",
                     "edu_komm2018","arts_komm2018","social_komm2018","business_komm2018","natural_sci_komm2018","ict_komm2018","engineering_komm2018","agri_komm2018","health_komm2018","services_komm2018")
broad_N<-rep(125016,30)

munge(files=broad_files, hm3 = "w_hm3.noMHC.snplist",trait.names=broad_trait_names,N=broad_N)

# ----------------------------------------------
traits<-c("edu.sumstats.gz","arts.sumstats.gz","social.sumstats.gz","business.sumstats.gz","natural_sci.sumstats.gz","ict.sumstats.gz","engineering.sumstats.gz","agri.sumstats.gz","health.sumstats.gz","services.sumstats.gz",
          "edu_komm.sumstats.gz","arts_komm.sumstats.gz","social_komm.sumstats.gz","business_komm.sumstats.gz","natural_sci_komm.sumstats.gz","ict_komm.sumstats.gz","engineering_komm.sumstats.gz","agri_komm.sumstats.gz","health_komm.sumstats.gz","services_komm.sumstats.gz",
            "edu_komm_par.sumstats.gz","arts_komm_par.sumstats.gz","social_komm_par.sumstats.gz","business_komm_par.sumstats.gz","natural_sci_komm_par.sumstats.gz","ict_komm_par.sumstats.gz","engineering_komm_par.sumstats.gz","agri_komm_par.sumstats.gz","health_komm_par.sumstats.gz","services_komm_par.sumstats.gz",
            "edu_komm2018.sumstats.gz","arts_komm2018.sumstats.gz","social_komm2018.sumstats.gz","business_komm2018.sumstats.gz","natural_sci_komm2018.sumstats.gz","ict_komm2018.sumstats.gz","engineering_komm2018.sumstats.gz","agri_komm2018.sumstats.gz","health_komm2018.sumstats.gz","services_komm2018.sumstats.gz")

names<-c("edu","arts","social","business","natural_sci","ict","engineering","agri","health","services",
                 "edu_komm","arts_komm","social_komm","business_komm","natural_sci_komm","ict_komm","engineering_komm","agri_komm","health_komm","services_komm",
                 "edu_komm_par","arts_komm_par","social_komm_par","business_komm_par","natural_sci_komm_par","ict_komm_par","engineering_komm_par","agri_komm_par","health_komm_par","services_komm_par",
                 "edu_komm2018","arts_komm2018","social_komm2018","business_komm2018","natural_sci_komm2018","ict_komm2018","engineering_komm2018","agri_komm2018","health_komm2018","services_komm2018")

sample.prev <- c(0.121  ,0.046  ,0.067  ,0.130  ,0.029  ,0.024  ,0.157  ,0.024  ,0.194  ,0.067  ,
                 0.121  ,0.046  ,0.067  ,0.130  ,0.029  ,0.024  ,0.157  ,0.024  ,0.194  ,0.067  ,
                 0.121  ,0.046  ,0.067  ,0.130  ,0.029  ,0.024  ,0.157  ,0.024  ,0.194  ,0.067  ,
                 0.121  ,0.046  ,0.067  ,0.130  ,0.029  ,0.024  ,0.157  ,0.024  ,0.194  ,0.067  )

population.prev <- c(0.070  ,0.043  ,0.042  ,0.120  ,0.020  ,0.017  ,0.174  ,0.023  ,0.126  ,0.067  ,
                     0.070  ,0.043  ,0.042  ,0.120  ,0.020  ,0.017  ,0.174  ,0.023  ,0.126  ,0.067  ,
                     0.070  ,0.043  ,0.042  ,0.120  ,0.020  ,0.017  ,0.174  ,0.023  ,0.126  ,0.067  ,
                     0.070  ,0.043  ,0.042  ,0.120  ,0.020  ,0.017  ,0.174  ,0.023  ,0.126  ,0.067  )

ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"

Geo_parEA_LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld,names)

save(Geo_parEA_LDSCoutput, file="Geo_par_LDSCoutput.RData")


# .......................................................................................................

setwd("~/Dropbox/PROMENTA/Choice of educational fields/summary_stats/cohorts_original/moba_raw")
load("Geo_par_LDSCoutput.RData")
# .......................................................................................................
# .......................................................................................................

# .......................................................................................................
# List of subjects
subjects <- c("edu","arts","social","business","natural_sci","ict","engineering","agri","health","services")
# List of subjects

# Placeholder for results
results_list <- list()

# Loop through each subject and apply the modeling process
for (subject in subjects) {
  
  # Construct the model string dynamically for each subject
  model <- paste0("
    ", subject, "_komm ~~ b*", subject, "_komm
    ", subject, "_komm_par ~~ c*", subject, "_komm_par
    ", subject, " ~~ d*", subject, "
    ", subject, "_komm_par ~~ ", subject, "_komm + ", subject, "
    ", subject, "_komm ~~ ", subject, "
    
    ", subject, "_komm_par.ratio :=  c/d
    ", subject, "_komm.ratio := b/d
  ")
  
  # Evaluate the model
  model_eval <- usermodel(covstruc = Geo_parEA_LDSCoutput, estimation = "DWLS", model = model)
  
  # Extract results
  significance_komm <- model_eval$results[8,]
  significance_komm_par <- model_eval$results[7,]
  
  # Compute p-value for whether it is different from 1
  significance_komm$p.diff1 <- pchisq(((1 - significance_komm$Unstand_Est) / significance_komm$Unstand_SE)^2, 1, lower = FALSE)
  significance_komm_par$p.diff1 <- pchisq(((1 - significance_komm_par$Unstand_Est) / significance_komm_par$Unstand_SE)^2, 1, lower = FALSE)
  
  # Adjust p-values for multiple comparisons using FDR
  significance_komm$h2ratio.p.fdr <- p.adjust(significance_komm$p.diff1, method = "fdr")
  significance_komm_par$h2ratio.p.fdr <- p.adjust(significance_komm_par$p.diff1, method = "fdr")
  
  # Combine results into a single data frame for this subject
  komm_results <- data.frame(
    subject = subject,
    type = "komm",
    estimate = significance_komm$Unstand_Est,
    se = significance_komm$Unstand_SE,
    p_diff1 = significance_komm$p.diff1,
    h2ratio_p_fdr = significance_komm$h2ratio.p.fdr
  )
  
  komm_par_results <- data.frame(
    subject = subject,
    type = "komm_par",
    estimate = significance_komm_par$Unstand_Est,
    se = significance_komm_par$Unstand_SE,
    p_diff1 = significance_komm_par$p.diff1,
    h2ratio_p_fdr = significance_komm_par$h2ratio.p.fdr
  )
  
  # Append the combined results to the results list
  results_list[[subject]] <- rbind(komm_results, komm_par_results)
}

# Combine all subjects' results into a single data frame
results_df <- do.call(rbind, results_list)

# Display the results
print(results_df)



