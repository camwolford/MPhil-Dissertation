### Perform MR ###
# This script is used to performed MR using each MR method on the harmonised data.
library(tidyverse)
library(MendelianRandomization)

# Load the data
metabolite_files <- list.files("Harmonised_T2DM_IVs", full.names = TRUE)
# How many contain ~
metabolite_files <- metabolite_files[!grepl("~", metabolite_files)]

# Initialise a counter for the number of rows equal to the number of metabolite_files
number_of_metabolites <- length(metabolite_files)

results_df <- data.frame(
  Metabolite = rep(NA, number_of_metabolites),
  Number_of_IVs = rep(NA, number_of_metabolites),
  Number_of_Proxies = rep(NA, number_of_metabolites),
  Weighted_Mode_Estimate = rep(NA, number_of_metabolites),
  Weighted_Mode_SE = rep(NA, number_of_metabolites),
  Weighted_Mode_Pval = rep(NA, number_of_metabolites),
  Weighted_Median_Estimate = rep(NA, number_of_metabolites),
  Weighted_Median_SE = rep(NA, number_of_metabolites),
  Weighted_Median_Pval = rep(NA, number_of_metabolites),
  Random_IVW_Estimate = rep(NA, number_of_metabolites),
  Random_IVW_SE = rep(NA, number_of_metabolites),
  Random_IVW_Pval = rep(NA, number_of_metabolites),
  Random_IVW_RSE = rep(NA, number_of_metabolites),
  Random_IVW_HetStat = rep(NA, number_of_metabolites),
  Random_IVW_HetStat_P = rep(NA, number_of_metabolites),
  Random_IVW_FStat = rep(NA, number_of_metabolites),
  Fixed_IVW_Estimate = rep(NA, number_of_metabolites),
  Fixed_IVW_SE = rep(NA, number_of_metabolites),
  Fixed_IVW_Pval = rep(NA, number_of_metabolites),
  Fixed_IVW_RSE = rep(NA, number_of_metabolites),
  Fixed_IVW_HetStat = rep(NA, number_of_metabolites),
  Fixed_IVW_HetStat_P = rep(NA, number_of_metabolites),
  Fixed_IVW_FStat = rep(NA, number_of_metabolites),
  Egger_Estimate = rep(NA, number_of_metabolites),
  Egger_SE = rep(NA, number_of_metabolites),
  Egger_Pval = rep(NA, number_of_metabolites),
  Egger_Intercept = rep(NA, number_of_metabolites),
  Egger_Intercept_SE = rep(NA, number_of_metabolites),
  Egger_Intercept_Pval = rep(NA, number_of_metabolites),
  Egger_RSE = rep(NA, number_of_metabolites),
  Egger_HetStat = rep(NA, number_of_metabolites),
  Egger_HetStat_P = rep(NA, number_of_metabolites),
  Egger_Isq = rep(NA, number_of_metabolites)
)

# Initialise a metabolite counter
metabolite_counter <- 0

for (metabolite_file in metabolite_files) {
  metabolite_counter <- metabolite_counter + 1
  metabolite_data <- readr::read_tsv(metabolite_file)
  
  # Store the metabolite name
  results_df$Metabolite[metabolite_counter] <- metabolite_data$Metabolite[1]
  # Store the number of IVs
  results_df$Number_of_IVs[metabolite_counter] <- nrow(metabolite_data)
  # Count the number of proxies by the number of rows where the Proxy column is not NA
  number_of_proxies <- sum(!is.na(metabolite_data$proxy))
  # Store the number of proxies
  results_df$Number_of_Proxies[metabolite_counter] <- number_of_proxies
  
  # Check the number of rows in the data
  if (nrow(metabolite_data) > 2) {
    MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$t2dm_Beta, byse = metabolite_data$t2dm_SE)
    MR_weighted_mode_out <- MendelianRandomization::mr_mbe(MRObject, weighting = "weighted")
    MR_weighted_median_out <- MendelianRandomization::mr_median(MRObject, weighting = "weighted")
    MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
    MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")
    MR_egger_out <- MendelianRandomization::mr_egger(MRObject)
    
    # Store the results
    results_df$Weighted_Mode_Estimate[metabolite_counter] <- MR_weighted_mode_out$Estimate
    results_df$Weighted_Mode_SE[metabolite_counter] <- MR_weighted_mode_out$StdError
    results_df$Weighted_Mode_Pval[metabolite_counter] <- MR_weighted_mode_out$Pvalue
    results_df$Weighted_Median_Estimate[metabolite_counter] <- MR_weighted_median_out$Estimate
    results_df$Weighted_Median_SE[metabolite_counter] <- MR_weighted_median_out$StdError
    results_df$Weighted_Median_Pval[metabolite_counter] <- MR_weighted_median_out$Pvalue
    results_df$Random_IVW_Estimate[metabolite_counter] <- MR_random_IVW_out$Estimate
    results_df$Random_IVW_SE[metabolite_counter] <- MR_random_IVW_out$StdError
    results_df$Random_IVW_Pval[metabolite_counter] <- MR_random_IVW_out$Pvalue
    results_df$Random_IVW_RSE[metabolite_counter] <- MR_random_IVW_out$RSE
    results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[1]
    results_df$Random_IVW_HetStat_P[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[2]
    results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
    results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
    results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
    results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
    results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
    results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[1]
    results_df$Fixed_IVW_HetStat_P[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[2]
    results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
    results_df$Egger_Estimate[metabolite_counter] <- MR_egger_out$Estimate
    results_df$Egger_SE[metabolite_counter] <- MR_egger_out$StdError.Est
    results_df$Egger_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Est
    results_df$Egger_Intercept[metabolite_counter] <- MR_egger_out$Intercept
    results_df$Egger_Intercept_SE[metabolite_counter] <- MR_egger_out$StdError.Int
    results_df$Egger_Intercept_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Int
    results_df$Egger_RSE[metabolite_counter] <- MR_egger_out$RSE
    results_df$Egger_HetStat[metabolite_counter] <- MR_egger_out$Heter.Stat[1]
    results_df$Egger_HetStat_P[metabolite_counter] <- MR_egger_out$Heter.Stat[2]
    results_df$Egger_Isq[metabolite_counter] <- MR_egger_out$I.sq
  }
  
  # If there are two IVs, perform the two IV test
  if (nrow(metabolite_data) == 2) {
    MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$t2dm_Beta, byse = metabolite_data$t2dm_SE)
    MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
    MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")

    # Store the results
    results_df$Random_IVW_Estimate[metabolite_counter] <- MR_random_IVW_out$Estimate
    results_df$Random_IVW_SE[metabolite_counter] <- MR_random_IVW_out$StdError
    results_df$Random_IVW_Pval[metabolite_counter] <- MR_random_IVW_out$Pvalue
    results_df$Random_IVW_RSE[metabolite_counter] <- MR_random_IVW_out$RSE
    results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[1]
    results_df$Random_IVW_HetStat_P[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[2]
    results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
    results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
    results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
    results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
    results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
    results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[1]
    results_df$Fixed_IVW_HetStat_P[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[2]
    results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
  }
  
  # If there is only one row in the data, perform the one IV test
  if (nrow(metabolite_data) == 1) {
    MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$t2dm_Beta, byse = metabolite_data$t2dm_SE)
    MR_WR_out <- MendelianRandomization::mr_ivw(MRObject)
    results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_WR_out$Estimate
    results_df$Fixed_IVW_SE[metabolite_counter] <- MR_WR_out$StdError
    results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_WR_out$Pvalue
    results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_WR_out$Fstat
  }
  # Log the progress
  print(paste("Metabolite", metabolite_counter, "completed"))
}

# Save the results as a tsv file
write.table(results_df, file = "T2DM_MR_Results.tsv", sep = "\t", row.names = FALSE)



### Do the same for Fasting Glucose ###  
# Load the data
metabolite_files <- list.files("Harmonised_FG_IVs", full.names = TRUE)
# How many contain ~
metabolite_files <- metabolite_files[!grepl("~", metabolite_files)]

# Initialise a counter for the number of rows equal to the number of metabolite_files
number_of_metabolites <- length(metabolite_files)

results_df <- data.frame(
  Metabolite = rep(NA, number_of_metabolites),
  Number_of_IVs = rep(NA, number_of_metabolites),
  Number_of_Proxies = rep(NA, number_of_metabolites),
  Weighted_Mode_Estimate = rep(NA, number_of_metabolites),
  Weighted_Mode_SE = rep(NA, number_of_metabolites),
  Weighted_Mode_Pval = rep(NA, number_of_metabolites),
  Weighted_Median_Estimate = rep(NA, number_of_metabolites),
  Weighted_Median_SE = rep(NA, number_of_metabolites),
  Weighted_Median_Pval = rep(NA, number_of_metabolites),
  Random_IVW_Estimate = rep(NA, number_of_metabolites),
  Random_IVW_SE = rep(NA, number_of_metabolites),
  Random_IVW_Pval = rep(NA, number_of_metabolites),
  Random_IVW_RSE = rep(NA, number_of_metabolites),
  Random_IVW_HetStat = rep(NA, number_of_metabolites),
  Random_IVW_HetStat_P = rep(NA, number_of_metabolites),
  Random_IVW_FStat = rep(NA, number_of_metabolites),
  Fixed_IVW_Estimate = rep(NA, number_of_metabolites),
  Fixed_IVW_SE = rep(NA, number_of_metabolites),
  Fixed_IVW_Pval = rep(NA, number_of_metabolites),
  Fixed_IVW_RSE = rep(NA, number_of_metabolites),
  Fixed_IVW_HetStat = rep(NA, number_of_metabolites),
  Fixed_IVW_HetStat_P = rep(NA, number_of_metabolites),
  Fixed_IVW_FStat = rep(NA, number_of_metabolites),
  Egger_Estimate = rep(NA, number_of_metabolites),
  Egger_SE = rep(NA, number_of_metabolites),
  Egger_Pval = rep(NA, number_of_metabolites),
  Egger_Intercept = rep(NA, number_of_metabolites),
  Egger_Intercept_SE = rep(NA, number_of_metabolites),
  Egger_Intercept_Pval = rep(NA, number_of_metabolites),
  Egger_RSE = rep(NA, number_of_metabolites),
  Egger_HetStat = rep(NA, number_of_metabolites),
  Egger_HetStat_P = rep(NA, number_of_metabolites),
  Egger_Isq = rep(NA, number_of_metabolites)
)

# Initialise a metabolite counter
metabolite_counter <- 0

for (metabolite_file in metabolite_files) {
  metabolite_counter <- metabolite_counter + 1
  metabolite_data <- readr::read_tsv(metabolite_file)
  
  # Store the metabolite name
  results_df$Metabolite[metabolite_counter] <- metabolite_data$Metabolite[1]
  # Store the number of IVs
  results_df$Number_of_IVs[metabolite_counter] <- nrow(metabolite_data)
  # Count the number of proxies by the number of rows where the Proxy column is not NA
  number_of_proxies <- sum(!is.na(metabolite_data$proxy))
  # Store the number of proxies
  results_df$Number_of_Proxies[metabolite_counter] <- number_of_proxies
  
  # Check the number of rows in the data
  if (nrow(metabolite_data) > 2) {
    MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$fg_Beta, byse = metabolite_data$fg_SE)
    MR_weighted_mode_out <- MendelianRandomization::mr_mbe(MRObject, weighting = "weighted")
    MR_weighted_median_out <- MendelianRandomization::mr_median(MRObject, weighting = "weighted")
    MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
    MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")
    MR_egger_out <- MendelianRandomization::mr_egger(MRObject)
    
    # Store the results
    results_df$Weighted_Mode_Estimate[metabolite_counter] <- MR_weighted_mode_out$Estimate
    results_df$Weighted_Mode_SE[metabolite_counter] <- MR_weighted_mode_out$StdError
    results_df$Weighted_Mode_Pval[metabolite_counter] <- MR_weighted_mode_out$Pvalue
    results_df$Weighted_Median_Estimate[metabolite_counter] <- MR_weighted_median_out$Estimate
    results_df$Weighted_Median_SE[metabolite_counter] <- MR_weighted_median_out$StdError
    results_df$Weighted_Median_Pval[metabolite_counter] <- MR_weighted_median_out$Pvalue
    results_df$Random_IVW_Estimate[metabolite_counter] <- MR_random_IVW_out$Estimate
    results_df$Random_IVW_SE[metabolite_counter] <- MR_random_IVW_out$StdError
    results_df$Random_IVW_Pval[metabolite_counter] <- MR_random_IVW_out$Pvalue
    results_df$Random_IVW_RSE[metabolite_counter] <- MR_random_IVW_out$RSE
    results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[1]
    results_df$Random_IVW_HetStat_P[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[2]
    results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
    results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
    results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
    results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
    results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
    results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[1]
    results_df$Fixed_IVW_HetStat_P[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[2]
    results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
    results_df$Egger_Estimate[metabolite_counter] <- MR_egger_out$Estimate
    results_df$Egger_SE[metabolite_counter] <- MR_egger_out$StdError.Est
    results_df$Egger_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Est
    results_df$Egger_Intercept[metabolite_counter] <- MR_egger_out$Intercept
    results_df$Egger_Intercept_SE[metabolite_counter] <- MR_egger_out$StdError.Int
    results_df$Egger_Intercept_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Int
    results_df$Egger_RSE[metabolite_counter] <- MR_egger_out$RSE
    results_df$Egger_HetStat[metabolite_counter] <- MR_egger_out$Heter.Stat[1]
    results_df$Egger_HetStat_P[metabolite_counter] <- MR_egger_out$Heter.Stat[2]
    results_df$Egger_Isq[metabolite_counter] <- MR_egger_out$I.sq
  }
  
  # If there are two IVs, perform the two IV test
  if (nrow(metabolite_data) == 2) {
    MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$fg_Beta, byse = metabolite_data$fg_SE)
    MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
    MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")
    
    # Store the results
    results_df$Random_IVW_Estimate[metabolite_counter] <- MR_random_IVW_out$Estimate
    results_df$Random_IVW_SE[metabolite_counter] <- MR_random_IVW_out$StdError
    results_df$Random_IVW_Pval[metabolite_counter] <- MR_random_IVW_out$Pvalue
    results_df$Random_IVW_RSE[metabolite_counter] <- MR_random_IVW_out$RSE
    results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[1]
    results_df$Random_IVW_HetStat_P[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[2]
    results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
    results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
    results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
    results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
    results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
    results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[1]
    results_df$Fixed_IVW_HetStat_P[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[2]
    results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
  }
  
  # If there is only one row in the data, perform the one IV test
  if (nrow(metabolite_data) == 1) {
    MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$fg_Beta, byse = metabolite_data$fg_SE)
    MR_WR_out <- MendelianRandomization::mr_ivw(MRObject)
    results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_WR_out$Estimate
    results_df$Fixed_IVW_SE[metabolite_counter] <- MR_WR_out$StdError
    results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_WR_out$Pvalue
    results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_WR_out$Fstat
  }
  # Log the progress
  print(paste("Metabolite", metabolite_counter, "completed"))
}

# Save the results as a tsv file
write.table(results_df, file = "FG_MR_Results.tsv", sep = "\t", row.names = FALSE)



### Do the same for HbA1c ###  
# Load the data
metabolite_files <- list.files("Harmonised_HBA1C_IVs", full.names = TRUE)
# How many contain ~
metabolite_files <- metabolite_files[!grepl("~", metabolite_files)]

# Initialise a counter for the number of rows equal to the number of metabolite_files
number_of_metabolites <- length(metabolite_files)

results_df <- data.frame(
  Metabolite = rep(NA, number_of_metabolites),
  Number_of_IVs = rep(NA, number_of_metabolites),
  Number_of_Proxies = rep(NA, number_of_metabolites),
  Weighted_Mode_Estimate = rep(NA, number_of_metabolites),
  Weighted_Mode_SE = rep(NA, number_of_metabolites),
  Weighted_Mode_Pval = rep(NA, number_of_metabolites),
  Weighted_Median_Estimate = rep(NA, number_of_metabolites),
  Weighted_Median_SE = rep(NA, number_of_metabolites),
  Weighted_Median_Pval = rep(NA, number_of_metabolites),
  Random_IVW_Estimate = rep(NA, number_of_metabolites),
  Random_IVW_SE = rep(NA, number_of_metabolites),
  Random_IVW_Pval = rep(NA, number_of_metabolites),
  Random_IVW_RSE = rep(NA, number_of_metabolites),
  Random_IVW_HetStat = rep(NA, number_of_metabolites),
  Random_IVW_HetStat_P = rep(NA, number_of_metabolites),
  Random_IVW_FStat = rep(NA, number_of_metabolites),
  Fixed_IVW_Estimate = rep(NA, number_of_metabolites),
  Fixed_IVW_SE = rep(NA, number_of_metabolites),
  Fixed_IVW_Pval = rep(NA, number_of_metabolites),
  Fixed_IVW_RSE = rep(NA, number_of_metabolites),
  Fixed_IVW_HetStat = rep(NA, number_of_metabolites),
  Fixed_IVW_HetStat_P = rep(NA, number_of_metabolites),
  Fixed_IVW_FStat = rep(NA, number_of_metabolites),
  Egger_Estimate = rep(NA, number_of_metabolites),
  Egger_SE = rep(NA, number_of_metabolites),
  Egger_Pval = rep(NA, number_of_metabolites),
  Egger_Intercept = rep(NA, number_of_metabolites),
  Egger_Intercept_SE = rep(NA, number_of_metabolites),
  Egger_Intercept_Pval = rep(NA, number_of_metabolites),
  Egger_RSE = rep(NA, number_of_metabolites),
  Egger_HetStat = rep(NA, number_of_metabolites),
  Egger_HetStat_P = rep(NA, number_of_metabolites),
  Egger_Isq = rep(NA, number_of_metabolites)
)

# Initialise a metabolite counter
metabolite_counter <- 0

for (metabolite_file in metabolite_files) {
  metabolite_counter <- metabolite_counter + 1
  metabolite_data <- readr::read_tsv(metabolite_file)
  
  # Store the metabolite name
  results_df$Metabolite[metabolite_counter] <- metabolite_data$Metabolite[1]
  # Store the number of IVs
  results_df$Number_of_IVs[metabolite_counter] <- nrow(metabolite_data)
  # Count the number of proxies by the number of rows where the Proxy column is not NA
  number_of_proxies <- sum(!is.na(metabolite_data$proxy))
  # Store the number of proxies
  results_df$Number_of_Proxies[metabolite_counter] <- number_of_proxies
  
  # Check the number of rows in the data
  if (nrow(metabolite_data) > 2) {
    MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$hba1c_Beta, byse = metabolite_data$hba1c_SE)
    MR_weighted_mode_out <- MendelianRandomization::mr_mbe(MRObject, weighting = "weighted")
    MR_weighted_median_out <- MendelianRandomization::mr_median(MRObject, weighting = "weighted")
    MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
    MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")
    MR_egger_out <- MendelianRandomization::mr_egger(MRObject)
    
    # Store the results
    results_df$Weighted_Mode_Estimate[metabolite_counter] <- MR_weighted_mode_out$Estimate
    results_df$Weighted_Mode_SE[metabolite_counter] <- MR_weighted_mode_out$StdError
    results_df$Weighted_Mode_Pval[metabolite_counter] <- MR_weighted_mode_out$Pvalue
    results_df$Weighted_Median_Estimate[metabolite_counter] <- MR_weighted_median_out$Estimate
    results_df$Weighted_Median_SE[metabolite_counter] <- MR_weighted_median_out$StdError
    results_df$Weighted_Median_Pval[metabolite_counter] <- MR_weighted_median_out$Pvalue
    results_df$Random_IVW_Estimate[metabolite_counter] <- MR_random_IVW_out$Estimate
    results_df$Random_IVW_SE[metabolite_counter] <- MR_random_IVW_out$StdError
    results_df$Random_IVW_Pval[metabolite_counter] <- MR_random_IVW_out$Pvalue
    results_df$Random_IVW_RSE[metabolite_counter] <- MR_random_IVW_out$RSE
    results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[1]
    results_df$Random_IVW_HetStat_P[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[2]
    results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
    results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
    results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
    results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
    results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
    results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[1]
    results_df$Fixed_IVW_HetStat_P[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[2]
    results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
    results_df$Egger_Estimate[metabolite_counter] <- MR_egger_out$Estimate
    results_df$Egger_SE[metabolite_counter] <- MR_egger_out$StdError.Est
    results_df$Egger_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Est
    results_df$Egger_Intercept[metabolite_counter] <- MR_egger_out$Intercept
    results_df$Egger_Intercept_SE[metabolite_counter] <- MR_egger_out$StdError.Int
    results_df$Egger_Intercept_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Int
    results_df$Egger_RSE[metabolite_counter] <- MR_egger_out$RSE
    results_df$Egger_HetStat[metabolite_counter] <- MR_egger_out$Heter.Stat[1]
    results_df$Egger_HetStat_P[metabolite_counter] <- MR_egger_out$Heter.Stat[2]
    results_df$Egger_Isq[metabolite_counter] <- MR_egger_out$I.sq
  }
  
  # If there are two IVs, perform the two IV test
  if (nrow(metabolite_data) == 2) {
    MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$hba1c_Beta, byse = metabolite_data$hba1c_SE)
    MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
    MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")
    
    # Store the results
    results_df$Random_IVW_Estimate[metabolite_counter] <- MR_random_IVW_out$Estimate
    results_df$Random_IVW_SE[metabolite_counter] <- MR_random_IVW_out$StdError
    results_df$Random_IVW_Pval[metabolite_counter] <- MR_random_IVW_out$Pvalue
    results_df$Random_IVW_RSE[metabolite_counter] <- MR_random_IVW_out$RSE
    results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[1]
    results_df$Random_IVW_HetStat_P[metabolite_counter] <- MR_random_IVW_out$Heter.Stat[2]
    results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
    results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
    results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
    results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
    results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
    results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[1]
    results_df$Fixed_IVW_HetStat_P[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat[2]
    results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
  }
  
  # If there is only one row in the data, perform the one IV test
  if (nrow(metabolite_data) == 1) {
    MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$hba1c_Beta, byse = metabolite_data$hba1c_SE)
    MR_WR_out <- MendelianRandomization::mr_ivw(MRObject)
    results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_WR_out$Estimate
    results_df$Fixed_IVW_SE[metabolite_counter] <- MR_WR_out$StdError
    results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_WR_out$Pvalue
    results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_WR_out$Fstat
  }
  # Log the progress
  print(paste("Metabolite", metabolite_counter, "completed"))
}

# Save the results as a tsv file
write.table(results_df, file = "HBA1C_MR_Results.tsv", sep = "\t", row.names = FALSE)



# Merge all of the results into a single dataframe
# Load the results
T2DM_results <- readr::read_tsv("T2DM_MR_Results.tsv")
FG_results <- readr::read_tsv("FG_MR_Results.tsv")
HBA1C_results <- readr::read_tsv("HBA1C_MR_Results.tsv")

# Make a new dataframe to store the merged results
# First find the number of unique metabolites between the three datasets
unique_metabolites <- unique(c(T2DM_results$Metabolite, FG_results$Metabolite, HBA1C_results$Metabolite))
print(length(unique_metabolites))
number_of_metabolites <- length(unique_metabolites)
# Make a dataframe with the first column as the unique metabolites
full_results_df <- data.frame(Metabolite = unique_metabolites)

# Add the T2DM results to the dataframe
full_results_df <- merge(full_results_df, T2DM_results, by = "Metabolite", all.x = TRUE)
# Add the FG results to the dataframe
full_results_df <- merge(full_results_df, FG_results, by = "Metabolite", all.x = TRUE)
# Add the HBA1C results to the dataframe
full_results_df <- merge(full_results_df, HBA1C_results, by = "Metabolite", all.x = TRUE)
# Rename the columns, any column ending in .x will be the T2DM results, any column ending in .y will be the FG results
colnames(full_results_df) <- gsub("\\.x", "_T2DM", colnames(full_results_df))
colnames(full_results_df) <- gsub("\\.y", "_FG", colnames(full_results_df))
# Add _HBA1C to the end of columns 60 onwards
colnames(full_results_df)[60:ncol(full_results_df)] <- paste(colnames(full_results_df)[60:ncol(full_results_df)], "_HBA1C", sep = "")

# Save the results as a tsv file
write.table(full_results_df, file = "Full_MR_Results.tsv", sep = "\t", row.names = FALSE)
