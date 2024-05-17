library(MendelianRandomization)
library(dplyr)

### MR Analysis using the MendelianRandomization package ###

# Load the data
T2DM_metabolite_files <- list.files("Reverse_MR/Harmonised_Reverse_T2DM_IVs", full.names = TRUE)

# Initialise a counter for the number of rows equal to the number of T2DM_metabolite_files
number_of_metabolites <- length(T2DM_metabolite_files)

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
  Random_IVW_FStat = rep(NA, number_of_metabolites),
  Fixed_IVW_Estimate = rep(NA, number_of_metabolites),
  Fixed_IVW_SE = rep(NA, number_of_metabolites),
  Fixed_IVW_Pval = rep(NA, number_of_metabolites),
  Fixed_IVW_RSE = rep(NA, number_of_metabolites),
  Fixed_IVW_HetStat = rep(NA, number_of_metabolites),
  Fixed_IVW_FStat = rep(NA, number_of_metabolites),
  Egger_Estimate = rep(NA, number_of_metabolites),
  Egger_SE = rep(NA, number_of_metabolites),
  Egger_Pval = rep(NA, number_of_metabolites),
  Egger_Intercept = rep(NA, number_of_metabolites),
  Egger_Intercept_SE = rep(NA, number_of_metabolites),
  Egger_Intercept_Pval = rep(NA, number_of_metabolites),
  Egger_RSE = rep(NA, number_of_metabolites),
  Egger_HetStat = rep(NA, number_of_metabolites),
  Egger_Isq = rep(NA, number_of_metabolites)
)

# Initialise a metabolite counter
metabolite_counter <- 0

for (metabolite_file in T2DM_metabolite_files) {
  metabolite_counter <- metabolite_counter + 1
  metabolite_data <- readr::read_tsv(metabolite_file)
  
  # Store the metabolite name
  results_df$Metabolite[metabolite_counter] <- metabolite_data$outcome
  # Store the number of IVs
  results_df$Number_of_IVs[metabolite_counter] <- nrow(metabolite_data)
  # Count the number of proxies by the number of rows where the Proxy column is not NA
  number_of_proxies <- sum(!is.na(metabolite_data$proxy))
  # Store the number of proxies
  results_df$Number_of_Proxies[metabolite_counter] <- number_of_proxies
  
  MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$out_Beta, byse = metabolite_data$out_SE)
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
  results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat
  results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
  results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
  results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
  results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
  results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
  results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat
  results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
  results_df$Egger_Estimate[metabolite_counter] <- MR_egger_out$Estimate
  results_df$Egger_SE[metabolite_counter] <- MR_egger_out$StdError.Est
  results_df$Egger_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Est
  results_df$Egger_Intercept[metabolite_counter] <- MR_egger_out$Intercept
  results_df$Egger_Intercept_SE[metabolite_counter] <- MR_egger_out$StdError.Int
  results_df$Egger_Intercept_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Int
  results_df$Egger_RSE[metabolite_counter] <- MR_egger_out$RSE
  results_df$Egger_HetStat[metabolite_counter] <- MR_egger_out$Heter.Stat
  results_df$Egger_Isq[metabolite_counter] <- MR_egger_out$I.sq
  # Log the progress
  print(paste("Metabolite", metabolite_counter, "completed"))
}

# Save the results as a tsv file
write.table(results_df, file = "Reverse_MR/T2DM_Reverse_MR_Results.tsv", sep = "\t", row.names = FALSE)



# Now for FG-metabolite pairs
FG_metabolite_files <- list.files("Reverse_MR/Harmonised_Reverse_FG_IVs", full.names = TRUE)

# Initialise a counter for the number of rows equal to the number of FG_metabolite_files
number_of_metabolites <- length(FG_metabolite_files)

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
  Random_IVW_FStat = rep(NA, number_of_metabolites),
  Fixed_IVW_Estimate = rep(NA, number_of_metabolites),
  Fixed_IVW_SE = rep(NA, number_of_metabolites),
  Fixed_IVW_Pval = rep(NA, number_of_metabolites),
  Fixed_IVW_RSE = rep(NA, number_of_metabolites),
  Fixed_IVW_HetStat = rep(NA, number_of_metabolites),
  Fixed_IVW_FStat = rep(NA, number_of_metabolites),
  Egger_Estimate = rep(NA, number_of_metabolites),
  Egger_SE = rep(NA, number_of_metabolites),
  Egger_Pval = rep(NA, number_of_metabolites),
  Egger_Intercept = rep(NA, number_of_metabolites),
  Egger_Intercept_SE = rep(NA, number_of_metabolites),
  Egger_Intercept_Pval = rep(NA, number_of_metabolites),
  Egger_RSE = rep(NA, number_of_metabolites),
  Egger_HetStat = rep(NA, number_of_metabolites),
  Egger_Isq = rep(NA, number_of_metabolites)
)

# Initialise a metabolite counter
metabolite_counter <- 0

for (metabolite_file in FG_metabolite_files) {
  metabolite_counter <- metabolite_counter + 1
  metabolite_data <- readr::read_tsv(metabolite_file)
  
  # Store the metabolite name
  results_df$Metabolite[metabolite_counter] <- metabolite_data$outcome
  # Store the number of IVs
  results_df$Number_of_IVs[metabolite_counter] <- nrow(metabolite_data)
  # Count the number of proxies by the number of rows where the Proxy column is not NA
  number_of_proxies <- sum(!is.na(metabolite_data$proxy))
  # Store the number of proxies
  results_df$Number_of_Proxies[metabolite_counter] <- number_of_proxies
  
  MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$out_Beta, byse = metabolite_data$out_SE)
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
  results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat
  results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
  results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
  results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
  results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
  results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
  results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat
  results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
  results_df$Egger_Estimate[metabolite_counter] <- MR_egger_out$Estimate
  results_df$Egger_SE[metabolite_counter] <- MR_egger_out$StdError.Est
  results_df$Egger_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Est
  results_df$Egger_Intercept[metabolite_counter] <- MR_egger_out$Intercept
  results_df$Egger_Intercept_SE[metabolite_counter] <- MR_egger_out$StdError.Int
  results_df$Egger_Intercept_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Int
  results_df$Egger_RSE[metabolite_counter] <- MR_egger_out$RSE
  results_df$Egger_HetStat[metabolite_counter] <- MR_egger_out$Heter.Stat
  results_df$Egger_Isq[metabolite_counter] <- MR_egger_out$I.sq
  # Log the progress
  print(paste("Metabolite", metabolite_counter, "completed"))
}

# Save the results as a tsv file
write.table(results_df, file = "Reverse_MR/FG_Reverse_MR_Results.tsv", sep = "\t", row.names = FALSE)








# Now for HBA1C-metabolite pairs
HBA1C_metabolite_files <- list.files("Reverse_MR/Harmonised_Reverse_HBA1C_IVs", full.names = TRUE)

# Initialise a counter for the number of rows equal to the number of HBA1C_metabolite_files
number_of_metabolites <- length(HBA1C_metabolite_files)

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
  Random_IVW_FStat = rep(NA, number_of_metabolites),
  Fixed_IVW_Estimate = rep(NA, number_of_metabolites),
  Fixed_IVW_SE = rep(NA, number_of_metabolites),
  Fixed_IVW_Pval = rep(NA, number_of_metabolites),
  Fixed_IVW_RSE = rep(NA, number_of_metabolites),
  Fixed_IVW_HetStat = rep(NA, number_of_metabolites),
  Fixed_IVW_FStat = rep(NA, number_of_metabolites),
  Egger_Estimate = rep(NA, number_of_metabolites),
  Egger_SE = rep(NA, number_of_metabolites),
  Egger_Pval = rep(NA, number_of_metabolites),
  Egger_Intercept = rep(NA, number_of_metabolites),
  Egger_Intercept_SE = rep(NA, number_of_metabolites),
  Egger_Intercept_Pval = rep(NA, number_of_metabolites),
  Egger_RSE = rep(NA, number_of_metabolites),
  Egger_HetStat = rep(NA, number_of_metabolites),
  Egger_Isq = rep(NA, number_of_metabolites)
)

# Initialise a metabolite counter
metabolite_counter <- 0

for (metabolite_file in HBA1C_metabolite_files) {
  metabolite_counter <- metabolite_counter + 1
  metabolite_data <- readr::read_tsv(metabolite_file)
  
  # Store the metabolite name
  results_df$Metabolite[metabolite_counter] <- metabolite_data$outcome
  # Store the number of IVs
  results_df$Number_of_IVs[metabolite_counter] <- nrow(metabolite_data)
  # Count the number of proxies by the number of rows where the Proxy column is not NA
  number_of_proxies <- sum(!is.na(metabolite_data$proxy))
  # Store the number of proxies
  results_df$Number_of_Proxies[metabolite_counter] <- number_of_proxies
  
  MRObject = mr_input(bx = metabolite_data$Beta, bxse = metabolite_data$SE, by = metabolite_data$out_Beta, byse = metabolite_data$out_SE)
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
  results_df$Random_IVW_HetStat[metabolite_counter] <- MR_random_IVW_out$Heter.Stat
  results_df$Random_IVW_FStat[metabolite_counter] <- MR_random_IVW_out$Fstat
  results_df$Fixed_IVW_Estimate[metabolite_counter] <- MR_fixed_IVW_out$Estimate
  results_df$Fixed_IVW_SE[metabolite_counter] <- MR_fixed_IVW_out$StdError
  results_df$Fixed_IVW_Pval[metabolite_counter] <- MR_fixed_IVW_out$Pvalue
  results_df$Fixed_IVW_RSE[metabolite_counter] <- MR_fixed_IVW_out$RSE
  results_df$Fixed_IVW_HetStat[metabolite_counter] <- MR_fixed_IVW_out$Heter.Stat
  results_df$Fixed_IVW_FStat[metabolite_counter] <- MR_fixed_IVW_out$Fstat
  results_df$Egger_Estimate[metabolite_counter] <- MR_egger_out$Estimate
  results_df$Egger_SE[metabolite_counter] <- MR_egger_out$StdError.Est
  results_df$Egger_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Est
  results_df$Egger_Intercept[metabolite_counter] <- MR_egger_out$Intercept
  results_df$Egger_Intercept_SE[metabolite_counter] <- MR_egger_out$StdError.Int
  results_df$Egger_Intercept_Pval[metabolite_counter] <- MR_egger_out$Pvalue.Int
  results_df$Egger_RSE[metabolite_counter] <- MR_egger_out$RSE
  results_df$Egger_HetStat[metabolite_counter] <- MR_egger_out$Heter.Stat
  results_df$Egger_Isq[metabolite_counter] <- MR_egger_out$I.sq
  # Log the progress
  print(paste("Metabolite", metabolite_counter, "completed"))
}

# Save the results as a tsv file
write.table(results_df, file = "Reverse_MR/HBA1C_Reverse_MR_Results.tsv", sep = "\t", row.names = FALSE)

