### Steiger Filtering ###
# This script is used to perform the Steiger filtering on the significant results from the MR analysis.
library(tidyverse)
library(TwoSampleMR)

# Load in the full significant results
sig_results <- readr::read_tsv("Full_Significant_Results_diss.tsv")
sig_metabolites <- sig_results$Metabolite

# Load the metabolite_gwas_associations_cleaned.tsv
metabolite_gwas <- readr::read_tsv("metabolite_gwas_associations_cleaned.tsv")
# Load the t2dm_gwas_cleaned.tsv
t2dm_gwas <- readr::read_tsv("t2dm_gwas_cleaned.tsv")
# Add a marker name column to t2dm gwas
t2dm_gwas$MarkerName <- paste0("chr", t2dm_gwas$Chromosome, "_", t2dm_gwas$Position)
# Load in the T2DM_files
T2DM_files <- list.files("Harmonised_T2DM_IVs", full.names = TRUE)

# Make a data frame for the sensitivity analysis results
metabolite_T2DM_Sensitivity_Analysis <- data.frame(Metabolite = sig_metabolites, R2_Exposure = NA, R2_Outcome = NA, Direction_Flag = NA, Steiger_Pval = NA, Liberal_Flag = NA)
# Initialise a counter
counter <- 0

for (metabolite in sig_results$Metabolite) {
  # Find the file in T2DM_files that contains the metabolite name
  metabolite_results <- sig_results[sig_results$Metabolite == metabolite,]
  # Check if there is a value in the Fixed_IVW_Pval_T2DM column
  if (!is.na(metabolite_results$Number_of_IVs_T2DM)) {
    # Select the metabolite file from the files by looking for any file that contains "/metabolite nameT2DM"
    file_name <- paste0("Harmonised_T2DM_IVs/", metabolite, "T2DMT2DM_Harmonised_IVs.tsv")
    metabolite_file <- T2DM_files[grep(file_name, T2DM_files, fixed = TRUE)]
    # Load in the metabolite file
    metabolite_data <- readr::read_tsv(metabolite_file)
    # Add a SNP_Pval column to the cholesterol data and match the SNP column to the metabolite_gwas data to get the SNP_Pval from the metabolite_gwas$Pval
    metabolite_data$SNP_metab_Pval <- metabolite_gwas$Pval[match(metabolite_data$SNP, metabolite_gwas$SNP)]
    # If the SNP_Pval is 0 then set it to 1e-300
    metabolite_data$SNP_metab_Pval[metabolite_data$SNP_metab_Pval == 0] <- 1e-300
    # Add a metab_samplesize column to the metabolite data
    metabolite_data$metab_samplesize <- metabolite_gwas$Samplesize[match(metabolite_data$SNP, metabolite_gwas$SNP)]
    
    metabolite_data$tempMarkerName <- paste0("chr", metabolite_data$t2dm_Chromosome, "_", metabolite_data$t2dm_Position)
    # match using marker name
    metabolite_data$ncase <- t2dm_gwas$Ncases[match(metabolite_data$tempMarkerName, t2dm_gwas$MarkerName)]
    metabolite_data$ncontrol <- t2dm_gwas$Ncontrols[match(metabolite_data$tempMarkerName, t2dm_gwas$MarkerName)]
    metabolite_data$neff <- t2dm_gwas$Neff[match(metabolite_data$tempMarkerName, t2dm_gwas$MarkerName)]
    
    metabolite_data$lor <- metabolite_data$t2dm_Beta
    metabolite_data$prevalence <- 0.063
    
    # Make an exposure_R2 and outcome_R2 column
    metabolite_data$r.exposure <- NA
    metabolite_data$r.outcome <- NA
    
    for (iv in 1:nrow(metabolite_data)) {
      # calculate the exposure R2
      metabolite_data$r.exposure[iv] <- get_r_from_pn(metabolite_data$SNP_metab_Pval[iv], metabolite_data$metab_samplesize[iv])
      # calculate the outcome R2
      metabolite_data$r.outcome[iv] <- get_r_from_lor(metabolite_data$lor[iv], metabolite_data$t2dm_EAF[iv], metabolite_data$ncase[iv],
                                                      metabolite_data$ncontrol[iv], metabolite_data$prevalence[iv])
    }
    
    metabolite_data$id.exposure <- "metabolite"
    metabolite_data$id.outcome <- "T2DM"
    metabolite_data$pval.exposure <- metabolite_data$Pval
    metabolite_data$pval.outcome <- metabolite_data$t2dm_Pval
    metabolite_data$samplesize.exposure <- metabolite_data$metab_samplesize
    metabolite_data$samplesize.outcome <- metabolite_data$neff
    metabolite_data$exposure <- "metabolite"
    metabolite_data$outcome <- "T2DM"
    
    direction_results <- directionality_test(metabolite_data)
    # Add the results to the metabolite_SNPs data frame to the correct metabolite
    metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "R2_Exposure"] <- direction_results$snp_r2.exposure
    metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "R2_Outcome"] <- direction_results$snp_r2.outcome
    metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "Direction_Flag"] <- direction_results$correct_causal_direction
    metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "Steiger_Pval"] <- direction_results$steiger_pval
    
    counter <- counter + 1
    print(paste0("Processed ", counter, " metabolite-T2DM associations"))
  }
}
# Remove NAs in the Direction_Flag column
metabolite_T2DM_Sensitivity_Analysis <- metabolite_T2DM_Sensitivity_Analysis[!is.na(metabolite_T2DM_Sensitivity_Analysis$Direction_Flag),]

# Save the results
write_tsv(metabolite_T2DM_Sensitivity_Analysis, "Metabolite_T2DM_Sensitivity_Analysis_diss.tsv")


### Do for fasting glucose ###
# Load in the full significant results
sig_results <- readr::read_tsv("Full_Significant_Results_diss.tsv")
sig_metabolites <- sig_results$Metabolite

# Load the metabolite_gwas_associations_cleaned.tsv
metabolite_gwas <- readr::read_tsv("metabolite_gwas_associations_cleaned.tsv")

# Load the fg_gwas_cleaned.tsv
fg_gwas <- readr::read_tsv("fasting_glucose_gwas_cleaned.tsv")
# Add a marker name column to fg gwas
fg_gwas$MarkerName <- paste0("chr", fg_gwas$Chromosome, "_", fg_gwas$Position)
# Load in the FG_files
FG_files <- list.files("Harmonised_FG_IVs", full.names = TRUE)

# Make a data frame for the sensitivity analysis results
metabolite_FG_Sensitivity_Analysis <- data.frame(Metabolite = sig_metabolites, R2_Exposure = NA, R2_Outcome = NA, Direction_Flag = NA, Steiger_Pval = NA, Liberal_Flag = NA)
# Initialise a counter
counter <- 0

for (metabolite in sig_results$Metabolite) {
  # Find the file in FG_files that contains the metabolite name
  metabolite_results <- sig_results[sig_results$Metabolite == metabolite,]
  # Check if there is a value in the Fixed_IVW_Pval_FG column
  if (!is.na(metabolite_results$Number_of_IVs_FG)) {
    # Select the metabolite file from the files by looking for any file that contains "/metabolite nameFG"
    file_name <- paste0("Harmonised_FG_IVs/", metabolite, "FGFG_Harmonised_IVs.tsv")
    metabolite_file <- FG_files[grep(file_name, FG_files, fixed = TRUE)]
    # Load in the metabolite file
    metabolite_data <- readr::read_tsv(metabolite_file)
    # Add a SNP_Pval column to the cholesterol data and match the SNP column to the metabolite_gwas data to get the SNP_Pval from the metabolite_gwas$Pval
    metabolite_data$SNP_metab_Pval <- metabolite_gwas$Pval[match(metabolite_data$SNP, metabolite_gwas$SNP)]
    # If the SNP_Pval is 0 then set it to 1e-300
    metabolite_data$SNP_metab_Pval[metabolite_data$SNP_metab_Pval == 0] <- 1e-300
    # Add a metab_samplesize column to the metabolite data
    metabolite_data$metab_samplesize <- metabolite_gwas$Samplesize[match(metabolite_data$SNP, metabolite_gwas$SNP)]
    
    metabolite_data$tempMarkerName <- paste0("chr", metabolite_data$fg_Chromosome, "_", metabolite_data$fg_Position)
    # match using marker name
    metabolite_data$fg_samplesize <- fg_gwas$SampleSize[match(metabolite_data$tempMarkerName, fg_gwas$MarkerName)]
    metabolite_data$outcome_Pval <- fg_gwas$Pval[match(metabolite_data$tempMarkerName, fg_gwas$MarkerName)]

    # Make an exposure_R2 and outcome_R2 column
    metabolite_data$r.exposure <- NA
    metabolite_data$r.outcome <- NA
    
    for (iv in 1:nrow(metabolite_data)) {
      # calculate the exposure R2
      metabolite_data$r.exposure[iv] <- get_r_from_pn(metabolite_data$SNP_metab_Pval[iv], metabolite_data$metab_samplesize[iv])
      # calculate the outcome R2
      metabolite_data$r.outcome[iv] <- get_r_from_pn(metabolite_data$outcome_Pval[iv], metabolite_data$fg_samplesize[iv])
    }
    
    metabolite_data$id.exposure <- "metabolite"
    metabolite_data$id.outcome <- "FG"
    metabolite_data$pval.exposure <- metabolite_data$Pval
    metabolite_data$pval.outcome <- metabolite_data$fg_Pval
    metabolite_data$samplesize.exposure <- metabolite_data$metab_samplesize
    metabolite_data$samplesize.outcome <- metabolite_data$fg_samplesize
    metabolite_data$exposure <- "metabolite"
    metabolite_data$outcome <- "FG"
    
    direction_results <- directionality_test(metabolite_data)
    # Add the results to the metabolite_SNPs data frame to the correct metabolite
    metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "R2_Exposure"] <- direction_results$snp_r2.exposure
    metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "R2_Outcome"] <- direction_results$snp_r2.outcome
    metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "Direction_Flag"] <- direction_results$correct_causal_direction
    metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "Steiger_Pval"] <- direction_results$steiger_pval

    counter <- counter + 1
    print(paste0("Processed ", counter, " metabolite-FG associations"))
  }
}
# Remove NAs in the Direction_Flag column
metabolite_FG_Sensitivity_Analysis <- metabolite_FG_Sensitivity_Analysis[!is.na(metabolite_FG_Sensitivity_Analysis$Direction_Flag),]

# Save the results
write_tsv(metabolite_FG_Sensitivity_Analysis, "Metabolite_FG_Sensitivity_Analysis_diss.tsv")


### Do for hba1c ###
# Load in the full significant results
sig_results <- readr::read_tsv("Full_Significant_Results_diss.tsv")
sig_metabolites <- sig_results$Metabolite

metabolite_gwas <- readr::read_tsv("metabolite_gwas_associations_cleaned.tsv")

# Load the hba1c_gwas_cleaned.tsv
hba1c_gwas <- readr::read_tsv("HBA1C_gwas_cleaned.tsv")
# Add a marker name column to hba1c gwas
hba1c_gwas$MarkerName <- paste0("chr", hba1c_gwas$Chromosome, "_", hba1c_gwas$Position)
# Load in the hba1c_files
hba1c_files <- list.files("Harmonised_HBA1C_IVs", full.names = TRUE)

# Make a data frame for the sensitivity analysis results
metabolite_HBA1C_Sensitivity_Analysis <- data.frame(Metabolite = sig_metabolites, R2_Exposure = NA, R2_Outcome = NA, Direction_Flag = NA, Steiger_Pval = NA, Liberal_Flag = NA)
# Initialise a counter
counter <- 0

for (metabolite in sig_results$Metabolite) {
  # Find the file in hba1c_files that contains the metabolite name
  metabolite_results <- sig_results[sig_results$Metabolite == metabolite,]
  # Check if there is a value in the Fixed_IVW_Pval_hba1c column
  if (!is.na(metabolite_results$Number_of_IVs_HBA1C)) {
    # Select the metabolite file from the hba1c files by looking for any file that contains "/metabolite nameHBA1C"
    file_name <- paste0("Harmonised_HBA1C_IVs/", metabolite, "HBA1CHBA1C_Harmonised_IVs.tsv")
    metabolite_file <- hba1c_files[grep(file_name, hba1c_files, fixed = TRUE)]
    # Load in the metabolite file
    metabolite_data <- readr::read_tsv(metabolite_file)
    # Add a SNP_Pval column to the cholesterol data and match the SNP column to the metabolite_gwas data to get the SNP_Pval from the metabolite_gwas$Pval
    metabolite_data$SNP_metab_Pval <- metabolite_gwas$Pval[match(metabolite_data$SNP, metabolite_gwas$SNP)]
    # If the SNP_Pval is 0 then set it to 1e-300
    metabolite_data$SNP_metab_Pval[metabolite_data$SNP_metab_Pval == 0] <- 1e-300
    # Add a metab_samplesize column to the metabolite data
    metabolite_data$metab_samplesize <- metabolite_gwas$Samplesize[match(metabolite_data$SNP, metabolite_gwas$SNP)]
    
    metabolite_data$tempMarkerName <- paste0("chr", metabolite_data$hba1c_Chromosome, "_", metabolite_data$hba1c_Position)
    # match using marker name
    metabolite_data$hba1c_samplesize <- hba1c_gwas$SampleSize[match(metabolite_data$tempMarkerName, hba1c_gwas$MarkerName)]
    metabolite_data$outcome_Pval <- hba1c_gwas$Pval[match(metabolite_data$tempMarkerName, hba1c_gwas$MarkerName)]
    
    # Make an exposure_R2 and outcome_R2 column
    metabolite_data$r.exposure <- NA
    metabolite_data$r.outcome <- NA
    
    for (iv in 1:nrow(metabolite_data)) {
      # calculate the exposure R2
      metabolite_data$r.exposure[iv] <- get_r_from_pn(metabolite_data$SNP_metab_Pval[iv], metabolite_data$metab_samplesize[iv])
      # calculate the outcome R2
      metabolite_data$r.outcome[iv] <- get_r_from_pn(metabolite_data$outcome_Pval[iv], metabolite_data$hba1c_samplesize[iv])
    }
    
    metabolite_data$id.exposure <- "metabolite"
    metabolite_data$id.outcome <- "HBA1C"
    metabolite_data$pval.exposure <- metabolite_data$Pval
    metabolite_data$pval.outcome <- metabolite_data$hba1c_Pval
    metabolite_data$samplesize.exposure <- metabolite_data$metab_samplesize
    metabolite_data$samplesize.outcome <- metabolite_data$hba1c_samplesize
    metabolite_data$exposure <- "metabolite"
    metabolite_data$outcome <- "HBA1C"
    
    direction_results <- directionality_test(metabolite_data)
    # Add the results to the metabolite_SNPs data frame to the correct metabolite
    metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "R2_Exposure"] <- direction_results$snp_r2.exposure
    metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "R2_Outcome"] <- direction_results$snp_r2.outcome
    metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "Direction_Flag"] <- direction_results$correct_causal_direction
    metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "Steiger_Pval"] <- direction_results$steiger_pval

    counter <- counter + 1
    print(paste0("Processed ", counter, " metabolite-HBA1C associations"))
  }
}
# Remove NAs in the Direction_Flag column
metabolite_HBA1C_Sensitivity_Analysis <- metabolite_HBA1C_Sensitivity_Analysis[!is.na(metabolite_HBA1C_Sensitivity_Analysis$Direction_Flag),]

# Save the results
write_tsv(metabolite_HBA1C_Sensitivity_Analysis, "Metabolite_HBA1C_Sensitivity_Analysis_diss.tsv")