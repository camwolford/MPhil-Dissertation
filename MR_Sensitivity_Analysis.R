library(tidyverse)
library(TwoSampleMR)
library(MRPRESSO)
library(mr.raps)

# Load in the full significant results
sig_results <- readr::read_tsv("Full_Significant_Results.tsv")
sig_metabolites <- sig_results$Metabolite

# Load the metabolite_gwas_associations_cleaned.tsv
metabolite_gwas <- readr::read_tsv("metabolite_gwas_associations_cleaned.tsv")
# Load the t2dm_gwas_cleaned.tsv
t2dm_gwas <- readr::read_tsv("t2dm_gwas_cleaned.tsv")
# Add a marker name column to t2dm gwas
t2dm_gwas$MarkerName <- paste0("chr", t2dm_gwas$Chromosome, "_", t2dm_gwas$Position)
# Load in the T2DM_files
T2DM_files <- list.files("Harmonised_T2DM_IVs", full.names = TRUE)
# Load the liberal results
Liberal_results <- readr::read_tsv("Liberal_Analysis/Full_MR_Results_Liberal.tsv")
Liberal_T2DM_files <- list.files("Liberal_Analysis/Harmonised_T2DM_IVs_Liberal", full.names = TRUE)

# Make a data frame for the sensitivity analysis results
metabolite_T2DM_Sensitivity_Analysis <- data.frame(Metabolite = sig_metabolites, R2_Exposure = NA, R2_Outcome = NA, Direction_Flag = NA, Steiger_Pval = NA, Liberal_Flag = NA,
                              MR_PRESSO_Global_Pval = NA, MR_PRESSO_Outlier_Indices = NA, MR_PRESSO_Outlier_Pvals = NA,
                              MR_RAPS_Beta = NA, MR_RAPS_SE = NA, MR_RAPS_Pval = NA)
# Initialise a counter
counter <- 0

for (metabolite in sig_results$Metabolite) {
  # Find the file in T2DM_files that contains the metabolite name
  metabolite_results <- sig_results[sig_results$Metabolite == metabolite,]
  # Check if there is a value in the Fixed_IVW_Pval_T2DM column
  if (!is.na(metabolite_results$Number_of_IVs_T2DM)) {
    # Select the metabolite file from the T2DM files by looking for any file that contains the metabolite name
    metabolite_file <- T2DM_files[str_detect(T2DM_files, fixed(metabolite, ignore_case = TRUE))]
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
    
    #metabolite_data$lor <- exp(metabolite_data$t2dm_Beta)
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

    if (metabolite_results$Number_of_IVs_T2DM > 3) {
      # Perform MR PRESSO
      mr_presso_data <- data.frame(BetaOutcome = metabolite_data$t2dm_Beta, BetaExposure = metabolite_data$Beta, SdOutcome = metabolite_data$t2dm_SE, 
                                   SdExposure = metabolite_data$SE)
      tryCatch(
        {
        mr_presso_results <- mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure", SdOutcome = "SdOutcome", 
                                       SdExposure = "SdExposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                       data = mr_presso_data, NbDistribution = 1000,  SignifThreshold = 0.05)
        metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Global_Pval"] <- mr_presso_results$`MR-PRESSO results`$`Global Test`$Pvalue
        # Make a character string of the outlier indices
        metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Indices"] <- paste(mr_presso_results$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, collapse = ",")
        metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Pvals"] <- mr_presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue
        }, error = function(e) {
          # pass
        })
      # Perform MR RAPS
      mr_raps_results <-mr.raps(mr_presso_data$BetaExposure, mr_presso_data$BetaOutcome, mr_presso_data$SdExposure, mr_presso_data$SdOutcome)
      metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Beta"] <- mr_raps_results$beta.hat
      metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_SE"] <- mr_raps_results$beta.se
      metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Pval"] <- mr_raps_results$beta.p.value
    } else {
      # Check the number of IVs in the liberal results for the metabolite
      tryCatch({
        if (Liberal_results[Liberal_results$Metabolite == metabolite,]$Number_of_IVs_T2DM > 3) {
          liberal_metabolite_file <- Liberal_T2DM_files[str_detect(Liberal_T2DM_files, fixed(metabolite, ignore_case = TRUE))]
          liberal_metabolite_data <- readr::read_tsv(liberal_metabolite_file)
          # Perform MR PRESSO
          mr_presso_data <- data.frame(BetaOutcome = metabolite_data$t2dm_Beta, BetaExposure = metabolite_data$Beta, SdOutcome = metabolite_data$t2dm_SE, 
                                     SdExposure = metabolite_data$SE)
          tryCatch(
            {
            mr_presso_results <- mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure", SdOutcome = "SdOutcome", 
                                         SdExposure = "SdExposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                         data = mr_presso_data, NbDistribution = 1000,  SignifThreshold = 0.05)
            metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Global_Pval"] <- mr_presso_results$`MR-PRESSO results`$`Global Test`$Pvalue
            # Make a character string of the outlier indices
            metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Indices"] <- paste(mr_presso_results$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, collapse = ",")
            metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Pvals"] <- mr_presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue
            }, error = function(e) {
              # pass
            })
          # Perform MR RAPS
          mr_raps_results <-mr.raps(mr_presso_data$BetaExposure, mr_presso_data$BetaOutcome, mr_presso_data$SdExposure, mr_presso_data$SdOutcome)
          metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Beta"] <- mr_raps_results$beta.hat
          metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_SE"] <- mr_raps_results$beta.se
          metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Pval"] <- mr_raps_results$beta.p.value
          metabolite_T2DM_Sensitivity_Analysis[metabolite_T2DM_Sensitivity_Analysis$Metabolite == metabolite, "Liberal_Flag"] <- TRUE
        }
        }, error = function(e) {
          # pass
        })
      }
    counter <- counter + 1
    print(paste0("Processed ", counter, " metabolite-T2DM associations"))
  }
}
# Remove NAs in the Direction_Flag column
metabolite_T2DM_Sensitivity_Analysis <- metabolite_T2DM_Sensitivity_Analysis[!is.na(metabolite_T2DM_Sensitivity_Analysis$Direction_Flag),]

# Save the results
write_tsv(metabolite_T2DM_Sensitivity_Analysis, "Metabolite_T2DM_Sensitivity_Analysis.tsv")





### Do for fasting glucose ###
# Load in the full significant results
sig_results <- readr::read_tsv("Full_Significant_Results.tsv")
sig_metabolites <- sig_results$Metabolite

# Load the metabolite_gwas_associations_cleaned.tsv
metabolite_gwas <- readr::read_tsv("metabolite_gwas_associations_cleaned.tsv")

# Load the fg_gwas_cleaned.tsv
fg_gwas <- readr::read_tsv("fasting_glucose_gwas_cleaned.tsv")
# Add a marker name column to fg gwas
fg_gwas$MarkerName <- paste0("chr", fg_gwas$Chromosome, "_", fg_gwas$Position)
# Load in the FG_files
FG_files <- list.files("Harmonised_FG_IVs", full.names = TRUE)
# Load the liberal results
Liberal_results <- readr::read_tsv("Liberal_Analysis/Full_MR_Results_Liberal.tsv")
Liberal_FG_files <- list.files("Liberal_Analysis/Harmonised_FG_IVs_Liberal", full.names = TRUE)

# Make a data frame for the sensitivity analysis results
metabolite_FG_Sensitivity_Analysis <- data.frame(Metabolite = sig_metabolites, R2_Exposure = NA, R2_Outcome = NA, Direction_Flag = NA, Steiger_Pval = NA, Liberal_Flag = NA,
                                                   MR_PRESSO_Global_Pval = NA, MR_PRESSO_Outlier_Indices = NA, MR_PRESSO_Outlier_Pvals = NA,
                                                   MR_RAPS_Beta = NA, MR_RAPS_SE = NA, MR_RAPS_Pval = NA)
# Initialise a counter
counter <- 0

for (metabolite in sig_results$Metabolite) {
  # Find the file in FG_files that contains the metabolite name
  metabolite_results <- sig_results[sig_results$Metabolite == metabolite,]
  # Check if there is a value in the Fixed_IVW_Pval_FG column
  if (!is.na(metabolite_results$Number_of_IVs_FG)) {
    # Select the metabolite file from the FG files by looking for any file that contains the metabolite name
    metabolite_file <- FG_files[str_detect(FG_files, fixed(metabolite, ignore_case = TRUE))]
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
    
    #metabolite_data$lor <- exp(metabolite_data$fg_Beta)
    metabolite_data$lor <- metabolite_data$fg_Beta
    metabolite_data$prevalence <- 0.063
    
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
    
    if (metabolite_results$Number_of_IVs_FG > 3) {
      # Perform MR PRESSO
      mr_presso_data <- data.frame(BetaOutcome = metabolite_data$fg_Beta, BetaExposure = metabolite_data$Beta, SdOutcome = metabolite_data$fg_SE, 
                                   SdExposure = metabolite_data$SE)
      tryCatch(
        {
          mr_presso_results <- mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure", SdOutcome = "SdOutcome", 
                                         SdExposure = "SdExposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                         data = mr_presso_data, NbDistribution = 1000,  SignifThreshold = 0.05)
          metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Global_Pval"] <- mr_presso_results$`MR-PRESSO results`$`Global Test`$Pvalue
          # Make a character string of the outlier indices
          metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Indices"] <- paste(mr_presso_results$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, collapse = ",")
          metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Pvals"] <- mr_presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue
        }, error = function(e) {
          # pass
        })
      # Perform MR RAPS
      mr_raps_results <-mr.raps(mr_presso_data$BetaExposure, mr_presso_data$BetaOutcome, mr_presso_data$SdExposure, mr_presso_data$SdOutcome)
      metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Beta"] <- mr_raps_results$beta.hat
      metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_SE"] <- mr_raps_results$beta.se
      metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Pval"] <- mr_raps_results$beta.p.value
    } else {
      # Check the number of IVs in the liberal results for the metabolite
      tryCatch({
        if (Liberal_results[Liberal_results$Metabolite == metabolite,]$Number_of_IVs_FG > 3) {
          liberal_metabolite_file <- Liberal_FG_files[str_detect(Liberal_FG_files, fixed(metabolite, ignore_case = TRUE))]
          liberal_metabolite_data <- readr::read_tsv(liberal_metabolite_file)
          # Perform MR PRESSO
          mr_presso_data <- data.frame(BetaOutcome = metabolite_data$fg_Beta, BetaExposure = metabolite_data$Beta, SdOutcome = metabolite_data$fg_SE, 
                                       SdExposure = metabolite_data$SE)
          tryCatch(
            {
              mr_presso_results <- mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure", SdOutcome = "SdOutcome", 
                                             SdExposure = "SdExposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                             data = mr_presso_data, NbDistribution = 1000,  SignifThreshold = 0.05)
              metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Global_Pval"] <- mr_presso_results$`MR-PRESSO results`$`Global Test`$Pvalue
              # Make a character string of the outlier indices
              metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Indices"] <- paste(mr_presso_results$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, collapse = ",")
              metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Pvals"] <- mr_presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue
            }, error = function(e) {
              # pass
            })
          # Perform MR RAPS
          mr_raps_results <-mr.raps(mr_presso_data$BetaExposure, mr_presso_data$BetaOutcome, mr_presso_data$SdExposure, mr_presso_data$SdOutcome)
          metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Beta"] <- mr_raps_results$beta.hat
          metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_SE"] <- mr_raps_results$beta.se
          metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Pval"] <- mr_raps_results$beta.p.value
          metabolite_FG_Sensitivity_Analysis[metabolite_FG_Sensitivity_Analysis$Metabolite == metabolite, "Liberal_Flag"] <- TRUE
        }
      }, error = function(e) {
        # pass
      })
    }
    counter <- counter + 1
    print(paste0("Processed ", counter, " metabolite-FG associations"))
  }
}
# Remove NAs in the Direction_Flag column
metabolite_FG_Sensitivity_Analysis <- metabolite_FG_Sensitivity_Analysis[!is.na(metabolite_FG_Sensitivity_Analysis$Direction_Flag),]

# Save the results
write_tsv(metabolite_FG_Sensitivity_Analysis, "Metabolite_FG_Sensitivity_Analysis.tsv")











### Do for hba1c ###
# Load in the full significant results
sig_results <- readr::read_tsv("Full_Significant_Results.tsv")
sig_metabolites <- sig_results$Metabolite

metabolite_gwas <- readr::read_tsv("metabolite_gwas_associations_cleaned.tsv")

# Load the hba1c_gwas_cleaned.tsv
hba1c_gwas <- readr::read_tsv("HBA1C_gwas_cleaned.tsv")
# Add a marker name column to hba1c gwas
hba1c_gwas$MarkerName <- paste0("chr", hba1c_gwas$Chromosome, "_", hba1c_gwas$Position)
# Load in the hba1c_files
hba1c_files <- list.files("Harmonised_HBA1C_IVs", full.names = TRUE)
# Load the liberal results
Liberal_results <- readr::read_tsv("Liberal_Analysis/Full_MR_Results_Liberal.tsv")
Liberal_hba1c_files <- list.files("Liberal_Analysis/Harmonised_HBA1C_IVs_Liberal", full.names = TRUE)

# Make a data frame for the sensitivity analysis results
metabolite_HBA1C_Sensitivity_Analysis <- data.frame(Metabolite = sig_metabolites, R2_Exposure = NA, R2_Outcome = NA, Direction_Flag = NA, Steiger_Pval = NA, Liberal_Flag = NA,
                                                 MR_PRESSO_Global_Pval = NA, MR_PRESSO_Outlier_Indices = NA, MR_PRESSO_Outlier_Pvals = NA,
                                                 MR_RAPS_Beta = NA, MR_RAPS_SE = NA, MR_RAPS_Pval = NA)
# Initialise a counter
counter <- 0

for (metabolite in sig_results$Metabolite) {
  # Find the file in hba1c_files that contains the metabolite name
  metabolite_results <- sig_results[sig_results$Metabolite == metabolite,]
  # Check if there is a value in the Fixed_IVW_Pval_hba1c column
  if (!is.na(metabolite_results$Number_of_IVs_HBA1C)) {
    # Select the metabolite file from the hba1c files by looking for any file that contains the metabolite name
    metabolite_file <- hba1c_files[str_detect(hba1c_files, fixed(metabolite, ignore_case = TRUE))]
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
    
    #metabolite_data$lor <- exp(metabolite_data$hba1c_Beta)
    metabolite_data$lor <- metabolite_data$hba1c_Beta
    metabolite_data$prevalence <- 0.063
    
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
    
    if (metabolite_results$Number_of_IVs_HBA1C > 3) {
      # Perform MR PRESSO
      mr_presso_data <- data.frame(BetaOutcome = metabolite_data$hba1c_Beta, BetaExposure = metabolite_data$Beta, SdOutcome = metabolite_data$hba1c_SE, 
                                   SdExposure = metabolite_data$SE)
      tryCatch(
        {
          mr_presso_results <- mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure", SdOutcome = "SdOutcome", 
                                         SdExposure = "SdExposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                         data = mr_presso_data, NbDistribution = 1000,  SignifThreshold = 0.05)
          metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Global_Pval"] <- mr_presso_results$`MR-PRESSO results`$`Global Test`$Pvalue
          # Make a character string of the outlier indices
          metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Indices"] <- paste(mr_presso_results$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, collapse = ",")
          metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Pvals"] <- mr_presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue
        }, error = function(e) {
          # pass
        })
      # Perform MR RAPS
      mr_raps_results <-mr.raps(mr_presso_data$BetaExposure, mr_presso_data$BetaOutcome, mr_presso_data$SdExposure, mr_presso_data$SdOutcome)
      metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Beta"] <- mr_raps_results$beta.hat
      metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_SE"] <- mr_raps_results$beta.se
      metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Pval"] <- mr_raps_results$beta.p.value
    } else {
      # Check the number of IVs in the liberal results for the metabolite
      tryCatch({
        if (Liberal_results[Liberal_results$Metabolite == metabolite,]$Number_of_IVs_HBA1C > 3) {
          liberal_metabolite_file <- Liberal_hba1c_files[str_detect(Liberal_hba1c_files, fixed(metabolite, ignore_case = TRUE))]
          liberal_metabolite_data <- readr::read_tsv(liberal_metabolite_file)
          # Perform MR PRESSO
          mr_presso_data <- data.frame(BetaOutcome = metabolite_data$hba1c_Beta, BetaExposure = metabolite_data$Beta, SdOutcome = metabolite_data$hba1c_SE, 
                                       SdExposure = metabolite_data$SE)
          tryCatch(
            {
              mr_presso_results <- mr_presso(BetaOutcome = "BetaOutcome", BetaExposure = "BetaExposure", SdOutcome = "SdOutcome", 
                                             SdExposure = "SdExposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, 
                                             data = mr_presso_data, NbDistribution = 1000,  SignifThreshold = 0.05)
              metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Global_Pval"] <- mr_presso_results$`MR-PRESSO results`$`Global Test`$Pvalue
              # Make a character string of the outlier indices
              metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Indices"] <- paste(mr_presso_results$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, collapse = ",")
              metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_PRESSO_Outlier_Pvals"] <- mr_presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue
            }, error = function(e) {
              # pass
            })
          # Perform MR RAPS
          mr_raps_results <-mr.raps(mr_presso_data$BetaExposure, mr_presso_data$BetaOutcome, mr_presso_data$SdExposure, mr_presso_data$SdOutcome)
          metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Beta"] <- mr_raps_results$beta.hat
          metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_SE"] <- mr_raps_results$beta.se
          metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "MR_RAPS_Pval"] <- mr_raps_results$beta.p.value
          metabolite_HBA1C_Sensitivity_Analysis[metabolite_HBA1C_Sensitivity_Analysis$Metabolite == metabolite, "Liberal_Flag"] <- TRUE
        }
      }, error = function(e) {
        # pass
      })
    }
    counter <- counter + 1
    print(paste0("Processed ", counter, " metabolite-HBA1C associations"))
  }
}
# Remove NAs in the Direction_Flag column
metabolite_HBA1C_Sensitivity_Analysis <- metabolite_HBA1C_Sensitivity_Analysis[!is.na(metabolite_HBA1C_Sensitivity_Analysis$Direction_Flag),]

# Save the results
write_tsv(metabolite_HBA1C_Sensitivity_Analysis, "Metabolite_HBA1C_Sensitivity_Analysis.tsv")










