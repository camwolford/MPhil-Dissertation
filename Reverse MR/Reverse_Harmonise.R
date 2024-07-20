### Harmonise IVs ###
# This script is used to finalise the clumped, matched IVs for each metabolite by filtering by harmonising alleles.
library(tidyverse)
library(TwoSampleMR)

# Load the data
T2DM_metabolite_files <- list.files("Reverse_MR/Matched_Reverse_MR_T2DM_IVs", full.names = TRUE)
FG_metabolite_files <- list.files("Reverse_MR/Matched_Reverse_MR_FG_IVs", full.names = TRUE)
HBA1C_metabolite_files <- list.files("Reverse_MR/Matched_Reverse_MR_HBA1C_IVs", full.names = TRUE)
# Combine the metabolite files
metabolite_files <- c(T2DM_metabolite_files, FG_metabolite_files, HBA1C_metabolite_files)

# Initialise a metabolite counter
metabolite_counter <- 0
# Initialise a remove counter
remove_counter <- 0
# Initialise a palindromic counter
palindromic_counter <- 0
# Initialise an ambiguous counter
ambiguous_counter <- 0
# Initialise an MR-keep counter
mr_keep_counter <- 0
# Initialise a multiple IVs counter
multiple_IVs <- 0

for (metabolite_file in metabolite_files) {
  # Load the metabolite data
  file_path <- metabolite_file
  metabolite_data <- readr::read_tsv(file_path)
  
  metabolite_data$EffectAllele <- ifelse(metabolite_data$EffectAllele == TRUE, "T", metabolite_data$EffectAllele)
  metabolite_data$NonEffectAllele <- ifelse(metabolite_data$NonEffectAllele == TRUE, "T", metabolite_data$NonEffectAllele)
  metabolite_data$out_EffectAllele <- ifelse(metabolite_data$out_EffectAllele == TRUE, "T", metabolite_data$out_EffectAllele)
  metabolite_data$out_NonEffectAllele <- ifelse(metabolite_data$out_NonEffectAllele == TRUE, "T", metabolite_data$out_NonEffectAllele)

  # If there are NAs in any of the rows, remove the row
  metabolite_data <- metabolite_data[!is.na(metabolite_data$EAF), ]
  
  # Make new dataframes with the column names for harmonising alleles using TwoSampleMR
  exposure_data <- metabolite_data %>% select(SNP, Beta, SE, EffectAllele, NonEffectAllele, EAF)
  outcome_data <- metabolite_data %>% select(SNP, out_Beta, out_SE, out_EffectAllele, out_NonEffectAllele, out_EAF)

  # Rename the columns
  colnames(exposure_data) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
  colnames(outcome_data) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")
  # Add other required columns
  exposure_data$exposure <- "exposure"
  exposure_data$id.exposure <- "exposure"
  outcome_data$outcome <- "outcome"
  outcome_data$id.outcome <- "outcome"
  
  # Make harmonised_IVs an empty dataframe
  harmonised_IVs <- data.frame()
  
  harmonised_IVs <- harmonise_data(exposure_data, outcome_data, action = 2)
    
  # Check if any of the SNPs should be removed
  palindromic_counter <- palindromic_counter + sum(harmonised_IVs$palindromic)
  ambiguous_counter <- ambiguous_counter + sum(harmonised_IVs$ambiguous)
  mr_keep_counter <- mr_keep_counter + sum(harmonised_IVs$mr_keep)
    
  # Remove all rows with both palindromic and ambiguous SNPs
  remove_counter <- remove_counter + sum(harmonised_IVs$palindromic & harmonised_IVs$ambiguous)
  rsids <- harmonised_IVs$SNP[harmonised_IVs$palindromic & harmonised_IVs$ambiguous]
  metabolite_data <- metabolite_data %>% filter(!SNP %in% rsids)
  exposure_data <- exposure_data %>% filter(!SNP %in% rsids)
  outcome_data <- outcome_data %>% filter(!SNP %in% rsids)
  harmonised_IVs <- harmonised_IVs %>% filter(!SNP %in% rsids)

  for (i in 1:nrow(harmonised_IVs)) {
    rsid <- harmonised_IVs$SNP[i]
    rsid_data <- harmonised_IVs[i, ]
    metabolite_data <- metabolite_data %>% mutate(Beta = ifelse(SNP == rsid, rsid_data$beta.exposure, Beta),
                                                  SE = ifelse(SNP == rsid, rsid_data$se.exposure, SE),
                                                  EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.exposure, EffectAllele),
                                                  NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.exposure, NonEffectAllele),
                                                  EAF = ifelse(SNP == rsid, rsid_data$eaf.exposure, EAF),
                                                  out_Beta = ifelse(SNP == rsid, rsid_data$beta.outcome, out_Beta),
                                                  out_SE = ifelse(SNP == rsid, rsid_data$se.outcome, out_SE),
                                                  out_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.outcome, out_EffectAllele),
                                                  out_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.outcome, out_NonEffectAllele),
                                                  out_EAF = ifelse(SNP == rsid, rsid_data$eaf.outcome, out_EAF))
  }
  multiple_IVs <- multiple_IVs + 1
  # Save the metabolite data
  # Extract the part of the file name before the "_matched"
  metabolite_name <- strsplit(basename(metabolite_file), "matched")[[1]][1]
  # Add "_harmonised_IVs" to the file name
  metabolite_name <- paste0(metabolite_name, "harmonised_IVs")
  # Save the data
  # If the metabolite names contains "T2DM", save the data in the T2DM folder
  if (grepl("T2DM", metabolite_name)) {
    write.table(metabolite_data, paste0("Reverse_MR/Harmonised_Reverse_T2DM_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
  }
  # If the metabolite names contains "FG", save the data in the FG folder
  else if (grepl("FG", metabolite_name)) {
    write.table(metabolite_data, paste0("Reverse_MR/Harmonised_Reverse_FG_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
  }
  # If the metabolite names contrain "HBA1C", save the data in the HBA1C folder
  else if (grepl("HBA1C", metabolite_name)) {
    write.table(metabolite_data, paste0("Reverse_MR/Harmonised_Reverse_HBA1C_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
  }
  metabolite_counter <- metabolite_counter + 1
  
  # Print the progress
  print(paste0("Metabolite ", metabolite_counter, " processed"))
  # Print the number of removed SNPs
  print(paste0("Number of removed SNPs: ", remove_counter))
  # Print the number of palindromic SNPs
  print(paste0("Number of palindromic SNPs: ", palindromic_counter))
  # Print the number of ambiguous SNPs
  print(paste0("Number of ambiguous SNPs: ", ambiguous_counter))
  # Print the number of MR-keep SNPs
  print(paste0("Number of MR-keep SNPs: ", mr_keep_counter))
  print(paste0("Number of Metabolites with multiple IVs: ", multiple_IVs))
}