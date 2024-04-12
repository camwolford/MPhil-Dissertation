### Harmonise IVs ###
# This script is used to finalise the clumped, matched IVs for each metabolite by filtering by harmonising alleles.

library(tidyverse)

# Test this on the first metabolite
# Load the data
metabolite_1 <- readr::read_tsv("Matched_IVs/1-(1-enyl-oleoyl)-GPC (P-18_1)*T2DM_Matched_IVs.tsv")

# Make new dataframes with the column names for harmonising alleles using TwoSampleMR
exposure_data <- metabolite_1 %>% select(SNP, Beta, SE, EffectAllele, NonEffectAllele, EAF)
outcome_data <- metabolite_1 %>% select(SNP, t2dm_Beta, t2dm_SE, t2dm_EffectAllele, t2dm_NonEffectAllele, t2dm_EAF)
# Rename the columns
colnames(exposure_data) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
colnames(outcome_data) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")
# Add other required columns
exposure_data$exposure <- "1-(1-enyl-oleoyl)-GPC (P-18:1)"
exposure_data$id.exposure <- "1-(1-enyl-oleoyl)-GPC (P-18:1)"
outcome_data$outcome <- "T2DM"
outcome_data$id.outcome <- "T2DM"

library(TwoSampleMR)

# Harmonise alleles
harmonised_IVs <- harmonise_data(exposure_data, outcome_data, action = 2)
# For each IV in the harmonised data, print the remove, palindromic, ambigous and mr_keep values
for (i in 1:nrow(harmonised_IVs)) {
  print(paste("Remove:", harmonised_IVs$remove[i]))
  print(paste("Palindromic:", harmonised_IVs$palindromic[i]))
  print(paste("Ambiguous:", harmonised_IVs$ambiguous[i]))
  print(paste("MR-keep:", harmonised_IVs$mr_keep[i]))
}
# Replace the original data with the harmonised data by matching on SNP
for (i in 1:nrow(harmonised_IVs)) {
  rsid <- harmonised_IVs$SNP[i]
  rsid_data <- harmonised_IVs[i, ]
  metabolite_1 <- metabolite_1 %>% mutate(Beta = ifelse(SNP == rsid, rsid_data$beta.exposure, Beta),
                                          SE = ifelse(SNP == rsid, rsid_data$se.exposure, SE),
                                          EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.exposure, EffectAllele),
                                          NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.exposure, NonEffectAllele),
                                          EAF = ifelse(SNP == rsid, rsid_data$eaf.exposure, EAF),
                                          t2dm_Beta = ifelse(SNP == rsid, rsid_data$beta.outcome, t2dm_Beta),
                                          t2dm_SE = ifelse(SNP == rsid, rsid_data$se.outcome, t2dm_SE),
                                          t2dm_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.outcome, t2dm_EffectAllele),
                                          t2dm_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.outcome, t2dm_NonEffectAllele),
                                          t2dm_EAF = ifelse(SNP == rsid, rsid_data$eaf.outcome, t2dm_EAF))
}


### Do this for all metabolites ###
# Load the data
metabolite_files <- list.files("Matched_IVs", full.names = TRUE)
# Remove files the contain "FG_Matched_IVs" and "HBA1C_Matched_IVs"
metabolite_files <- metabolite_files[!grepl("FG_Matched_IVs", metabolite_files)]
metabolite_files <- metabolite_files[!grepl("HBA1C_Matched_IVs", metabolite_files)]

# Initialise a metabolite counter
metabolite_counter <- 0
# Initialise a counter for the number of proxies used
proxy_counter <- 0
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
# Initialise a one IV counter
one_IV <- 0
# Initialise a no IVs counter
no_IVs <- 0
# Initialise an ambiguous proxy counter
ambiguous_proxy <- 0

for (metabolite_file in metabolite_files) {
  # Load the metabolite data
  file_path <- metabolite_file
  metabolite_data <- readr::read_tsv(file_path)
  
  # If the value for EffectAllele, NonEffectAllele, t2dm_EffectAllele or t2dm_NonEffectAllele is TRUE, change it to T as a character
  metabolite_data$EffectAllele <- ifelse(metabolite_data$EffectAllele == TRUE, "T", metabolite_data$EffectAllele)
  metabolite_data$NonEffectAllele <- ifelse(metabolite_data$NonEffectAllele == TRUE, "T", metabolite_data$NonEffectAllele)
  metabolite_data$t2dm_EffectAllele <- ifelse(metabolite_data$t2dm_EffectAllele == TRUE, "T", metabolite_data$t2dm_EffectAllele)
  metabolite_data$t2dm_NonEffectAllele <- ifelse(metabolite_data$t2dm_NonEffectAllele == TRUE, "T", metabolite_data$t2dm_NonEffectAllele)
  metabolite_data$proxy_EffectAllele <- ifelse(metabolite_data$proxy_EffectAllele == TRUE, "T", metabolite_data$proxy_EffectAllele)
  metabolite_data$proxy_NonEffectAllele <- ifelse(metabolite_data$proxy_NonEffectAllele == TRUE, "T", metabolite_data$proxy_NonEffectAllele)
  
  # If there are NAs in any of the columns, remove the row
  metabolite_data <- metabolite_data[!is.na(metabolite_data$EAF), ]
  
  # Make new dataframes with the column names for harmonising alleles using TwoSampleMR
  exposure_data <- metabolite_data %>% select(SNP, Beta, SE, EffectAllele, NonEffectAllele, EAF)
  outcome_data <- metabolite_data %>% select(SNP, t2dm_Beta, t2dm_SE, t2dm_EffectAllele, t2dm_NonEffectAllele, t2dm_EAF)
  
  tolerance = 0.08
  
  # Check the value of the proxy column for each SNP
  for (i in 1:nrow(metabolite_data)) {
    # If the proxy is not NA, proceed
    if (!is.na(metabolite_data$proxy[i])) {
      exposure_data$EAF[i] <- metabolite_data$proxy_EAF[i]
      exposure_data$EffectAllele[i] <- metabolite_data$proxy_EffectAllele[i]
      exposure_data$NonEffectAllele[i] <- metabolite_data$proxy_NonEffectAllele[i]
      proxy_counter <- proxy_counter + 1
      
      # Check if the proxy EAF is similar to the EAF or 1 - (EAF) - tolerance is 0.08
      upper_tolerance = metabolite_data$EAF[i] + tolerance
      lower_tolerance = metabolite_data$EAF[i] - tolerance
      # If proxy_EAF is not between the tolerances check if it is 1 - (EAF)
      if (metabolite_data$proxy_EAF[i] < lower_tolerance | metabolite_data$proxy_EAF[i] > upper_tolerance) {
        inverse_EAF = 1 - metabolite_data$EAF[i]
        upper_inverse_tolerance = inverse_EAF + tolerance
        lower_inverse_tolerance = inverse_EAF - tolerance
        # If proxy_EAF is between the tolerances, change the exposure_data$Effect_Allele to the proxy_NonEffectAllele...
        if (metabolite_data$proxy_EAF[i] > lower_inverse_tolerance & metabolite_data$proxy_EAF[i] < upper_inverse_tolerance) {
          exposure_data$EffectAllele[i] <- metabolite_data$proxy_NonEffectAllele[i]
          exposure_data$NonEffectAllele[i] <- metabolite_data$proxy_EffectAllele[i]
          exposure_data$EAF[i] <- metabolite_data$EAF[i]
        } else {
          # Remove the IV from exposure_data, outcome_data and metabolite data
          metabolite_data <- metabolite_data[metabolite_data$SNP != metabolite_data$SNP[i], ]
          exposure_data <- exposure_data[exposure_data$SNP != exposure_data$SNP[i], ]
          outcome_data <- outcome_data[outcome_data$SNP != outcome_data$SNP[i], ]
          print("Removed IV due to ambigous proxies")
          ambiguous_proxy <- ambiguous_proxy + 1
        }
      }
    }
  }
  # Rename the columns
  colnames(exposure_data) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
  colnames(outcome_data) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")
  # Add other required columns
  exposure_data$exposure <- "exposure"
  exposure_data$id.exposure <- "exposure"
  outcome_data$outcome <- "outcome"
  outcome_data$id.outcome <- "outcome"
  
  # Check if any rows remain before harmonsing
  # If a row is all NA, remove the row
  metabolite_data <- metabolite_data[rowSums(is.na(metabolite_data)) != ncol(metabolite_data), ]
  exposure_data <- exposure_data[rowSums(is.na(exposure_data)) != ncol(exposure_data), ]
  outcome_data <- outcome_data[rowSums(is.na(outcome_data)) != ncol(outcome_data), ]
  
  # Make harmonised_IVs an empty dataframe
  harmonised_IVs <- data.frame()
  
  # Harmonise alleles if metabolite_data, exposure_data and outcome_data all have more than 0 rows
  if (nrow(metabolite_data) > 0 & nrow(exposure_data) > 0 & nrow(outcome_data) > 0) {
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
    # If harmonised_IVs has more than 0 rows, proceed
    if (nrow(harmonised_IVs) > 0) {
      for (i in 1:nrow(harmonised_IVs)) {
        rsid <- harmonised_IVs$SNP[i]
        rsid_data <- harmonised_IVs[i, ]
        proxy_value <- metabolite_data$proxy[metabolite_data$SNP == rsid]
        # If the proxy is NA, proceed
        if (is.na(proxy_value)) {
          metabolite_data <- metabolite_data %>% mutate(Beta = ifelse(SNP == rsid, rsid_data$beta.exposure, Beta),
                                                        SE = ifelse(SNP == rsid, rsid_data$se.exposure, SE),
                                                        EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.exposure, EffectAllele),
                                                        NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.exposure, NonEffectAllele),
                                                        EAF = ifelse(SNP == rsid, rsid_data$eaf.exposure, EAF),
                                                        t2dm_Beta = ifelse(SNP == rsid, rsid_data$beta.outcome, t2dm_Beta),
                                                        t2dm_SE = ifelse(SNP == rsid, rsid_data$se.outcome, t2dm_SE),
                                                        t2dm_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.outcome, t2dm_EffectAllele),
                                                        t2dm_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.outcome, t2dm_NonEffectAllele),
                                                        t2dm_EAF = ifelse(SNP == rsid, rsid_data$eaf.outcome, t2dm_EAF))
        }
        # If the proxy is not NA, proceed
        if (!is.na(proxy_value)) {
          metabolite_data <- metabolite_data %>% mutate(Beta = ifelse(SNP == rsid, rsid_data$beta.exposure, Beta),
                                                        SE = ifelse(SNP == rsid, rsid_data$se.exposure, SE),
                                                        proxy_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.exposure, proxy_EffectAllele),
                                                        proxy_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.exposure, proxy_NonEffectAllele),
                                                        proxy_EAF = ifelse(SNP == rsid, rsid_data$eaf.exposure, EAF),
                                                        t2dm_Beta = ifelse(SNP == rsid, rsid_data$beta.outcome, t2dm_Beta),
                                                        t2dm_SE = ifelse(SNP == rsid, rsid_data$se.outcome, t2dm_SE),
                                                        t2dm_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.outcome, t2dm_EffectAllele),
                                                        t2dm_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.outcome, t2dm_NonEffectAllele),
                                                        t2dm_EAF = ifelse(SNP == rsid, rsid_data$eaf.outcome, t2dm_EAF))
        }
      }
    }
  }
  # Check the number of IVs left
  if (nrow(metabolite_data) > 1) {
    multiple_IVs <- multiple_IVs + 1
    # Save the metabolite data
    # Extract the part of the file name before the "_Matched_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Matched_IVs")[[1]][1]
    # Add "T2DM_Harmonised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "T2DM_Harmonised_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Harmonised_T2DM_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 1) {
    one_IV <- one_IV + 1
    # Save the metabolite data
    # Extract the part of the file name before the "_Matched_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Matched_IVs")[[1]][1]
    # Add "T2DM_Harmonised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "T2DM_Harmonised_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Harmonised_T2DM_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 0) {
    no_IVs <- no_IVs + 1
    metabolite_counter <- metabolite_counter + 1
    # Do not save the data
  }
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
  print(paste0("Number of ambiguous proxies: ", ambiguous_proxy))
  print(paste0("Number of Metabolites with multiple IVs: ", multiple_IVs))
  print(paste0("Number of Metabolites with one IV: ", one_IV))
  print(paste0("Number of Metabolites removed: ", no_IVs))
}







### Now for fasting glucose ###
# Load the data
metabolite_files <- list.files("Matched_IVs", full.names = TRUE)
# Remove files the contain "FG_Matched_IVs" and "HBA1C_Matched_IVs"
metabolite_files <- metabolite_files[!grepl("T2DM_Matched_IVs", metabolite_files)]
metabolite_files <- metabolite_files[!grepl("HBA1C_Matched_IVs", metabolite_files)]

# Initialise a metabolite counter
metabolite_counter <- 0
# Initialise a counter for the number of proxies used
proxy_counter <- 0
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
# Initialise a one IV counter
one_IV <- 0
# Initialise a no IVs counter
no_IVs <- 0
# Initialise an ambiguous proxy counter
ambiguous_proxy <- 0

for (metabolite_file in metabolite_files) {
  # Load the metabolite data
  file_path <- metabolite_file
  metabolite_data <- readr::read_tsv(file_path)
  
  # If the value for EffectAllele, NonEffectAllele, fg_EffectAllele or fg_NonEffectAllele is TRUE, change it to T as a character
  metabolite_data$EffectAllele <- ifelse(metabolite_data$EffectAllele == TRUE, "T", metabolite_data$EffectAllele)
  metabolite_data$NonEffectAllele <- ifelse(metabolite_data$NonEffectAllele == TRUE, "T", metabolite_data$NonEffectAllele)
  metabolite_data$fg_EffectAllele <- ifelse(metabolite_data$fg_EffectAllele == TRUE, "T", metabolite_data$fg_EffectAllele)
  metabolite_data$fg_NonEffectAllele <- ifelse(metabolite_data$fg_NonEffectAllele == TRUE, "T", metabolite_data$fg_NonEffectAllele)
  metabolite_data$proxy_EffectAllele <- ifelse(metabolite_data$proxy_EffectAllele == TRUE, "T", metabolite_data$proxy_EffectAllele)
  metabolite_data$proxy_NonEffectAllele <- ifelse(metabolite_data$proxy_NonEffectAllele == TRUE, "T", metabolite_data$proxy_NonEffectAllele)
  
  # If there are NAs in any of the columns, remove the row
  metabolite_data <- metabolite_data[!is.na(metabolite_data$EAF), ]
  
  # Make new dataframes with the column names for harmonising alleles using TwoSampleMR
  exposure_data <- metabolite_data %>% select(SNP, Beta, SE, EffectAllele, NonEffectAllele, EAF)
  outcome_data <- metabolite_data %>% select(SNP, fg_Beta, fg_SE, fg_EffectAllele, fg_NonEffectAllele, fg_EAF)
  
  tolerance = 0.08
  
  # Check the value of the proxy column for each SNP
  for (i in 1:nrow(metabolite_data)) {
    # If the proxy is not NA, proceed
    if (!is.na(metabolite_data$proxy[i])) {
      exposure_data$EAF[i] <- metabolite_data$proxy_EAF[i]
      exposure_data$EffectAllele[i] <- metabolite_data$proxy_EffectAllele[i]
      exposure_data$NonEffectAllele[i] <- metabolite_data$proxy_NonEffectAllele[i]
      proxy_counter <- proxy_counter + 1
      
      # Check if the proxy EAF is similar to the EAF or 1 - (EAF) - tolerance is 0.08
      upper_tolerance = metabolite_data$EAF[i] + tolerance
      lower_tolerance = metabolite_data$EAF[i] - tolerance
      # If proxy_EAF is not between the tolerances check if it is 1 - (EAF)
      if (metabolite_data$proxy_EAF[i] < lower_tolerance | metabolite_data$proxy_EAF[i] > upper_tolerance) {
        inverse_EAF = 1 - metabolite_data$EAF[i]
        upper_inverse_tolerance = inverse_EAF + tolerance
        lower_inverse_tolerance = inverse_EAF - tolerance
        # If proxy_EAF is between the tolerances, change the exposure_data$Effect_Allele to the proxy_NonEffectAllele...
        if (metabolite_data$proxy_EAF[i] > lower_inverse_tolerance & metabolite_data$proxy_EAF[i] < upper_inverse_tolerance) {
          exposure_data$EffectAllele[i] <- metabolite_data$proxy_NonEffectAllele[i]
          exposure_data$NonEffectAllele[i] <- metabolite_data$proxy_EffectAllele[i]
          exposure_data$EAF[i] <- metabolite_data$EAF[i]
        } else {
          # Remove the IV from exposure_data, outcome_data and metabolite data
          metabolite_data <- metabolite_data[metabolite_data$SNP != metabolite_data$SNP[i], ]
          exposure_data <- exposure_data[exposure_data$SNP != exposure_data$SNP[i], ]
          outcome_data <- outcome_data[outcome_data$SNP != outcome_data$SNP[i], ]
          print("Removed IV due to ambigous proxies")
          ambiguous_proxy <- ambiguous_proxy + 1
        }
      }
    }
  }
  # Rename the columns
  colnames(exposure_data) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
  colnames(outcome_data) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")
  # Add other required columns
  exposure_data$exposure <- "exposure"
  exposure_data$id.exposure <- "exposure"
  outcome_data$outcome <- "outcome"
  outcome_data$id.outcome <- "outcome"
  
  # Check if any rows remain before harmonsing
  # If a row is all NA, remove the row
  metabolite_data <- metabolite_data[rowSums(is.na(metabolite_data)) != ncol(metabolite_data), ]
  exposure_data <- exposure_data[rowSums(is.na(exposure_data)) != ncol(exposure_data), ]
  outcome_data <- outcome_data[rowSums(is.na(outcome_data)) != ncol(outcome_data), ]
  
  # Make harmonised_IVs an empty dataframe
  harmonised_IVs <- data.frame()
  
  # Harmonise alleles if metabolite_data, exposure_data and outcome_data all have more than 0 rows
  if (nrow(metabolite_data) > 0 & nrow(exposure_data) > 0 & nrow(outcome_data) > 0) {
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
    # If harmonised_IVs has more than 0 rows, proceed
    if (nrow(harmonised_IVs) > 0) {
      for (i in 1:nrow(harmonised_IVs)) {
        rsid <- harmonised_IVs$SNP[i]
        rsid_data <- harmonised_IVs[i, ]
        proxy_value <- metabolite_data$proxy[metabolite_data$SNP == rsid]
        # If the proxy is NA, proceed
        if (is.na(proxy_value)) {
          metabolite_data <- metabolite_data %>% mutate(Beta = ifelse(SNP == rsid, rsid_data$beta.exposure, Beta),
                                                        SE = ifelse(SNP == rsid, rsid_data$se.exposure, SE),
                                                        EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.exposure, EffectAllele),
                                                        NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.exposure, NonEffectAllele),
                                                        EAF = ifelse(SNP == rsid, rsid_data$eaf.exposure, EAF),
                                                        fg_Beta = ifelse(SNP == rsid, rsid_data$beta.outcome, fg_Beta),
                                                        fg_SE = ifelse(SNP == rsid, rsid_data$se.outcome, fg_SE),
                                                        fg_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.outcome, fg_EffectAllele),
                                                        fg_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.outcome, fg_NonEffectAllele),
                                                        fg_EAF = ifelse(SNP == rsid, rsid_data$eaf.outcome, fg_EAF))
        }
        # If the proxy is not NA, proceed
        if (!is.na(proxy_value)) {
          metabolite_data <- metabolite_data %>% mutate(Beta = ifelse(SNP == rsid, rsid_data$beta.exposure, Beta),
                                                        SE = ifelse(SNP == rsid, rsid_data$se.exposure, SE),
                                                        proxy_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.exposure, proxy_EffectAllele),
                                                        proxy_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.exposure, proxy_NonEffectAllele),
                                                        proxy_EAF = ifelse(SNP == rsid, rsid_data$eaf.exposure, EAF),
                                                        fg_Beta = ifelse(SNP == rsid, rsid_data$beta.outcome, fg_Beta),
                                                        fg_SE = ifelse(SNP == rsid, rsid_data$se.outcome, fg_SE),
                                                        fg_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.outcome, fg_EffectAllele),
                                                        fg_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.outcome, fg_NonEffectAllele),
                                                        fg_EAF = ifelse(SNP == rsid, rsid_data$eaf.outcome, fg_EAF))
        }
      }
    }
  }
  # Check the number of IVs left
  if (nrow(metabolite_data) > 1) {
    multiple_IVs <- multiple_IVs + 1
    # Save the metabolite data
    # Extract the part of the file name before the "_Matched_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Matched_IVs")[[1]][1]
    # Add "FG_Harmonised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "FG_Harmonised_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Harmonised_FG_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 1) {
    one_IV <- one_IV + 1
    # Save the metabolite data
    # Extract the part of the file name before the "_Matched_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Matched_IVs")[[1]][1]
    # Add "FG_Harmonised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "FG_Harmonised_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Harmonised_FG_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 0) {
    no_IVs <- no_IVs + 1
    metabolite_counter <- metabolite_counter + 1
    # Do not save the data
  }
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
  print(paste0("Number of ambiguous proxies: ", ambiguous_proxy))
  print(paste0("Number of Metabolites with multiple IVs: ", multiple_IVs))
  print(paste0("Number of Metabolites with one IV: ", one_IV))
  print(paste0("Number of Metabolites removed: ", no_IVs))
}







### Now for hbA1c ###
# Load the data
metabolite_files <- list.files("Matched_IVs", full.names = TRUE)
# Remove files the contain "FG_Matched_IVs" and "HBA1C_Matched_IVs"
metabolite_files <- metabolite_files[!grepl("FG_Matched_IVs", metabolite_files)]
metabolite_files <- metabolite_files[!grepl("T2DM_Matched_IVs", metabolite_files)]

# Initialise a metabolite counter
metabolite_counter <- 0
# Initialise a counter for the number of proxies used
proxy_counter <- 0
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
# Initialise a one IV counter
one_IV <- 0
# Initialise a no IVs counter
no_IVs <- 0
# Initialise an ambiguous proxy counter
ambiguous_proxy <- 0

for (metabolite_file in metabolite_files) {
  # Load the metabolite data
  file_path <- metabolite_file
  metabolite_data <- readr::read_tsv(file_path)
  
  # If the value for EffectAllele, NonEffectAllele, hba1c_EffectAllele or hba1c_NonEffectAllele is TRUE, change it to T as a character
  metabolite_data$EffectAllele <- ifelse(metabolite_data$EffectAllele == TRUE, "T", metabolite_data$EffectAllele)
  metabolite_data$NonEffectAllele <- ifelse(metabolite_data$NonEffectAllele == TRUE, "T", metabolite_data$NonEffectAllele)
  metabolite_data$hba1c_EffectAllele <- ifelse(metabolite_data$hba1c_EffectAllele == TRUE, "T", metabolite_data$hba1c_EffectAllele)
  metabolite_data$hba1c_NonEffectAllele <- ifelse(metabolite_data$hba1c_NonEffectAllele == TRUE, "T", metabolite_data$hba1c_NonEffectAllele)
  metabolite_data$proxy_EffectAllele <- ifelse(metabolite_data$proxy_EffectAllele == TRUE, "T", metabolite_data$proxy_EffectAllele)
  metabolite_data$proxy_NonEffectAllele <- ifelse(metabolite_data$proxy_NonEffectAllele == TRUE, "T", metabolite_data$proxy_NonEffectAllele)
  
  # If there are NAs in any of the columns, remove the row
  metabolite_data <- metabolite_data[!is.na(metabolite_data$EAF), ]
  
  # Make new dataframes with the column names for harmonising alleles using TwoSampleMR
  exposure_data <- metabolite_data %>% select(SNP, Beta, SE, EffectAllele, NonEffectAllele, EAF)
  outcome_data <- metabolite_data %>% select(SNP, hba1c_Beta, hba1c_SE, hba1c_EffectAllele, hba1c_NonEffectAllele, hba1c_EAF)
  
  tolerance = 0.08
  
  # Check the value of the proxy column for each SNP
  for (i in 1:nrow(metabolite_data)) {
    # If the proxy is not NA, proceed
    if (!is.na(metabolite_data$proxy[i])) {
      exposure_data$EAF[i] <- metabolite_data$proxy_EAF[i]
      exposure_data$EffectAllele[i] <- metabolite_data$proxy_EffectAllele[i]
      exposure_data$NonEffectAllele[i] <- metabolite_data$proxy_NonEffectAllele[i]
      proxy_counter <- proxy_counter + 1
      
      # Check if the proxy EAF is similar to the EAF or 1 - (EAF) - tolerance is 0.08
      upper_tolerance = metabolite_data$EAF[i] + tolerance
      lower_tolerance = metabolite_data$EAF[i] - tolerance
      # If proxy_EAF is not between the tolerances check if it is 1 - (EAF)
      if (metabolite_data$proxy_EAF[i] < lower_tolerance | metabolite_data$proxy_EAF[i] > upper_tolerance) {
        inverse_EAF = 1 - metabolite_data$EAF[i]
        upper_inverse_tolerance = inverse_EAF + tolerance
        lower_inverse_tolerance = inverse_EAF - tolerance
        # If proxy_EAF is between the tolerances, change the exposure_data$Effect_Allele to the proxy_NonEffectAllele...
        if (metabolite_data$proxy_EAF[i] > lower_inverse_tolerance & metabolite_data$proxy_EAF[i] < upper_inverse_tolerance) {
          exposure_data$EffectAllele[i] <- metabolite_data$proxy_NonEffectAllele[i]
          exposure_data$NonEffectAllele[i] <- metabolite_data$proxy_EffectAllele[i]
          exposure_data$EAF[i] <- metabolite_data$EAF[i]
        } else {
          # Remove the IV from exposure_data, outcome_data and metabolite data
          metabolite_data <- metabolite_data[metabolite_data$SNP != metabolite_data$SNP[i], ]
          exposure_data <- exposure_data[exposure_data$SNP != exposure_data$SNP[i], ]
          outcome_data <- outcome_data[outcome_data$SNP != outcome_data$SNP[i], ]
          print("Removed IV due to ambigous proxies")
          ambiguous_proxy <- ambiguous_proxy + 1
        }
      }
    }
  }
  # Rename the columns
  colnames(exposure_data) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure")
  colnames(outcome_data) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome")
  # Add other required columns
  exposure_data$exposure <- "exposure"
  exposure_data$id.exposure <- "exposure"
  outcome_data$outcome <- "outcome"
  outcome_data$id.outcome <- "outcome"
  
  # Check if any rows remain before harmonsing
  # If a row is all NA, remove the row
  metabolite_data <- metabolite_data[rowSums(is.na(metabolite_data)) != ncol(metabolite_data), ]
  exposure_data <- exposure_data[rowSums(is.na(exposure_data)) != ncol(exposure_data), ]
  outcome_data <- outcome_data[rowSums(is.na(outcome_data)) != ncol(outcome_data), ]
  
  # Make harmonised_IVs an empty dataframe
  harmonised_IVs <- data.frame()
  
  # Harmonise alleles if metabolite_data, exposure_data and outcome_data all have more than 0 rows
  if (nrow(metabolite_data) > 0 & nrow(exposure_data) > 0 & nrow(outcome_data) > 0) {
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
    # If harmonised_IVs has more than 0 rows, proceed
    if (nrow(harmonised_IVs) > 0) {
      for (i in 1:nrow(harmonised_IVs)) {
        rsid <- harmonised_IVs$SNP[i]
        rsid_data <- harmonised_IVs[i, ]
        proxy_value <- metabolite_data$proxy[metabolite_data$SNP == rsid]
        # If the proxy is NA, proceed
        if (is.na(proxy_value)) {
          metabolite_data <- metabolite_data %>% mutate(Beta = ifelse(SNP == rsid, rsid_data$beta.exposure, Beta),
                                                        SE = ifelse(SNP == rsid, rsid_data$se.exposure, SE),
                                                        EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.exposure, EffectAllele),
                                                        NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.exposure, NonEffectAllele),
                                                        EAF = ifelse(SNP == rsid, rsid_data$eaf.exposure, EAF),
                                                        hba1c_Beta = ifelse(SNP == rsid, rsid_data$beta.outcome, hba1c_Beta),
                                                        hba1c_SE = ifelse(SNP == rsid, rsid_data$se.outcome, hba1c_SE),
                                                        hba1c_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.outcome, hba1c_EffectAllele),
                                                        hba1c_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.outcome, hba1c_NonEffectAllele),
                                                        hba1c_EAF = ifelse(SNP == rsid, rsid_data$eaf.outcome, hba1c_EAF))
        }
        # If the proxy is not NA, proceed
        if (!is.na(proxy_value)) {
          metabolite_data <- metabolite_data %>% mutate(Beta = ifelse(SNP == rsid, rsid_data$beta.exposure, Beta),
                                                        SE = ifelse(SNP == rsid, rsid_data$se.exposure, SE),
                                                        proxy_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.exposure, proxy_EffectAllele),
                                                        proxy_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.exposure, proxy_NonEffectAllele),
                                                        proxy_EAF = ifelse(SNP == rsid, rsid_data$eaf.exposure, EAF),
                                                        hba1c_Beta = ifelse(SNP == rsid, rsid_data$beta.outcome, hba1c_Beta),
                                                        hba1c_SE = ifelse(SNP == rsid, rsid_data$se.outcome, hba1c_SE),
                                                        hba1c_EffectAllele = ifelse(SNP == rsid, rsid_data$effect_allele.outcome, hba1c_EffectAllele),
                                                        hba1c_NonEffectAllele = ifelse(SNP == rsid, rsid_data$other_allele.outcome, hba1c_NonEffectAllele),
                                                        hba1c_EAF = ifelse(SNP == rsid, rsid_data$eaf.outcome, hba1c_EAF))
        }
      }
    }
  }
  # Check the number of IVs left
  if (nrow(metabolite_data) > 1) {
    multiple_IVs <- multiple_IVs + 1
    # Save the metabolite data
    # Extract the part of the file name before the "_Matched_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Matched_IVs")[[1]][1]
    # Add "HBA1C_Harmonised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "HBA1C_Harmonised_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Harmonised_HBA1C_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 1) {
    one_IV <- one_IV + 1
    # Save the metabolite data
    # Extract the part of the file name before the "_Matched_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Matched_IVs")[[1]][1]
    # Add "HBA1C_Harmonised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "HBA1C_Harmonised_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Harmonised_HBA1C_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 0) {
    no_IVs <- no_IVs + 1
    metabolite_counter <- metabolite_counter + 1
    # Do not save the data
  }
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
  print(paste0("Number of ambiguous proxies: ", ambiguous_proxy))
  print(paste0("Number of Metabolites with multiple IVs: ", multiple_IVs))
  print(paste0("Number of Metabolites with one IV: ", one_IV))
  print(paste0("Number of Metabolites removed: ", no_IVs))
}



























