library(tidyverse)
library(stringr)
library(MendelianRandomization)

# Load metabolite_t2dm_coloc_pass.tsv
metabolite_t2dm_coloc_pass <- read_tsv("Colocalisation/metabolite_t2dm_coloc_pass.tsv")
# Select the metabolites that do not pass any of the conditions
metabolite_t2dm_coloc_fail <- metabolite_t2dm_coloc_pass %>% filter(!unconditional_pass & !conditional_pass)

# Add the following columns to metabolite_t2dm_coloc_fail: num_IVs_v2, IVW_Estimate, IVW_SE, IVW_Pval, IVW_Fstat
metabolite_t2dm_coloc_fail <- metabolite_t2dm_coloc_fail %>% 
  mutate(num_IVs_v2 = NA_real_, IVW_Estimate = NA_real_, IVW_SE = NA_real_, IVW_Pval = NA_real_, IVW_Fstat = NA_real_)

# Load pwcoco_t2dm_coloc_out.coloc
pwcoco_t2dm_coloc_out <- read_delim("Colocalisation/pwcoco_out_t2dm.coloc")

# For each metabolite in metabolite_t2dm_coloc_fail
for (metabolite in metabolite_t2dm_coloc_fail$Metabolite) {
  # Select the rows in pwcoco_t2dm_coloc_out that contain the metabolite
  pwcoco_t2dm_coloc_out_metabolite <- pwcoco_t2dm_coloc_out %>% 
    filter(str_detect(Dataset1, fixed(metabolite)))
  
  # Remove the _metabolite_snps_pwcoco.tsv from the unique_Dataset1_values
  pwcoco_t2dm_coloc_out_metabolite$Dataset1 <- str_replace_all(pwcoco_t2dm_coloc_out_metabolite$Dataset1, "_metabolite_snps_pwcoco.tsv", "")
  # Remove all characters before "rs" in the unique_Dataset1_values
  pwcoco_t2dm_coloc_out_metabolite$Dataset1 <- str_replace_all(pwcoco_t2dm_coloc_out_metabolite$Dataset1, ".*rs", "rs")
  
  # Extract the unique rsids from the remaining Dataset1 values
  rsids <- pwcoco_t2dm_coloc_out_metabolite$Dataset1 %>% unique()
  
  # Check if any of the H4 values are greater than or equal to 0.8
  for (rsid in rsids) {
    if (any(pwcoco_t2dm_coloc_out_metabolite$H4[pwcoco_t2dm_coloc_out_metabolite$Dataset1 == rsid] >= 0.8)) {
      # Get the rsid of that row
      rsid_high_h4 <- rsid
      # Remove all rows with that rsid_high_h4
      pwcoco_t2dm_coloc_out_metabolite <- pwcoco_t2dm_coloc_out_metabolite %>% filter(Dataset1 != rsid_high_h4)
    }
  }
  
  # Load the harmonsied IVs to re-run Mendelian Randomisation
  harm_IVs <- read_tsv(paste0("Harmonised_T2DM_IVs/", metabolite, "T2DMT2DM_Harmonised_IVs.tsv"))
  
  # Remove the rows from harm_IVs if they have a harm_IVs$SNP that is in pwcoco_t2dm_coloc_out_metabolite$Dataset1
  harm_IVs <- harm_IVs %>% filter(!SNP %in% pwcoco_t2dm_coloc_out_metabolite$Dataset1)
  
  num_IVs_v2 <- nrow(harm_IVs)
  metabolite_t2dm_coloc_fail$num_IVs_v2[metabolite_t2dm_coloc_fail$Metabolite == metabolite] <- num_IVs_v2
  
  # If there are no IVs left, skip the metabolite
  if (num_IVs_v2 == 0) {
    next
  }
  
  MRObject = mr_input(bx = harm_IVs$Beta, bxse = harm_IVs$SE, by = harm_IVs$t2dm_Beta, byse = harm_IVs$t2dm_SE)
  
  # If there are 3 or more IVs use a random effects model
  if (num_IVs_v2 >= 3) {
    MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
    metabolite_t2dm_coloc_fail$IVW_Estimate[metabolite_t2dm_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$Estimate
    metabolite_t2dm_coloc_fail$IVW_SE[metabolite_t2dm_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$StdError
    metabolite_t2dm_coloc_fail$IVW_Pval[metabolite_t2dm_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$Pvalue
    metabolite_t2dm_coloc_fail$IVW_Fstat[metabolite_t2dm_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$Fstat
  } else if (num_IVs_v2 < 3) {
    # Use a fixed effects model
    MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")
    metabolite_t2dm_coloc_fail$IVW_Estimate[metabolite_t2dm_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$Estimate
    metabolite_t2dm_coloc_fail$IVW_SE[metabolite_t2dm_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$StdError
    metabolite_t2dm_coloc_fail$IVW_Pval[metabolite_t2dm_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$Pvalue
    metabolite_t2dm_coloc_fail$IVW_Fstat[metabolite_t2dm_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$Fstat
  }
}

# Store metabolites with 0 IVs in a separate dataframe
metabolite_t2dm_coloc_fail_remove <- metabolite_t2dm_coloc_fail %>% filter(num_IVs_v2 == 0)
# Add metabolites with IVW_Pval > 3.006614e-5 to metabolite_t2dm_coloc_fail_remove
metabolite_t2dm_coloc_fail_remove <- rbind(metabolite_t2dm_coloc_fail_remove, metabolite_t2dm_coloc_fail %>% filter(IVW_Pval > 3.006614e-5))

# Load Filtered_T2DM_Results.tsv
filtered_t2dm_IVs <- read_tsv("Filtered_T2DM_Results.tsv")

# Add metabolites with opposite signs for the metabolite_t2dm_coloc_fail$IVW_Estimate and filtered_t2dm_IVs$IVW_Estimate
for (metabolite in metabolite_t2dm_coloc_fail$Metabolite) {
  if (metabolite_t2dm_coloc_fail$num_IVs_v2[metabolite_t2dm_coloc_fail$Metabolite == metabolite] > 0) {
    rerun_MR_effect <- metabolite_t2dm_coloc_fail$IVW_Estimate[metabolite_t2dm_coloc_fail$Metabolite == metabolite]
    original_MR_effect <- filtered_t2dm_IVs$Fixed_IVW_Estimate[filtered_t2dm_IVs$Metabolite == metabolite]
    if (sign(rerun_MR_effect) != sign(original_MR_effect)) {
      metabolite_t2dm_coloc_fail_remove <- rbind(metabolite_t2dm_coloc_fail_remove, metabolite_t2dm_coloc_fail %>% filter(Metabolite == metabolite))
    }
  }
}

# Add a column to metabolite_t2dm_coloc_pass to indicate if the metabolite passed or failed Coloc_MR
metabolite_t2dm_coloc_pass$Coloc_MR_Pass <- TRUE
# Change it to FALSE if the metabolite is in metabolite_t2dm_coloc_fail_remove
metabolite_t2dm_coloc_pass$Coloc_MR_Pass[metabolite_t2dm_coloc_pass$Metabolite %in% metabolite_t2dm_coloc_fail_remove$Metabolite] <- FALSE

# Save the results
write_tsv(metabolite_t2dm_coloc_pass, "Colocalisation/Final_Metabolite_T2DM_Coloc_Results.tsv")
write_tsv(metabolite_t2dm_coloc_fail, "Colocalisation/Metabolite_T2DM_Coloc_Fail_MR_Results.tsv")




# Perform for FG
# Load metabolite_fg_coloc_pass.tsv
metabolite_fg_coloc_pass <- read_tsv("Colocalisation/metabolite_fg_coloc_pass.tsv")
# Select the metabolites that do not pass any of the conditions
metabolite_fg_coloc_fail <- metabolite_fg_coloc_pass %>% filter(!unconditional_pass & !conditional_pass)

# Add the following columns to metabolite_fg_coloc_fail: num_IVs_v2, IVW_Estimate, IVW_SE, IVW_Pval, IVW_Fstat
metabolite_fg_coloc_fail <- metabolite_fg_coloc_fail %>% 
  mutate(num_IVs_v2 = NA_real_, IVW_Estimate = NA_real_, IVW_SE = NA_real_, IVW_Pval = NA_real_, IVW_Fstat = NA_real_)

# Load pwcoco_fg_coloc_out.coloc
pwcoco_fg_coloc_out <- read_delim("Colocalisation/pwcoco_out_fg.coloc")

# For each metabolite in metabolite_fg_coloc_fail
for (metabolite in metabolite_fg_coloc_fail$Metabolite) {
  # Select the rows in pwcoco_fg_coloc_out that contain the metabolite
  pwcoco_fg_coloc_out_metabolite <- pwcoco_fg_coloc_out %>% 
    filter(str_detect(Dataset1, fixed(metabolite)))
  
  # Remove the _metabolite_snps_pwcoco.tsv from the unique_Dataset1_values
  pwcoco_fg_coloc_out_metabolite$Dataset1 <- str_replace_all(pwcoco_fg_coloc_out_metabolite$Dataset1, "_metabolite_snps_pwcoco.tsv", "")
  # Remove all characters before "rs" in the unique_Dataset1_values
  pwcoco_fg_coloc_out_metabolite$Dataset1 <- str_replace_all(pwcoco_fg_coloc_out_metabolite$Dataset1, ".*rs", "rs")
  
  # Extract the unique rsids from the remaining Dataset1 values
  rsids <- pwcoco_fg_coloc_out_metabolite$Dataset1 %>% unique()
  
  # Check if any of the H4 values are greater than or equal to 0.8
  for (rsid in rsids) {
    if (any(pwcoco_fg_coloc_out_metabolite$H4[pwcoco_fg_coloc_out_metabolite$Dataset1 == rsid] >= 0.8)) {
      # Get the rsid of that row
      rsid_high_h4 <- rsid
      # Remove all rows with that rsid_high_h4
      pwcoco_fg_coloc_out_metabolite <- pwcoco_fg_coloc_out_metabolite %>% filter(Dataset1 != rsid_high_h4)
    }
  }
  
  # Load the harmonsied IVs to re-run Mendelian Randomisation
  harm_IVs <- read_tsv(paste0("Harmonised_FG_IVs/", metabolite, "FGFG_Harmonised_IVs.tsv"))
  
  # Remove the rows from harm_IVs if they have a harm_IVs$SNP that is in pwcoco_fg_coloc_out_metabolite$Dataset1
  harm_IVs <- harm_IVs %>% filter(!SNP %in% pwcoco_fg_coloc_out_metabolite$Dataset1)
  
  num_IVs_v2 <- nrow(harm_IVs)
  metabolite_fg_coloc_fail$num_IVs_v2[metabolite_fg_coloc_fail$Metabolite == metabolite] <- num_IVs_v2
  
  # If there are no IVs left, skip the metabolite
  if (num_IVs_v2 == 0) {
    next
  }
  
  MRObject = mr_input(bx = harm_IVs$Beta, bxse = harm_IVs$SE, by = harm_IVs$fg_Beta, byse = harm_IVs$fg_SE)
  
  # If there are 3 or more IVs use a random effects model
  if (num_IVs_v2 >= 3) {
    MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
    metabolite_fg_coloc_fail$IVW_Estimate[metabolite_fg_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$Estimate
    metabolite_fg_coloc_fail$IVW_SE[metabolite_fg_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$StdError
    metabolite_fg_coloc_fail$IVW_Pval[metabolite_fg_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$Pvalue
    metabolite_fg_coloc_fail$IVW_Fstat[metabolite_fg_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$Fstat
  } else if (num_IVs_v2 < 3) {
    # Use a fixed effects model
    MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")
    metabolite_fg_coloc_fail$IVW_Estimate[metabolite_fg_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$Estimate
    metabolite_fg_coloc_fail$IVW_SE[metabolite_fg_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$StdError
    metabolite_fg_coloc_fail$IVW_Pval[metabolite_fg_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$Pvalue
    metabolite_fg_coloc_fail$IVW_Fstat[metabolite_fg_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$Fstat
  }
}

# Store metabolites with 0 IVs in a separate dataframe
metabolite_fg_coloc_fail_remove <- metabolite_fg_coloc_fail %>% filter(num_IVs_v2 == 0)
# Add metabolites with IVW_Pval > 3.006614e-5 to metabolite_fg_coloc_fail_remove
metabolite_fg_coloc_fail_remove <- rbind(metabolite_fg_coloc_fail_remove, metabolite_fg_coloc_fail %>% filter(IVW_Pval > 3.006614e-5))

# Load Filtered_FG_Results.tsv
filtered_fg_IVs <- read_tsv("Filtered_FG_Results.tsv")

# Add metabolites with opposite signs for the metabolite_fg_coloc_fail$IVW_Estimate and filtered_fg_IVs$IVW_Estimate
for (metabolite in metabolite_fg_coloc_fail$Metabolite) {
  if (metabolite_fg_coloc_fail$num_IVs_v2[metabolite_fg_coloc_fail$Metabolite == metabolite] > 0) {
    rerun_MR_effect <- metabolite_fg_coloc_fail$IVW_Estimate[metabolite_fg_coloc_fail$Metabolite == metabolite]
    original_MR_effect <- filtered_fg_IVs$Fixed_IVW_Estimate[filtered_fg_IVs$Metabolite == metabolite]
    if (sign(rerun_MR_effect) != sign(original_MR_effect)) {
      metabolite_fg_coloc_fail_remove <- rbind(metabolite_fg_coloc_fail_remove, metabolite_fg_coloc_fail %>% filter(Metabolite == metabolite))
    }
  }
}

# Add a column to metabolite_fg_coloc_pass to indicate if the metabolite passed or failed Coloc_MR
metabolite_fg_coloc_pass$Coloc_MR_Pass <- TRUE
# Change it to FALSE if the metabolite is in metabolite_fg_coloc_fail_remove
metabolite_fg_coloc_pass$Coloc_MR_Pass[metabolite_fg_coloc_pass$Metabolite %in% metabolite_fg_coloc_fail_remove$Metabolite] <- FALSE

# Save the results
write_tsv(metabolite_fg_coloc_pass, "Colocalisation/Final_Metabolite_FG_Coloc_Results.tsv")
write_tsv(metabolite_fg_coloc_fail, "Colocalisation/Metabolite_FG_Coloc_Fail_MR_Results.tsv")



# Perform for HBA1C
# Load metabolite_hba1c_coloc_pass.tsv
metabolite_hba1c_coloc_pass <- read_tsv("Colocalisation/metabolite_hba1c_coloc_pass.tsv")
# Select the metabolites that do not pass any of the conditions
metabolite_hba1c_coloc_fail <- metabolite_hba1c_coloc_pass %>% filter(!unconditional_pass & !conditional_pass)

# Add the following columns to metabolite_hba1c_coloc_fail: num_IVs_v2, IVW_Estimate, IVW_SE, IVW_Pval, IVW_Fstat
metabolite_hba1c_coloc_fail <- metabolite_hba1c_coloc_fail %>% 
  mutate(num_IVs_v2 = NA_real_, IVW_Estimate = NA_real_, IVW_SE = NA_real_, IVW_Pval = NA_real_, IVW_Fstat = NA_real_)

# Load pwcoco_hba1c_coloc_out.coloc
pwcoco_hba1c_coloc_out <- read_delim("Colocalisation/pwcoco_out_hba1c.coloc")

# For each metabolite in metabolite_hba1c_coloc_fail
for (metabolite in metabolite_hba1c_coloc_fail$Metabolite) {
  # Select the rows in pwcoco_hba1c_coloc_out that contain the metabolite
  pwcoco_hba1c_coloc_out_metabolite <- pwcoco_hba1c_coloc_out %>% 
    filter(str_detect(Dataset1, fixed(metabolite)))
  
  # Remove the _metabolite_snps_pwcoco.tsv from the unique_Dataset1_values
  pwcoco_hba1c_coloc_out_metabolite$Dataset1 <- str_replace_all(pwcoco_hba1c_coloc_out_metabolite$Dataset1, "_metabolite_snps_pwcoco.tsv", "")
  # Remove all characters before "rs" in the unique_Dataset1_values
  pwcoco_hba1c_coloc_out_metabolite$Dataset1 <- str_replace_all(pwcoco_hba1c_coloc_out_metabolite$Dataset1, ".*rs", "rs")
  
  # Extract the unique rsids from the remaining Dataset1 values
  rsids <- pwcoco_hba1c_coloc_out_metabolite$Dataset1 %>% unique()
  
  # Check if any of the H4 values are greater than or equal to 0.8
  for (rsid in rsids) {
    if (any(pwcoco_hba1c_coloc_out_metabolite$H4[pwcoco_hba1c_coloc_out_metabolite$Dataset1 == rsid] >= 0.8)) {
      # Get the rsid of that row
      rsid_high_h4 <- rsid
      # Remove all rows with that rsid_high_h4
      pwcoco_hba1c_coloc_out_metabolite <- pwcoco_hba1c_coloc_out_metabolite %>% filter(Dataset1 != rsid_high_h4)
    }
  }
  
  # Load the harmonsied IVs to re-run Mendelian Randomisation
  harm_IVs <- read_tsv(paste0("Harmonised_HBA1C_IVs/", metabolite, "HBA1CHBA1C_Harmonised_IVs.tsv"))
  
  # Remove the rows from harm_IVs if they have a harm_IVs$SNP that is in pwcoco_hba1c_coloc_out_metabolite$Dataset1
  harm_IVs <- harm_IVs %>% filter(!SNP %in% pwcoco_hba1c_coloc_out_metabolite$Dataset1)
  
  num_IVs_v2 <- nrow(harm_IVs)
  metabolite_hba1c_coloc_fail$num_IVs_v2[metabolite_hba1c_coloc_fail$Metabolite == metabolite] <- num_IVs_v2
  
  # If there are no IVs left, skip the metabolite
  if (num_IVs_v2 == 0) {
    next
  }
  
  MRObject = mr_input(bx = harm_IVs$Beta, bxse = harm_IVs$SE, by = harm_IVs$hba1c_Beta, byse = harm_IVs$hba1c_SE)
  
  # If there are 3 or more IVs use a random effects model
  if (num_IVs_v2 >= 3) {
    MR_random_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "random")
    metabolite_hba1c_coloc_fail$IVW_Estimate[metabolite_hba1c_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$Estimate
    metabolite_hba1c_coloc_fail$IVW_SE[metabolite_hba1c_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$StdError
    metabolite_hba1c_coloc_fail$IVW_Pval[metabolite_hba1c_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$Pvalue
    metabolite_hba1c_coloc_fail$IVW_Fstat[metabolite_hba1c_coloc_fail$Metabolite == metabolite] <- MR_random_IVW_out$Fstat
  } else if (num_IVs_v2 < 3) {
    # Use a fixed effects model
    MR_fixed_IVW_out <- MendelianRandomization::mr_ivw(MRObject, model = "fixed")
    metabolite_hba1c_coloc_fail$IVW_Estimate[metabolite_hba1c_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$Estimate
    metabolite_hba1c_coloc_fail$IVW_SE[metabolite_hba1c_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$StdError
    metabolite_hba1c_coloc_fail$IVW_Pval[metabolite_hba1c_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$Pvalue
    metabolite_hba1c_coloc_fail$IVW_Fstat[metabolite_hba1c_coloc_fail$Metabolite == metabolite] <- MR_fixed_IVW_out$Fstat
  }
}

# Store metabolites with 0 IVs in a separate dataframe
metabolite_hba1c_coloc_fail_remove <- metabolite_hba1c_coloc_fail %>% filter(num_IVs_v2 == 0)
# Add metabolites with IVW_Pval > 3.006614e-5 to metabolite_hba1c_coloc_fail_remove
metabolite_hba1c_coloc_fail_remove <- rbind(metabolite_hba1c_coloc_fail_remove, metabolite_hba1c_coloc_fail %>% filter(IVW_Pval > 3.006614e-5))

# Load Filtered_HBA1C_Results.tsv
filtered_hba1c_IVs <- read_tsv("Filtered_HBA1C_Results.tsv")

# Add metabolites with opposite signs for the metabolite_hba1c_coloc_fail$IVW_Estimate and filtered_hba1c_IVs$IVW_Estimate
for (metabolite in metabolite_hba1c_coloc_fail$Metabolite) {
  if (metabolite_hba1c_coloc_fail$num_IVs_v2[metabolite_hba1c_coloc_fail$Metabolite == metabolite] > 0) {
    rerun_MR_effect <- metabolite_hba1c_coloc_fail$IVW_Estimate[metabolite_hba1c_coloc_fail$Metabolite == metabolite]
    original_MR_effect <- filtered_hba1c_IVs$Fixed_IVW_Estimate[filtered_hba1c_IVs$Metabolite == metabolite]
    if (sign(rerun_MR_effect) != sign(original_MR_effect)) {
      metabolite_hba1c_coloc_fail_remove <- rbind(metabolite_hba1c_coloc_fail_remove, metabolite_hba1c_coloc_fail %>% filter(Metabolite == metabolite))
    }
  }
}

# Add a column to metabolite_hba1c_coloc_pass to indicate if the metabolite passed or failed Coloc_MR
metabolite_hba1c_coloc_pass$Coloc_MR_Pass <- TRUE
# Change it to FALSE if the metabolite is in metabolite_hba1c_coloc_fail_remove
metabolite_hba1c_coloc_pass$Coloc_MR_Pass[metabolite_hba1c_coloc_pass$Metabolite %in% metabolite_hba1c_coloc_fail_remove$Metabolite] <- FALSE

# Save the results
write_tsv(metabolite_hba1c_coloc_pass, "Colocalisation/Final_Metabolite_HBA1C_Coloc_Results.tsv")
write_tsv(metabolite_hba1c_coloc_fail, "Colocalisation/Metabolite_HBA1C_Coloc_Fail_MR_Results.tsv")
