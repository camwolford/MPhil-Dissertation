### Sensitivity Analysis Filtering and Steiger Filtering ###
# This script is used to filter out metabolites that do not pass the sensitivity analysis filtering criteria, before performing the filtering of Steiger filtering.
library(tidyverse)

# Load in the significant_T2DM_results.tsv file
significant_t2dm_results <- readr::read_tsv("significant_T2DM_results_diss.tsv")



# F Stat
# Keep rows with Fixed_IVW_FStat > 10
significant_t2dm_results <- significant_t2dm_results[which(significant_t2dm_results$Fixed_IVW_FStat > 10),]

# Q Stat
# Identify which metabolites have a Fixed_IVW_HetStat_P < 0.05
remove_t2dm_het <- significant_t2dm_results[which(significant_t2dm_results$Fixed_IVW_HetStat_P <= 0.05),]

# Pleiotropy P Value
# Select the metabolites from the significant T2DM results that have an Egger_Intercept_Pval < 0.05
remove_t2dm_ple <- significant_t2dm_results %>% 
  # Select the remaining metabolites from the metabolites_to_remove that have an Egger_Intercept_Pval < 0.05
  dplyr::filter(Egger_Intercept_Pval <= 0.05) %>% 
  dplyr::select(Metabolite)

# Remove the metabolites from the significant T2DM results data frame
significant_t2dm_results <- significant_t2dm_results[which(!significant_t2dm_results$Metabolite %in% remove_t2dm_het$Metabolite),]
# Remove the metabolites from the significant T2DM results data frame
significant_t2dm_results <- significant_t2dm_results[which(!significant_t2dm_results$Metabolite %in% remove_t2dm_ple$Metabolite),]
# Removed three metabolites 62 to 62 to 35 to 35


# Load in the significant_FG_results.tsv file
significant_fg_results <- readr::read_tsv("significant_FG_results_diss.tsv")

# F Stat
# Keep rows with Fixed_IVW_FStat > 10
significant_fg_results <- significant_fg_results[which(significant_fg_results$Fixed_IVW_FStat > 10),]

# Q Stat
# Identify which metabolites hae a Fixed_IVW_HetStat_P < 0.05
remove_fg_het <- significant_fg_results[which(significant_fg_results$Fixed_IVW_HetStat_P <= 0.05),]

# Pleiotropy P Value
# Select the metabolites from the significant FG results that have an Egger_Intercept_Pval < 0.05
remove_fg_ple <- significant_fg_results %>% 
  # Select the remaining metabolites from the metabolites_to_remove that have an Egger_Intercept_Pval < 0.05
  dplyr::filter(Egger_Intercept_Pval <= 0.05) %>% 
  dplyr::select(Metabolite)

# Remove the metabolites from the significant FG results data frame
significant_fg_results <- significant_fg_results[which(!significant_fg_results$Metabolite %in% remove_fg_het$Metabolite),]
# Remove the metabolites from the significant FG results data frame
significant_fg_results <- significant_fg_results[which(!significant_fg_results$Metabolite %in% remove_fg_ple$Metabolite),]
# Removed three metabolites 32 to 32 to 14 to 14


# Load in the significant_HBA1C_results.tsv file
significant_hba1c_results <- readr::read_tsv("significant_HBA1C_results_diss.tsv")

# F Stat
# Keep rows with Fixed_IVW_FStat > 10
significant_hba1c_results <- significant_hba1c_results[which(significant_hba1c_results$Fixed_IVW_FStat > 10),]

# Q Stat
# Identify which metabolites hae a Fixed_IVW_HetStat_P < 0.05
remove_hba1c_het <- significant_hba1c_results[which(significant_hba1c_results$Fixed_IVW_HetStat_P <= 0.05),]

# Pleiotropy P Value
# Select the metabolites from the significant HbA1C results that have an Egger_Intercept_Pval < 0.05
remove_hba1c_ple <- significant_hba1c_results %>% 
  # Select the remaining metabolites from the metabolites_to_remove that have an Egger_Intercept_Pval < 0.05
  dplyr::filter(Egger_Intercept_Pval <= 0.05) %>% 
  dplyr::select(Metabolite)

# Remove the metabolites from the significant HbA1C results data frame
significant_hba1c_results <- significant_hba1c_results[which(!significant_hba1c_results$Metabolite %in% remove_hba1c_het$Metabolite),]
# Remove the metabolites from the significant HbA1C results data frame
significant_hba1c_results <- significant_hba1c_results[which(!significant_hba1c_results$Metabolite %in% remove_hba1c_ple$Metabolite),]
# Removed three metabolites 31 to 31 to 16 to 16



### Steiger Filtering ###
# Load the Metabolite_T2DM_Sensitivity_Analysis.tsv file
metabolite_t2dm_sensitivity_analysis <- readr::read_tsv("Metabolite_T2DM_Sensitivity_Analysis_diss.tsv")

# Keep only the metabolites in metabolite_t2dm_sensitivity_analysis that are still in significant_t2dm_results
metabolite_t2dm_sensitivity_analysis <- metabolite_t2dm_sensitivity_analysis[which(metabolite_t2dm_sensitivity_analysis$Metabolite %in% significant_t2dm_results$Metabolite),]

remove_t2dm_steiger <- metabolite_t2dm_sensitivity_analysis %>% 
  dplyr::filter(Direction_Flag == FALSE) %>%
  dplyr::select(Metabolite) 
# Remove the metabolites from the significant T2DM results data frame
significant_t2dm_results <- significant_t2dm_results[which(!significant_t2dm_results$Metabolite %in% remove_t2dm_steiger$Metabolite),]
# 0 Removed 

# Load the Metabolite_FG_Sensitivity_Analysis.tsv file
metabolite_fg_sensitivity_analysis <- readr::read_tsv("Metabolite_FG_Sensitivity_Analysis_diss.tsv")

# Keep only the metabolites in metabolite_fg_sensitivity_analysis that are still in significant_fg_results
metabolite_fg_sensitivity_analysis <- metabolite_fg_sensitivity_analysis[which(metabolite_fg_sensitivity_analysis$Metabolite %in% significant_fg_results$Metabolite),]

remove_fg_steiger <- metabolite_fg_sensitivity_analysis %>% 
  dplyr::filter(Direction_Flag == FALSE) %>%
  dplyr::select(Metabolite)
# Remove the metabolites from the significant FG results data frame
significant_fg_results <- significant_fg_results[which(!significant_fg_results$Metabolite %in% remove_fg_steiger$Metabolite),]
# 0 Removed 14 to 14

# Load the Metabolite_HBA1C_Sensitivity_Analysis.tsv file
metabolite_hba1c_sensitivity_analysis <- readr::read_tsv("Metabolite_HBA1C_Sensitivity_Analysis_diss.tsv")

# Keep only the metabolites in metabolite_hba1c_sensitivity_analysis that are still in significant_hba1c_results
metabolite_hba1c_sensitivity_analysis <- metabolite_hba1c_sensitivity_analysis[which(metabolite_hba1c_sensitivity_analysis$Metabolite %in% significant_hba1c_results$Metabolite),]

remove_hba1c_steiger <- metabolite_hba1c_sensitivity_analysis %>% 
  dplyr::filter(Direction_Flag == FALSE) %>%
  dplyr::select(Metabolite)
# Remove the metabolites from the significant HbA1C results data frame
significant_hba1c_results <- significant_hba1c_results[which(!significant_hba1c_results$Metabolite %in% remove_hba1c_steiger$Metabolite),]
# 1 Removed 16 to 15


# Save the significant results data frames to a file
readr::write_tsv(significant_t2dm_results, "Filtered_T2DM_Results_diss.tsv")
readr::write_tsv(significant_fg_results, "Filtered_FG_Results_diss.tsv")
readr::write_tsv(significant_hba1c_results, "Filtered_HbA1C_Results_diss.tsv")

# Load in the Filtered_T2DM_Results.tsv file
filtered_t2dm_results <- readr::read_tsv("Filtered_T2DM_Results_diss.tsv")
# Load in the Filtered_FG_Results.tsv file
filtered_fg_results <- readr::read_tsv("Filtered_FG_Results_diss.tsv")
# Load in the Filtered_HbA1C_Results.tsv file
filtered_hba1c_results <- readr::read_tsv("Filtered_HbA1C_Results_diss.tsv")
