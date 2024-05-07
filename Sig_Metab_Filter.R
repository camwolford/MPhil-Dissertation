# Load the Metabolite_T2DM_Sensitivity_Analysis.tsv file
metabolite_t2dm_sensitivity_analysis <- readr::read_tsv("Metabolite_T2DM_Sensitivity_Analysis.tsv")

# Load in the significant_T2DM_results.tsv file
significant_t2dm_results <- readr::read_tsv("significant_T2DM_results.tsv")

# Select the metabolites that have a MR_PRESSO_Outlier_Pvals < 0.05
metabolites_to_remove <- metabolite_t2dm_sensitivity_analysis %>% 
  dplyr::filter(MR_PRESSO_Outlier_Pvals < 0.05)
# Check the value of the Liberal_Flag
metabolites_to_remove$Liberal_Flag # NA
# Add metabolites with a Direction_Flag of FALSE to the metabolites_to_remove data frame keeping the current metabolites in the data frame
metabolites_to_remove <- metabolites_to_remove %>% 
  dplyr::filter(Direction_Flag == FALSE) %>% 
  dplyr::select(Metabolite) %>% 
  dplyr::bind_rows(metabolites_to_remove)
# Add metabolites with a MR_RAPS_Pval > 0.05 to the metabolites_to_remove data frame keeping the current metabolites in the data frame
metabolites_to_remove <- metabolites_to_remove %>% 
  dplyr::filter(MR_RAPS_Pval > 0.05) %>% 
  dplyr::select(Metabolite) %>% 
  dplyr::bind_rows(metabolites_to_remove)

significant_t2dm_results[significant_t2dm_results$Metabolite %in% metabolites_to_remove$Metabolite,] <- NA

# One metabolite removed 57 to 57


# Load the Metabolite_FG_Sensitivity_Analysis.tsv file
metabolite_fg_sensitivity_analysis <- readr::read_tsv("Metabolite_FG_Sensitivity_Analysis.tsv")

# Load in the significant_FG_results.tsv file
significant_fg_results <- readr::read_tsv("significant_FG_results.tsv")

# Select the metabolites that have a MR_PRESSO_Outlier_Pvals < 0.05
metabolites_to_remove <- metabolite_fg_sensitivity_analysis %>% 
  dplyr::filter(MR_PRESSO_Outlier_Pvals < 0.05)
# Check the value of the Liberal_Flag
metabolites_to_remove$Liberal_Flag # NA
# Add metabolites with a Direction_Flag of FALSE to the metabolites_to_remove data frame keeping the current metabolites in the data frame
metabolites_to_remove <- metabolites_to_remove %>% 
  dplyr::filter(Direction_Flag == FALSE) %>% 
  dplyr::select(Metabolite) %>% 
  dplyr::bind_rows(metabolites_to_remove)
# Add metabolites with a MR_RAPS_Pval > 0.05 to the metabolites_to_remove data frame keeping the current metabolites in the data frame
metabolites_to_remove <- metabolites_to_remove %>% 
  dplyr::filter(MR_RAPS_Pval > 0.05) %>% 
  dplyr::select(Metabolite) %>% 
  dplyr::bind_rows(metabolites_to_remove)

significant_fg_results[significant_fg_results$Metabolite %in% metabolites_to_remove$Metabolite,] <- NA

# Two metabolites removed 29 to 29


# Load the Metabolite_HBA1C_Sensitivity_Analysis.tsv file
metabolite_hba1c_sensitivity_analysis <- readr::read_tsv("Metabolite_HBA1C_Sensitivity_Analysis.tsv")

# Load in the significant_HBA1C_results.tsv file
significant_hba1c_results <- readr::read_tsv("significant_HBA1C_results.tsv")

# Select the metabolites that have a MR_PRESSO_Outlier_Pvals < 0.05
metabolites_to_remove <- metabolite_hba1c_sensitivity_analysis %>% 
  dplyr::filter(MR_PRESSO_Outlier_Pvals < 0.05)
# Check the value of the Liberal_Flag
metabolites_to_remove$Liberal_Flag # NA
# Add metabolites with a Direction_Flag of FALSE to the metabolites_to_remove data frame keeping the current metabolites in the data frame
metabolites_to_remove <- metabolites_to_remove %>% 
  dplyr::filter(Direction_Flag == FALSE) %>% 
  dplyr::select(Metabolite) %>% 
  dplyr::bind_rows(metabolites_to_remove)
# Add metabolites with a MR_RAPS_Pval > 0.05 to the metabolites_to_remove data frame keeping the current metabolites in the data frame
metabolites_to_remove <- metabolites_to_remove %>% 
  dplyr::filter(MR_RAPS_Pval > 0.05) %>% 
  dplyr::select(Metabolite) %>% 
  dplyr::bind_rows(metabolites_to_remove)

significant_hba1c_results[significant_hba1c_results$Metabolite %in% metabolites_to_remove$Metabolite,] <- NA

# One metabolite removed 30 to 30

# Remove the rows with NA values in the Metabolite column for each of the significant results data frames
significant_t2dm_results <- significant_t2dm_results %>% 
  dplyr::filter(!is.na(Metabolite))
significant_fg_results <- significant_fg_results %>%
  dplyr::filter(!is.na(Metabolite))
significant_hba1c_results <- significant_hba1c_results %>%
  dplyr::filter(!is.na(Metabolite))


# Select the metabolites from the significant T2DM results that have an Egger_Pval < 0.05
metabolites_to_remove <- significant_t2dm_results %>% 
# Select the remaining metabolites from the metabolites_to_remove that have an Egger_Intercept_Pval < 0.05
  dplyr::filter(Egger_Intercept_Pval < 0.05) %>% 
  dplyr::select(Metabolite)
# Remove the metabolites from the significant T2DM results data frame
significant_t2dm_results[significant_t2dm_results$Metabolite %in% metabolites_to_remove$Metabolite,] <- NA
# Remove the rows with NA values in the Metabolite column for the significant T2DM results data frame
significant_t2dm_results <- significant_t2dm_results %>% 
  dplyr::filter(!is.na(Metabolite))
# Removed four metabolites 57 to 53

# Select the metabolites from the significant FG results that have an Egger_Pval < 0.05
metabolites_to_remove <- significant_fg_results %>% 
# Select the remaining metabolites from the metabolites_to_remove that have an Egger_Intercept_Pval < 0.05
  dplyr::filter(Egger_Intercept_Pval < 0.05) %>% 
  dplyr::select(Metabolite)
# Remove the metabolites from the significant FG results data frame
significant_fg_results[significant_fg_results$Metabolite %in% metabolites_to_remove$Metabolite,] <- NA
# Remove the rows with NA values in the Metabolite column for the significant FG results data frame
significant_fg_results <- significant_fg_results %>% 
  dplyr::filter(!is.na(Metabolite))
# Removed one metabolite 29 to 28

# Select the metabolites from the significant HbA1C results that have an Egger_Pval < 0.05
metabolites_to_remove <- significant_hba1c_results %>% 
# Select the remaining metabolites from the metabolites_to_remove that have an Egger_Intercept_Pval < 0.05
  dplyr::filter(Egger_Intercept_Pval < 0.05) %>% 
  dplyr::select(Metabolite)
# Remove the metabolites from the significant HbA1C results data frame
significant_hba1c_results[significant_hba1c_results$Metabolite %in% metabolites_to_remove$Metabolite,] <- NA
# Remove the rows with NA values in the Metabolite column for the significant HbA1C results data frame
significant_hba1c_results <- significant_hba1c_results %>% 
  dplyr::filter(!is.na(Metabolite))
# Removed two metabolites 30 to 28


# Save the significant results data frames to a file
readr::write_tsv(significant_t2dm_results, "Filtered_T2DM_Results.tsv")
readr::write_tsv(significant_fg_results, "Filtered_FG_Results.tsv")
readr::write_tsv(significant_hba1c_results, "Filtered_HbA1C_Results.tsv")

# Load in the Filtered_T2DM_Results.tsv file
filtered_t2dm_results <- readr::read_tsv("Filtered_T2DM_Results.tsv")
# Load in the Filtered_FG_Results.tsv file
filtered_fg_results <- readr::read_tsv("Filtered_FG_Results.tsv")
# Load in the Filtered_HbA1C_Results.tsv file
filtered_hba1c_results <- readr::read_tsv("Filtered_HbA1C_Results.tsv")

# Generate a Full_Filtered_Results.tsv file
### Merge the significant results into one table ###
# Make a new dataframe to store the merged results
# First find the number of unique metabolites between the three datasets
unique_metabolites <- unique(c(filtered_t2dm_results$Metabolite, filtered_fg_results$Metabolite, filtered_hba1c_results$Metabolite))
print(length(unique_metabolites))
number_of_metabolites <- length(unique_metabolites)
# Make a dataframe with the first column as the unique metabolites
Full_Filtered_Results_df <- data.frame(Metabolite = unique_metabolites)

# Add the T2DM results to the dataframe
Full_Filtered_Results_df <- merge(Full_Filtered_Results_df, filtered_t2dm_results, by = "Metabolite", all.x = TRUE)
# Add the FG results to the dataframe
Full_Filtered_Results_df <- merge(Full_Filtered_Results_df, filtered_fg_results, by = "Metabolite", all.x = TRUE)
# Add the HBA1C results to the dataframe
Full_Filtered_Results_df <- merge(Full_Filtered_Results_df, filtered_hba1c_results, by = "Metabolite", all.x = TRUE)
# Rename the columns, any column ending in .x will be the T2DM results, any column ending in .y will be the FG results
colnames(Full_Filtered_Results_df) <- gsub("\\.x", "_T2DM", colnames(Full_Filtered_Results_df))
colnames(Full_Filtered_Results_df) <- gsub("\\.y", "_FG", colnames(Full_Filtered_Results_df))
# Add _HBA1C to the end of columns 60 onwards
colnames(Full_Filtered_Results_df)[60:ncol(Full_Filtered_Results_df)] <- paste(colnames(Full_Filtered_Results_df)[60:ncol(Full_Filtered_Results_df)], "_HBA1C", sep = "")

# Save the Full_Filtered_Results_df to a file
readr::write_tsv(Full_Filtered_Results_df, "Full_Filtered_Results.tsv")

# Read in responses_combined.tsv
responses_combined <- readr::read_tsv("responses_combined.tsv")

# Make a new dataframe to store the Metabolite Names and compids for the significant results
Sig_Metabolites <- data.frame(Metabolite = Full_Filtered_Results_df$Metabolite)
# Attempt to match the Metabolite names to responses_combined$name using fuzzy matching
Sig_Metabolites$Index <- stringdist::amatch(Sig_Metabolites$Metabolite, responses_combined$name, maxDist = 3, method = "jw")
Sig_Metabolites$name <- responses_combined$name[Sig_Metabolites$Index]
Sig_Metabolites$compid <- responses_combined$compid[Sig_Metabolites$Index]
# Remove the Index column and Metabolite column
Sig_Metabolites <- Sig_Metabolites %>% dplyr::select(-Index, -Metabolite)
# Save the Sig_Metabolites dataframe to a file
readr::write_tsv(Sig_Metabolites, "Sig_Metabolites_Compid.tsv")





