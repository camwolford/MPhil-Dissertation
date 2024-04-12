# Load in the T2DM_MR_Results file
T2DM_results <- readr::read_tsv("T2DM_MR_Results.tsv")

# Make a new table with only the rows that have a Fixed_IVW_Pval < 3.06e-5
significant_T2DM_results <- T2DM_results[which(T2DM_results$Fixed_IVW_Pval < 3.06e-5),]

# Keep rows with Fixed_IVW_FStat > 10
significant_T2DM_results <- significant_T2DM_results[which(significant_T2DM_results$Fixed_IVW_FStat > 10),]

# Keep rows with Random_IVW_Hetstat and Fixed_IVW_Hetstat > Number_of_IVs - 1
holder <- significant_T2DM_results[which(significant_T2DM_results$Fixed_IVW_HetStat > (significant_T2DM_results$Number_of_IVs - 1)),]
significant_T2DM_results <- significant_T2DM_results[which(is.na(significant_T2DM_results$Random_IVW_HetStat)),]
significant_T2DM_results <- rbind(significant_T2DM_results, holder)

# Keep rows with Egger_Pval < 0.05
holder <- significant_T2DM_results[which(significant_T2DM_results$Egger_Pval < 0.05),]
significant_T2DM_results <- significant_T2DM_results[which(is.na(significant_T2DM_results$Egger_Pval)),]
significant_T2DM_results <- rbind(significant_T2DM_results, holder)

# Keep rows with Egger_Intercept_Pval > 0.05
holder <- significant_T2DM_results[which(significant_T2DM_results$Egger_Intercept_Pval > 0.05),]
significant_T2DM_results <- significant_T2DM_results[which(is.na(significant_T2DM_results$Egger_Intercept_Pval)),]
significant_T2DM_results <- rbind(significant_T2DM_results, holder)
# 524 to 77 to 46


### Now for fasting glucose ###
# Load in the FG_MR_Results file
FG_results <- readr::read_tsv("FG_MR_Results.tsv")

# Make a new table with only the rows that have a Fixed_IVW_Pval < 3.06e-5
significant_FG_results <- FG_results[which(FG_results$Fixed_IVW_Pval < 3.06e-5),]

# Keep rows with Fixed_IVW_FStat > 10
significant_FG_results <- significant_FG_results[which(significant_FG_results$Fixed_IVW_FStat > 10),]

# Keep rows with Random_IVW_Hetstat and Fixed_IVW_Hetstat > Number_of_IVs - 1
holder <- significant_FG_results[which(significant_FG_results$Fixed_IVW_HetStat > (significant_FG_results$Number_of_IVs - 1)),]
significant_FG_results <- significant_FG_results[which(is.na(significant_FG_results$Random_IVW_HetStat)),]
significant_FG_results <- rbind(significant_FG_results, holder)

# Keep rows with Egger_Pval < 0.05
holder <- significant_FG_results[which(significant_FG_results$Egger_Pval < 0.05),]
significant_FG_results <- significant_FG_results[which(is.na(significant_FG_results$Egger_Pval)),]
significant_FG_results <- rbind(significant_FG_results, holder)

# Keep rows with Egger_Intercept_Pval > 0.05
holder <- significant_FG_results[which(significant_FG_results$Egger_Intercept_Pval > 0.05),]
significant_FG_results <- significant_FG_results[which(is.na(significant_FG_results$Egger_Intercept_Pval)),]
significant_FG_results <- rbind(significant_FG_results, holder)
# 553 to 46 to 25



### Now for hbA1c ###
# Load in the HBA1C_MR_Results file
HBA1C_results <- readr::read_tsv("HBA1C_MR_Results.tsv")

# Make a new table with only the rows that have a Fixed_IVW_Pval < 3.06e-5
significant_HBA1C_results <- HBA1C_results[which(HBA1C_results$Fixed_IVW_Pval < 3.06e-5),]

# Keep rows with Fixed_IVW_FStat > 10
significant_HBA1C_results <- significant_HBA1C_results[which(significant_HBA1C_results$Fixed_IVW_FStat > 10),]

# Keep rows with Random_IVW_Hetstat and Fixed_IVW_Hetstat > Number_of_IVs - 1
holder <- significant_HBA1C_results[which(significant_HBA1C_results$Fixed_IVW_HetStat > (significant_HBA1C_results$Number_of_IVs - 1)),]
significant_HBA1C_results <- significant_HBA1C_results[which(is.na(significant_HBA1C_results$Random_IVW_HetStat)),]
significant_HBA1C_results <- rbind(significant_HBA1C_results, holder)

# Keep rows with Egger_Pval < 0.05
holder <- significant_HBA1C_results[which(significant_HBA1C_results$Egger_Pval < 0.05),]
significant_HBA1C_results <- significant_HBA1C_results[which(is.na(significant_HBA1C_results$Egger_Pval)),]
significant_HBA1C_results <- rbind(significant_HBA1C_results, holder)

# Keep rows with Egger_Intercept_Pval > 0.05
holder <- significant_HBA1C_results[which(significant_HBA1C_results$Egger_Intercept_Pval > 0.05),]
significant_HBA1C_results <- significant_HBA1C_results[which(is.na(significant_HBA1C_results$Egger_Intercept_Pval)),]
significant_HBA1C_results <- rbind(significant_HBA1C_results, holder)
# 553 to 43 to 26


### Merge the significant results into one table ###
# Make a new dataframe to store the merged results
# First find the number of unique metabolites between the three datasets
unique_metabolites <- unique(c(significant_T2DM_results$Metabolite, significant_FG_results$Metabolite, significant_HBA1C_results$Metabolite))
print(length(unique_metabolites))
number_of_metabolites <- length(unique_metabolites)
# Make a dataframe with the first column as the unique metabolites
full_significant_results_df <- data.frame(Metabolite = unique_metabolites)

# Add the T2DM results to the dataframe
full_significant_results_df <- merge(full_significant_results_df, significant_T2DM_results, by = "Metabolite", all.x = TRUE)
# Add the FG results to the dataframe
full_significant_results_df <- merge(full_significant_results_df, significant_FG_results, by = "Metabolite", all.x = TRUE)
# Add the HBA1C results to the dataframe
full_significant_results_df <- merge(full_significant_results_df, significant_HBA1C_results, by = "Metabolite", all.x = TRUE)
# Rename the columns, any column ending in .x will be the T2DM results, any column ending in .y will be the FG results
colnames(full_significant_results_df) <- gsub("\\.x", "_T2DM", colnames(full_significant_results_df))
colnames(full_significant_results_df) <- gsub("\\.y", "_FG", colnames(full_significant_results_df))
# Add _HBA1C to the end of columns 60 onwards
colnames(full_significant_results_df)[60:ncol(full_significant_results_df)] <- paste(colnames(full_significant_results_df)[60:ncol(full_significant_results_df)], "_HBA1C", sep = "")

# Save the results as a tsv file
write.table(full_significant_results_df, file = "Full_Significant_Results.tsv", sep = "\t", row.names = FALSE)

