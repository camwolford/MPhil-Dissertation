### Analyse MR Results ###
# This script is used to filter the MR results for significant associations.
library(tidyverse)

# Load in the T2DM_MR_Results file
T2DM_results <- readr::read_tsv("T2DM_MR_Results.tsv")

# Make a new table with only the rows that have a Fixed_IVW_Pval < 3.004807e-5 if the metabolite had 3 or less IVs
significant_T2DM_results <- T2DM_results[which(T2DM_results$Fixed_IVW_Pval < 3.004807e-5 & T2DM_results$Number_of_IVs <= 3),]
# Add in the metabolites that have more than 3 IVs and have a Random_IVW_Pval < 3.004807e-5
significant_T2DM_results <- rbind(significant_T2DM_results, T2DM_results[which(T2DM_results$Random_IVW_Pval < 3.004807e-5 & T2DM_results$Number_of_IVs > 3),])

# Check if any of the metabolties in significant_T2DM_results have a Egger_Pval, Weighted_Median_Pval, and Weighted_Mode_Pval > 0.05
check <- significant_T2DM_results[which(significant_T2DM_results$Egger_Pval > 0.05 & significant_T2DM_results$Weighted_Median_Pval > 0.05 & significant_T2DM_results$Weighted_Mode_Pval > 0.05),]

# Keep rows with Egger_Pval or Weighted_Median_Pval < 0.05 or Simple_Median_Pval < 0.05
holder <- significant_T2DM_results[which(significant_T2DM_results$Egger_Pval < 0.05 | significant_T2DM_results$Weighted_Median_Pval < 0.05 | significant_T2DM_results$Weighted_Mode_Pval < 0.05),]
# Check which ones did not pass the below check
holder[which(ifelse(holder$Egger_Pval < 0.05 & sign(holder$Egger_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) == FALSE |
               ifelse(holder$Weighted_Median_Pval < 0.05 & sign(holder$Weighted_Median_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) == FALSE |
               ifelse(holder$Weighted_Mode_Pval < 0.05 & sign(holder$Weighted_Mode_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) == FALSE),]
# Check that if the metabolite is significant in Egger, Weighted_Median, or Weighted_Mode, the effect is in the same direction as the Fixed_IVW
holder <- holder[which(ifelse(holder$Egger_Pval < 0.05 & sign(holder$Egger_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) &
                         ifelse(holder$Weighted_Median_Pval < 0.05 & sign(holder$Weighted_Median_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) &
                         ifelse(holder$Weighted_Mode_Pval < 0.05 & sign(holder$Weighted_Mode_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE)),]
significant_T2DM_results <- significant_T2DM_results[which(is.na(significant_T2DM_results$Egger_Pval)),]
significant_T2DM_results <- rbind(significant_T2DM_results, holder)
# Save the significant_T2DM_results table to a new file
write.table(significant_T2DM_results, "significant_T2DM_results_diss.tsv", sep = "\t", row.names = FALSE)
# 545 to 62


### Now for fasting glucose ###
# Load in the FG_MR_Results file
FG_results <- readr::read_tsv("FG_MR_Results.tsv")

# Make a new table with only the rows that have a Fixed_IVW_Pval < 3.004807e-5 if the metabolite had 3 or less IVs
significant_FG_results <- FG_results[which(FG_results$Fixed_IVW_Pval < 3.004807e-5 & FG_results$Number_of_IVs <= 3),]
# Add in the metabolites that have more than 3 IVs and have a Random_IVW_Pval < 3.004807e-5
significant_FG_results <- rbind(significant_FG_results, FG_results[which(FG_results$Random_IVW_Pval < 3.004807e-5 & FG_results$Number_of_IVs > 3),])

# Check if any of the metabolties in significant_FG_results have a Egger_Pval, Weighted_Median_Pval, and Weighted_Mode_Pval > 0.05
check <- significant_FG_results[which(significant_FG_results$Egger_Pval > 0.05 & significant_FG_results$Weighted_Median_Pval > 0.05 & significant_FG_results$Weighted_Mode_Pval > 0.05),]

# Keep rows with Egger_Pval < 0.05 or Weighted_Median_Pval < 0.05 or Simple_Median_Pval < 0.05
holder <- significant_FG_results[which(significant_FG_results$Egger_Pval < 0.05 | significant_FG_results$Weighted_Median_Pval < 0.05 | significant_FG_results$Weighted_Mode_Pval < 0.05),]
# Check which ones did not pass the below check
holder[which(ifelse(holder$Egger_Pval < 0.05 & sign(holder$Egger_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) == FALSE |
               ifelse(holder$Weighted_Median_Pval < 0.05 & sign(holder$Weighted_Median_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) == FALSE |
               ifelse(holder$Weighted_Mode_Pval < 0.05 & sign(holder$Weighted_Mode_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) == FALSE),]
# Check that if the metabolite is significant in Egger, Weighted_Median, or Weighted_Mode, the effect is in the same direction as the Fixed_IVW
holder <- holder[which(ifelse(holder$Egger_Pval < 0.05 & sign(holder$Egger_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) &
                         ifelse(holder$Weighted_Median_Pval < 0.05 & sign(holder$Weighted_Median_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) &
                         ifelse(holder$Weighted_Mode_Pval < 0.05 & sign(holder$Weighted_Mode_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE)),]
significant_FG_results <- significant_FG_results[which(is.na(significant_FG_results$Egger_Pval)),]
significant_FG_results <- rbind(significant_FG_results, holder)
# Save the significant_FG_results table to a new file
write.table(significant_FG_results, "significant_FG_results_diss.tsv", sep = "\t", row.names = FALSE)
# 559 to 32


### Now for hbA1c ###
# Load in the HBA1C_MR_Results file
HBA1C_results <- readr::read_tsv("HBA1C_MR_Results.tsv")

# Make a new table with only the rows that have a Fixed_IVW_Pval < 3.004807e-5 if the metabolite had 3 or less IVs
significant_HBA1C_results <- HBA1C_results[which(HBA1C_results$Fixed_IVW_Pval < 3.004807e-5 & HBA1C_results$Number_of_IVs <= 3),]
# Add in the metabolites that have more than 3 IVs and have a Random_IVW_Pval < 3.004807e-5
significant_HBA1C_results <- rbind(significant_HBA1C_results, HBA1C_results[which(HBA1C_results$Random_IVW_Pval < 3.004807e-5 & HBA1C_results$Number_of_IVs > 3),])

# Check if any of the metabolties in significant_HBA1C_results have a Egger_Pval, Weighted_Median_Pval, and Weighted_Mode_Pval > 0.05
check <- significant_HBA1C_results[which(significant_HBA1C_results$Egger_Pval > 0.05 & significant_HBA1C_results$Weighted_Median_Pval > 0.05 & significant_HBA1C_results$Weighted_Mode_Pval > 0.05),]

# Keep rows with Egger_Pval < 0.05 or Weighted_Median_Pval < 0.05 or Simple_Median_Pval < 0.05
holder <- significant_HBA1C_results[which(significant_HBA1C_results$Egger_Pval < 0.05 | significant_HBA1C_results$Weighted_Median_Pval < 0.05 | significant_HBA1C_results$Weighted_Mode_Pval < 0.05),]
# Check which ones did not pass the below check
holder[which(ifelse(holder$Egger_Pval < 0.05 & sign(holder$Egger_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) == FALSE |
               ifelse(holder$Weighted_Median_Pval < 0.05 & sign(holder$Weighted_Median_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) == FALSE |
               ifelse(holder$Weighted_Mode_Pval < 0.05 & sign(holder$Weighted_Mode_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) == FALSE),]
# Check that if the metabolite is significant in Egger, Weighted_Median, or Weighted_Mode, the effect is in the same direction as the Fixed_IVW
holder <- holder[which(ifelse(holder$Egger_Pval < 0.05 & sign(holder$Egger_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) &
                         ifelse(holder$Weighted_Median_Pval < 0.05 & sign(holder$Weighted_Median_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE) &
                         ifelse(holder$Weighted_Mode_Pval < 0.05 & sign(holder$Weighted_Mode_Estimate) != sign(holder$Fixed_IVW_Estimate), FALSE, TRUE)),]
significant_HBA1C_results <- significant_HBA1C_results[which(is.na(significant_HBA1C_results$Egger_Pval)),]
significant_HBA1C_results <- rbind(significant_HBA1C_results, holder)
# Save the significant_HBA1C_results table to a new file
write.table(significant_HBA1C_results, "significant_HBA1C_results_diss.tsv", sep = "\t", row.names = FALSE)
# 559 to 31
