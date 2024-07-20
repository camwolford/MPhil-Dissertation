### Analyse the Reverse MR Results ###
# This script was used to analyse the reverse MR results.
library(tidyverse)

# Load the data
reverse_t2dm <- readr::read_tsv("Reverse_MR/T2DM_Reverse_MR_Results.tsv")
reverse_fg <- readr::read_tsv("Reverse_MR/FG_Reverse_MR_Results.tsv")
reverse_hba1c <- readr::read_tsv("Reverse_MR/HBA1C_Reverse_MR_Results.tsv")
filtered_t2dm <- readr::read_tsv("Filtered_T2DM_Results_diss.tsv")
filtered_fg <- readr::read_tsv("Filtered_FG_Results_diss.tsv")
filtered_hba1c <- readr::read_tsv("Filtered_HbA1C_Results_diss.tsv")
metab_compids <- readr::read_tsv("Sig_Metabolites_Compid_diss.tsv")

# Convert the ":" and "/" in metab_compids$name to "_"
metab_compids$name <- gsub(":", "_", metab_compids$name)
metab_compids$name <- gsub("/", "_", metab_compids$name)

# Select the metabolites from reverse_t2dm that have a value of less than 0.05 in any of the following columns: 
# Weighted_Mode_Pval, Weighted_Median_Pval, Random_IVW_Pval, Egger_Pval
reverse_sig_t2dm <- reverse_t2dm %>% 
  filter(Random_IVW_Pval < 0.05 | Egger_Pval < 0.05)
# Do the same for reverse_fg and reverse_hba1c
reverse_sig_fg <- reverse_fg %>% 
  filter(Random_IVW_Pval < 0.05 | Egger_Pval < 0.05)
reverse_sig_hba1c <- reverse_hba1c %>%
  filter(Random_IVW_Pval < 0.05 | Egger_Pval < 0.05)

length(intersect(intersect(reverse_sig_t2dm$Metabolite, reverse_sig_fg$Metabolite), reverse_sig_hba1c$Metabolite))
length(intersect(reverse_sig_t2dm$Metabolite, reverse_sig_fg$Metabolite))
length(intersect(reverse_sig_t2dm$Metabolite, reverse_sig_hba1c$Metabolite))
length(intersect(reverse_sig_fg$Metabolite, reverse_sig_hba1c$Metabolite))


# Save a list of the unique metabolites from the combination of the three reverse dataframes
reverse_sig_metabs <- unique(c(reverse_sig_t2dm$Metabolite, reverse_sig_fg$Metabolite, reverse_sig_hba1c$Metabolite))

# Remove all but the first column from the filtered_t2dm dataframe
filtered_t2dm <- filtered_t2dm %>% 
  select(Metabolite)
# Remove all but the first column from the filtered_fg dataframe
filtered_fg <- filtered_fg %>% 
  select(Metabolite)
# Remove all but the first column from the filtered_hba1c dataframe
filtered_hba1c <- filtered_hba1c %>% 
  select(Metabolite)
# Add the compids from metab_compids to the filtered_t2dm dataframe
filtered_t2dm <- merge(filtered_t2dm, metab_compids, by.x = "Metabolite", by.y = "name", all.x = TRUE)
# Add the compids from metab_compids to the filtered_fg dataframe
filtered_fg <- merge(filtered_fg, metab_compids, by.x = "Metabolite", by.y = "name", all.x = TRUE)
# Add the compids from metab_compids to the filtered_hba1c dataframe
filtered_hba1c <- merge(filtered_hba1c, metab_compids, by.x = "Metabolite", by.y = "name", all.x = TRUE)
# Remove the Metabolites from filtered_t2dm that are in reverse_sig_metabs list
filtered_t2dm <- filtered_t2dm[!filtered_t2dm$compid %in% reverse_sig_metabs,]
# Remove the Metabolites from filtered_fg that are in reverse_sig_metabs list
filtered_fg <- filtered_fg[!filtered_fg$compid %in% reverse_sig_metabs,]
# Remove the Metabolites from filtered_hba1c that are in reverse_sig_metabs list
filtered_hba1c <- filtered_hba1c[!filtered_hba1c$compid %in% reverse_sig_metabs,]


# Save the reverse_sig_t2dm dataframe to a tsv file
write_tsv(reverse_sig_t2dm, "Reverse_MR/Reverse_Sig_T2DM_Results.tsv")
# Save the RV_filtered_t2dm dataframe to a tsv file
write_tsv(filtered_t2dm, "RV_Filtered_T2DM_Results.tsv")
# Save the reverse_sig_fg dataframe to a tsv file
write_tsv(reverse_sig_fg, "Reverse_MR/Reverse_Sig_FG_Results.tsv")
# Save the RV_filtered_fg dataframe to a tsv file
write_tsv(filtered_fg, "RV_Filtered_FG_Results.tsv")
# Save the reverse_sig_hba1c dataframe to a tsv file
write_tsv(reverse_sig_hba1c, "Reverse_MR/Reverse_Sig_HbA1C_Results.tsv")
# Save the RV_filtered_hba1c dataframe to a tsv file
write_tsv(filtered_hba1c, "RV_Filtered_HbA1C_Results.tsv")

# Check which metabolites are overlapping between the three filtered dataframes
overlap_metabs <- Reduce(intersect, list(filtered_t2dm$Metabolite, filtered_fg$Metabolite, filtered_hba1c$Metabolite))
overlap_metabs
total_metabs <- c(filtered_t2dm$Metabolite, filtered_fg$Metabolite, filtered_hba1c$Metabolite)
unique_metabs <- unique(total_metabs)

total_reverse_metabs <- c(reverse_sig_t2dm$Metabolite, reverse_sig_fg$Metabolite, reverse_sig_hba1c$Metabolite)
unique_reverse_metabs <- unique(total_reverse_metabs)