library(tidyverse)

# Load in the full significant results
sig_results <- readr::read_tsv("Full_Significant_Results.tsv")
# Load in the responses_combined.tsv file
responses <- readr::read_tsv("responses_combined.tsv")

# Convert any ":" or "/" in response names to "_"
responses$name <- gsub(":", "_", responses$name)
responses$name <- gsub("/", "_", responses$name)

sig_metabolites <- sig_results$Metabolite

sig_metabolites_data <- responses %>% filter(name %in% sig_metabolites)

# Make a table of the unique value counts of the superpathway
superpathway_counts <- sig_metabolites_data %>% 
  group_by(superpathway) %>% 
  summarise(n = n_distinct(name)) %>% 
  arrange(desc(n))
# Make a table of the unique value counts of the subpathway
subpathway_counts <- sig_metabolites_data %>% 
  group_by(subpathway) %>% 
  summarise(n = n_distinct(name)) %>% 
  arrange(desc(n))

# Make a dataframe with metabolite names and SNPs
metabolite_SNPs <- data.frame(name = sig_metabolites, T2DM_SNPs = NA, FG_SNPs = NA, HBA1C_SNPs = NA, stringsAsFactors = FALSE)

# For each sig_metabolite, retrieve the snps used
T2DM_files <- list.files("Harmonised_T2DM_IVs", full.names = TRUE)
# Remove files with ~
T2DM_files <- T2DM_files[!grepl("~", T2DM_files)]
FG_files <- list.files("Harmonised_FG_IVs", full.names = TRUE)
# Remove files with ~
FG_files <- FG_files[!grepl("~", FG_files)]
HBA1C_files <- list.files("Harmonised_HBA1C_IVs", full.names = TRUE)
# Remove files with ~
HBA1C_files <- HBA1C_files[!grepl("~", HBA1C_files)]

for (metabolite_file in T2DM_files) {
  # Read in the file 
  metabolite_data <- readr::read_tsv(metabolite_file)
  metabolite_name <- strsplit(basename(metabolite_file), "_Harmonised_IVs")[[1]][1]
  metabolite_name <- gsub("T2DM", "", metabolite_name)
  T2DM_IVs <- metabolite_data$SNP
  metabolite_SNPs$T2DM_SNPs[metabolite_SNPs$name == metabolite_name] <- paste(T2DM_IVs, collapse = ", ")
}
for (metabolite_file in FG_files) {
  # Read in the file 
  metabolite_data <- readr::read_tsv(metabolite_file)
  metabolite_name <- strsplit(basename(metabolite_file), "_Harmonised_IVs")[[1]][1]
  metabolite_name <- gsub("FG", "", metabolite_name)
  FG_IVs <- metabolite_data$SNP
  metabolite_SNPs$FG_SNPs[metabolite_SNPs$name == metabolite_name] <- paste(FG_IVs, collapse = ", ")
}
for (metabolite_file in HBA1C_files) {
  # Read in the file 
  metabolite_data <- readr::read_tsv(metabolite_file)
  metabolite_name <- strsplit(basename(metabolite_file), "_Harmonised_IVs")[[1]][1]
  metabolite_name <- gsub("HBA1C", "", metabolite_name)
  HBA1C_IVs <- metabolite_data$SNP
  metabolite_SNPs$HBA1C_SNPs[metabolite_SNPs$name == metabolite_name] <- paste(HBA1C_IVs, collapse = ", ")
}

# Select only the sig_metabolites from the metabolite_SNPs dataframe
sig_metabolites_data$T2DM_SNPs <- metabolite_SNPs$T2DM_SNPs[match(sig_metabolites_data$name, metabolite_SNPs$name)]
sig_metabolites_data$FG_SNPs <- metabolite_SNPs$FG_SNPs[match(sig_metabolites_data$name, metabolite_SNPs$name)]
sig_metabolites_data$HBA1C_SNPs <- metabolite_SNPs$HBA1C_SNPs[match(sig_metabolites_data$name, metabolite_SNPs$name)]

# Process T2DM SNPs
T2DM_expanded_snps <- sig_metabolites_data %>%
  filter(!is.na(T2DM_SNPs)) %>%
  mutate(T2DM_SNPs = strsplit(as.character(T2DM_SNPs), ",\\s*")) %>%
  unnest(T2DM_SNPs)
T2DM_snp_counts <- T2DM_expanded_snps %>%
  group_by(T2DM_SNPs) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(desc(n))
# For the Snps with more than one associated metabolite T2DM make a new column with a list of the superpathways and subpathways
T2DM_snp_counts$superpathway <- NA
for (snp in T2DM_snp_counts$T2DM_SNPs) {
  metabolites <- sig_metabolites_data %>%
    filter(!is.na(T2DM_SNPs)) %>%
    filter(str_detect(T2DM_SNPs, snp)) %>%
    pull(superpathway)
  T2DM_snp_counts$superpathway[T2DM_snp_counts$T2DM_SNPs == snp] <- paste(unique(metabolites), collapse = ", ")
}
T2DM_snp_counts$subpathway <- NA
for (snp in T2DM_snp_counts$T2DM_SNPs) {
  metabolites <- sig_metabolites_data %>%
    filter(!is.na(T2DM_SNPs)) %>%
    filter(str_detect(T2DM_SNPs, snp)) %>%
    pull(subpathway)
  T2DM_snp_counts$subpathway[T2DM_snp_counts$T2DM_SNPs == snp] <- paste(unique(metabolites), collapse = ", ")
}


# Process FG SNPs
FG_expanded_snps <- sig_metabolites_data %>%
  filter(!is.na(FG_SNPs)) %>%
  mutate(FG_SNPs = strsplit(as.character(FG_SNPs), ",\\s*")) %>%
  unnest(FG_SNPs)
FG_snp_counts <- FG_expanded_snps %>%
  group_by(FG_SNPs) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(desc(n))
# For the Snps with more than one associated metabolite FG make a new column with a list of the superpathways and subpathways
FG_snp_counts$superpathway <- NA
for (snp in FG_snp_counts$FG_SNPs) {
  metabolites <- sig_metabolites_data %>%
    filter(!is.na(FG_SNPs)) %>%
    filter(str_detect(FG_SNPs, snp)) %>%
    pull(superpathway)
  FG_snp_counts$superpathway[FG_snp_counts$FG_SNPs == snp] <- paste(unique(metabolites), collapse = ", ")
}
FG_snp_counts$subpathway <- NA
for (snp in FG_snp_counts$FG_SNPs) {
  metabolites <- sig_metabolites_data %>%
    filter(!is.na(FG_SNPs)) %>%
    filter(str_detect(FG_SNPs, snp)) %>%
    pull(subpathway)
  FG_snp_counts$subpathway[FG_snp_counts$FG_SNPs == snp] <- paste(unique(metabolites), collapse = ", ")
}

# Process HBA1C SNPs
HBA1C_expanded_snps <- sig_metabolites_data %>%
  filter(!is.na(HBA1C_SNPs)) %>%
  mutate(HBA1C_SNPs = strsplit(as.character(HBA1C_SNPs), ",\\s*")) %>%
  unnest(HBA1C_SNPs)
HBA1C_snp_counts <- HBA1C_expanded_snps %>%
  group_by(HBA1C_SNPs) %>%
  summarise(n = n(), .groups = 'drop') %>%
  arrange(desc(n))
# For the Snps with more than one associated metabolite HBA1C make a new column with a list of the superpathways and subpathways
HBA1C_snp_counts$superpathway <- NA
for (snp in HBA1C_snp_counts$HBA1C_SNPs) {
  metabolites <- sig_metabolites_data %>%
    filter(!is.na(HBA1C_SNPs)) %>%
    filter(str_detect(HBA1C_SNPs, snp)) %>%
    pull(superpathway)
  HBA1C_snp_counts$superpathway[HBA1C_snp_counts$HBA1C_SNPs == snp] <- paste(unique(metabolites), collapse = ", ")
}
HBA1C_snp_counts$subpathway <- NA
for (snp in HBA1C_snp_counts$HBA1C_SNPs) {
  metabolites <- sig_metabolites_data %>%
    filter(!is.na(HBA1C_SNPs)) %>%
    filter(str_detect(HBA1C_SNPs, snp)) %>%
    pull(subpathway)
  HBA1C_snp_counts$subpathway[HBA1C_snp_counts$HBA1C_SNPs == snp] <- paste(unique(metabolites), collapse = ", ")
}


# Save the superpathway_counts, subpathway_counts, T2DM_snp_counts, FG_snp_counts, HBA1C_snp_counts and significant_metabolites_data dataframes as TSV files
write_tsv(superpathway_counts, "superpathway_counts.tsv")
write_tsv(subpathway_counts, "subpathway_counts.tsv")
write_tsv(T2DM_snp_counts, "T2DM_snp_counts.tsv")
write_tsv(FG_snp_counts, "FG_snp_counts.tsv")
write_tsv(HBA1C_snp_counts, "HBA1C_snp_counts.tsv")
write_tsv(sig_metabolites_data, "significant_metabolites_data.tsv")

# Plot the distribution of the number of metabolites associated with each SNP
ggplot(T2DM_snp_counts, aes(n)) +
  geom_histogram(binwidth = 1) +
  labs(title = "Distribution of the number of metabolites associated with each T2DM SNP",
       x = "Number of metabolites associated with SNP",
       y = "Number of SNPs") +
  theme_minimal()

# Calculate the number of T2DM SNP associations that have more than one associated metabolite for each superpathway
superpathway_T2DM_counts <- T2DM_snp_counts %>%
  filter(!is.na(superpathway)) %>%
  group_by(superpathway) %>%
  summarise(n = sum(n > 1), .groups = 'drop') %>%
  arrange(desc(n))

# Remove the superpathways with no SNPs with more than one associated metabolite
superpathway_T2DM_counts <- superpathway_T2DM_counts %>%
  filter(n > 0)

# Plot the number of T2DM SNP associations that have more than one associated metabolite for each superpathway 
ggplot(superpathway_T2DM_counts, aes(reorder(superpathway, n), n)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of T2DM SNP associations with more than one associated metabolite \n for each superpathway",
       x = "Superpathway",
       y = "Number of SNPs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Adjust horizontal justification for x-axis labels
        plot.title = element_text(hjust = 0.5))  # Center the title

# Calculate the number of T2DM SNP associations that have more than one associated metabolite for each subpathway
subpathway_T2DM_counts <- T2DM_snp_counts %>%
  filter(!is.na(subpathway)) %>%
  group_by(subpathway) %>%
  summarise(n = sum(n > 1), .groups = 'drop') %>%
  arrange(desc(n))

# Remove the subpathways with no SNPs with more than one associated metabolite
subpathway_T2DM_counts <- subpathway_T2DM_counts %>%
  filter(n > 0)

# Remove any rows with more than one subpathway
subpathway_T2DM_counts <- subpathway_T2DM_counts %>%
  filter(!str_detect(subpathway, ", "))

# Plot the number of T2DM SNP associations that have more than one associated metabolite for each subpathway
ggplot(subpathway_T2DM_counts, aes(reorder(subpathway, n), n)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of T2DM SNP associations with more than one associated metabolite \n for each subpathway",
       x = "Subpathway",
       y = "Number of SNPs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Adjust horizontal justification for x-axis labels
        plot.title = element_text(hjust = 0.5))  # Center the title
  
