### Filter IVs ###
# This script filters the IVs based on the following criteria:
# 1. Exclude SNPs with a p-value > 5e-8 # Check if this is the right threshold
# 2. Exclude SNPs with a minor allele frequency (MAF) < 0.01 # Check if this is the right threshold
# 3. Exclude SNPs that are shared with more than 5 metabolites]
library(tidyverse)

### Perform the filtering on all metabolites ###
# Load in each metabolite data in /Individual_Metabolite_GWAS
metabolite_files <- list.files("Individual_Metabolite_GWAS", full.names = TRUE)

# Initialise a counter for the number of metabolites
metabolite_counter <- 0
metabolites_with_more_than_0_IVs <- 0

for (metabolite_file in metabolite_files) {
  # Load the metabolite data
  file_path <- metabolite_file
  metabolite_data <- readr::read_tsv(file_path)
  # Filter the SNPs
  filtered_SNPs <- metabolite_data %>%
    filter(Pval < p_value_threshold) %>%
    filter(EAF > MAF_threshold) %>%
    filter(!SNP %in% shared_SNPs$SNP)
  # Check if there are any IVs left
  if (nrow(filtered_SNPs) > 0) {
    metabolites_with_more_than_0_IVs <- metabolites_with_more_than_0_IVs + 1
    # Extract the part of the file name before the "_GWAS"
    metabolite_name <- strsplit(basename(metabolite_file), "_GWAS")[[1]][1]
    # Add "_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "_IVs")
    # Save the filtered SNPs
    write_tsv(filtered_SNPs, paste0("Filtered_IVs/", metabolite_name, ".tsv"))
  }
  # Increment the counter
  metabolite_counter <- metabolite_counter + 1
  print(paste("Processed", metabolite_counter, "metabolites"))
  print(paste("Metabolites with more than 0 IVs:", metabolites_with_more_than_0_IVs))
}