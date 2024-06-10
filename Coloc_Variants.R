# Load the coloc_filtered_metabolites_data.tsv file
metabolites_data_frame <- readr::read_tsv("coloc_filtered_metabolites_data.tsv")
# Group the metabolites by subpathway
subpathways <- unique(metabolites_data_frame$subpathway)

# Load pwcoco_t2dm_coloc_out.coloc
pwcoco_t2dm_coloc_out <- read_delim("Colocalisation/pwcoco_out_t2dm.coloc")
# Load the pwcoco_fg_coloc_out.coloc
pwcoco_fgwas_t2dm_coloc_out <- read_delim("Colocalisation/pwcoco_out_fg.coloc")
# Load the pwcoco_hba1c_coloc_out.coloc
pwcoco_hba1c_coloc_out <- read_delim("Colocalisation/pwcoco_out_hba1c.coloc")
# Combine the pwcoco_t2dm_coloc_out, pwcoco_fg_coloc_out, and pwcoco_hba1c_coloc_out data frames
pwcoco_coloc_out <- bind_rows(pwcoco_t2dm_coloc_out, pwcoco_fgwas_t2dm_coloc_out, pwcoco_hba1c_coloc_out)

# Load in the T2DM_files
T2DM_files <- list.files("Harmonised_T2DM_IVs", full.names = TRUE)

# Initialise a dataframe with columns MarkerName and SNP
coloc_metab_variants <- data.frame(MarkerName = character(), SNP = character())

for (metabolite in metabolites_data_frame$name) {
  # Select the rows in pwcoco_coloc_out that contain the metabolite
  pwcoco_coloc_out_metabolite <- pwcoco_coloc_out %>% 
    filter(str_detect(Dataset1, fixed(metabolite)))
    
  # Remove the _metabolite_snps_pwcoco.tsv from the unique_Dataset1_values
  pwcoco_coloc_out_metabolite$Dataset1 <- str_replace_all(pwcoco_coloc_out_metabolite$Dataset1, "_metabolite_snps_pwcoco.tsv", "")
  # Remove all characters before "rs" in the unique_Dataset1_values
  pwcoco_coloc_out_metabolite$Dataset1 <- str_replace_all(pwcoco_coloc_out_metabolite$Dataset1, ".*rs", "rs")
    
  # Extract the unique rsids from the remaining Dataset1 values
  rsids <- pwcoco_coloc_out_metabolite$Dataset1 %>% unique()
    
  # Check if any of the H4 values are greater than or equal to 0.8
  for (rsid in rsids) {
    if (any(pwcoco_coloc_out_metabolite$H4[pwcoco_coloc_out_metabolite$Dataset1 == rsid] >= 0.8)) {
      # Get the rsid of that row
      rsid_high_h4 <- rsid
      # Remove all rows with that rsid_high_h4
      pwcoco_coloc_out_metabolite <- pwcoco_coloc_out_metabolite %>% filter(Dataset1 != rsid_high_h4)
    }
  }
    
  variants_to_remove <- unique(pwcoco_coloc_out_metabolite$Dataset1)
    
  # Get the IVs of each metabolite
  metabolite_file <- T2DM_files[str_detect(T2DM_files, fixed(metabolite, ignore_case = TRUE))]
  metabolite_data <- readr::read_tsv(metabolite_file)
    
  # Remove the variants in variants_to_remove from the metabolite_data
  metabolite_data <- metabolite_data %>% filter(!SNP %in% variants_to_remove)
    
  # Change _ to : in the MarkerName column
  metabolite_data$MarkerName <- str_replace_all(metabolite_data$MarkerName, "_", ":")
  # Remove "chr" from the MarkerName column
  metabolite_data$MarkerName <- str_replace_all(metabolite_data$MarkerName, "chr", "")
  
  # Remove columns 1,2,5-25
  metabolite_data <- metabolite_data[, -c(1, 2, 5:25)]
    
  # Add the metabolite_data to the coloc_metab_variants dataframe
  coloc_metab_variants <- bind_rows(coloc_metab_variants, metabolite_data)
}
# Distinct the coloc_metab_variants dataframe
coloc_metab_variants <- distinct(coloc_metab_variants)
# Save the tsv
write_tsv(coloc_metab_variants, "coloc_metab_variants.tsv")



# Pefrom PheWAS
library(tidyverse)
library(dplyr)
library(ieugwasr)
api_status()

coloc_metab_variants <- read_tsv("coloc_metab_variants.tsv")
coloc_metab_variants$traits <- NA

# For each SNP in the coloc_metab_variants dataframe, perform a PheWAS
for (i in 1:nrow(coloc_metab_variants)) {
  # Get the SNP
  SNP <- coloc_metab_variants$SNP[i]
  # Perform the PheWAS
  phewas_results <- ieugwasr::phewas(variants=SNP, pval=1e-5)
  # Make a list of the phewas_results$trait
  traits <- phewas_results$trait
  # Unique the traits
  unique_traits <- unique(traits)
  # Add the unique traits to the coloc_metab_variants dataframe
  coloc_metab_variants$traits[i] <- paste(unique_traits, collapse = ", ")
}

# Make a list of all the values in the traits column, split by ", "
traits_list <- coloc_metab_variants$traits %>% str_split(", ")
traits_list <- unlist(traits_list)
# Identify what the most common traits are
traits_table <- table(traits_list)
traits_table <- as.data.frame(traits_table)
# Save the tsv
write_tsv(traits_table, "traits_table.tsv")

# Read the traits_table.tsv
traits_table <- read_tsv("traits_table.tsv")
