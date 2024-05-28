library(tidyverse)

# Load the t2dm_gwas_cleaned.tsv
t2dm_gwas_cleaned <- read_tsv("t2dm_gwas_cleaned.tsv")

# Remove rows with EAF of NA
t2dm_gwas_cleaned <- t2dm_gwas_cleaned %>% filter(!is.na(EAF))

# Load the RV_Filtered_T2DM_Results.tsv
rv_filtered_t2dm_results <- read_tsv("RV_Filtered_T2DM_Results.tsv")

# Initialise lists and counters
t2dm_coloc <- list()
metabolite_count <- 0

# For each metabolite in rv_filtered_t2dm_results, load the corresponding file from Harmonised_T2DM_IVs
for (metabolite in rv_filtered_t2dm_results$Metabolite) {
  harm_IVs <- read_tsv(paste0("Harmonised_T2DM_IVs/", metabolite, "T2DMT2DM_Harmonised_IVs.tsv"))
  
  # Add the metabolite to the list
  t2dm_coloc[[metabolite]] <- list()
  
  # Remove any rows with a value in the proxy column, keep if the value is NA
  harm_IVs <- harm_IVs %>% filter(is.na(proxy))
  
  # Add the number of rows to the list
  t2dm_coloc[[metabolite]]$num_IVs <- nrow(harm_IVs)
  # If there are no remaining IVs, print a message and continue to the next metabolite
  if (t2dm_coloc[[metabolite]]$num_IVs == 0) {
    message(paste("No IVs for", metabolite))
    next
  }
  print(paste("Metabolite has", t2dm_coloc[[metabolite]]$num_IVs, "IVs"))
  
  # Get the compid for the metabolite from rv_filtered_t2dm_results$compid
  compid <- rv_filtered_t2dm_results$compid[rv_filtered_t2dm_results$Metabolite == metabolite]
  
  # If the compid is less than 5 digits, add leading zeros
  if (nchar(compid) < 5) {
    compid <- paste0(strrep("0", 5 - nchar(compid)), compid)
  }
  
  # If the compid does not contain "M", add it to the beginning
  if (!grepl("M", compid)) {
    compid <- paste0("M", compid)
  }
  
  # For each IV in harm_IVs, do the following
  for (i in 1:nrow(harm_IVs)) {
    # Store the chromosome and position of the IV
    chromosome <- harm_IVs$Chromosome[i]
    position <- harm_IVs$Position[i]
    
    iv_chromosome <- paste0("chr", "_", chromosome)
    
    raw_gz_file <- paste0("Reverse_MR/Raw_Metab_GWAS_Data/", compid, "/", "INTERVAL_MRC-Epi_",  compid, "_sorted_", iv_chromosome, ".tbl.gz")
    R.utils::gunzip(raw_gz_file, overwrite = TRUE)

    raw_tbl_file <- paste0("Reverse_MR/Raw_Metab_GWAS_Data/", compid, "/","INTERVAL_MRC-Epi_", compid, "_sorted_", iv_chromosome, ".tbl")
    combined_data <- readr::read_tsv(raw_tbl_file)
    # Re-zip the file
    R.utils::gzip(raw_tbl_file, overwrite = TRUE)
    
    # Create a window +- 0.5 Mb around the IV
    window_05 <- c(chromosome, position - 500000, position + 500000)
    
    # Select the SNPs from t2dm_gwas_cleaned that are within the windows
    t2dm_snps_05 <- t2dm_gwas_cleaned %>% filter(Chromosome == window_05[1], Position >= window_05[2], Position <= window_05[3])

    # Select the SNPs from combined_data that are within the windows
    metabolite_snps_05 <- combined_data %>% filter(chrom == window_05[1], chromStart >= window_05[2], chromStart <= window_05[3])
    
    # Save the dataframes in Colocalisation/Preped_Regions/T2DM_Regions
    write_tsv(t2dm_snps_05, paste0("Colocalisation/Preped_Regions/T2DM_Regions/", metabolite, "_", harm_IVs$SNP[i], "_t2dm_snps.tsv"))
    write_tsv(metabolite_snps_05, paste0("Colocalisation/Preped_Regions/T2DM_Regions/", metabolite, "_", harm_IVs$SNP[i], "_metabolite_snps.tsv"))
  }
  print(paste("Finished", metabolite))
  metabolite_count <- metabolite_count + 1
  print(paste("Metabolite count:", metabolite_count))
}
# Make the list of lists into a dataframe
t2dm_coloc_df <- bind_rows(t2dm_coloc)
# Add the metabolite names as a column
t2dm_coloc_df$Metabolite <- names(t2dm_coloc)
# Switch the columns around
t2dm_coloc_df <- t2dm_coloc_df %>% select(Metabolite, num_IVs)
# Save the dataframe as a tsv
write_tsv(t2dm_coloc_df, "Colocalisation/T2DM_Coloc_IVs.tsv")








# Load the fasting_glucose_gwas_cleaned.tsv
fg_gwas_cleaned <- read_tsv("fasting_glucose_gwas_cleaned.tsv")

# Remove rows with EAF of NA
fg_gwas_cleaned <- fg_gwas_cleaned %>% filter(!is.na(EAF))

# Load the RV_Filtered_FG_Results.tsv
rv_filtered_fg_results <- read_tsv("RV_Filtered_FG_Results.tsv")

# Initialise lists and counters
fg_coloc <- list()
metabolite_count <- 0

# For each metabolite in rv_filtered_fg_results, load the corresponding file from Harmonised_FG_IVs
for (metabolite in rv_filtered_fg_results$Metabolite) {
  harm_IVs <- read_tsv(paste0("Harmonised_FG_IVs/", metabolite, "FGFG_Harmonised_IVs.tsv"))
  
  # Add the metabolite to the list
  fg_coloc[[metabolite]] <- list()
  
  # Remove any rows with a value in the proxy column, keep if the value is NA
  harm_IVs <- harm_IVs %>% filter(is.na(proxy))
  
  # Add the number of rows to the list
  fg_coloc[[metabolite]]$num_IVs <- nrow(harm_IVs)
  # If there are no remaining IVs, print a message and continue to the next metabolite
  if (fg_coloc[[metabolite]]$num_IVs == 0) {
    message(paste("No IVs for", metabolite))
    next
  }
  print(paste("Metabolite has", fg_coloc[[metabolite]]$num_IVs, "IVs"))
  
  # Get the compid for the metabolite from rv_filtered_fg_results$compid
  compid <- rv_filtered_fg_results$compid[rv_filtered_fg_results$Metabolite == metabolite]
  
  # If the compid is less than 5 digits, add leading zeros
  if (nchar(compid) < 5) {
    compid <- paste0(strrep("0", 5 - nchar(compid)), compid)
  }
  
  # If the compid does not contain "M", add it to the beginning
  if (!grepl("M", compid)) {
    compid <- paste0("M", compid)
  }
  
  # For each IV in harm_IVs, do the following
  for (i in 1:nrow(harm_IVs)) {
    # Store the chromosome and position of the IV
    chromosome <- harm_IVs$Chromosome[i]
    position <- harm_IVs$Position[i]
    
    iv_chromosome <- paste0("chr", "_", chromosome)
    
    raw_gz_file <- paste0("Reverse_MR/Raw_Metab_GWAS_Data/", compid, "/", "INTERVAL_MRC-Epi_",  compid, "_sorted_", iv_chromosome, ".tbl.gz")
    R.utils::gunzip(raw_gz_file, overwrite = TRUE)
    
    raw_tbl_file <- paste0("Reverse_MR/Raw_Metab_GWAS_Data/", compid, "/","INTERVAL_MRC-Epi_", compid, "_sorted_", iv_chromosome, ".tbl")
    combined_data <- readr::read_tsv(raw_tbl_file)
    # Re-zip the file
    R.utils::gzip(raw_tbl_file, overwrite = TRUE)
    
    # Create a window +- 0.5 Mb around the IV
    window_05 <- c(chromosome, position - 500000, position + 500000)
    
    # Select the SNPs from fg_gwas_cleaned that are within the windows
    fg_snps_05 <- fg_gwas_cleaned %>% filter(Chromosome == window_05[1], Position >= window_05[2], Position <= window_05[3])
    
    # Select the SNPs from combined_data that are within the windows
    metabolite_snps_05 <- combined_data %>% filter(chrom == window_05[1], chromStart >= window_05[2], chromStart <= window_05[3])
    
    # Save the dataframes in Colocalisation/Preped_Regions/FG_Regions
    write_tsv(fg_snps_05, paste0("Colocalisation/Preped_Regions/FG_Regions/", metabolite, "_", harm_IVs$SNP[i], "_fg_snps.tsv"))
    write_tsv(metabolite_snps_05, paste0("Colocalisation/Preped_Regions/FG_Regions/", metabolite, "_", harm_IVs$SNP[i], "_metabolite_snps.tsv"))
  }
  print(paste("Finished", metabolite))
  metabolite_count <- metabolite_count + 1
  print(paste("Metabolite count:", metabolite_count))
}
# Make the list of lists into a dataframe
fg_coloc_df <- bind_rows(fg_coloc)
# Add the metabolite names as a column
fg_coloc_df$Metabolite <- names(fg_coloc)
# Switch the columns around
fg_coloc_df <- fg_coloc_df %>% select(Metabolite, num_IVs)
# Save the dataframe as a tsv
write_tsv(fg_coloc_df, "Colocalisation/FG_Coloc_IVs.tsv")











# Load the hbA1c_gwas_cleaned.tsv
hba1c_gwas_cleaned <- read_tsv("hbA1c_gwas_cleaned.tsv")

# Remove rows with EAF of NA
hba1c_gwas_cleaned <- hba1c_gwas_cleaned %>% filter(!is.na(EAF))

# Load the RV_Filtered_HbA1C_Results.tsv
rv_filtered_hba1c_results <- read_tsv("RV_Filtered_HbA1C_Results.tsv")

# Initialise lists and counters
hba1c_coloc <- list()
metabolite_count <- 0

# For each metabolite in rv_filtered_hba1c_results, load the corresponding file from Harmonised_HBA1C_IVs
for (metabolite in rv_filtered_hba1c_results$Metabolite) {
  harm_IVs <- read_tsv(paste0("Harmonised_HBA1C_IVs/", metabolite, "HBA1CHBA1C_Harmonised_IVs.tsv"))
  
  # Add the metabolite to the list
  hba1c_coloc[[metabolite]] <- list()
  
  # Remove any rows with a value in the proxy column, keep if the value is NA
  harm_IVs <- harm_IVs %>% filter(is.na(proxy))
  
  # Add the number of rows to the list
  hba1c_coloc[[metabolite]]$num_IVs <- nrow(harm_IVs)
  # If there are no remaining IVs, print a message and continue to the next metabolite
  if (hba1c_coloc[[metabolite]]$num_IVs == 0) {
    message(paste("No IVs for", metabolite))
    next
  }
  print(paste("Metabolite has", hba1c_coloc[[metabolite]]$num_IVs, "IVs"))
  
  # Get the compid for the metabolite from rv_filtered_hba1c_results$compid
  compid <- rv_filtered_hba1c_results$compid[rv_filtered_hba1c_results$Metabolite == metabolite]
  
  # If the compid is less than 5 digits, add leading zeros
  if (nchar(compid) < 5) {
    compid <- paste0(strrep("0", 5 - nchar(compid)), compid)
  }
  
  # If the compid does not contain "M", add it to the beginning
  if (!grepl("M", compid)) {
    compid <- paste0("M", compid)
  }
  
  # For each IV in harm_IVs, do the following
  for (i in 1:nrow(harm_IVs)) {
    # Store the chromosome and position of the IV
    chromosome <- harm_IVs$Chromosome[i]
    position <- harm_IVs$Position[i]
    
    iv_chromosome <- paste0("chr", "_", chromosome)
    
    raw_gz_file <- paste0("Reverse_MR/Raw_Metab_GWAS_Data/", compid, "/", "INTERVAL_MRC-Epi_",  compid, "_sorted_", iv_chromosome, ".tbl.gz")
    R.utils::gunzip(raw_gz_file, overwrite = TRUE)
    
    raw_tbl_file <- paste0("Reverse_MR/Raw_Metab_GWAS_Data/", compid, "/","INTERVAL_MRC-Epi_", compid, "_sorted_", iv_chromosome, ".tbl")
    combined_data <- readr::read_tsv(raw_tbl_file)
    # Re-zip the file
    R.utils::gzip(raw_tbl_file, overwrite = TRUE)
    
    # Create a window +- 0.5 Mb around the IV
    window_05 <- c(chromosome, position - 500000, position + 500000)
    
    # Select the SNPs from hba1c_gwas_cleaned that are within the windows
    hba1c_snps_05 <- hba1c_gwas_cleaned %>% filter(Chromosome == window_05[1], Position >= window_05[2], Position <= window_05[3])
    
    # Select the SNPs from combined_data that are within the windows
    metabolite_snps_05 <- combined_data %>% filter(chrom == window_05[1], chromStart >= window_05[2], chromStart <= window_05[3])
    
    # Save the dataframes in Colocalisation/Preped_Regions/HBA1C_Regions
    write_tsv(hba1c_snps_05, paste0("Colocalisation/Preped_Regions/HBA1C_Regions/", metabolite, "_", harm_IVs$SNP[i], "_hba1c_snps.tsv"))
    write_tsv(metabolite_snps_05, paste0("Colocalisation/Preped_Regions/HBA1C_Regions/", metabolite, "_", harm_IVs$SNP[i], "_metabolite_snps.tsv"))
  }
  print(paste("Finished", metabolite))
  metabolite_count <- metabolite_count + 1
  print(paste("Metabolite count:", metabolite_count))
}
# Make the list of lists into a dataframe
hba1c_coloc_df <- bind_rows(hba1c_coloc)
# Add the metabolite names as a column
hba1c_coloc_df$Metabolite <- names(hba1c_coloc)
# Switch the columns around
hba1c_coloc_df <- hba1c_coloc_df %>% select(Metabolite, num_IVs)
# Save the dataframe as a tsv
write_tsv(hba1c_coloc_df, "Colocalisation/HBA1C_Coloc_IVs.tsv")

