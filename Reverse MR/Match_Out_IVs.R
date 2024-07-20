### Match Outcome IVs to Metabolite Data ###
# This script was used to match the outcome IVs to the metabolite data. The outcome IVs' F-statistics were calculated. 
# The metabolite data was then read in and matched to the outcome IVs based on the chromosome, position, effect allele, and other allele.
library(tidyverse)
library(R.utils)

# Read in clumped_t2dm_ivs
clumped_t2dm_ivs <- readr::read_tsv("Reverse_MR/clumped_t2dm_ivs.tsv")
# Read in clumped_fg_ivs
clumped_fg_ivs <- readr::read_tsv("Reverse_MR/clumped_fg_ivs.tsv")
# Read in clumped_hba1c_ivs
clumped_hba1c_ivs <- readr::read_tsv("Reverse_MR/clumped_hba1c_ivs.tsv")

# Make MarkerName column in clumped_t2dm_ivs, clumped_fg_ivs, and clumped_hba1c_ivs
clumped_t2dm_ivs$MarkerName <- paste0("chr", clumped_t2dm_ivs$Chromosome, ":", clumped_t2dm_ivs$`Position     (bp, b37)`,
                                      ":", clumped_t2dm_ivs$Risk, ":", clumped_t2dm_ivs$Other)
clumped_fg_ivs$MarkerName <- paste0("chr", clumped_fg_ivs$Chr, ":", clumped_fg_ivs$`Pos (bp)`,
                                    ":", clumped_fg_ivs$`Effect Allele`, ":", clumped_fg_ivs$`Other Allele`)
clumped_hba1c_ivs$MarkerName <- paste0("chr", clumped_hba1c_ivs$Chr, ":", clumped_hba1c_ivs$`Pos (bp)`,
                                       ":", clumped_hba1c_ivs$`Effect Allele`, ":", clumped_hba1c_ivs$`Other Allele`)
# Make reverseMarkerName column in clumped_t2dm_ivs, clumped_fg_ivs, and clumped_hba1c_ivs
clumped_t2dm_ivs$reverseMarkerName <- paste0("chr", clumped_t2dm_ivs$Chromosome, ":", clumped_t2dm_ivs$`Position     (bp, b37)`,
                                      ":", clumped_t2dm_ivs$Other, ":", clumped_t2dm_ivs$Risk)
clumped_fg_ivs$reverseMarkerName <- paste0("chr", clumped_fg_ivs$Chr, ":", clumped_fg_ivs$`Pos (bp)`,
                                    ":", clumped_fg_ivs$`Other Allele`, ":", clumped_fg_ivs$`Effect Allele`)
clumped_hba1c_ivs$reverseMarkerName <- paste0("chr", clumped_hba1c_ivs$Chr, ":", clumped_hba1c_ivs$`Pos (bp)`,
                                       ":", clumped_hba1c_ivs$`Other Allele`, ":", clumped_hba1c_ivs$`Effect Allele`)

# Calculate the Fstat column in clumped_t2dm_ivs, clumped_fg_ivs, and clumped_hba1c_ivs
clumped_t2dm_ivs$Fstat <- (clumped_t2dm_ivs$`Log-OR` / clumped_t2dm_ivs$SE)^2
clumped_fg_ivs$Fstat <- (clumped_fg_ivs$Effect / clumped_fg_ivs$SE)^2
clumped_hba1c_ivs$Fstat <- (clumped_hba1c_ivs$Effect / clumped_hba1c_ivs$SE)^2
# Check for any Fstat values that are less than 10
clumped_t2dm_ivs[clumped_t2dm_ivs$Fstat < 10, ]
clumped_fg_ivs[clumped_fg_ivs$Fstat < 10, ]
clumped_hba1c_ivs[clumped_hba1c_ivs$Fstat < 10, ]
# No weak instruments

# Load in the file names in the Reverse_MR/Raw_Metab_GWAS_Data folder
file_names <- list.files("Reverse_MR/Raw_Metab_GWAS_Data", full.names = TRUE)

# Initialise counters
metabolite_counter <- 0
chromosome_counter <- 0
chromosome_counter_list <- c()
matched_T2DM_IVs <- 0
matched_FG_IVs <- 0
matched_HBA1C_IVs <- 0
T2DM_removed_IVs <- 0
FG_removed_IVs <- 0
HBA1C_removed_IVs <- 0
compid_IVs_list <- c()


for (i in 1:length(file_names)) {
  # Set the base directory for files
  base_dir <- "Reverse_MR/Raw_Metab_GWAS_Data"
  
  # Get list of .tbl.gz files in the directory
  raw_gz_files <- list.files(file_names[i], pattern = "\\.gz$", full.names = TRUE)
  
  # Decompress .gz files that contain 'chr' in their names
  for (file_path in raw_gz_files) {
    if (grepl("chr", file_path)) {
      R.utils::gunzip(file_path, overwrite = TRUE)
    }
  }
  
  # Get list of .tbl files now available
  raw_tbl_files <- list.files(file_names[i], pattern = "\\.tbl$", full.names = TRUE)
  
  # Initialize an empty list to store data frames
  all_data <- list()
  
  # Read each .tbl file, assuming they are TSVs
  for (file_path in raw_tbl_files) {
    if (grepl("chr", file_path)) {
      tbl_data <- readr::read_tsv(file_path)
      all_data[[basename(file_path)]] <- tbl_data
      chromosome_counter <- chromosome_counter + 1
      # Re-zip the file
      R.utils::gzip(file_path, overwrite = TRUE)
    }
  }
  
  # Combine all data frames into one
  combined_data <- dplyr::bind_rows(all_data, .id = "source_file")
  metabolite_counter <- metabolite_counter + 1
  chromosome_counter_list <- c(chromosome_counter_list, chromosome_counter)
  chromosome_counter <- 0
  # Separate the metabolite compid from the file_names[i] as it is the part after the final /
  metabolite_compid <- strsplit(file_names[i], "/")[[1]][3]
  
  
  ### Match IVs ###
  # Initialise a match flag
  match_flag <- FALSE
  # Make an empty data frame to store the matched IVs with the same number of rows as clumped_t2dm_ivs
  matched_ivs <- as.data.frame(matrix(NA, nrow = nrow(clumped_t2dm_ivs)))
  # Initialise columns
  matched_ivs$Chromosome <- clumped_t2dm_ivs$Chromosome
  matched_ivs$Position <- clumped_t2dm_ivs$`Position     (bp, b37)`
  matched_ivs$MarkerName <- clumped_t2dm_ivs$MarkerName
  matched_ivs$reverseMarkerName <- clumped_t2dm_ivs$reverseMarkerName
  matched_ivs$SNP <- clumped_t2dm_ivs$rsid
  matched_ivs$EffectAllele <- clumped_t2dm_ivs$Risk
  matched_ivs$NonEffectAllele <- clumped_t2dm_ivs$Other
  matched_ivs$EAF <- clumped_t2dm_ivs$EAF
  matched_ivs$Beta <- clumped_t2dm_ivs$`Log-OR`
  matched_ivs$SE <- clumped_t2dm_ivs$SE
  matched_ivs$Pval <- clumped_t2dm_ivs$pval
  matched_ivs$Fstat <- clumped_t2dm_ivs$Fstat
  matched_ivs$exposure <- "T2DM"
  matched_ivs$outcome <- metabolite_compid
  matched_ivs$out_Chromosome <- NA
  matched_ivs$out_Position <- NA
  matched_ivs$out_EffectAllele <- NA
  matched_ivs$out_NonEffectAllele <- NA
  matched_ivs$out_Beta <- NA
  matched_ivs$out_SE <- NA
  matched_ivs$out_EAF <- NA
  matched_ivs$out_Pval <- NA
  # Remove the first column
  matched_ivs <- matched_ivs[, -1]
  # For each outcome
  # Use the MarkerName and reverseMarkerName columns to match the IVs
  for (k in 1:nrow(matched_ivs)) {
    if (matched_ivs$MarkerName[k] %in% combined_data$MarkerName) {
      # Select the row from combined_data where the MarkerName matches the IV
      out_values <- combined_data[combined_data$MarkerName == matched_ivs$MarkerName[k], ]
      matched_ivs$out_Chromosome[k] <- out_values$chrom
      matched_ivs$out_Position[k] <- out_values$chromStart
      matched_ivs$out_EffectAllele[k] <- out_values$Allele1
      matched_ivs$out_NonEffectAllele[k] <- out_values$Allele2
      matched_ivs$out_Beta[k] <- out_values$Effect
      matched_ivs$out_SE[k] <- out_values$StdErr
      matched_ivs$out_EAF[k] <- out_values$Freq1
      matched_ivs$out_Pval[k] <- out_values$Pvalue
      match_flag <- TRUE
      matched_T2DM_IVs <- matched_T2DM_IVs + 1
    } else if (matched_ivs$reverseMarkerName[k] %in% combined_data$MarkerName) {
      # Select the row from combined_data where the reverseMarkerName matches the IV
      out_values <- combined_data[combined_data$MarkerName == matched_ivs$reverseMarkerName[k], ]
      matched_ivs$out_Chromosome[k] <- out_values$chrom
      matched_ivs$out_Position[k] <- out_values$chromStart
      matched_ivs$out_EffectAllele[k] <- out_values$Allele1
      matched_ivs$out_NonEffectAllele[k] <- out_values$Allele2
      matched_ivs$out_Beta[k] <- out_values$Effect
      matched_ivs$out_SE[k] <- out_values$StdErr
      matched_ivs$out_EAF[k] <- out_values$Freq1
      matched_ivs$out_Pval[k] <- out_values$Pvalue
      match_flag <- TRUE
      matched_T2DM_IVs <- matched_T2DM_IVs + 1
    }
    # If the match_flag is still FALSE, remove the IV
    if (match_flag == FALSE) {
      T2DM_removed_IVs <- T2DM_removed_IVs + 1
    }
    match_flag <- FALSE
  }
  # Make the out_EffectAllele and out_NonEffectAllele columns capitalised
  matched_ivs$out_EffectAllele <- toupper(matched_ivs$out_EffectAllele)
  matched_ivs$out_NonEffectAllele <- toupper(matched_ivs$out_NonEffectAllele)
  # Remove rows with NA in the out_Beta column
  matched_ivs <- matched_ivs[!is.na(matched_ivs$out_Beta), ]
  # Add the metabolite_compid to the compid_IVs_list
  compid_IVs_list <- c(compid_IVs_list, metabolite_compid)
  # Add the number of IVs to the compid_IVs_list
  compid_IVs_list <- c(compid_IVs_list, nrow(matched_ivs))
  # Save the matched IVs to a new folder (make sure this folder exists or add code to create it)
  output_dir <- "Reverse_MR/Matched_Reverse_MR_T2DM_IVs"
  # Select the part after the last / in the file path
  output_file_name <- paste0(basename(file_names[i]), "_T2DM_matched.tsv")
  write.table(matched_ivs, file.path(output_dir, output_file_name), sep = "\t", row.names = FALSE)
  
  
  ## Do the same for fg_ivs ##
  # Initialise a match flag
  match_flag <- FALSE
  # Make an empty data frame to store the matched IVs
  matched_ivs <- as.data.frame(matrix(NA, nrow = nrow(clumped_fg_ivs)))
  # Initialise columns
  matched_ivs$Chromosome <- clumped_fg_ivs$Chr
  matched_ivs$Position <- clumped_fg_ivs$`Pos (bp)`
  matched_ivs$MarkerName <- clumped_fg_ivs$MarkerName
  matched_ivs$reverseMarkerName <- clumped_fg_ivs$reverseMarkerName
  matched_ivs$SNP <- clumped_fg_ivs$rsid
  matched_ivs$EffectAllele <- clumped_fg_ivs$`Effect Allele`
  matched_ivs$NonEffectAllele <- clumped_fg_ivs$`Other Allele`
  matched_ivs$EAF <- clumped_fg_ivs$EAF
  matched_ivs$Beta <- clumped_fg_ivs$Effect
  matched_ivs$SE <- clumped_fg_ivs$SE
  matched_ivs$Pval <- clumped_fg_ivs$pval
  matched_ivs$Fstat <- clumped_fg_ivs$Fstat
  matched_ivs$exposure <- "FG"
  matched_ivs$outcome <- metabolite_compid
  matched_ivs$out_Chromosome <- NA
  matched_ivs$out_Position <- NA
  matched_ivs$out_EffectAllele <- NA
  matched_ivs$out_NonEffectAllele <- NA
  matched_ivs$out_Beta <- NA
  matched_ivs$out_SE <- NA
  matched_ivs$out_EAF <- NA
  matched_ivs$out_Pval <- NA
  # Remove the first column
  matched_ivs <- matched_ivs[, -1]
  # For each outcome
  # Use the MarkerName and reverseMarkerName columns to match the IVs
  for (k in 1:nrow(matched_ivs)) {
    if (matched_ivs$MarkerName[k] %in% combined_data$MarkerName) {
      # Select the row from combined_data where the MarkerName matches the IV
      out_values <- combined_data[combined_data$MarkerName == matched_ivs$MarkerName[k], ]
      matched_ivs$out_Chromosome[k] <- out_values$chrom
      matched_ivs$out_Position[k] <- out_values$chromStart
      matched_ivs$out_EffectAllele[k] <- out_values$Allele1
      matched_ivs$out_NonEffectAllele[k] <- out_values$Allele2
      matched_ivs$out_Beta[k] <- out_values$Effect
      matched_ivs$out_SE[k] <- out_values$StdErr
      matched_ivs$out_EAF[k] <- out_values$Freq1
      matched_ivs$out_Pval[k] <- out_values$Pvalue
      match_flag <- TRUE
      matched_FG_IVs <- matched_FG_IVs + 1
    } else if (matched_ivs$reverseMarkerName[k] %in% combined_data$MarkerName) {
      # Select the row from combined_data where the reverseMarkerName matches the IV
      out_values <- combined_data[combined_data$MarkerName == matched_ivs$reverseMarkerName[k], ]
      matched_ivs$out_Chromosome[k] <- out_values$chrom
      matched_ivs$out_Position[k] <- out_values$chromStart
      matched_ivs$out_EffectAllele[k] <- out_values$Allele1
      matched_ivs$out_NonEffectAllele[k] <- out_values$Allele2
      matched_ivs$out_Beta[k] <- out_values$Effect
      matched_ivs$out_SE[k] <- out_values$StdErr
      matched_ivs$out_EAF[k] <- out_values$Freq1
      matched_ivs$out_Pval[k] <- out_values$Pvalue
      match_flag <- TRUE
      matched_FG_IVs <- matched_FG_IVs + 1
    }
    # If the match_flag is still FALSE, remove the IV
    if (match_flag == FALSE) {
      FG_removed_IVs <- FG_removed_IVs + 1
    }
    match_flag <- FALSE
  }
  # Make the out_EffectAllele and out_NonEffectAllele columns capitalised
  matched_ivs$out_EffectAllele <- toupper(matched_ivs$out_EffectAllele)
  matched_ivs$out_NonEffectAllele <- toupper(matched_ivs$out_NonEffectAllele)
  # Remove rows with NA in the out_Beta column
  matched_ivs <- matched_ivs[!is.na(matched_ivs$out_Beta), ]
  # Add the number of IVs to the compid_IVs_list
  compid_IVs_list <- c(compid_IVs_list, nrow(matched_ivs))
  # Save the matched IVs to a new folder (make sure this folder exists or add code to create it)
  output_dir <- "Reverse_MR/Matched_Reverse_MR_FG_IVs"
  # Select the part after the last / in the file path
  output_file_name <- paste0(basename(file_names[i]), "_FG_matched.tsv")
  write.table(matched_ivs, file.path(output_dir, output_file_name), sep = "\t", row.names = FALSE)
  
  
  ## Do the same for hba1c_ivs ##
  # Initialise a match flag
  match_flag <- FALSE
  # Make an empty data frame to store the matched IVs
  matched_ivs <- as.data.frame(matrix(NA, nrow = nrow(clumped_hba1c_ivs)))
  # Initialise columns
  matched_ivs$Chromosome <- clumped_hba1c_ivs$Chr
  matched_ivs$Position <- clumped_hba1c_ivs$`Pos (bp)`
  matched_ivs$MarkerName <- clumped_hba1c_ivs$MarkerName
  matched_ivs$reverseMarkerName <- clumped_hba1c_ivs$reverseMarkerName
  matched_ivs$SNP <- clumped_hba1c_ivs$rsid
  matched_ivs$EffectAllele <- clumped_hba1c_ivs$`Effect Allele`
  matched_ivs$NonEffectAllele <- clumped_hba1c_ivs$`Other Allele`
  matched_ivs$EAF <- clumped_hba1c_ivs$EAF
  matched_ivs$Beta <- clumped_hba1c_ivs$Effect
  matched_ivs$SE <- clumped_hba1c_ivs$SE
  matched_ivs$Pval <- clumped_hba1c_ivs$pval
  matched_ivs$Fstat <- clumped_hba1c_ivs$Fstat
  matched_ivs$exposure <- "HbA1c"
  matched_ivs$outcome <- metabolite_compid
  matched_ivs$out_Chromosome <- NA
  matched_ivs$out_Position <- NA
  matched_ivs$out_EffectAllele <- NA
  matched_ivs$out_NonEffectAllele <- NA
  matched_ivs$out_Beta <- NA
  matched_ivs$out_SE <- NA
  matched_ivs$out_EAF <- NA
  matched_ivs$out_Pval <- NA
  # Remove the first column
  matched_ivs <- matched_ivs[, -1]
  # For each outcome
  # Use the MarkerName and reverseMarkerName columns to match the IVs
  for (k in 1:nrow(matched_ivs)) {
    if (matched_ivs$MarkerName[k] %in% combined_data$MarkerName) {
      # Select the row from combined_data where the MarkerName matches the IV
      out_values <- combined_data[combined_data$MarkerName == matched_ivs$MarkerName[k], ]
      matched_ivs$out_Chromosome[k] <- out_values$chrom
      matched_ivs$out_Position[k] <- out_values$chromStart
      matched_ivs$out_EffectAllele[k] <- out_values$Allele1
      matched_ivs$out_NonEffectAllele[k] <- out_values$Allele2
      matched_ivs$out_Beta[k] <- out_values$Effect
      matched_ivs$out_SE[k] <- out_values$StdErr
      matched_ivs$out_EAF[k] <- out_values$Freq1
      matched_ivs$out_Pval[k] <- out_values$Pvalue
      matched_HBA1C_IVs <- matched_HBA1C_IVs + 1
      match_flag <- TRUE
    } else if (matched_ivs$reverseMarkerName[k] %in% combined_data$MarkerName) {
      # Select the row from combined_data where the reverseMarkerName matches the IV
      out_values <- combined_data[combined_data$MarkerName == matched_ivs$reverseMarkerName[k], ]
      matched_ivs$out_Chromosome[k] <- out_values$chrom
      matched_ivs$out_Position[k] <- out_values$chromStart
      matched_ivs$out_EffectAllele[k] <- out_values$Allele1
      matched_ivs$out_NonEffectAllele[k] <- out_values$Allele2
      matched_ivs$out_Beta[k] <- out_values$Effect
      matched_ivs$out_SE[k] <- out_values$StdErr
      matched_ivs$out_EAF[k] <- out_values$Freq1
      matched_ivs$out_Pval[k] <- out_values$Pvalue
      matched_HBA1C_IVs <- matched_HBA1C_IVs + 1
      match_flag <- TRUE
    }
    # If the match_flag is still FALSE, remove the IV
    if (match_flag == FALSE) {
      HBA1C_removed_IVs <- HBA1C_removed_IVs + 1
    }
    match_flag <- FALSE
  }
  # Make the out_EffectAllele and out_NonEffectAllele columns capitalised
  matched_ivs$out_EffectAllele <- toupper(matched_ivs$out_EffectAllele)
  matched_ivs$out_NonEffectAllele <- toupper(matched_ivs$out_NonEffectAllele)
  # Remove rows with NA in the out_Beta column
  matched_ivs <- matched_ivs[!is.na(matched_ivs$out_Beta), ]
  # Add the number of IVs to the compid_IVs_list
  compid_IVs_list <- c(compid_IVs_list, nrow(matched_ivs))
  # Save the matched IVs to a new folder (make sure this folder exists or add code to create it)
  output_dir <- "Reverse_MR/Matched_Reverse_MR_HBA1C_IVs"
  # Select the part after the last / in the file path
  output_file_name <- paste0(basename(file_names[i]), "_HBA1C_matched.tsv")
  write.table(matched_ivs, file.path(output_dir, output_file_name), sep = "\t", row.names = FALSE)
  
  # Print the counters
  print(paste("Processed", metabolite_counter, "metabolites"))
  print(paste("Matched T2DM IVs:", matched_T2DM_IVs))
  print(paste("Matched T2DM IVs:", matched_FG_IVs))
  print(paste("Matched T2DM IVs:", matched_HBA1C_IVs))
  print(paste("T2DM Removed IVs:", T2DM_removed_IVs))
  print(paste("FG Removed IVs:", FG_removed_IVs))
  print(paste("HbA1c Removed IVs:", HBA1C_removed_IVs))
}