### Match IVs ###
# This script is used to match the clumped IVs for each metabolite by filtering by F statistics, matching variants and finding proxies.
library(tidyverse)

# Load in each metabolite data in /Clumped_IVs
metabolite_files <- list.files("Clumped_IVs", full.names = TRUE)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
remotes::install_github("mrcieu/gwasvcf")
gwasvcf::set_plink(genetics.binaRies::get_plink_binary())

# Load in the outcome GWAS data
t2dm_data <- readr::read_tsv("t2dm_gwas_cleaned.tsv")

# For each of the outcome GWAS data, make a new column called "MarkerName"
t2dm_data$MarkerName <- paste("chr", t2dm_data$Chromosome, "_", t2dm_data$Position, "_", t2dm_data$EffectAllele, "_", t2dm_data$NonEffectAllele, sep = "")

#Initialise a metabolite counter
metabolite_counter <- 0
# Initialise a counter for the number of metabolites that have weak IVs
weak_IVs <- 0
# Initialise a counter for the number of metabolites that have a match
matched_IVs <- 0
# Initialise a counter for the number of metabolites that need a proxy
need_proxy <- 0
# Initialise a counter for the number of metabolites that had an IV removed
removed_IVs <- 0
# Initialise a match found flag
match_found <- FALSE
# Initialise a proxy found flag
proxy_found <- FALSE
# Initialise a counter for the number of metabolites that have more than one IV
multiple_IVs <- 0
# Initialise a counter for the number of metabolites that have a one IV
one_IV <- 0
# Initialise a counter for the number of metabolites that have no remaining IVs
no_IVs <- 0

for (metabolite_file in metabolite_files) {
  # Load the metabolite data
  file_path <- metabolite_file
  metabolite_data <- readr::read_tsv(file_path)
  # Calculate the F statistics
  metabolite_data$Fstat <- (metabolite_data$Beta / metabolite_data$SE)^2
  # Check if any of the IVs have a low F statistic
  if (any(metabolite_data$Fstat < 10)) {
    weak_IVs <- weak_IVs + 1
  }
  # If the value for EffectAllele, NonEffectAllele, t2dm_EffectAllele or t2dm_NonEffectAllele is TRUE, change it to T as a character
  metabolite_data$EffectAllele <- ifelse(metabolite_data$EffectAllele == TRUE, "T", metabolite_data$EffectAllele)
  metabolite_data$NonEffectAllele <- ifelse(metabolite_data$NonEffectAllele == TRUE, "T", metabolite_data$NonEffectAllele)
  # Make a tempMarkerName column
  metabolite_data$tempMarkerName <- paste("chr", metabolite_data$Chromosome, "_", metabolite_data$Position, "_", metabolite_data$EffectAllele, "_", metabolite_data$NonEffectAllele, sep = "")
  # Make a reverseMarkerName column
  metabolite_data$reverseMarkerName <- paste("chr", metabolite_data$Chromosome, "_", metabolite_data$Position, "_", metabolite_data$NonEffectAllele, "_", metabolite_data$EffectAllele, sep = "")
  # Initialise the proxy and outcome data columns
  metabolite_data$proxy <- NA
  metabolite_data$proxy_EAF <- NA
  metabolite_data$proxy_R2 <- NA
  metabolite_data$proxy_EffectAllele <- NA
  metabolite_data$proxy_NonEffectAllele <- NA
  metabolite_data$t2dm_Chromosome <- NA
  metabolite_data$t2dm_Position <- NA
  metabolite_data$t2dm_EffectAllele <- NA
  metabolite_data$t2dm_NonEffectAllele <- NA
  metabolite_data$t2dm_Beta <- NA
  metabolite_data$t2dm_SE <- NA
  metabolite_data$t2dm_EAF <- NA
  metabolite_data$t2dm_Pval <- NA
  # For each IV, check if it is in the outcome GWAS data
  for (i in 1:nrow(metabolite_data)) {
    if (metabolite_data$tempMarkerName[i] %in% t2dm_data$MarkerName) {
      # Select the row from t2dm_data where the tempMarkerName matches the IV
      t2dm_values <- t2dm_data[t2dm_data$MarkerName == metabolite_data$tempMarkerName[i], ]
      metabolite_data$proxy[i] <- NA
      metabolite_data$proxy_EAF[i] <- NA
      metabolite_data$proxy_R2[i] <- NA
      metabolite_data$proxy_EffectAllele[i] <- NA
      metabolite_data$proxy_NonEffectAllele[i] <- NA
      metabolite_data$t2dm_Chromosome[i] <- t2dm_values$Chromosome
      metabolite_data$t2dm_Position[i] <- t2dm_values$Position
      metabolite_data$t2dm_EffectAllele[i] <- t2dm_values$EffectAllele
      metabolite_data$t2dm_NonEffectAllele[i] <- t2dm_values$NonEffectAllele
      metabolite_data$t2dm_Beta[i] <- t2dm_values$Beta
      metabolite_data$t2dm_SE[i] <- t2dm_values$SE
      metabolite_data$t2dm_EAF[i] <- t2dm_values$EAF
      metabolite_data$t2dm_Pval[i] <- t2dm_values$Pval
      match_found <- TRUE
      matched_IVs <- matched_IVs + 1
      # Print SNP found 
      print("Found SNP")
    } 
    if (metabolite_data$reverseMarkerName[i] %in% t2dm_data$MarkerName) {
      t2dm_values <- t2dm_data[t2dm_data$MarkerName == metabolite_data$reverseMarkerName[i], ]
      metabolite_data$proxy[i] <- NA
      metabolite_data$proxy_EAF[i] <- NA
      metabolite_data$proxy_R2[i] <- NA
      metabolite_data$proxy_EffectAllele[i] <- NA
      metabolite_data$proxy_NonEffectAllele[i] <- NA
      metabolite_data$t2dm_Chromosome[i] <- t2dm_values$Chromosome
      metabolite_data$t2dm_Position[i] <- t2dm_values$Position
      metabolite_data$t2dm_EffectAllele[i] <- t2dm_values$EffectAllele
      metabolite_data$t2dm_NonEffectAllele[i] <- t2dm_values$NonEffectAllele
      metabolite_data$t2dm_Beta[i] <- t2dm_values$Beta
      metabolite_data$t2dm_SE[i] <- t2dm_values$SE
      metabolite_data$t2dm_EAF[i] <- t2dm_values$EAF
      metabolite_data$t2dm_Pval[i] <- t2dm_values$Pval
      match_found <- TRUE
      matched_IVs <- matched_IVs + 1
      # Print SNP found 
      print("Found SNP")
    } 
    # If the iv is not in the outcome GWAS data, find a proxy
    if (match_found == FALSE) {
      # Find a proxy for the IV
      ld_proxies <- gwasvcf::get_ld_proxies(
        rsid = metabolite_data$SNP[i],
        bfile = "1kg.v3/EUR",
        searchspace = NULL,
        tag_kb = 5000,
        tag_nsnp = 5000,
        tag_r2 = 0.9,
        threads = 1,
        out = tempfile())
      # Remove proxies with ld_proxies$R < 0.9
      ld_proxies <- ld_proxies[ld_proxies$R >= 0.9, ]
      # Check if any of the proxies are in the outcome GWAS data using a for loop, break if a proxy is found
      for (j in 1:nrow(ld_proxies)) {
        if (any(t2dm_data$Chromosome == ld_proxies$CHR_B[j] & t2dm_data$Position == ld_proxies$BP_B[j] & t2dm_data$EffectAllele == ld_proxies$B1[j] & 
                t2dm_data$NonEffectAllele == ld_proxies$B2[j])) {
          # Select the row from t2dm_data where the tempMarkerName matches the proxy
          t2dm_values <- t2dm_data[t2dm_data$Chromosome == ld_proxies$CHR_B[j] & t2dm_data$Position == ld_proxies$BP_B[j] & t2dm_data$EffectAllele == ld_proxies$B1[j] & 
                                     t2dm_data$NonEffectAllele == ld_proxies$B2[j], ]
          metabolite_data$proxy[i] <- ld_proxies$SNP_B[j]
          metabolite_data$proxy_EAF[i] <- ld_proxies$MAF_B[j]
          metabolite_data$proxy_R2[i] <- ld_proxies$R[j]
          metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B1[j]
          metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B2[j]
          metabolite_data$t2dm_Chromosome[i] <- t2dm_values$Chromosome
          metabolite_data$t2dm_Position[i] <- t2dm_values$Position
          metabolite_data$t2dm_EffectAllele[i] <- t2dm_values$EffectAllele
          metabolite_data$t2dm_NonEffectAllele[i] <- t2dm_values$NonEffectAllele
          metabolite_data$t2dm_Beta[i] <- t2dm_values$Beta
          metabolite_data$t2dm_SE[i] <- t2dm_values$SE
          metabolite_data$t2dm_EAF[i] <- t2dm_values$EAF
          metabolite_data$t2dm_Pval[i] <- t2dm_values$Pval
          if (metabolite_data$EffectAllele[i] == ld_proxies$A1[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A2[j]) {
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
          if (metabolite_data$EffectAllele[i] == ld_proxies$A2[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A1[j]) {
            metabolite_data$proxy_EAF[i] <- metabolite_data$EAF[i]
            # Switch the proxy alleles
            metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B2[j]
            metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B1[j]
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
        }
        if (any(t2dm_data$Chromosome == ld_proxies$CHR_B[j] & t2dm_data$Position == ld_proxies$BP_B[j] &
                t2dm_data$EffectAllele == ld_proxies$B2[j] & t2dm_data$NonEffectAllele == ld_proxies$B1[j])) {
          t2dm_values <- t2dm_data[t2dm_data$Chromosome == ld_proxies$CHR_B[j] & t2dm_data$Position == ld_proxies$BP_B[j] & t2dm_data$EffectAllele == ld_proxies$B2[j] & 
                                     t2dm_data$NonEffectAllele == ld_proxies$B1[j], ]
          metabolite_data$proxy[i] <- ld_proxies$SNP_B[j]
          metabolite_data$proxy_EAF[i] <- ld_proxies$MAF_B[j]
          metabolite_data$proxy_R2[i] <- ld_proxies$R[j]
          metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B1[j]
          metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B2[j]
          metabolite_data$t2dm_Chromosome[i] <- t2dm_values$Chromosome
          metabolite_data$t2dm_Position[i] <- t2dm_values$Position
          metabolite_data$t2dm_EffectAllele[i] <- t2dm_values$EffectAllele
          metabolite_data$t2dm_NonEffectAllele[i] <- t2dm_values$NonEffectAllele
          metabolite_data$t2dm_Beta[i] <- t2dm_values$Beta
          metabolite_data$t2dm_SE[i] <- t2dm_values$SE
          metabolite_data$t2dm_EAF[i] <- t2dm_values$EAF
          metabolite_data$t2dm_Pval[i] <- t2dm_values$Pval
          if (metabolite_data$EffectAllele[i] == ld_proxies$A1[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A2[j]) {
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
          if (metabolite_data$EffectAllele[i] == ld_proxies$A2[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A1[j]) {
            metabolite_data$proxy_EAF[i] <- metabolite_data$EAF[i]
            # Switch the proxy alleles
            metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B2[j]
            metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B1[j]
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
        }
      }
    }
    # If the match_found and proxy_found are still FALSE, remove the IV
    if ((match_found == FALSE) & (proxy_found == FALSE)) {
        removed_IVs <- removed_IVs + 1
    }
    match_found <- FALSE
    proxy_found <- FALSE
  }
  # If a row is all NA, remove the row
  metabolite_data <- metabolite_data[rowSums(is.na(metabolite_data)) != ncol(metabolite_data), ]
  # Remove any rows with an empty metabolite_data$t2dm_Beta
  metabolite_data <- metabolite_data[!is.na(metabolite_data$t2dm_Beta), ]
  # Remove the tempMarkerName column
  metabolite_data$tempMarkerName <- NULL
  # Remove the reverseMarkerName column
  metabolite_data$reverseMarkerName <- NULL
  # Check the number of IVs left
  if (nrow(metabolite_data) > 1) {
    multiple_IVs <- multiple_IVs + 1
    # Save the metabolite data with the IV information
    # Extract the part of the file name before the "_Clumped_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Clumped_IVs")[[1]][1]
    # Add "T2DM_Finalised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "T2DM_Matched_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Matched_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 1) {
    one_IV <- one_IV + 1
    # Save the metabolite data with the IV information
    # Extract the part of the file name before the "_Clumped_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Clumped_IVs")[[1]][1]
    # Add "T2DM_Finalised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "T2DM_Matched_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Matched_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 0) {
    no_IVs <- no_IVs + 1
    metabolite_counter <- metabolite_counter + 1
    # Do not save the data
  }
  print(paste("Processed", metabolite_counter, "metabolites"))
  print(paste("Metabolites with multiple IVs:", multiple_IVs))
  print(paste("Metabolites with one IV:", one_IV))
  print(paste("Metabolites with no IVs:", no_IVs))
  print(paste("Metabolites with weak IVs:", weak_IVs))
  print(paste("Matched IVs:", matched_IVs))
  print(paste("Proxies:", need_proxy))
  print(paste("Removed IVs:", removed_IVs))
}


### Do the same for fasting glucose ###

fg_data <- readr::read_tsv("fasting_glucose_gwas_cleaned.tsv")
# Remove rows with NAs
fg_data <- fg_data[!is.na(fg_data$EAF), ]

fg_data$MarkerName <- paste("chr", fg_data$Chromosome, "_", fg_data$Position, "_", fg_data$EffectAllele, "_", fg_data$NonEffectAllele, sep = "")
fg_data$MarkerName[1]

# Are there any repeats of the MarkerNames?
sum(duplicated(fg_data$MarkerName))

### Perform the IV matching on all metabolites ###
# Load in each metabolite data in /Clumped_IVs
metabolite_files <- list.files("Clumped_IVs", full.names = TRUE)

gwasvcf::set_plink(genetics.binaRies::get_plink_binary())

#Initialise a metabolite counter
metabolite_counter <- 0
# Initialise a counter for the number of metabolites that have weak IVs
weak_IVs <- 0
# Initialise a counter for the number of metabolites that have a match
matched_IVs <- 0
# Initialise a counter for the number of metabolites that need a proxy
need_proxy <- 0
# Initialise a counter for the number of metabolites that had an IV removed
removed_IVs <- 0
# Initialise a match found flag
match_found <- FALSE
# Initialise a proxy found flag
proxy_found <- FALSE
# Initialise a counter for the number of metabolites that have more than one IV
multiple_IVs <- 0
# Initialise a counter for the number of metabolites that have a one IV
one_IV <- 0
# Initialise a counter for the number of metabolites that have no remaining IVs
no_IVs <- 0

for (metabolite_file in metabolite_files) {
  # Load the metabolite data
  file_path <- metabolite_file
  metabolite_data <- readr::read_tsv(file_path)
  # Calculate the F statistics
  metabolite_data$Fstat <- (metabolite_data$Beta / metabolite_data$SE)^2
  # Check if any of the IVs have a low F statistic
  if (any(metabolite_data$Fstat < 10)) {
    weak_IVs <- weak_IVs + 1
  }
  # If the value for EffectAllele, NonEffectAllele, fg_EffectAllele or fg_NonEffectAllele is TRUE, change it to T as a character
  metabolite_data$EffectAllele <- ifelse(metabolite_data$EffectAllele == TRUE, "T", metabolite_data$EffectAllele)
  metabolite_data$NonEffectAllele <- ifelse(metabolite_data$NonEffectAllele == TRUE, "T", metabolite_data$NonEffectAllele)
  # Make a tempMarkerName column
  metabolite_data$tempMarkerName <- paste("chr", metabolite_data$Chromosome, "_", metabolite_data$Position, "_", metabolite_data$EffectAllele, "_", metabolite_data$NonEffectAllele, sep = "")
  # Make a reverseMarkerName column
  metabolite_data$reverseMarkerName <- paste("chr", metabolite_data$Chromosome, "_", metabolite_data$Position, "_", metabolite_data$NonEffectAllele, "_", metabolite_data$EffectAllele, sep = "")
  # Initialise the proxy and outcome data columns
  metabolite_data$proxy <- NA
  metabolite_data$proxy_EAF <- NA
  metabolite_data$proxy_R2 <- NA
  metabolite_data$proxy_EffectAllele <- NA
  metabolite_data$proxy_NonEffectAllele <- NA
  metabolite_data$fg_Chromosome <- NA
  metabolite_data$fg_Position <- NA
  metabolite_data$fg_EffectAllele <- NA
  metabolite_data$fg_NonEffectAllele <- NA
  metabolite_data$fg_Beta <- NA
  metabolite_data$fg_SE <- NA
  metabolite_data$fg_EAF <- NA
  metabolite_data$fg_Pval <- NA
  # For each IV, check if it is in the outcome GWAS data
  for (i in 1:nrow(metabolite_data)) {
    if (metabolite_data$tempMarkerName[i] %in% fg_data$MarkerName) {
      # Select the row from fg_data where the tempMarkerName matches the IV
      fg_values <- fg_data[fg_data$MarkerName == metabolite_data$tempMarkerName[i], ]
      metabolite_data$proxy[i] <- NA
      metabolite_data$proxy_EAF[i] <- NA
      metabolite_data$proxy_R2[i] <- NA
      metabolite_data$proxy_EffectAllele[i] <- NA
      metabolite_data$proxy_NonEffectAllele[i] <- NA
      metabolite_data$fg_Chromosome[i] <- fg_values$Chromosome
      metabolite_data$fg_Position[i] <- fg_values$Position
      metabolite_data$fg_EffectAllele[i] <- fg_values$EffectAllele
      metabolite_data$fg_NonEffectAllele[i] <- fg_values$NonEffectAllele
      metabolite_data$fg_Beta[i] <- fg_values$Beta
      metabolite_data$fg_SE[i] <- fg_values$SE
      metabolite_data$fg_EAF[i] <- fg_values$EAF
      metabolite_data$fg_Pval[i] <- fg_values$Pval
      match_found <- TRUE
      matched_IVs <- matched_IVs + 1
      # Print SNP found 
      print("Found SNP")
    } 
    if (metabolite_data$reverseMarkerName[i] %in% fg_data$MarkerName) {
      fg_values <- fg_data[fg_data$MarkerName == metabolite_data$reverseMarkerName[i], ]
      metabolite_data$proxy[i] <- NA
      metabolite_data$proxy_EAF[i] <- NA
      metabolite_data$proxy_R2[i] <- NA
      metabolite_data$proxy_EffectAllele[i] <- NA
      metabolite_data$proxy_NonEffectAllele[i] <- NA
      metabolite_data$fg_Chromosome[i] <- fg_values$Chromosome
      metabolite_data$fg_Position[i] <- fg_values$Position
      metabolite_data$fg_EffectAllele[i] <- fg_values$EffectAllele
      metabolite_data$fg_NonEffectAllele[i] <- fg_values$NonEffectAllele
      metabolite_data$fg_Beta[i] <- fg_values$Beta
      metabolite_data$fg_SE[i] <- fg_values$SE
      metabolite_data$fg_EAF[i] <- fg_values$EAF
      metabolite_data$fg_Pval[i] <- fg_values$Pval
      match_found <- TRUE
      matched_IVs <- matched_IVs + 1
      # Print SNP found 
      print("Found SNP")
    } 
    # If the iv is not in the outcome GWAS data, find a proxy
    if (match_found == FALSE) {
      # Find a proxy for the IV
      ld_proxies <- gwasvcf::get_ld_proxies(
        rsid = metabolite_data$SNP[i],
        bfile = "1kg.v3/EUR",
        searchspace = NULL,
        tag_kb = 5000,
        tag_nsnp = 5000,
        tag_r2 = 0.9,
        threads = 1,
        out = tempfile())
      # Remove proxies with ld_proxies$R < 0.9
      ld_proxies <- ld_proxies[ld_proxies$R >= 0.9, ]
      # Check if any of the proxies are in the outcome GWAS data using a for loop, break if a proxy is found
      for (j in 1:nrow(ld_proxies)) {
        if (any(fg_data$Chromosome == ld_proxies$CHR_B[j] & fg_data$Position == ld_proxies$BP_B[j] & fg_data$EffectAllele == ld_proxies$B1[j] & 
                fg_data$NonEffectAllele == ld_proxies$B2[j])) {
          # Select the row from fg_data where the tempMarkerName matches the proxy
          fg_values <- fg_data[fg_data$Chromosome == ld_proxies$CHR_B[j] & fg_data$Position == ld_proxies$BP_B[j] & fg_data$EffectAllele == ld_proxies$B1[j] & 
                                     fg_data$NonEffectAllele == ld_proxies$B2[j], ]
          metabolite_data$proxy[i] <- ld_proxies$SNP_B[j]
          metabolite_data$proxy_EAF[i] <- ld_proxies$MAF_B[j]
          metabolite_data$proxy_R2[i] <- ld_proxies$R[j]
          metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B1[j]
          metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B2[j]
          metabolite_data$fg_Chromosome[i] <- fg_values$Chromosome
          metabolite_data$fg_Position[i] <- fg_values$Position
          metabolite_data$fg_EffectAllele[i] <- fg_values$EffectAllele
          metabolite_data$fg_NonEffectAllele[i] <- fg_values$NonEffectAllele
          metabolite_data$fg_Beta[i] <- fg_values$Beta
          metabolite_data$fg_SE[i] <- fg_values$SE
          metabolite_data$fg_EAF[i] <- fg_values$EAF
          metabolite_data$fg_Pval[i] <- fg_values$Pval
          need_proxy <- need_proxy + 1
          if (metabolite_data$EffectAllele[i] == ld_proxies$A1[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A2[j]) {
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
          if (metabolite_data$EffectAllele[i] == ld_proxies$A2[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A1[j]) {
            metabolite_data$proxy_EAF[i] <- metabolite_data$EAF[i]
            # Switch the proxy alleles
            metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B2[j]
            metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B1[j]
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
        }
        if (any(fg_data$Chromosome == ld_proxies$CHR_B[j] & fg_data$Position == ld_proxies$BP_B[j] &
                fg_data$EffectAllele == ld_proxies$B2[j] & fg_data$NonEffectAllele == ld_proxies$B1[j])) {
          fg_values <- fg_data[fg_data$Chromosome == ld_proxies$CHR_B[j] & fg_data$Position == ld_proxies$BP_B[j] & fg_data$EffectAllele == ld_proxies$B2[j] & 
                                     fg_data$NonEffectAllele == ld_proxies$B1[j], ]
          metabolite_data$proxy[i] <- ld_proxies$SNP_B[j]
          metabolite_data$proxy_EAF[i] <- ld_proxies$MAF_B[j]
          metabolite_data$proxy_R2[i] <- ld_proxies$R[j]
          metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B1[j]
          metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B2[j]
          metabolite_data$fg_Chromosome[i] <- fg_values$Chromosome
          metabolite_data$fg_Position[i] <- fg_values$Position
          metabolite_data$fg_EffectAllele[i] <- fg_values$EffectAllele
          metabolite_data$fg_NonEffectAllele[i] <- fg_values$NonEffectAllele
          metabolite_data$fg_Beta[i] <- fg_values$Beta
          metabolite_data$fg_SE[i] <- fg_values$SE
          metabolite_data$fg_EAF[i] <- fg_values$EAF
          metabolite_data$fg_Pval[i] <- fg_values$Pval
          need_proxy <- need_proxy + 1
          if (metabolite_data$EffectAllele[i] == ld_proxies$A1[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A2[j]) {
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
          if (metabolite_data$EffectAllele[i] == ld_proxies$A2[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A1[j]) {
            metabolite_data$proxy_EAF[i] <- metabolite_data$EAF[i]
            # Switch the proxy alleles
            metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B2[j]
            metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B1[j]
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
        }
      }
    }
    # If the match_found and proxy_found are still FALSE, remove the IV
    if ((match_found == FALSE) & (proxy_found == FALSE)) {
      removed_IVs <- removed_IVs + 1
    }
    match_found <- FALSE
    proxy_found <- FALSE
  }
  # If a row is all NA, remove the row
  metabolite_data <- metabolite_data[rowSums(is.na(metabolite_data)) != ncol(metabolite_data), ]
  # Remove any rows with an empty metabolite_data$fg_Beta
  metabolite_data <- metabolite_data[!is.na(metabolite_data$fg_Beta), ]
  # Remove the tempMarkerName column
  metabolite_data$tempMarkerName <- NULL
  # Remove the reverseMarkerName column
  metabolite_data$reverseMarkerName <- NULL
  # Check the number of IVs left
  if (nrow(metabolite_data) > 1) {
    multiple_IVs <- multiple_IVs + 1
    # Save the metabolite data with the IV information
    # Extract the part of the file name before the "_Clumped_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Clumped_IVs")[[1]][1]
    # Add "FG_Finalised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "FG_Matched_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Matched_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 1) {
    one_IV <- one_IV + 1
    # Save the metabolite data with the IV information
    # Extract the part of the file name before the "_Clumped_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Clumped_IVs")[[1]][1]
    # Add "FG_Finalised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "FG_Matched_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Matched_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 0) {
    no_IVs <- no_IVs + 1
    metabolite_counter <- metabolite_counter + 1
    # Do not save the data
  }
  print(paste("Processed", metabolite_counter, "metabolites"))
  print(paste("Metabolites with multiple IVs:", multiple_IVs))
  print(paste("Metabolites with one IV:", one_IV))
  print(paste("Metabolites with no IVs:", no_IVs))
  print(paste("Metabolites with weak IVs:", weak_IVs))
  print(paste("Matched IVs:", matched_IVs))
  print(paste("Proxies:", need_proxy))
  print(paste("Removed IVs:", removed_IVs))
}


### Do the same for HbA1c ###

hba1c_data <- readr::read_tsv("hba1c_gwas_cleaned.tsv")
# Remove rows with NAs
hba1c_data <- hba1c_data[!is.na(hba1c_data$EAF), ]

hba1c_data$MarkerName <- paste("chr", hba1c_data$Chromosome, "_", hba1c_data$Position, "_", hba1c_data$EffectAllele, "_", hba1c_data$NonEffectAllele, sep = "")
hba1c_data$MarkerName[1]

# Are there any repeats of the MarkerNames?
sum(duplicated(hba1c_data$MarkerName))

### Perform the IV matching on all metabolites ###
# Load in each metabolite data in /Clumped_IVs
metabolite_files <- list.files("Clumped_IVs", full.names = TRUE)

gwasvcf::set_plink(genetics.binaRies::get_plink_binary())

#Initialise a metabolite counter
metabolite_counter <- 0
# Initialise a counter for the number of metabolites that have weak IVs
weak_IVs <- 0
# Initialise a counter for the number of metabolites that have a match
matched_IVs <- 0
# Initialise a counter for the number of metabolites that need a proxy
need_proxy <- 0
# Initialise a counter for the number of metabolites that had an IV removed
removed_IVs <- 0
# Initialise a match found flag
match_found <- FALSE
# Initialise a proxy found flag
proxy_found <- FALSE
# Initialise a counter for the number of metabolites that have more than one IV
multiple_IVs <- 0
# Initialise a counter for the number of metabolites that have a one IV
one_IV <- 0
# Initialise a counter for the number of metabolites that have no remaining IVs
no_IVs <- 0

for (metabolite_file in metabolite_files) {
  # Load the metabolite data
  file_path <- metabolite_file
  metabolite_data <- readr::read_tsv(file_path)
  # Calculate the F statistics
  metabolite_data$Fstat <- (metabolite_data$Beta / metabolite_data$SE)^2
  # Check if any of the IVs have a low F statistic
  if (any(metabolite_data$Fstat < 10)) {
    weak_IVs <- weak_IVs + 1
  }
  # If the value for EffectAllele, NonEffectAllele, hbA1c_EffectAllele or hbA1c_NonEffectAllele is TRUE, change it to T as a character
  metabolite_data$EffectAllele <- ifelse(metabolite_data$EffectAllele == TRUE, "T", metabolite_data$EffectAllele)
  metabolite_data$NonEffectAllele <- ifelse(metabolite_data$NonEffectAllele == TRUE, "T", metabolite_data$NonEffectAllele)
  # Make a tempMarkerName column
  metabolite_data$tempMarkerName <- paste("chr", metabolite_data$Chromosome, "_", metabolite_data$Position, "_", metabolite_data$EffectAllele, "_", metabolite_data$NonEffectAllele, sep = "")
  # Make a reverseMarkerName column
  metabolite_data$reverseMarkerName <- paste("chr", metabolite_data$Chromosome, "_", metabolite_data$Position, "_", metabolite_data$NonEffectAllele, "_", metabolite_data$EffectAllele, sep = "")
  # Initialise the proxy and outcome data columns
  metabolite_data$proxy <- NA
  metabolite_data$proxy_EAF <- NA
  metabolite_data$proxy_R2 <- NA
  metabolite_data$proxy_EffectAllele <- NA
  metabolite_data$proxy_NonEffectAllele <- NA
  metabolite_data$hba1c_Chromosome <- NA
  metabolite_data$hba1c_Position <- NA
  metabolite_data$hba1c_EffectAllele <- NA
  metabolite_data$hba1c_NonEffectAllele <- NA
  metabolite_data$hba1c_Beta <- NA
  metabolite_data$hba1c_SE <- NA
  metabolite_data$hba1c_EAF <- NA
  metabolite_data$hba1c_Pval <- NA
  # For each IV, check if it is in the outcome GWAS data
  for (i in 1:nrow(metabolite_data)) {
    if (metabolite_data$tempMarkerName[i] %in% hba1c_data$MarkerName) {
      # Select the row from hba1c_data where the tempMarkerName matches the IV
      hba1c_values <- hba1c_data[hba1c_data$MarkerName == metabolite_data$tempMarkerName[i], ]
      metabolite_data$proxy[i] <- NA
      metabolite_data$proxy_EAF[i] <- NA
      metabolite_data$proxy_R2[i] <- NA
      metabolite_data$proxy_EffectAllele[i] <- NA
      metabolite_data$proxy_NonEffectAllele[i] <- NA
      metabolite_data$hba1c_Chromosome[i] <- hba1c_values$Chromosome
      metabolite_data$hba1c_Position[i] <- hba1c_values$Position
      metabolite_data$hba1c_EffectAllele[i] <- hba1c_values$EffectAllele
      metabolite_data$hba1c_NonEffectAllele[i] <- hba1c_values$NonEffectAllele
      metabolite_data$hba1c_Beta[i] <- hba1c_values$Beta
      metabolite_data$hba1c_SE[i] <- hba1c_values$SE
      metabolite_data$hba1c_EAF[i] <- hba1c_values$EAF
      metabolite_data$hba1c_Pval[i] <- hba1c_values$Pval
      match_found <- TRUE
      matched_IVs <- matched_IVs + 1
      # Print SNP found 
      print("Found SNP")
    } 
    if (metabolite_data$reverseMarkerName[i] %in% hba1c_data$MarkerName) {
      hba1c_values <- hba1c_data[hba1c_data$MarkerName == metabolite_data$reverseMarkerName[i], ]
      metabolite_data$proxy[i] <- NA
      metabolite_data$proxy_EAF[i] <- NA
      metabolite_data$proxy_R2[i] <- NA
      metabolite_data$proxy_EffectAllele[i] <- NA
      metabolite_data$proxy_NonEffectAllele[i] <- NA
      metabolite_data$hba1c_Chromosome[i] <- hba1c_values$Chromosome
      metabolite_data$hba1c_Position[i] <- hba1c_values$Position
      metabolite_data$hba1c_EffectAllele[i] <- hba1c_values$EffectAllele
      metabolite_data$hba1c_NonEffectAllele[i] <- hba1c_values$NonEffectAllele
      metabolite_data$hba1c_Beta[i] <- hba1c_values$Beta
      metabolite_data$hba1c_SE[i] <- hba1c_values$SE
      metabolite_data$hba1c_EAF[i] <- hba1c_values$EAF
      metabolite_data$hba1c_Pval[i] <- hba1c_values$Pval
      match_found <- TRUE
      matched_IVs <- matched_IVs + 1
      # Print SNP found 
      print("Found SNP")
    } 
    # If the iv is not in the outcome GWAS data, find a proxy
    if (match_found == FALSE) {
      # Find a proxy for the IV
      ld_proxies <- gwasvcf::get_ld_proxies(
        rsid = metabolite_data$SNP[i],
        bfile = "1kg.v3/EUR",
        searchspace = NULL,
        tag_kb = 5000,
        tag_nsnp = 5000,
        tag_r2 = 0.9,
        threads = 1,
        out = tempfile())
      # Remove proxies with ld_proxies$R < 0.9
      ld_proxies <- ld_proxies[ld_proxies$R >= 0.9, ]
      # Check if any of the proxies are in the outcome GWAS data using a for loop, break if a proxy is found
      for (j in 1:nrow(ld_proxies)) {
        if (any(hba1c_data$Chromosome == ld_proxies$CHR_B[j] & hba1c_data$Position == ld_proxies$BP_B[j] & hba1c_data$EffectAllele == ld_proxies$B1[j] & 
                hba1c_data$NonEffectAllele == ld_proxies$B2[j])) {
          # Select the row from hba1c_data where the tempMarkerName matches the proxy
          hba1c_values <- hba1c_data[hba1c_data$Chromosome == ld_proxies$CHR_B[j] & hba1c_data$Position == ld_proxies$BP_B[j] & hba1c_data$EffectAllele == ld_proxies$B1[j] & 
                                       hba1c_data$NonEffectAllele == ld_proxies$B2[j], ]
          metabolite_data$proxy[i] <- ld_proxies$SNP_B[j]
          metabolite_data$proxy_EAF[i] <- ld_proxies$MAF_B[j]
          metabolite_data$proxy_R2[i] <- ld_proxies$R[j]
          metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B1[j]
          metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B2[j]
          metabolite_data$hba1c_Chromosome[i] <- hba1c_values$Chromosome
          metabolite_data$hba1c_Position[i] <- hba1c_values$Position
          metabolite_data$hba1c_EffectAllele[i] <- hba1c_values$EffectAllele
          metabolite_data$hba1c_NonEffectAllele[i] <- hba1c_values$NonEffectAllele
          metabolite_data$hba1c_Beta[i] <- hba1c_values$Beta
          metabolite_data$hba1c_SE[i] <- hba1c_values$SE
          metabolite_data$hba1c_EAF[i] <- hba1c_values$EAF
          metabolite_data$hba1c_Pval[i] <- hba1c_values$Pval
          if (metabolite_data$EffectAllele[i] == ld_proxies$A1[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A2[j]) {
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
          if (metabolite_data$EffectAllele[i] == ld_proxies$A2[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A1[j]) {
            metabolite_data$proxy_EAF[i] <- metabolite_data$EAF[i]
            # Switch the proxy alleles
            metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B2[j]
            metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B1[j]
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
        }
        if (any(hba1c_data$Chromosome == ld_proxies$CHR_B[j] & hba1c_data$Position == ld_proxies$BP_B[j] &
                hba1c_data$EffectAllele == ld_proxies$B2[j] & hba1c_data$NonEffectAllele == ld_proxies$B1[j])) {
          hba1c_values <- hba1c_data[hba1c_data$Chromosome == ld_proxies$CHR_B[j] & hba1c_data$Position == ld_proxies$BP_B[j] & hba1c_data$EffectAllele == ld_proxies$B2[j] & 
                                       hba1c_data$NonEffectAllele == ld_proxies$B1[j], ]
          metabolite_data$proxy[i] <- ld_proxies$SNP_B[j]
          metabolite_data$proxy_EAF[i] <- ld_proxies$MAF_B[j]
          metabolite_data$proxy_R2[i] <- ld_proxies$R[j]
          metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B1[j]
          metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B2[j]
          metabolite_data$hba1c_Chromosome[i] <- hba1c_values$Chromosome
          metabolite_data$hba1c_Position[i] <- hba1c_values$Position
          metabolite_data$hba1c_EffectAllele[i] <- hba1c_values$EffectAllele
          metabolite_data$hba1c_NonEffectAllele[i] <- hba1c_values$NonEffectAllele
          metabolite_data$hba1c_Beta[i] <- hba1c_values$Beta
          metabolite_data$hba1c_SE[i] <- hba1c_values$SE
          metabolite_data$hba1c_EAF[i] <- hba1c_values$EAF
          metabolite_data$hba1c_Pval[i] <- hba1c_values$Pval
          if (metabolite_data$EffectAllele[i] == ld_proxies$A1[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A2[j]) {
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
          if (metabolite_data$EffectAllele[i] == ld_proxies$A2[j] & metabolite_data$NonEffectAllele[i] == ld_proxies$A1[j]) {
            metabolite_data$proxy_EAF[i] <- metabolite_data$EAF[i]
            # Switch the proxy alleles
            metabolite_data$proxy_EffectAllele[i] <- ld_proxies$B2[j]
            metabolite_data$proxy_NonEffectAllele[i] <- ld_proxies$B1[j]
            need_proxy <- need_proxy + 1
            print("Proxy SNP")
            proxy_found <- TRUE
            break
          }
        }
      }
    }
    # If the match_found and proxy_found are still FALSE, remove the IV
    if ((match_found == FALSE) & (proxy_found == FALSE)) {
      removed_IVs <- removed_IVs + 1
    }
    match_found <- FALSE
    proxy_found <- FALSE
  }
  # If a row is all NA, remove the row
  metabolite_data <- metabolite_data[rowSums(is.na(metabolite_data)) != ncol(metabolite_data), ]
  # Remove any rows with an empty metabolite_data$hba1c_Beta
  metabolite_data <- metabolite_data[!is.na(metabolite_data$hba1c_Beta), ]
  # Remove the tempMarkerName column
  metabolite_data$tempMarkerName <- NULL
  # Remove the reverseMarkerName column
  metabolite_data$reverseMarkerName <- NULL
  # Check the number of IVs left
  if (nrow(metabolite_data) > 1) {
    multiple_IVs <- multiple_IVs + 1
    # Save the metabolite data with the IV information
    # Extract the part of the file name before the "_Clumped_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Clumped_IVs")[[1]][1]
    # Add "HBA1C_Finalised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "HBA1C_Matched_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Matched_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 1) {
    one_IV <- one_IV + 1
    # Save the metabolite data with the IV information
    # Extract the part of the file name before the "_Clumped_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_Clumped_IVs")[[1]][1]
    # Add "HBA1C_Finalised_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "HBA1C_Matched_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Matched_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  if (nrow(metabolite_data) == 0) {
    no_IVs <- no_IVs + 1
    metabolite_counter <- metabolite_counter + 1
    # Do not save the data
  }
  print(paste("Processed", metabolite_counter, "metabolites"))
  print(paste("Metabolites with multiple IVs:", multiple_IVs))
  print(paste("Metabolites with one IV:", one_IV))
  print(paste("Metabolites with no IVs:", no_IVs))
  print(paste("Metabolites with weak IVs:", weak_IVs))
  print(paste("Matched IVs:", matched_IVs))
  print(paste("Proxies:", need_proxy))
  print(paste("Removed IVs:", removed_IVs))
}
