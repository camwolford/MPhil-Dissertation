### Run the PwCoCo Software ###
# This script was used to automate the running of the PwCoCo software for the colocalisation analysis.
library(tidyverse)
library(readr)

# Load in T2DM_Coloc_IVs.tsv
T2DM_Coloc_IVs <- read_tsv("Colocalisation/T2DM_Coloc_IVs.tsv")

# For metabolites with num_IVs > 0
for (metabolite in T2DM_Coloc_IVs$Metabolite){
  if (T2DM_Coloc_IVs$num_IVs[metabolite == T2DM_Coloc_IVs$Metabolite] > 0){
    harm_IVs <- read_tsv(paste0("Harmonised_T2DM_IVs/", metabolite, "T2DMT2DM_Harmonised_IVs.tsv"))
    # Remove rows with a non-NA value in the proxy column
    harm_IVs <- harm_IVs %>% filter(is.na(proxy))
    for (i in 1:nrow(harm_IVs)){
      SNP <- harm_IVs$SNP[i]
      Chromosome <- harm_IVs$Chromosome[i]
      
      # Construct the metabolite_filename
      metabolite_filename <- paste0(metabolite, "_", SNP, "_metabolite_snps_pwcoco.tsv")
      # Construct the t2dm_filename
      t2dm_filename <- paste0(metabolite, "_", SNP, "_t2dm_snps_pwcoco.tsv")
      
      # Use shQuote to handle filenames with special characters
      metabolite_filename <- shQuote(metabolite_filename)
      t2dm_filename <- shQuote(t2dm_filename)
      
      # Add the file paths to the filenames
      metabolite_filenamepath <- paste0("Colocalisation/PwCoCo_Ready_T2DM/", metabolite_filename)
      t2dm_filenamepath <- paste0("Colocalisation/PwCoCo_Ready_T2DM/", t2dm_filename)
      
      # Construct the command to run PwCoCo
      command <- paste(
        "Colocalisation/pwcoco/build/pwcoco",
        "--bfile 1kg.v3/interval_merged",
        "--sum_stats1", metabolite_filenamepath,
        "--sum_stats2", t2dm_filenamepath,
        "--out Colocalisation/pwcoco_out_t2dm",
        "--chr ", Chromosome,
        "--maf 0.01"
      )
      
      # Run the command
      system(command)
    }
  }
}


# Load in FG_Coloc_IVs.tsv
FG_Coloc_IVs <- read_tsv("Colocalisation/FG_Coloc_IVs.tsv")

# For metabolites with num_IVs > 0
for (metabolite in FG_Coloc_IVs$Metabolite){
  if (FG_Coloc_IVs$num_IVs[metabolite == FG_Coloc_IVs$Metabolite] > 0){
    harm_IVs <- read_tsv(paste0("Harmonised_FG_IVs/", metabolite, "FGFG_Harmonised_IVs.tsv"))
    # Remove rows with a non-NA value in the proxy column
    harm_IVs <- harm_IVs %>% filter(is.na(proxy))
    for (i in 1:nrow(harm_IVs)){
      SNP <- harm_IVs$SNP[i]
      Chromosome <- harm_IVs$Chromosome[i]
      
      # Construct the metabolite_filename
      metabolite_filename <- paste0(metabolite, "_", SNP, "_metabolite_snps_pwcoco.tsv")
      # Construct the fg_filename
      fg_filename <- paste0(metabolite, "_", SNP, "_fg_snps_pwcoco.tsv")
      
      # Use shQuote to handle filenames with special characters
      metabolite_filename <- shQuote(metabolite_filename)
      fg_filename <- shQuote(fg_filename)
      
      # Add the file paths to the filenames
      metabolite_filenamepath <- paste0("Colocalisation/PwCoCo_Ready_FG/", metabolite_filename)
      fg_filenamepath <- paste0("Colocalisation/PwCoCo_Ready_FG/", fg_filename)
      
      # Construct the command to run PwCoCo
      command <- paste(
        "Colocalisation/pwcoco/build/pwcoco",
        "--bfile 1kg.v3/interval_merged",
        "--sum_stats1", metabolite_filenamepath,
        "--sum_stats2", fg_filenamepath,
        "--out Colocalisation/pwcoco_out_fg",
        "--chr ", Chromosome,
        "--maf 0.01"
      )
      
      # Run the command
      system(command)
    }
  }
}


# Load in HBA1C_Coloc_IVs.tsv
HBA1C_Coloc_IVs <- read_tsv("Colocalisation/HBA1C_Coloc_IVs.tsv")

# For metabolites with num_IVs > 0
for (metabolite in HBA1C_Coloc_IVs$Metabolite){
  if (HBA1C_Coloc_IVs$num_IVs[metabolite == HBA1C_Coloc_IVs$Metabolite] > 0){
    harm_IVs <- read_tsv(paste0("Harmonised_HBA1C_IVs/", metabolite, "HBA1CHBA1C_Harmonised_IVs.tsv"))
    # Remove rows with a non-NA value in the proxy column
    harm_IVs <- harm_IVs %>% filter(is.na(proxy))
    for (i in 1:nrow(harm_IVs)){
      SNP <- harm_IVs$SNP[i]
      Chromosome <- harm_IVs$Chromosome[i]
      
      # Construct the metabolite_filename
      metabolite_filename <- paste0(metabolite, "_", SNP, "_metabolite_snps_pwcoco.tsv")
      # Construct the hba1c_filename
      hba1c_filename <- paste0(metabolite, "_", SNP, "_hba1c_snps_pwcoco.tsv")
      
      # Use shQuote to handle filenames with special characters
      metabolite_filename <- shQuote(metabolite_filename)
      hba1c_filename <- shQuote(hba1c_filename)
      
      # Add the file paths to the filenames
      metabolite_filenamepath <- paste0("Colocalisation/PwCoCo_Ready_HBA1C/", metabolite_filename)
      hba1c_filenamepath <- paste0("Colocalisation/PwCoCo_Ready_HBA1C/", hba1c_filename)
      
      # Construct the command to run PwCoCo
      command <- paste(
        "Colocalisation/pwcoco/build/pwcoco",
        "--bfile 1kg.v3/interval_merged",
        "--sum_stats1", metabolite_filenamepath,
        "--sum_stats2", hba1c_filenamepath,
        "--out Colocalisation/pwcoco_out_hba1c",
        "--chr ", Chromosome,
        "--maf 0.01"
      )
      
      # Run the command
      system(command)
    }
  }
}