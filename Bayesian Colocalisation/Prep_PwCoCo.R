### Prepare the wrangled data for the PwCoCo analysis ###
# This script was used to prepare the data for PwCoCo analysis by converting it to the necessary format.
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
      snp <- harm_IVs$SNP[i]
      # Load in metabolite SNPs data
      metabolite_snps <- read_tsv(paste0("Colocalisation/Preped_Regions/T2DM_Regions/", metabolite, "_", snp, "_metabolite_snps.tsv"))
      # Load in T2DM SNPs data
      t2dm_snps <- read_tsv(paste0("Colocalisation/Preped_Regions/T2DM_Regions/", metabolite, "_", snp, "_t2dm_snps.tsv"))

      # Remove any rows in metabolite_snps with "d" or "i" in the Allele1 column
      metabolite_snps <- metabolite_snps %>% filter(!grepl("d", Allele1) & !grepl("i", Allele1))

      # Remove columns 3, 8, 9, 10, 14, 15, 16, 17, 18, 20 from metabolite_snps
      metabolite_snps <- metabolite_snps %>% select(-c(3, 8, 9, 10, 14, 15, 16, 17, 18, 20))

      # Make the Allele1 and Allele2 columns in metabolite_snps capitalised
      metabolite_snps$Allele1 <- toupper(metabolite_snps$Allele1)
      metabolite_snps$Allele2 <- toupper(metabolite_snps$Allele2)

      # Make a MarkerName column in the t2dm_snps data
      t2dm_snps$MarkerName <- paste0("chr", t2dm_snps$Chromosome, ":", t2dm_snps$Position, ":", t2dm_snps$EffectAllele, ":", t2dm_snps$NonEffectAllele)

      # Select the following columns from metabolite_snps: MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, Pvalue, N
      metabolite_snps <- metabolite_snps %>% select(MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, Pvalue, N)
      # Select the following columns from t2dm_snps: MarkerName, EffectAllele, NonEffectAllele, EAF, Beta, SE, Pval, Ncontrols, Ncases
      t2dm_snps <- t2dm_snps %>% select(MarkerName, EffectAllele, NonEffectAllele, EAF, Beta, SE, Pval, Ncontrols, Ncases)

      # Make the Ncontrols column the sum of Ncases to the Ncontrols column
      t2dm_snps$Ncontrols <- t2dm_snps$Ncontrols + t2dm_snps$Ncases

      # Rename the columns in metabolite_snps to SNP, A1, A2, A1_freq, beta, se, p, n
      colnames(metabolite_snps) <- c("SNP", "A1", "A2", "A1_freq", "beta", "se", "p", "n")
      # Rename the columns in t2dm_snps to SNP, A1, A2, A1_freq, beta, se, p, n, ncase
      colnames(t2dm_snps) <- c("SNP", "A1", "A2", "A1_freq", "beta", "se", "p", "n", "ncase")

      # Save the prepared data
      write_tsv(metabolite_snps, paste0("Colocalisation/PwCoCo_Ready_T2DM/", metabolite, "_", snp, "_metabolite_snps_pwcoco.tsv"))
      write_tsv(t2dm_snps, paste0("Colocalisation/PwCoCo_Ready_T2DM/", metabolite, "_", snp, "_t2dm_snps_pwcoco.tsv"))
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
      snp <- harm_IVs$SNP[i]
      # Load in metabolite SNPs data
      metabolite_snps <- read_tsv(paste0("Colocalisation/Preped_Regions/FG_Regions/", metabolite, "_", snp, "_metabolite_snps.tsv"))
      # Load in FG SNPs data
      fg_snps <- read_tsv(paste0("Colocalisation/Preped_Regions/FG_Regions/", metabolite, "_", snp, "_fg_snps.tsv"))

      # Remove any rows in metabolite_snps with "d" or "i" in the Allele1 column
      metabolite_snps <- metabolite_snps %>% filter(!grepl("d", Allele1) & !grepl("i", Allele1))

      # Remove columns 3, 8, 9, 10, 14, 15, 16, 17, 18, 20 from metabolite_snps
      metabolite_snps <- metabolite_snps %>% select(-c(3, 8, 9, 10, 14, 15, 16, 17, 18, 20))

      # Make the Allele1 and Allele2 columns in metabolite_snps capitalised
      metabolite_snps$Allele1 <- toupper(metabolite_snps$Allele1)
      metabolite_snps$Allele2 <- toupper(metabolite_snps$Allele2)

      # Make a MarkerName column in the fg_snps data
      fg_snps$MarkerName <- paste0("chr", fg_snps$Chromosome, ":", fg_snps$Position, ":", fg_snps$EffectAllele, ":", fg_snps$NonEffectAllele)

      # Select the following columns from metabolite_snps: MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, Pvalue, N
      metabolite_snps <- metabolite_snps %>% select(MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, Pvalue, N)
      # Select the following columns from fg_snps: MarkerName, EffectAllele, NonEffectAllele, EAF, Beta, SE, Pval, N
      fg_snps <- fg_snps %>% select(MarkerName, EffectAllele, NonEffectAllele, EAF, Beta, SE, Pval, SampleSize)

      # Rename the columns in metabolite_snps to SNP, A1, A2, A1_freq, beta, se, p, n
      colnames(metabolite_snps) <- c("SNP", "A1", "A2", "A1_freq", "beta", "se", "p", "n")
      colnames(metabolite_snps_025) <- c("SNP", "A1", "A2", "A1_freq", "beta", "se", "p", "n")
      # Rename the columns in fg_snps to SNP, A1, A2, A1_freq, beta, se, p, n
      colnames(fg_snps) <- c("SNP", "A1", "A2", "A1_freq", "beta", "se", "p", "n")

      # Save the prepared data
      write_tsv(metabolite_snps, paste0("Colocalisation/PwCoCo_Ready_FG/", metabolite, "_", snp, "_metabolite_snps_pwcoco.tsv"))
      write_tsv(fg_snps, paste0("Colocalisation/PwCoCo_Ready_FG/", metabolite, "_", snp, "_fg_snps_pwcoco.tsv"))
      write_tsv(metabolite_snps_025, paste0("Colocalisation/PwCoCo_Ready_FG/", metabolite, "_", snp, "_metabolite_snps_025_pwcoco.tsv"))
      write_tsv(fg_snps_025, paste0("Colocalisation/PwCoCo_Ready_FG/", metabolite, "_", snp, "_fg_snps_025_pwcoco.tsv"))
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
      snp <- harm_IVs$SNP[i]
      # Load in metabolite SNPs data
      metabolite_snps <- read_tsv(paste0("Colocalisation/Preped_Regions/HBA1C_Regions/", metabolite, "_", snp, "_metabolite_snps.tsv"))
      # Load in HBA1C SNPs data
      hba1c_snps <- read_tsv(paste0("Colocalisation/Preped_Regions/HBA1C_Regions/", metabolite, "_", snp, "_hba1c_snps.tsv"))

      # Remove any rows in metabolite_snps with "d" or "i" in the Allele1 column
      metabolite_snps <- metabolite_snps %>% filter(!grepl("d", Allele1) & !grepl("i", Allele1))

      # Remove columns 3, 8, 9, 10, 14, 15, 16, 17, 18, 20 from metabolite_snps
      metabolite_snps <- metabolite_snps %>% select(-c(3, 8, 9, 10, 14, 15, 16, 17, 18, 20))

      # Make the Allele1 and Allele2 columns in metabolite_snps capitalised
      metabolite_snps$Allele1 <- toupper(metabolite_snps$Allele1)
      metabolite_snps$Allele2 <- toupper(metabolite_snps$Allele2)

      # Make a MarkerName column in the hba1c_snps data
      hba1c_snps$MarkerName <- paste0("chr", hba1c_snps$Chromosome, ":", hba1c_snps$Position, ":", hba1c_snps$EffectAllele, ":", hba1c_snps$NonEffectAllele)

      # Select the following columns from metabolite_snps: MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, Pvalue, N
      metabolite_snps <- metabolite_snps %>% select(MarkerName, Allele1, Allele2, Freq1, Effect, StdErr, Pvalue, N)
      # Select the following columns from hba1c_snps: MarkerName, EffectAllele, NonEffectAllele, EAF, Beta, SE, Pval, N
      hba1c_snps <- hba1c_snps %>% select(MarkerName, EffectAllele, NonEffectAllele, EAF, Beta, SE, Pval, SampleSize)

      # Rename the columns in metabolite_snps to SNP, A1, A2, A1_freq, beta, se, p, n
      colnames(metabolite_snps) <- c("SNP", "A1", "A2", "A1_freq", "beta", "se", "p", "n")
      # Rename the columns in fg_snps to SNP, A1, A2, A1_freq, beta, se, p, n
      colnames(hba1c_snps) <- c("SNP", "A1", "A2", "A1_freq", "beta", "se", "p", "n")

      # Save the files
      write_tsv(metabolite_snps, paste0("Colocalisation/PwCoCo_Ready_HBA1C/", metabolite, "_", snp, "_metabolite_snps_pwcoco.tsv"))
      write_tsv(hba1c_snps, paste0("Colocalisation/PwCoCo_Ready_HBA1C/", metabolite, "_", snp, "_hba1c_snps_pwcoco.tsv"))
    }
  }
}