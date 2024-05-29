library(tidyverse)
devtools::install_github("myles-lewis/locuszoomr")
library(locuszoomr)
library(cowplot)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("AnnotationHub")
BiocManager::install("AnnotationFilter", force = TRUE)
library(ensembldb)
library(AnnotationFilter)
library(AnnotationHub)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v76
edb

### Example Implementation ###
# Load Colocalisation/t2dm_included/pwcoco_out_t2dm.1-(1-enyl-palmitoyl)-2-arachidonoyl-GPE (P-16_0_20_4)*_rs3115669_metabolite_snps_pwcoco.tsv.included
exposure_data_included <- read.delim("Colocalisation/t2dm_included/pwcoco_out_t2dm.1-(1-enyl-palmitoyl)-2-arachidonoyl-GPE (P-16_0_20_4)*_rs3115669_metabolite_snps_pwcoco.tsv.included")
# Load Colocalisation/t2dm_included/pwcoco_out_t2dm.1-(1-enyl-palmitoyl)-2-arachidonoyl-GPE (P-16_0_20_4)*_rs3115669_t2dm_snps_pwcoco.tsv.included
outcome_data_included <- read.delim("Colocalisation/t2dm_included/pwcoco_out_t2dm.1-(1-enyl-palmitoyl)-2-arachidonoyl-GPE (P-16_0_20_4)*_rs3115669_t2dm_snps_pwcoco.tsv.included")
# Load Colocalisation/Preped_Regions/T2DM_Regions/1-(1-enyl-palmitoyl)-2-arachidonoyl-GPE (P-16_0_20_4)*_rs3115669_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/T2DM_Regions/1-(1-enyl-palmitoyl)-2-arachidonoyl-GPE (P-16_0_20_4)*_rs3115669_metabolite_snps.tsv")

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data 
# after merging the two dataframes by MarkerName and SNP
exposure_data <- merge(exposure_data_raw, exposure_data_included, by.x = "MarkerName", by.y = "SNP")
exposure_data_plot <- exposure_data %>% dplyr::select(chrom, chromStart, MarkerName, Allele1, Allele2, B, SE, Pval, Freq)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(exposure_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data
# after merging the two dataframes by MarkerName and SNP
outcome_data <- merge(exposure_data_raw, outcome_data_included, by.x = "MarkerName", by.y = "SNP")
outcome_data_plot <- outcome_data %>% dplyr::select(chrom, chromStart, MarkerName, Allele1, Allele2, B, SE, Pval, Freq)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(outcome_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

loc_exposure <- locus(gene = "MSH5", data = exposure_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "MSH5", data = outcome_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_outcome)
loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("index"),
           label_x = c(4, -5), maxrows = 1)
locus_plot(loc_outcome, labels = c("index"),
           label_x = c(4, -5), maxrows = 1)


### cysteine sulfinic acid Implementation ###
# Load Colocalisation/PwCoCo_Ready_T2DM/cysteine sulfinic acid_rs2727271_metabolite_snps_pwcoco.tsv
exposure_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_T2DM/cysteine sulfinic acid_rs2727271_metabolite_snps_pwcoco.tsv")
# Load Colocalisation/PwCoCo_Ready_T2DM/cysteine sulfinic acid_rs2727271_t2dm_snps_pwcoco.tsv
outcome_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_T2DM/cysteine sulfinic acid_rs2727271_t2dm_snps_pwcoco.tsv")
# Load Colocalisation/Preped_Regions/T2DM_Regions/cysteine sulfinic acid_rs2727271_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/T2DM_Regions/cysteine sulfinic acid_rs2727271_metabolite_snps.tsv")
# Load Colocalisation/Preped_Regions/T2DM_Regions/cysteine sulfinic acid_rs2727271_t2dm_snps.tsv
outcome_data_raw <- read_tsv("Colocalisation/Preped_Regions/T2DM_Regions/cysteine sulfinic acid_rs2727271_t2dm_snps.tsv")
# Make a MarkerName column as chr{Chromosome}:{Position}:{EffectAllele}:{NonEffectAllele}
outcome_data_raw$MarkerName <- paste0("chr", outcome_data_raw$Chromosome, ":", outcome_data_raw$Position, ":", outcome_data_raw$EffectAllele, ":", outcome_data_raw$NonEffectAllele)

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data
# after merging the two dataframes by MarkerName and SNP
exposure_data <- merge(exposure_data_raw, exposure_data_included, by.x = "MarkerName", by.y = "SNP")
exposure_data_plot <- exposure_data %>% dplyr::select(chrom, chromStart, MarkerName, Allele1, Allele2, Effect, StdErr, Pvalue, Freq1)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(exposure_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data
# after merging the two dataframes by MarkerName and SNP
outcome_data <- merge(outcome_data_raw, outcome_data_included, by.x = "MarkerName", by.y = "SNP")
outcome_data_plot <- outcome_data %>% dplyr::select(Chromosome, Position, MarkerName, A1, A2, beta, se, p, A1_freq)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(outcome_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

# Only use the snps in both exposure_data_plot and outcome_data_plot
exposure_data_plot <- exposure_data_plot[exposure_data_plot$snp %in% outcome_data_plot$snp,]
outcome_data_plot <- outcome_data_plot[outcome_data_plot$snp %in% exposure_data_plot$snp,]

# In both exposure_data_plot and outcome_data_plot, change the value for snp:chr11:61603358:A:T to rs2727271
exposure_data_plot$snp[exposure_data_plot$snp == "chr11:61603358:A:T"] <- "rs2727271"
outcome_data_plot$snp[outcome_data_plot$snp == "chr11:61603358:A:T"] <- "rs2727271"

loc_exposure <- locus(gene = "FEN1", data = exposure_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "FEN1", data = outcome_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_outcome)
loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("rs2727271"),
           label_x = c(4, -5), maxrows = 3)
locus_plot(loc_outcome, labels = c("rs2727271"),
           label_x = c(4, -5), maxrows = 3)

### 5-acetylamino-6-amino-3-methyluracil_rs35246381 Implementation ###
# Load Colocalisation/PwCoCo_Ready_T2DM/5-acetylamino-6-amino-3-methyluracil_rs35246381_metabolite_snps_pwcoco.tsv
exposure_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_T2DM/5-acetylamino-6-amino-3-methyluracil_rs35246381_metabolite_snps_pwcoco.tsv")
# Load Colocalisation/PwCoCo_Ready_T2DM/5-acetylamino-6-amino-3-methyluracil_rs35246381_t2dm_snps_pwcoco.tsv
outcome_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_T2DM/5-acetylamino-6-amino-3-methyluracil_rs35246381_t2dm_snps_pwcoco.tsv")
# Load Colocalisation/Preped_Regions/T2DM_Regions/5-acetylamino-6-amino-3-methyluracil_rs35246381_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/T2DM_Regions/5-acetylamino-6-amino-3-methyluracil_rs35246381_metabolite_snps.tsv")
# Load Colocalisation/Preped_Regions/T2DM_Regions/5-acetylamino-6-amino-3-methyluracil_rs35246381_t2dm_snps.tsv
outcome_data_raw <- read_tsv("Colocalisation/Preped_Regions/T2DM_Regions/5-acetylamino-6-amino-3-methyluracil_rs35246381_t2dm_snps.tsv")
# Make a MarkerName column as chr{Chromosome}:{Position}:{EffectAllele}:{NonEffectAllele}
outcome_data_raw$MarkerName <- paste0("chr", outcome_data_raw$Chromosome, ":", outcome_data_raw$Position, ":", outcome_data_raw$EffectAllele, ":", outcome_data_raw$NonEffectAllele)

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data
# after merging the two dataframes by MarkerName and SNP
exposure_data <- merge(exposure_data_raw, exposure_data_included, by.x = "MarkerName", by.y = "SNP")
exposure_data_plot <- exposure_data %>% dplyr::select(chrom, chromStart, MarkerName, Allele1, Allele2, Effect, StdErr, Pvalue, Freq1)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(exposure_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data
# after merging the two dataframes by MarkerName and SNP
outcome_data <- merge(outcome_data_raw, outcome_data_included, by.x = "MarkerName", by.y = "SNP")
outcome_data_plot <- outcome_data %>% dplyr::select(Chromosome, Position, MarkerName, A1, A2, beta, se, p, A1_freq)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(outcome_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

# Only use the snps in both exposure_data_plot and outcome_data_plot
exposure_data_plot <- exposure_data_plot[exposure_data_plot$snp %in% outcome_data_plot$snp,]
outcome_data_plot <- outcome_data_plot[outcome_data_plot$snp %in% exposure_data_plot$snp,]

loc_exposure <- locus(gene = "NAT2", data = exposure_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "NAT2", data = outcome_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_outcome)

loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("index"),
           label_x = c(4, -5), maxrows = 3)
locus_plot(loc_outcome, labels = c("index", "chr8:18272466:A:G"),
           label_x = c(4, -5), maxrows = 3)        

### sphinganine_rs676457 Implementation ###
# Load Colocalisation/PwCoCo_Ready_HBA1c/sphinganine_rs676457_metabolite_snps_pwcoco.tsv
exposure_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_HBA1C/sphinganine_rs676457_metabolite_snps_pwcoco.tsv")
# Load Colocalisation/PwCoCo_Ready_HBA1C/sphinganine_rs676457_hba1c_snps_pwcoco.tsv
outcome_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_HBA1C/sphinganine_rs676457_hba1c_snps_pwcoco.tsv")
# Load Colocalisation/Preped_Regions/HBA1C_Regions/sphinganine_rs676457_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/HBA1C_Regions/sphinganine_rs676457_metabolite_snps.tsv")
# Load Colocalisation/Preped_Regions/HBA1C_Regions/sphinganine_rs676457_hba1c_snps.tsv
outcome_data_raw <- read_tsv("Colocalisation/Preped_Regions/HBA1C_Regions/sphinganine_rs676457_hba1c_snps.tsv")
# Make a MarkerName column as chr{Chromosome}:{Position}:{EffectAllele}:{NonEffectAllele}
outcome_data_raw$MarkerName <- paste0("chr", outcome_data_raw$Chromosome, ":", outcome_data_raw$Position, ":", outcome_data_raw$EffectAllele, ":", outcome_data_raw$NonEffectAllele)

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data
# after merging the two dataframes by MarkerName and SNP
exposure_data <- merge(exposure_data_raw, exposure_data_included, by.x = "MarkerName", by.y = "SNP")
exposure_data_plot <- exposure_data %>% dplyr::select(chrom, chromStart, MarkerName, Allele1, Allele2, Effect, StdErr, Pvalue, Freq1)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(exposure_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data
# after merging the two dataframes by MarkerName and SNP
outcome_data <- merge(outcome_data_raw, outcome_data_included, by.x = "MarkerName", by.y = "SNP")
outcome_data_plot <- outcome_data %>% dplyr::select(Chromosome, Position, MarkerName, A1, A2, beta, se, p, A1_freq)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(outcome_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

# Only use the snps in both exposure_data_plot and outcome_data_plot
exposure_data_plot <- exposure_data_plot[exposure_data_plot$snp %in% outcome_data_plot$snp,]
outcome_data_plot <- outcome_data_plot[outcome_data_plot$snp %in% exposure_data_plot$snp,]

loc_exposure <- locus(gene = "ABO", data = exposure_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "ABO", data = outcome_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_outcome)
loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("index", "chr9:136149399:A:G"),
           label_x = c(4, -5), maxrows = 3)
locus_plot(loc_outcome, labels = c("index", "chr9:136146227:A:T"),
           label_x = c(4, -5), maxrows = 3)



### 1-eicosapentaenoyl-GPC (20_5)*_rs174545 Implementation ###
# Load Colocalisation/PwCoCo_Ready_HBA1C/1-eicosapentaenoyl-GPC (20_5)*_rs174545_metabolite_snps_pwcoco.tsv
exposure_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_HBA1C/1-eicosapentaenoyl-GPC (20_5)*_rs174545_metabolite_snps_pwcoco.tsv")
# Load Colocalisation/PwCoCo_Ready_HBA1C/1-eicosapentaenoyl-GPC (20_5)*_rs174545_hba1c_snps_pwcoco.tsv
outcome_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_HBA1C/1-eicosapentaenoyl-GPC (20_5)*_rs174545_hba1c_snps_pwcoco.tsv")
# Load Colocalisation/Preped_Regions/HBA1C_Regions/1-eicosapentaenoyl-GPC (20_5)*_rs174545_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/HBA1C_Regions/1-eicosapentaenoyl-GPC (20_5)*_rs174545_metabolite_snps.tsv")
# Load Colocalisation/Preped_Regions/HBA1C_Regions/1-eicosapentaenoyl-GPC (20_5)*_rs174545_hba1c_snps.tsv
outcome_data_raw <- read_tsv("Colocalisation/Preped_Regions/HBA1C_Regions/1-eicosapentaenoyl-GPC (20_5)*_rs174545_hba1c_snps.tsv")
# Make a MarkerName column as chr{Chromosome}:{Position}:{EffectAllele}:{NonEffectAllele}
outcome_data_raw$MarkerName <- paste0("chr", outcome_data_raw$Chromosome, ":", outcome_data_raw$Position, ":", outcome_data_raw$EffectAllele, ":", outcome_data_raw$NonEffectAllele)

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data
# after merging the two dataframes by MarkerName and SNP
exposure_data <- merge(exposure_data_raw, exposure_data_included, by.x = "MarkerName", by.y = "SNP")
exposure_data_plot <- exposure_data %>% dplyr::select(chrom, chromStart, MarkerName, Allele1, Allele2, Effect, StdErr, Pvalue, Freq1)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(exposure_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

# Make a new dataframe including the chrom, chromStart, Allele1, Allele2 columns from the raw data and the SNP, B, SE, Pval, Freq columns from the included data
# after merging the two dataframes by MarkerName and SNP
outcome_data <- merge(outcome_data_raw, outcome_data_included, by.x = "MarkerName", by.y = "SNP")
outcome_data_plot <- outcome_data %>% dplyr::select(Chromosome, Position, MarkerName, A1, A2, beta, se, p, A1_freq)
# Rename the columns chrom, pos, snp, a1, a2, b, se, p, freq
colnames(outcome_data_plot) <- c("chrom", "pos", "snp", "a1", "a2", "b", "se", "p", "freq")

# Only use the snps in both exposure_data_plot and outcome_data_plot
exposure_data_plot <- exposure_data_plot[exposure_data_plot$snp %in% outcome_data_plot$snp,]
outcome_data_plot <- outcome_data_plot[outcome_data_plot$snp %in% exposure_data_plot$snp,]

loc_exposure <- locus(gene = "FADS2", data = exposure_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "FADS2", data = outcome_data_plot, fix_window = 1e6,
             ens_db = edb)
summary(loc_outcome)
loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("index", "chr11:61581656:A:G"),
           label_x = c(4, -5), maxrows = 3)
locus_plot(loc_outcome, labels = c("index", "chr11:61569306:C:G"),
           label_x = c(4, -5), maxrows = 3)





### 800 650 Crop Top ###
