### Select IVs from the outcomes ###
# This scriprt was used to select IVs from the outcomes.
library(tidyverse)
library(ieugwasr)
library(devtools)
devtools::install_github("explodecomputer/genetics.binaRies")
genetics.binaRies::get_plink_binary()

# Read in t2dm_sig_snps_ancestry.csv
t2dm_sig_snps_ancestry <- read_csv("Reverse_MR/t2dm_sig_snps_ancestry.csv", skip = 2)

# Rename columns
colnames(t2dm_sig_snps_ancestry)[6] <- "Residual"
colnames(t2dm_sig_snps_ancestry)[7] <- "Risk"
colnames(t2dm_sig_snps_ancestry)[8] <- "Other"
# Remove columns 9 to 18
t2dm_sig_snps_ancestry <- t2dm_sig_snps_ancestry[,-c(9:18)]
# Remove columns 14 to 28
t2dm_sig_snps_ancestry <- t2dm_sig_snps_ancestry[,-c(14:28)]
# Rename columns
colnames(t2dm_sig_snps_ancestry)[9] <- "Effective_sample_size"
colnames(t2dm_sig_snps_ancestry)[10] <- "EAF"
colnames(t2dm_sig_snps_ancestry)[11] <- "Log-OR"
colnames(t2dm_sig_snps_ancestry)[12] <- "SE"
colnames(t2dm_sig_snps_ancestry)[13] <- "Pval"
# Remove first row
t2dm_sig_snps_ancestry <- t2dm_sig_snps_ancestry[-1,]

# Values in the Chromosome column are NAs, if they are between two of the same values, fill them in
t2dm_sig_snps_ancestry$Chromosome <- zoo::na.locf(t2dm_sig_snps_ancestry$Chromosome)

### Manually edit the regions between Chromosomes
# All good 

t2dm_sig_snps_eur <- t2dm_sig_snps_ancestry
str(t2dm_sig_snps_eur)
# Convert the Pval column to numeric
t2dm_sig_snps_eur$Pval <- as.numeric(t2dm_sig_snps_eur$Pval) # Some NAs made - not huge
# Remove rows with NA in Pval
t2dm_sig_snps_eur <- t2dm_sig_snps_eur %>% filter(!is.na(Pval))
# Convert the EAF column to numeric
t2dm_sig_snps_eur$EAF <- as.numeric(t2dm_sig_snps_eur$EAF)

# Select IVs with Pval < 5e-8
t2dm_ivs <- t2dm_sig_snps_eur %>% filter(Pval < 5e-8)

# Select IVs with EAF > 0.01
t2dm_ivs <- t2dm_ivs %>% filter(EAF > 0.01)

# Remove the first column
t2dm_ivs <- t2dm_ivs[,-1]
# Check for the number of NA values in each column
colSums(is.na(t2dm_ivs))

colnames(t2dm_ivs)[2] = "rsid"
colnames(t2dm_ivs)[12] = "pval"
clumped_t2dm_ivs <- ieugwasr::ld_clump_local(t2dm_ivs,  clump_kb = 10000,
                         clump_r2 = 0.001,
                         clump_p = 1,
                         bfile = "1kg.v3/EUR",
                         plink_bin = genetics.binaRies::get_plink_binary())

# Save clumped_t2dm_ivs as a tsv file
write.table(clumped_t2dm_ivs, file = "Reverse_MR/clumped_t2dm_ivs.tsv", sep = "\t", row.names = FALSE)


# Read in fg_hba1c_sig_snps_ancestry.csv
fg_hba1c_sig_snps_ancestry <- read_csv("Reverse_MR/fg_hba1c_sig_snps_ancestry.csv")
# Check the unique values of the trait column
unique(fg_hba1c_sig_snps_ancestry$Trait)
# Select rows with Trait == "FG"
fg_sig_snps_ancestry <- fg_hba1c_sig_snps_ancestry %>% filter(Trait == "FG")
# Select rows with Trait == "HbA1c"
hba1c_sig_snps_ancestry <- fg_hba1c_sig_snps_ancestry %>% filter(Trait == "HbA1c")

str(fg_sig_snps_ancestry)
str(hba1c_sig_snps_ancestry)

# Convert the 'P-value (METAL)' column to numeric
fg_sig_snps_ancestry$`P-value (METAL)` <- as.numeric(fg_sig_snps_ancestry$`P-value (METAL)`) # Some values are "NA" not huge
hba1c_sig_snps_ancestry$`P-value (METAL)` <- as.numeric(hba1c_sig_snps_ancestry$`P-value (METAL)`) # Some values are "NA" not huge
# Remove rows with NA values in the 'P-value (METAL)' column
fg_sig_snps_ancestry <- fg_sig_snps_ancestry %>% filter(!is.na(`P-value (METAL)`))
hba1c_sig_snps_ancestry <- hba1c_sig_snps_ancestry %>% filter(!is.na(`P-value (METAL)`))
# Conver the EAF column to numeric
fg_sig_snps_ancestry$EAF <- as.numeric(fg_sig_snps_ancestry$EAF)
hba1c_sig_snps_ancestry$EAF <- as.numeric(hba1c_sig_snps_ancestry$EAF)

# Select IVs with Pval < 5e-8
fg_ivs <- fg_sig_snps_ancestry %>% filter(`P-value (METAL)` < 5e-8)
hba1c_ivs <- hba1c_sig_snps_ancestry %>% filter(`P-value (METAL)` < 5e-8)

# Slect IVs with EAF > 0.01
fg_ivs <- fg_ivs %>% filter(EAF > 0.01)
hba1c_ivs <- hba1c_ivs %>% filter(EAF > 0.01)


colnames(fg_ivs)[9] = "rsid"
colnames(fg_ivs)[18] = "pval"
clumped_fg_ivs <- ieugwasr::ld_clump_local(fg_ivs,  clump_kb = 10000,
                                             clump_r2 = 0.001,
                                             clump_p = 1,
                                             bfile = "1kg.v3/EUR",
                                             plink_bin = genetics.binaRies::get_plink_binary())

colnames(hba1c_ivs)[9] = "rsid"
colnames(hba1c_ivs)[18] = "pval"
clumped_hba1c_ivs <- ieugwasr::ld_clump_local(hba1c_ivs,  clump_kb = 10000,
                                           clump_r2 = 0.001,
                                           clump_p = 1,
                                           bfile = "1kg.v3/EUR",
                                           plink_bin = genetics.binaRies::get_plink_binary())

# Remove columns 1-5, 14, 19
clumped_fg_ivs <- clumped_fg_ivs[,-c(1:5, 14, 19)]
clumped_hba1c_ivs <- clumped_hba1c_ivs[,-c(1:5, 14, 19)]

# Save clumped_fg_ivs and clumped_hba1c_ivs as tsv files
write.table(clumped_fg_ivs, file = "Reverse_MR/clumped_fg_ivs.tsv", sep = "\t", row.names = FALSE)
write.table(clumped_hba1c_ivs, file = "Reverse_MR/clumped_hba1c_ivs.tsv", sep = "\t", row.names = FALSE)
