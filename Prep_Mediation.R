library(tidyverse)
library(stringr)
library(R.utils)
library(ieugwasr)
library(devtools)

devtools::install_github("explodecomputer/genetics.binaRies")
genetics.binaRies::get_plink_binary()

# Gunzip the variants.tsv.bgz file
gunzip("Mediation/variants.tsv.bgz", "Mediation/ukbb_variants_data.tsv")

# Read the ukbb_variants_data.tsv file
ukbb_variants_data <- readr::read_tsv("Mediation/ukbb_variants_data.tsv")

# Gunzip the 30020_irnt.gwas.imputed_v3.both_sexes.tsv.bgz file
gunzip("Mediation/30020_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "Mediation/haem_data.tsv")

# Read the haem_data.tsv file
haem_data <- readr::read_tsv("Mediation/haem_data.tsv")

# Use the haem_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
haem_ivs <- haem_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))
# Remove columns 7-25
haem_ivs <- haem_ivs[, -c(7:25)]
# Add the columns from haem_data to haem_ivs by matching the variant column
haem_ivs <- haem_ivs %>%
  dplyr::left_join(haem_data, by = c("variant" = "variant"))

# Remove haem_data
rm(haem_data)

# Select the IV variants 
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
haem_ivs <- haem_ivs %>%
  dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
haem_ivs <- haem_ivs[, -c(11, 12, 15)]

clumped_haem_ivs <- ieugwasr::ld_clump_local(haem_ivs,  clump_kb = 10000,
                                             clump_r2 = 0.001,
                                             clump_p = 1,
                                             bfile = "1kg.v3/EUR",
                                             plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_haem_ivs data frame to a file
write_tsv(clumped_haem_ivs, "Mediation/clumped_haem_ivs.tsv")


# Gunzip the 21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz file
gunzip("Mediation/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "Mediation/bmi_data.tsv")

# Read the bmi_data.tsv file
bmi_data <- readr::read_tsv("Mediation/bmi_data.tsv")

# Use the bmi_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
bmi_ivs <- bmi_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))

# Remove columns 7-25
bmi_ivs <- bmi_ivs[, -c(7:25)]

# Add the columns from bmi_data to bmi_ivs by matching the variant column
bmi_ivs <- bmi_ivs %>%
  dplyr::left_join(bmi_data, by = c("variant" = "variant"))

# Remove bmi_data
rm(bmi_data)

# Select the IV variants
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
bmi_ivs <- bmi_ivs %>%
  dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
bmi_ivs <- bmi_ivs[, -c(11, 12, 15)]

clumped_bmi_ivs <- ieugwasr::ld_clump_local(bmi_ivs,  clump_kb = 10000,
                                            clump_r2 = 0.001,
                                            clump_p = 1,
                                            bfile = "1kg.v3/EUR",
                                            plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_bmi_ivs data frame to a file
write_tsv(clumped_bmi_ivs, "Mediation/clumped_bmi_ivs.tsv")


# Gunzip the 48_irnt.gwas.imputed_v3.both_sexes.tsv.bgz file
gunzip("Mediation/48_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "Mediation/waist_data.tsv")

# Read the waist_data.tsv file
waist_data <- readr::read_tsv("Mediation/waist_data.tsv")

# Use the waist_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
waist_ivs <- waist_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))

# Remove columns 7-25
waist_ivs <- waist_ivs[, -c(7:25)]

# Add the columns from waist_data to waist_ivs by matching the variant column
waist_ivs <- waist_ivs %>%
  dplyr::left_join(waist_data, by = c("variant" = "variant"))

# Remove waist_data
rm(waist_data)

# Select the IV variants
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
waist_ivs <- waist_ivs %>%
  dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
waist_ivs <- waist_ivs[, -c(11, 12, 15)]

clumped_waist_ivs <- ieugwasr::ld_clump_local(waist_ivs,  clump_kb = 10000,
                                              clump_r2 = 0.001,
                                              clump_p = 1,
                                              bfile = "1kg.v3/EUR",
                                              plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_waist_ivs data frame to a file
write_tsv(clumped_waist_ivs, "Mediation/clumped_waist_ivs.tsv")


# Gunzip the 23099_irnt.gwas.imputed_v3.both_sexes.tsv.bgz file
gunzip("Mediation/23099_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "Mediation/body_fat_data.tsv")

# Read the body_fat_data.tsv file
body_fat_data <- readr::read_tsv("Mediation/body_fat_data.tsv")

# Use the body_fat_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
body_fat_ivs <- body_fat_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))

# Remove columns 7-25
body_fat_ivs <- body_fat_ivs[, -c(7:25)]

# Add the columns from body_fat_data to body_fat_ivs by matching the variant column
body_fat_ivs <- body_fat_ivs %>%
  dplyr::left_join(body_fat_data, by = c("variant" = "variant"))

# Remove body_fat_data
rm(body_fat_data)

# Select the IV variants
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
body_fat_ivs <- body_fat_ivs %>%
  dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
body_fat_ivs <- body_fat_ivs[, -c(11, 12, 15)]

clumped_body_fat_ivs <- ieugwasr::ld_clump_local(body_fat_ivs,  clump_kb = 10000,
                                                 clump_r2 = 0.001,
                                                 clump_p = 1,
                                                 bfile = "1kg.v3/EUR",
                                                 plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_body_fat_ivs data frame to a file
write_tsv(clumped_body_fat_ivs, "Mediation/clumped_body_fat_ivs.tsv")


# Gunzip the 1558.gwas.imputed_v3.both_sexes.tsv.bgz file
gunzip("Mediation/1558.gwas.imputed_v3.both_sexes.tsv.bgz", "Mediation/alcohol_data.tsv")

# Read the alcohol_data.tsv file
alcohol_data <- readr::read_tsv("Mediation/alcohol_data.tsv")

# Use the alcohol_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
alcohol_ivs <- alcohol_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))

# Remove columns 7-25
alcohol_ivs <- alcohol_ivs[, -c(7:25)]

# Add the columns from alcohol_data to alcohol_ivs by matching the variant column
alcohol_ivs <- alcohol_ivs %>%
  dplyr::left_join(alcohol_data, by = c("variant" = "variant"))

# Remove alcohol_data
rm(alcohol_data)

# Select the IV variants
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
alcohol_ivs <- alcohol_ivs %>%
  dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
alcohol_ivs <- alcohol_ivs[, -c(11, 12, 15)]

clumped_alcohol_ivs <- ieugwasr::ld_clump_local(alcohol_ivs,  clump_kb = 10000,
                                               clump_r2 = 0.001,
                                               clump_p = 1,
                                               bfile = "1kg.v3/EUR",
                                               plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_alcohol_ivs data frame to a file
write_tsv(clumped_alcohol_ivs, "Mediation/clumped_alcohol_ivs.tsv")


# Gunzip the 2887.gwas.imputed_v3.both_sexes.tsv.bgz
gunzip("Mediation/2887.gwas.imputed_v3.both_sexes.tsv.bgz", "Mediation/smoke_data.tsv")

# Read the smoke_data.tsv file
smoke_data <- readr::read_tsv("Mediation/smoke_data.tsv")

# Use the smoke_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
smoke_ivs <- smoke_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))

# Remove columns 7-25
smoke_ivs <- smoke_ivs[, -c(7:25)]

# Add the columns from smoke_data to smoke_ivs by matching the variant column
smoke_ivs <- smoke_ivs %>%
  dplyr::left_join(smoke_data, by = c("variant" = "variant"))

# Remove smoke_data
rm(smoke_data)

# Select the IV variants
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
smoke_ivs <- smoke_ivs %>%
  #dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
smoke_ivs <- smoke_ivs[, -c(11, 12, 15)]

clumped_smoke_ivs <- ieugwasr::ld_clump_local(smoke_ivs,  clump_kb = 10000,
                                              clump_r2 = 0.001,
                                              clump_p = 1,
                                              bfile = "1kg.v3/EUR",
                                              plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_smoke_ivs data frame to a file
write_tsv(clumped_smoke_ivs, "Mediation/clumped_smoke_ivs.tsv")


# Gunzip the 30100_irnt.gwas.imputed_v3.both_sexes.tsv.bgz file
gunzip("Mediation/30100_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "Mediation/mpv_data.tsv")

# Read the mpv_data.tsv file
mpv_data <- readr::read_tsv("Mediation/mpv_data.tsv")

# Use the mpv_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
mpv_ivs <- mpv_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))

# Remove columns 7-25
mpv_ivs <- mpv_ivs[, -c(7:25)]

# Add the columns from mpv_data to mpv_ivs by matching the variant column
mpv_ivs <- mpv_ivs %>%
  dplyr::left_join(mpv_data, by = c("variant" = "variant"))

# Remove mpv_data
rm(mpv_data)

# Select the IV variants
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
mpv_ivs <- mpv_ivs %>%
  dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
mpv_ivs <- mpv_ivs[, -c(11, 12, 15)]

clumped_mpv_ivs <- ieugwasr::ld_clump_local(mpv_ivs,  clump_kb = 10000,
                                            clump_r2 = 0.001,
                                            clump_p = 1,
                                            bfile = "1kg.v3/EUR",
                                            plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_mpv_ivs data frame to a file
write_tsv(clumped_mpv_ivs, "Mediation/clumped_mpv_ivs.tsv")

# Gunzip the 30640_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz file
gunzip("Mediation/30640_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz", "Mediation/ApoB_data.tsv")

# Read the ApoB_data.tsv file
ApoB_data <- readr::read_tsv("Mediation/ApoB_data.tsv")

# Use the ApoB_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
ApoB_ivs <- ApoB_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))

# Remove columns 7-25
ApoB_ivs <- ApoB_ivs[, -c(7:25)]

# Add the columns from ApoB_data to ApoB_ivs by matching the variant column
ApoB_ivs <- ApoB_ivs %>%
  dplyr::left_join(ApoB_data, by = c("variant" = "variant"))

# Remove ApoB_data
rm(ApoB_data)

# Select the IV variants
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
ApoB_ivs <- ApoB_ivs %>%
  dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
ApoB_ivs <- ApoB_ivs[, -c(11, 12, 15)]

clumped_ApoB_ivs <- ieugwasr::ld_clump_local(ApoB_ivs,  clump_kb = 10000,
                                            clump_r2 = 0.001,
                                            clump_p = 1,
                                            bfile = "1kg.v3/EUR",
                                            plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_ApoB_ivs data frame to a file
write_tsv(clumped_ApoB_ivs, "Mediation/clumped_ApoB_ivs.tsv")

# Gunzip the 30710_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz
gunzip("Mediation/30710_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz", "Mediation/crp_data.tsv")

# Read the crp_data.tsv file
crp_data <- readr::read_tsv("Mediation/crp_data.tsv")

# Use the crp_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
crp_ivs <- crp_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))

# Remove columns 7-25
crp_ivs <- crp_ivs[, -c(7:25)]

# Add the columns from crp_data to crp_ivs by matching the variant column
crp_ivs <- crp_ivs %>%
  dplyr::left_join(crp_data, by = c("variant" = "variant"))

# Remove crp_data
rm(crp_data)

# Select the IV variants
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
crp_ivs <- crp_ivs %>%
  dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
crp_ivs <- crp_ivs[, -c(11, 12, 15)]

clumped_crp_ivs <- ieugwasr::ld_clump_local(crp_ivs,  clump_kb = 10000,
                                            clump_r2 = 0.001,
                                            clump_p = 1,
                                            bfile = "1kg.v3/EUR",
                                            plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_crp_ivs data frame to a file
write_tsv(clumped_crp_ivs, "Mediation/clumped_crp_ivs.tsv")

# Gunzip the 30010_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
gunzip("Mediation/30010_irnt.gwas.imputed_v3.both_sexes.tsv.bgz", "Mediation/rbc_data.tsv")

# Read the rbc_data.tsv file
rbc_data <- readr::read_tsv("Mediation/rbc_data.tsv")

# Use the rbc_data$variant to get the rsid of the variants from the ukbb_variants_data$rsid by matching with ukbb_variants_data$variant
rbc_ivs <- rbc_data %>%
  dplyr::select(variant) %>%
  dplyr::left_join(ukbb_variants_data, by = c("variant" = "variant"))

# Remove columns 7-25
rbc_ivs <- rbc_ivs[, -c(7:25)]

# Add the columns from rbc_data to rbc_ivs by matching the variant column
rbc_ivs <- rbc_ivs %>%
  dplyr::left_join(rbc_data, by = c("variant" = "variant"))

# Remove rbc_data
rm(rbc_data)

# Select the IV variants
# Remove variants with low confidence, minor allele frequency < 0.01, and variants with a p-value > 5e-8
# ref or alt of more than one character
rbc_ivs <- rbc_ivs %>%
  dplyr::filter(low_confidence_variant == FALSE) %>%
  dplyr::filter(minor_AF > 0.01) %>%
  dplyr::filter(pval < 5e-8) %>%
  dplyr::filter(str_length(alt) == 1) %>%
  dplyr::filter(str_length(ref) == 1)

# Remove columns 11, 12 and 15
rbc_ivs <- rbc_ivs[, -c(11, 12, 15)]

clumped_rbc_ivs <- ieugwasr::ld_clump_local(rbc_ivs,  clump_kb = 10000,
                                            clump_r2 = 0.001,
                                            clump_p = 1,
                                            bfile = "1kg.v3/EUR",
                                            plink_bin = genetics.binaRies::get_plink_binary())

# Save the clumped_rbc_ivs data frame to a file
write_tsv(clumped_rbc_ivs, "Mediation/clumped_rbc_ivs.tsv")

# Load clumped_waist_ivs, clumped_smoke_ivs, clumped_mpv_ivs, clumped_haem_ivs, clumped_body_fat_ivs, clumped_bmi_ivs, clumped_alcohol_ivs tsv files
clumped_waist_ivs <- readr::read_tsv("Mediation/clumped_waist_ivs.tsv")
clumped_smoke_ivs <- readr::read_tsv("Mediation/clumped_smoke_ivs.tsv")
clumped_mpv_ivs <- readr::read_tsv("Mediation/clumped_mpv_ivs.tsv")
clumped_haem_ivs <- readr::read_tsv("Mediation/clumped_haem_ivs.tsv")
clumped_body_fat_ivs <- readr::read_tsv("Mediation/clumped_body_fat_ivs.tsv")
clumped_bmi_ivs <- readr::read_tsv("Mediation/clumped_bmi_ivs.tsv")
clumped_alcohol_ivs <- readr::read_tsv("Mediation/clumped_alcohol_ivs.tsv")
clumped_ApoB_ivs <- readr::read_tsv("Mediation/clumped_ApoB_ivs.tsv")
clumped_crp_ivs <- readr::read_tsv("Mediation/clumped_crp_ivs.tsv")
clumped_rbc_ivs <- readr::read_tsv("Mediation/clumped_rbc_ivs.tsv")



# Combine the clumped data frames
clumped_ivs <- bind_rows(clumped_waist_ivs, clumped_smoke_ivs, clumped_mpv_ivs, clumped_haem_ivs, clumped_body_fat_ivs, clumped_bmi_ivs, clumped_alcohol_ivs,
                         clumped_ApoB_ivs, clumped_crp_ivs, clumped_rbc_ivs)
# Keep only the rows with unique variants
clumped_ivs <- clumped_ivs %>% distinct(variant, .keep_all = TRUE)
# Remove columns 2,3,4,5,7-16
clumped_ivs <- clumped_ivs[, -c(2:5, 7:16)]

# Remove the ukbb_variants_data
rm(ukbb_variants_data)


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

all_possible_variants <- c()

# For each subpathway
for (subpathway in subpathways) {
  # Initialise a unique_variants dataframe for the subpathway
  assign(paste0("unique_variants_", subpathway), clumped_ivs)
  
  # Collect the metabolites in the subpathway
  metabolites <- metabolites_data_frame$name[subpathways == subpathway]
  
  for (metabolite in metabolites) {
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
    # Change the column names to be the same as clumped_ivs
    colnames(metabolite_data) <- colnames(clumped_ivs)

    # Get current unique variants dataframe
    current_variants <- get(paste0("unique_variants_", subpathway))
    
    # Add the IVs to the {pathway}_unique_variants
    new_variants <- rbind(current_variants, metabolite_data)
    
    # Assign back the updated unique variants
    assign(paste0("unique_variants_", subpathway), new_variants)
  }
  
  # Keep the distinct variants from {pathway}_unique_variants
  final_variants <- get(paste0("unique_variants_", subpathway)) %>% distinct()
  assign(paste0("unique_variants_", subpathway), final_variants)
  
  # Add the final_variants to all_possible_variants
  all_possible_variants <- rbind(all_possible_variants, final_variants)

  # Save the unique variants to a tsv file
  write_tsv(final_variants, paste0("Mediation/unique_variants_", subpathway, ".tsv"))
}
# Keep the distinct variants from all_possible_variants
all_possible_variants <- all_possible_variants %>% distinct()
# Save the all possible variants to a tsv file
write_tsv(data.frame(variant = all_possible_variants), "Mediation/all_possible_variants.tsv")
