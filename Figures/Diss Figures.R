### Dissertation Figure Creation ###
# This script was used to create the figures shown in the dissertation.
library(ggplot2)
library(TwoSampleMR)
library(ggrepel)
library(tidyverse)
library(conflicted)
devtools::install_github("myles-lewis/locuszoomr")
library(locuszoomr)
library(cowplot)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v75", force = TRUE)
BiocManager::install("AnnotationHub")
BiocManager::install("AnnotationFilter", force = TRUE)
library(ensembldb)
library(AnnotationFilter)
library(AnnotationHub)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75


# Read in the T2DM_MR_Results.tsv file
t2dm_mr_results <- read.table("T2DM_MR_Results.tsv", header = TRUE, sep = "\t")
# Read in the responses_combined.tsv file
responses_combined <- read.table("responses_combined.tsv", header = TRUE, sep = "\t")

# In a new data frame, collect the metabolite names, Fixed_IVW_Estimate, and Fixed_IVW_Pval columns
t2dm_mr_results_subset <- t2dm_mr_results[, c("Metabolite", "Fixed_IVW_Estimate", "Fixed_IVW_Pval", "Number_of_IVs", "Random_IVW_Pval")]

# Correcting the setup for top_significant calculation
t2dm_mr_results_subset$logp <- -log10(ifelse(t2dm_mr_results_subset$Number_of_IVs <= 3, 
                                             t2dm_mr_results_subset$Fixed_IVW_Pval, 
                                             t2dm_mr_results_subset$Random_IVW_Pval))
# Ensure the data is sorted by logp correctly and top_significant is calculated based on updated logp values
t2dm_mr_results_subset <- t2dm_mr_results_subset %>%
  arrange(desc(logp)) %>%
  mutate(top_significant = logp >= 8)
# Set the color based on significance and effect size
t2dm_mr_results_subset$color <- ifelse(t2dm_mr_results_subset$logp >= -log10(3.004807e-5) & t2dm_mr_results_subset$Fixed_IVW_Estimate > 0, "red", 
                                       ifelse(t2dm_mr_results_subset$logp >= -log10(3.004807e-5) & t2dm_mr_results_subset$Fixed_IVW_Estimate < 0, "blue", "grey"))

# In a new column make responses_combined$editmetabnames where "/" and ":" are replaced with "_"
responses_combined$editmetabnames <- gsub("/", "_", responses_combined$name)
responses_combined$editmetabnames <- gsub(":", "_", responses_combined$editmetabnames)

# Add the name column to the t2dm_mr_results_subset data frame by matching responses_combined$editmetabnames to t2dm_mr_results_subset$Metabolite
t2dm_mr_results_subset <- merge(t2dm_mr_results_subset, responses_combined, by.x = "Metabolite", by.y = "editmetabnames", all.x = TRUE)

# Exponentiate the Fixed_IVW_Estimate column to get the odds ratio
t2dm_mr_results_subset$Fixed_IVW_Estimate <- exp(t2dm_mr_results_subset$Fixed_IVW_Estimate)

# Build the volcano plot
volcano_plot <- ggplot(t2dm_mr_results_subset, aes(x = Fixed_IVW_Estimate, y = logp, color = color)) +
  geom_point(alpha = 0.8) +  
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
  coord_cartesian(xlim = c(0.35, 1.65), ylim = c(0, 35)) +
  theme_minimal() +
  labs(x = "Odds Ratio",
       y = "-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_hline(yintercept = -log10(3.004807e-5), linetype = "dashed", color = "black")
# Adding labels only to the top significant metabolites, corrected placement and usage of aesthetics
volcano_plot <- volcano_plot + 
  geom_text_repel(aes(label = ifelse(top_significant, as.character(name), "")),
                  size = 3.5, box.padding = 0.5, point.padding = 0.5, segment.color = 'grey', max.overlaps = 15)

# Print the plot
print(volcano_plot)


# Do the same for FG
# Read in the FG_MR_Results.tsv file
fg_mr_results <- read.table("FG_MR_Results.tsv", header = TRUE, sep = "\t")

# In a new data frame, collect the metabolite names, Fixed_IVW_Estimate, and Fixed_IVW_Pval columns
fg_mr_results_subset <- fg_mr_results[, c("Metabolite", "Fixed_IVW_Estimate", "Fixed_IVW_Pval", "Number_of_IVs", "Random_IVW_Pval")]

# Correcting the setup for top_significant calculation 
fg_mr_results_subset$logp <- -log10(ifelse(fg_mr_results_subset$Number_of_IVs <= 3, 
                                           fg_mr_results_subset$Fixed_IVW_Pval, 
                                           fg_mr_results_subset$Random_IVW_Pval))
# Ensure the data is sorted by logp correctly and top_significant is calculated based on updated logp values
fg_mr_results_subset <- fg_mr_results_subset %>%
  arrange(desc(logp)) %>%
  mutate(top_significant = logp >= 18)
# Set the color based on significance and effect size
fg_mr_results_subset$color <- ifelse(fg_mr_results_subset$logp >= -log10(3.004807e-5) & fg_mr_results_subset$Fixed_IVW_Estimate > 0, "red", 
                                     ifelse(fg_mr_results_subset$logp >= -log10(3.004807e-5) & fg_mr_results_subset$Fixed_IVW_Estimate < 0, "blue", "grey"))

# Add the name column to the fg_mr_results_subset data frame by matching responses_combined$editmetabnames to fg_mr_results_subset$Metabolite
fg_mr_results_subset <- merge(fg_mr_results_subset, responses_combined, by.x = "Metabolite", by.y = "editmetabnames", all.x = TRUE)

# Build the volcano plot
volcano_plot <- ggplot(fg_mr_results_subset, aes(x = Fixed_IVW_Estimate, y = logp, color = color)) +
  geom_point(alpha = 0.8) +  
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
  coord_cartesian(xlim = c(-0.4, 0.4), ylim = c(0, 40)) +
  theme_minimal() +
  labs(x = "Estimate",
       y = "-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_hline(yintercept = -log10(3.004807e-5), linetype = "dashed", color = "black")
# Adding labels only to the top significant metabolites, corrected placement and usage of aesthetics
volcano_plot <- volcano_plot + 
  geom_text_repel(aes(label = ifelse(top_significant, as.character(name), "")),
                  size = 3.5, box.padding = 0.5, point.padding = 0.5, segment.color = 'grey', max.overlaps = 8)

# Print the plot
print(volcano_plot)


# Do the same for HbA1c
# Read in the HbA1c_MR_Results.tsv file
hba1c_mr_results <- read.table("HbA1c_MR_Results.tsv", header = TRUE, sep = "\t")

# In a new data frame, collect the metabolite names, Fixed_IVW_Estimate, and Fixed_IVW_Pval columns
hba1c_mr_results_subset <- hba1c_mr_results[, c("Metabolite", "Fixed_IVW_Estimate", "Fixed_IVW_Pval", "Number_of_IVs", "Random_IVW_Pval")]

# Correcting the setup for top_significant calculation
hba1c_mr_results_subset$logp <- -log10(ifelse(hba1c_mr_results_subset$Number_of_IVs <= 3, 
                                              hba1c_mr_results_subset$Fixed_IVW_Pval, 
                                              hba1c_mr_results_subset$Random_IVW_Pval))
# Ensure the data is sorted by logp correctly and top_significant is calculated based on updated logp values
hba1c_mr_results_subset <- hba1c_mr_results_subset %>%
  arrange(desc(logp)) %>%
  mutate(top_significant = logp >= 7)
# Set the color based on significance and effect size
hba1c_mr_results_subset$color <- ifelse(hba1c_mr_results_subset$logp >= -log10(3.004807e-5) & hba1c_mr_results_subset$Fixed_IVW_Estimate > 0, "red", 
                                        ifelse(hba1c_mr_results_subset$logp >= -log10(3.004807e-5) & hba1c_mr_results_subset$Fixed_IVW_Estimate < 0, "blue", "grey"))

# Add the name column to the hba1c_mr_results_subset data frame by matching responses_combined$editmetabnames to hba1c_mr_results_subset$Metabolite
hba1c_mr_results_subset <- merge(hba1c_mr_results_subset, responses_combined, by.x = "Metabolite", by.y = "editmetabnames", all.x = TRUE)

# Build the volcano plot
volcano_plot <- ggplot(hba1c_mr_results_subset, aes(x = Fixed_IVW_Estimate, y = logp, color = color)) +
  geom_point(alpha = 0.8) +  
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(0, 310)) +
  theme_minimal() +
  labs(x = "Estimate",
       y = "-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_hline(yintercept = -log10(3.004807e-5), linetype = "dashed", color = "black")
# Adding labels only to the top significant metabolites, corrected placement and usage of aesthetics
volcano_plot <- volcano_plot + 
  geom_text_repel(aes(label = ifelse(top_significant, as.character(name), "")),
                  size = 3.5, box.padding = 0.5, point.padding = 0.5, segment.color = 'grey', max.overlaps = 4.5)

# Print the plot
print(volcano_plot)

# Build the volcano plot
volcano_plot <- ggplot(hba1c_mr_results_subset, aes(x = Fixed_IVW_Estimate, y = logp, color = color)) +
  geom_point(alpha = 0.8) +  
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
  coord_cartesian(xlim = c(-0.15, 0.15), ylim = c(0, 13)) +
  theme_minimal() +
  labs(x = "Estimate",
       y = "-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_hline(yintercept = -log10(3.004807e-5), linetype = "dashed", color = "black")
# Adding labels only to the top significant metabolites, corrected placement and usage of aesthetics
volcano_plot <- volcano_plot + 
  geom_text_repel(aes(label = ifelse(top_significant, as.character(name), "")),
                  size = 3.5, box.padding = 0.7, point.padding = 0.5, segment.color = 'grey', max.overlaps = 6)

# Print the plot
print(volcano_plot)



### Individual MR Plots ###

# Read in the T2DM_MR_Results.tsv file
t2dm_mr_results <- read.table("T2DM_MR_Results.tsv", header = TRUE, sep = "\t")
# Read in the FG_MR_Results.tsv file
fg_mr_results <- read.table("FG_MR_Results.tsv", header = TRUE, sep = "\t")
# Read in the HbA1c_MR_Results.tsv file
hba1c_mr_results <- read.table("HbA1c_MR_Results.tsv", header = TRUE, sep = "\t")

# Load the Harmonsied_T2DM_IVs files
t2dm_ivs_files <- list.files("Harmonised_T2DM_IVs", full.names = TRUE)
# Load the Harmonsied_FG_IVs files
fg_ivs_files <- list.files("Harmonised_FG_IVs", full.names = TRUE)
# Load the Harmonsied_HbA1c_IVs files
hba1c_ivs_files <- list.files("Harmonised_HbA1c_IVs", full.names = TRUE)


# Load the T2DM IVs for the metabolite: 1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)*
t2dm_metab_ivs <- read.table("Harmonised_T2DM_IVs/1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)*T2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Plot the MR results for cholesterol
dat <- data.frame(beta.exposure = t2dm_metab_ivs$Beta, se.exposure = t2dm_metab_ivs$SE,
                  beta.outcome = t2dm_metab_ivs$t2dm_Beta, se.outcome = t2dm_metab_ivs$t2dm_SE,
                  id.exposure = "1-palmitoyl-2-arachidonoyl-GPI (16:0/20:4)*", id.outcome = "T2DM", mr_keep = TRUE, exposure = "1-palmitoyl-2-arachidonoyl-GPI (16:0/20:4)*", outcome = "Type 2 Diabetes",
                  SNP = t2dm_metab_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_forest_plot(res_single, exponentiate = TRUE)


# Load the FG IVs for the metabolite: sphingomyelin (d18_1_18_1, d18_2_18_0) 
fg_metab_ivs <- read.table("Harmonised_FG_IVs/sphingomyelin (d18_1_18_1, d18_2_18_0)FGFG_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Plot the MR results for cholesterol
dat <- data.frame(beta.exposure = fg_metab_ivs$Beta, se.exposure = fg_metab_ivs$SE,
                  beta.outcome = fg_metab_ivs$fg_Beta, se.outcome = fg_metab_ivs$fg_SE,
                  id.exposure = "sphingomyelin (d18:1/18:1, d18:2/18:0)", id.outcome = "Fasting Glucose", mr_keep = TRUE, exposure = "sphingomyelin (d18:1/18:1, d18:2/18:0)", outcome = "Fasting Glucose",
                  SNP = fg_metab_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_forest_plot(res_single)


# Load the T2DM IVs for the metabolite: 1-myristoyl-2-linoleoyl-GPC (14_0_18_2)*
t2dm_metab_ivs <- read.table("Harmonised_T2DM_IVs/1-myristoyl-2-linoleoyl-GPC (14_0_18_2)*T2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Plot the MR results for cholesterol
dat <- data.frame(beta.exposure = t2dm_metab_ivs$Beta, se.exposure = t2dm_metab_ivs$SE,
                  beta.outcome = t2dm_metab_ivs$t2dm_Beta, se.outcome = t2dm_metab_ivs$t2dm_SE,
                  id.exposure = "1-myristoyl-2-linoleoyl-GPC (14:0/18:2)*", id.outcome = "T2DM", mr_keep = TRUE, exposure = "1-myristoyl-2-linoleoyl-GPC (14:0/18:2)*", outcome = "Type 2 Diabetes",
                  SNP = t2dm_metab_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_fe", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_forest_plot(res_single)
# height 700 width 550



### Reverse MR Plots ###

# Load significant_T2DM_results.tsv
t2dm_results <- read.table("significant_T2DM_results_diss.tsv", header = TRUE, sep = "\t")

# Load Sig_Metabolites_Compid_diss.tsv
responses <- read.table("Sig_Metabolites_Compid_diss.tsv", header = TRUE, sep = "\t")

# Load Reverse_MR/T2DM_Reverse_MR_Results.tsv
t2dm_reverse_mr_results <- read.table("Reverse_MR/Reverse_Sig_T2DM_Results.tsv", header = TRUE, sep = "\t")
# Add the name of the metabolite to t2dm_reverse_mr_results by matching the Metabolite column with the edit_name column in responses
t2dm_reverse_mr_results <- merge(t2dm_reverse_mr_results, responses, by.x = "Metabolite", by.y = "compid", all.x = TRUE)
# Change / and : to _ in the name column
t2dm_reverse_mr_results$name <- gsub("/", "_", t2dm_reverse_mr_results$name)
t2dm_reverse_mr_results$name <- gsub(":", "_", t2dm_reverse_mr_results$name)

# Keep only the metabolites from t2dm_results that are in t2dm_reverse_mr_results
t2dm_results <- t2dm_results[t2dm_results$Metabolite %in% t2dm_reverse_mr_results$name,]
t2dm_reverse_mr_results <- t2dm_reverse_mr_results[t2dm_reverse_mr_results$name %in% t2dm_results$Metabolite,]

# Merge the forward and reverse results by metabolite name
merged_results <- merge(t2dm_results, t2dm_reverse_mr_results, by.x = "Metabolite", by.y = "name")
# Merge responses with merged_results by Metabolite.y and compid
merged_results <- merge(merged_results, responses, by.x = "Metabolite.y", by.y = "compid")

# Calculate the CIs for the IVW estimates
merged_results$Fixed_IVW_CI_Lower.x <- merged_results$Fixed_IVW_Estimate.x - 1.96 * merged_results$Fixed_IVW_SE.x
merged_results$Fixed_IVW_CI_Upper.x <- merged_results$Fixed_IVW_Estimate.x + 1.96 * merged_results$Fixed_IVW_SE.x
merged_results$Fixed_IVW_CI_Lower.y <- merged_results$Fixed_IVW_Estimate.y - 1.96 * merged_results$Fixed_IVW_SE.y
merged_results$Fixed_IVW_CI_Upper.y <- merged_results$Fixed_IVW_Estimate.y + 1.96 * merged_results$Fixed_IVW_SE.y

# Create the scatter plot with error bars
plot <- ggplot(merged_results, aes(x = Fixed_IVW_Estimate.x, y = Fixed_IVW_Estimate.y)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_errorbar(aes(ymin = Fixed_IVW_CI_Lower.y, ymax = Fixed_IVW_CI_Upper.y), width = 0.02) +
  geom_errorbarh(aes(xmin = Fixed_IVW_CI_Lower.x, xmax = Fixed_IVW_CI_Upper.x), height = 0.02) +
  geom_text_repel(aes(label = name), size = 4, box.padding = 0.75, point.padding = 1, segment.color = 'grey50', max.overlaps = 7) +
  theme_minimal() +
  labs(
    x = "Forward MR Estimate (Metabolite to Type 2 Diabetes)",
    y = "Reverse MR Estimate (Type 2 Diabetes to Metabolite)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold")) +
  geom_hline(yintercept = 0, size = 0.5, color = "black") +
  geom_vline(xintercept = 0, size = 0.5, color = "black") +
  coord_cartesian(xlim = c(-0.9, 0.9), ylim = c(-0.9, 0.9))

# Print the plot
print(plot)


# Load the HBA1C IVs for the metabolite: M03147
hba1c_M03147_ivs <- read.table("Reverse_MR/Harmonised_Reverse_HbA1C_IVs/M03147_HBA1C_harmonised_IVs.tsv", header = TRUE, sep = "\t")
dat <- data.frame(beta.exposure = hba1c_M03147_ivs$Beta, se.exposure = hba1c_M03147_ivs$SE,
                  beta.outcome = hba1c_M03147_ivs$out_Beta, se.outcome = hba1c_M03147_ivs$out_SE,
                  id.exposure = "HbA1c", id.outcome = "Xanthine", mr_keep = TRUE, exposure = "HbA1c", outcome = "Xanthine",
                  SNP = hba1c_M03147_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
mr_forest_plot(res_single)
mr_funnel_plot(res_single)

# Load the FG IVs for the metabolite: M48047
fg_M48047_ivs <- read.table("Reverse_MR/Harmonised_Reverse_FG_IVs/M48047_FG_harmonised_IVs.tsv", header = TRUE, sep = "\t")
dat <- data.frame(beta.exposure = fg_M48047_ivs$Beta, se.exposure = fg_M48047_ivs$SE,
                  beta.outcome = fg_M48047_ivs$out_Beta, se.outcome = fg_M48047_ivs$out_SE,
                  id.exposure = "Fasting Glucose", id.outcome = "X - 18886", mr_keep = TRUE, exposure = "Fasting Glucose", outcome = "X - 18886",
                  SNP = fg_M48047_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
mr_forest_plot(res_single)
mr_funnel_plot(res_single)



### Colocalisation Plots

### High Unconditional Colocalisation ####
# Load Colocalisation/PwCoCo_Ready_T2DM/sphinganine_rs676457_metabolite_snps_pwcoco.tsv
exposure_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_T2DM/sphinganine_rs676457_metabolite_snps_pwcoco.tsv")
# Load Colocalisation/PwCoCo_Ready_T2DM/sphinganine_rs676457_t2dm_snps_pwcoco.tsv
outcome_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_T2DM/sphinganine_rs676457_t2dm_snps_pwcoco.tsv")
# Load Colocalisation/Preped_Regions/T2DM_Regions/sphinganine_rs676457_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/T2DM_Regions/sphinganine_rs676457_metabolite_snps.tsv")
# Load Colocalisation/Preped_Regions/T2DM_Regions/sphinganine_rs676457_t2dm_snps.tsv.tsv
outcome_data_raw <- read_tsv("Colocalisation/Preped_Regions/T2DM_Regions/sphinganine_rs676457_t2dm_snps.tsv")
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

# In both exposure_data_plot and outcome_data_plot, change the value for snp chr9:136146227:A:T to rs676457
exposure_data_plot$snp[exposure_data_plot$snp == "chr9:136146227:A:T"] <- "rs676457"
outcome_data_plot$snp[outcome_data_plot$snp == "chr9:136146227:A:T"] <- "rs676457"

loc_exposure <- locus(gene = "ABO", data = exposure_data_plot, fix_window = 1e6,
                      ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "ABO", data = outcome_data_plot, fix_window = 1e6,
                     ens_db = edb)
summary(loc_outcome)
loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("rs676457"),
           label_x = c(4, -5), maxrows = 3)
locus_plot(loc_outcome, labels = c("rs676457"),
           label_x = c(4, -5), maxrows = 3)


### Low Unconditional and Conditional Colocalisation ###
# Load Colocalisation/PwCoCo_Ready_HBA1C/sphinganine_rs676457_metabolite_snps_pwcoco.tsv
exposure_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_HBA1C/sphinganine_rs676457_metabolite_snps_pwcoco.tsv")
# Load Colocalisation/PwCoCo_Ready_HBA1C/sphinganine_rs676457_t2dm_snps_pwcoco.tsv
outcome_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_HBA1C/sphinganine_rs676457_hba1c_snps_pwcoco.tsv")
# Load Colocalisation/Preped_Regions/HBA1C_Regions/sphinganine_rs676457_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/HBA1C_Regions/sphinganine_rs676457_metabolite_snps.tsv")
# Load Colocalisation/Preped_Regions/HBA1C_Regions/sphinganine_rs676457_t2dm_snps.tsv
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

# In both exposure_data_plot and outcome_data_plot, change the value for snp chr9:136146227:A:T to rs676457
exposure_data_plot$snp[exposure_data_plot$snp == "chr9:136146227:A:T"] <- "rs676457"
outcome_data_plot$snp[outcome_data_plot$snp == "chr9:136146227:A:T"] <- "rs676457"

loc_exposure <- locus(gene = "ABO", data = exposure_data_plot, fix_window = 1e6,
                      ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "ABO", data = outcome_data_plot, fix_window = 1e6,
                     ens_db = edb)
summary(loc_outcome)
loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("rs676457"),
           label_x = c(4, -5), maxrows = 3)
locus_plot(loc_outcome, labels = c("rs676457"),
           label_x = c(7, -10), maxrows = 3)


### Low Unconditional and Conditional Colocalisation ###
# Load Colocalisation/PwCoCo_Ready_FG/sphinganine_rs676457_metabolite_snps_pwcoco.tsv
exposure_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_FG/sphinganine_rs676457_metabolite_snps_pwcoco.tsv")
# Load Colocalisation/PwCoCo_Ready_FG/sphinganine_rs676457_fg_snps_pwcoco.tsv
outcome_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_FG/sphinganine_rs676457_fg_snps_pwcoco.tsv")
# Load Colocalisation/Preped_Regions/FG_Regions/sphinganine_rs676457_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/FG_Regions/sphinganine_rs676457_metabolite_snps.tsv")
# Load Colocalisation/Preped_Regions/FG_Regions/sphinganine_rs676457_fg_snps.tsv
outcome_data_raw <- read_tsv("Colocalisation/Preped_Regions/FG_Regions/sphinganine_rs676457_fg_snps.tsv")
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

# In both exposure_data_plot and outcome_data_plot, change the value for snp chr9:136146227:A:T to rs676457
exposure_data_plot$snp[exposure_data_plot$snp == "chr9:136146227:A:T"] <- "rs676457"
outcome_data_plot$snp[outcome_data_plot$snp == "chr9:136146227:A:T"] <- "rs676457"

loc_exposure <- locus(gene = "ABO", data = exposure_data_plot, fix_window = 1e6,
                      ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "ABO", data = outcome_data_plot, fix_window = 1e6,
                     ens_db = edb)
summary(loc_outcome)
loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("rs676457"),
           label_x = c(4, -5), maxrows = 3)
locus_plot(loc_outcome, labels = c("rs676457"),
           label_x = c(7, -10), maxrows = 3)


### Low Unconditional but High Conditional Colocalisation ###
# Load Colocalisation/PwCoCo_Ready_T2DM/5-acetylamino-6-formylamino-3-methyluracil_rs4921913_metabolite_snps_pwcoco.tsv
exposure_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_T2DM/5-acetylamino-6-formylamino-3-methyluracil_rs4921913_metabolite_snps_pwcoco.tsv")
# Load Colocalisation/PwCoCo_Ready_T2DM/5-acetylamino-6-formylamino-3-methyluracil_rs4921913_t2dm_snps_pwcoco.tsv
outcome_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_T2DM/5-acetylamino-6-formylamino-3-methyluracil_rs4921913_t2dm_snps_pwcoco.tsv")
# Load Colocalisation/Preped_Regions/T2DM_Regions/5-acetylamino-6-formylamino-3-methyluracil_rs4921913_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/T2DM_Regions/5-acetylamino-6-formylamino-3-methyluracil_rs4921913_metabolite_snps.tsv")
# Load Colocalisation/Preped_Regions/T2DM_Regions/5-acetylamino-6-formylamino-3-methyluracil_rs4921913_t2dm_snps.tsv
outcome_data_raw <- read_tsv("Colocalisation/Preped_Regions/T2DM_Regions/5-acetylamino-6-formylamino-3-methyluracil_rs4921913_t2dm_snps.tsv")
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

# In both exposure_data_plot and outcome_data_plot, change the value for snp chr8:18272377:C:T to rs4921913
exposure_data_plot$snp[exposure_data_plot$snp == "chr8:18272377:C:T"] <- "rs4921913"
outcome_data_plot$snp[outcome_data_plot$snp == "chr8:18272377:T:C"] <- "rs4921913"

loc_exposure <- locus(gene = "NAT2", data = exposure_data_plot, fix_window = 1e6,
                      ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "NAT2", data = outcome_data_plot, fix_window = 1e6,
                     ens_db = edb)
summary(loc_outcome)

loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("rs4921913"),
           label_x = c(4, -5), maxrows = 3)
locus_plot(loc_outcome, labels = c("rs4921913"),
           label_x = c(4, -5), maxrows = 3)        


### Ambigous Colocalisation ###
# Load Colocalisation/PwCoCo_Ready_FG/sphingomyelin (d18_2_14_0, d18_1_14_1)*_rs174544_metabolite_snps_pwcoco.tsv
exposure_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_FG/sphingomyelin (d18_2_14_0, d18_1_14_1)*_rs174544_metabolite_snps_pwcoco.tsv")
# Load Colocalisation/PwCoCo_Ready_FG/sphingomyelin (d18_2_14_0, d18_1_14_1)*_rs174544_t2dm_snps_pwcoco.tsv
outcome_data_included <- read_tsv("Colocalisation/PwCoCo_Ready_FG/sphingomyelin (d18_2_14_0, d18_1_14_1)*_rs174544_fg_snps_pwcoco.tsv")
# Load Colocalisation/Preped_Regions/FG_Regions/sphingomyelin (d18_2_14_0, d18_1_14_1)*_rs174544_metabolite_snps.tsv
exposure_data_raw <- read_tsv("Colocalisation/Preped_Regions/FG_Regions/sphingomyelin (d18_2_14_0, d18_1_14_1)*_rs174544_metabolite_snps.tsv")
# Load Colocalisation/Preped_Regions/FG_Regions/sphingomyelin (d18_2_14_0, d18_1_14_1)*_rs174544_t2dm_snps.tsv
outcome_data_raw <- read_tsv("Colocalisation/Preped_Regions/FG_Regions/sphingomyelin (d18_2_14_0, d18_1_14_1)*_rs174544_fg_snps.tsv")
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

# In both exposure_data_plot and outcome_data_plot, change the value for snp chr11:61567753:A:C to rs174544
exposure_data_plot$snp[exposure_data_plot$snp == "chr11:61567753:A:C"] <- "rs174544"
outcome_data_plot$snp[outcome_data_plot$snp == "chr11:61567753:A:C"] <- "rs174544"

loc_exposure <- locus(gene = "FEN1", data = exposure_data_plot, fix_window = 1e6,
                      ens_db = edb)
summary(loc_exposure)
loc_exposure <- link_recomb(loc_exposure, genome = "hg19")

loc_outcome <- locus(gene = "FEN1", data = outcome_data_plot, fix_window = 1e6,
                     ens_db = edb)
summary(loc_outcome)
loc_outcome <- link_recomb(loc_outcome, genome = "hg19")

locus_plot(loc_exposure, labels = c("rs174544"),
           label_x = c(4, -5), maxrows = 3)
locus_plot(loc_outcome, labels = c("rs174544"),
           label_x = c(-7, 0), maxrows = 3)
