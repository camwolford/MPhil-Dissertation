library(ggplot2)
library(TwoSampleMR)
library(dplyr)

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
  mutate(top_significant = logp >= 12)
# Set the color based on significance and effect size
t2dm_mr_results_subset$color <- ifelse(t2dm_mr_results_subset$logp >= -log10(3.06e-5) & t2dm_mr_results_subset$Fixed_IVW_Estimate > 0, "red", 
                                       ifelse(t2dm_mr_results_subset$logp >= -log10(3.06e-5) & t2dm_mr_results_subset$Fixed_IVW_Estimate < 0, "blue", "grey"))
# Build the volcano plot
volcano_plot <- ggplot(t2dm_mr_results_subset, aes(x = Fixed_IVW_Estimate, y = logp, color = color)) +
  geom_point(alpha = 0.8) +  
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
  coord_cartesian(xlim = c(-0.8, 0.8), ylim = c(0, 35)) +
  theme_minimal() +
  labs(title = "Volcano Plot of MR Metabolite-T2DM Associations",
       x = "Effect Size (Beta)",
       y = "-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_hline(yintercept = -log10(3.06e-5), linetype = "dashed", color = "black")
# Adding labels only to the top significant metabolites, corrected placement and usage of aesthetics
volcano_plot <- volcano_plot + 
  geom_text(aes(label = ifelse(top_significant, as.character(t2dm_mr_results_subset$Metabolite), ""), color = color),
            vjust = 1.5, hjust = 0.5, size = 3.5, check_overlap = TRUE)
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
fg_mr_results_subset$color <- ifelse(fg_mr_results_subset$logp >= -log10(3.06e-5) & fg_mr_results_subset$Fixed_IVW_Estimate > 0, "red", 
                                     ifelse(fg_mr_results_subset$logp >= -log10(3.06e-5) & fg_mr_results_subset$Fixed_IVW_Estimate < 0, "blue", "grey"))
# Build the volcano plot
volcano_plot <- ggplot(fg_mr_results_subset, aes(x = Fixed_IVW_Estimate, y = logp, color = color)) +
  geom_point(alpha = 0.8) +  
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(0, 40)) +
  theme_minimal() +
  labs(title = "Volcano Plot of MR Metabolite-FG Associations",
       x = "Effect Size (Beta)",
       y = "-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_hline(yintercept = -log10(3.06e-5), linetype = "dashed", color = "black")
# Adding labels only to the top significant metabolites, corrected placement and usage of aesthetics
volcano_plot <- volcano_plot + 
  geom_text(data = fg_mr_results_subset[fg_mr_results_subset$top_significant, ],
            aes(label = Metabolite, color = color),
            vjust = 1.5, hjust = 0.5, size = 3.5, check_overlap = TRUE)
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
  mutate(top_significant = logp >= 9)
# Set the color based on significance and effect size
hba1c_mr_results_subset$color <- ifelse(hba1c_mr_results_subset$logp >= -log10(3.06e-5) & hba1c_mr_results_subset$Fixed_IVW_Estimate > 0, "red", 
                                       ifelse(hba1c_mr_results_subset$logp >= -log10(3.06e-5) & hba1c_mr_results_subset$Fixed_IVW_Estimate < 0, "blue", "grey"))
# Build the volcano plot
volcano_plot <- ggplot(hba1c_mr_results_subset, aes(x = Fixed_IVW_Estimate, y = logp, color = color)) +
  geom_point(alpha = 0.8) +  
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
  coord_cartesian(xlim = c(-0.15, 0.15), ylim = c(0, 15)) +
  theme_minimal() +
  labs(title = "Volcano Plot of MR Metabolite-HbA1c Associations",
       x = "Effect Size (Beta)",
       y = "-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  geom_hline(yintercept = -log10(3.06e-5), linetype = "dashed", color = "black")
# Adding labels only to the top significant metabolites, corrected placement and usage of aesthetics
volcano_plot <- volcano_plot + 
  geom_text(aes(label = ifelse(top_significant, as.character(hba1c_mr_results_subset$Metabolite), ""), color = color),
            vjust = 1.5, hjust = 0.5, size = 3.5, check_overlap = TRUE)
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

# Load the T2DM IVs for the metabolite: cholesterol
t2dm_cholesterol_ivs <- read.table("Harmonised_T2DM_IVs/cholesterolT2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Collect the previous MR results for cholesterol from the T2DM_MR_Results.tsv file
t2dm_mr_results[t2dm_mr_results$Metabolite == "cholesterol", ]
# Replicate previous MR analysis for glutamine using the TwoSampleMR package
mr_ivw_fe(t2dm_cholesterol_ivs$Beta, t2dm_cholesterol_ivs$t2dm_Beta, t2dm_cholesterol_ivs$SE, t2dm_cholesterol_ivs$t2dm_SE) # Replicated
mr_ivw_mre(t2dm_cholesterol_ivs$Beta, t2dm_cholesterol_ivs$t2dm_Beta, t2dm_cholesterol_ivs$SE, t2dm_cholesterol_ivs$t2dm_SE) # Replicated
mr_penalised_weighted_median(t2dm_cholesterol_ivs$Beta, t2dm_cholesterol_ivs$t2dm_Beta, t2dm_cholesterol_ivs$SE, t2dm_cholesterol_ivs$t2dm_SE) # ~ Replicated but used penalised version
mr_weighted_mode(t2dm_cholesterol_ivs$Beta, t2dm_cholesterol_ivs$t2dm_Beta, t2dm_cholesterol_ivs$SE, t2dm_cholesterol_ivs$t2dm_SE) # ~ Replicated
mr_egger_regression(t2dm_cholesterol_ivs$Beta, t2dm_cholesterol_ivs$t2dm_Beta, t2dm_cholesterol_ivs$SE, t2dm_cholesterol_ivs$t2dm_SE) # Replicated
# Plot the MR results for glutamine
dat <- data.frame(beta.exposure = t2dm_cholesterol_ivs$Beta, se.exposure = t2dm_cholesterol_ivs$SE,
                  beta.outcome = t2dm_cholesterol_ivs$t2dm_Beta, se.outcome = t2dm_cholesterol_ivs$t2dm_SE,
                  id.exposure = "Cholesterol", id.outcome = "T2DM", mr_keep = TRUE, exposure = "Cholesterol", outcome = "T2DM",
                  SNP = t2dm_cholesterol_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_fe", "mr_ivw_mre", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_forest_plot(res_single)
res_loo <- mr_leaveoneout(dat, method = mr_ivw_fe)
mr_leaveoneout_plot(res_loo)
mr_funnel_plot(res_single)


# Load the T2DM IVs for the metabolite: 1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)*
t2dm_16_0_20_4_ivs <- read.table("Harmonised_T2DM_IVs/1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)*T2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Collect the previous MR results for 1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)* from the T2DM_MR_Results.tsv file
t2dm_mr_results[t2dm_mr_results$Metabolite == "1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)*", ]
# Replicate previous MR analysis for 1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)* using the TwoSampleMR package
mr_ivw_fe(t2dm_16_0_20_4_ivs$Beta, t2dm_16_0_20_4_ivs$t2dm_Beta, t2dm_16_0_20_4_ivs$SE, t2dm_16_0_20_4_ivs$t2dm_SE) # Replicated
mr_ivw_mre(t2dm_16_0_20_4_ivs$Beta, t2dm_16_0_20_4_ivs$t2dm_Beta, t2dm_16_0_20_4_ivs$SE, t2dm_16_0_20_4_ivs$t2dm_SE) # Replicated
mr_penalised_weighted_median(t2dm_16_0_20_4_ivs$Beta, t2dm_16_0_20_4_ivs$t2dm_Beta, t2dm_16_0_20_4_ivs$SE, t2dm_16_0_20_4_ivs$t2dm_SE) # ~ Replicated but used penalised version
mr_weighted_mode(t2dm_16_0_20_4_ivs$Beta, t2dm_16_0_20_4_ivs$t2dm_Beta, t2dm_16_0_20_4_ivs$SE, t2dm_16_0_20_4_ivs$t2dm_SE) # ~ Replicated
mr_egger_regression(t2dm_16_0_20_4_ivs$Beta, t2dm_16_0_20_4_ivs$t2dm_Beta, t2dm_16_0_20_4_ivs$SE, t2dm_16_0_20_4_ivs$t2dm_SE) # Replicated
# Plot the MR results for 1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)*
dat <- data.frame(beta.exposure = t2dm_16_0_20_4_ivs$Beta, se.exposure = t2dm_16_0_20_4_ivs$SE,
                  beta.outcome = t2dm_16_0_20_4_ivs$t2dm_Beta, se.outcome = t2dm_16_0_20_4_ivs$t2dm_SE,
                  id.exposure = "1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)*", id.outcome = "T2DM", mr_keep = TRUE, exposure = "1-palmitoyl-2-arachidonoyl-GPI (16_0_20_4)*", outcome = "T2DM",
                  SNP = t2dm_16_0_20_4_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_fe", "mr_ivw_mre", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_forest_plot(res_single)
res_loo <- mr_leaveoneout(dat, method = mr_ivw_fe)
mr_leaveoneout_plot(res_loo)
mr_funnel_plot(res_single)

# Load the T2DM IVs for the metabolite: carnitine
t2dm_carnitine_ivs <- read.table("Harmonised_T2DM_IVs/carnitineT2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Collect the previous MR results for carnitine from the T2DM_MR_Results.tsv file
t2dm_mr_results[t2dm_mr_results$Metabolite == "carnitine", ]
# Replicate previous MR analysis for carnitine using the TwoSampleMR package
mr_ivw_fe(t2dm_carnitine_ivs$Beta, t2dm_carnitine_ivs$t2dm_Beta, t2dm_carnitine_ivs$SE, t2dm_carnitine_ivs$t2dm_SE) # Replicated
mr_ivw_mre(t2dm_carnitine_ivs$Beta, t2dm_carnitine_ivs$t2dm_Beta, t2dm_carnitine_ivs$SE, t2dm_carnitine_ivs$t2dm_SE) # Replicated
mr_penalised_weighted_median(t2dm_carnitine_ivs$Beta, t2dm_carnitine_ivs$t2dm_Beta, t2dm_carnitine_ivs$SE, t2dm_carnitine_ivs$t2dm_SE) # Replicated
mr_weighted_mode(t2dm_carnitine_ivs$Beta, t2dm_carnitine_ivs$t2dm_Beta, t2dm_carnitine_ivs$SE, t2dm_carnitine_ivs$t2dm_SE) # Replicated
mr_egger_regression(t2dm_carnitine_ivs$Beta, t2dm_carnitine_ivs$t2dm_Beta, t2dm_carnitine_ivs$SE, t2dm_carnitine_ivs$t2dm_SE) # Replicated
# Plot the MR results for carnitine
dat <- data.frame(beta.exposure = t2dm_carnitine_ivs$Beta, se.exposure = t2dm_carnitine_ivs$SE,
                  beta.outcome = t2dm_carnitine_ivs$t2dm_Beta, se.outcome = t2dm_carnitine_ivs$t2dm_SE,
                  id.exposure = "carnitine", id.outcome = "T2DM", mr_keep = TRUE, exposure = "carnitine", outcome = "T2DM",
                  SNP = t2dm_carnitine_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_fe", "mr_ivw_mre", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_forest_plot(res_single)

# Load the T2DM IVs for the metabolite: gamma-glutamyltyrosine
t2dm_gamma_glutamyltyrosine_ivs <- read.table("Harmonised_T2DM_IVs/gamma-glutamyltyrosineT2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Collect the previous MR results for gamma-glutamyltyrosine from the T2DM_MR_Results.tsv file
t2dm_mr_results[t2dm_mr_results$Metabolite == "gamma-glutamyltyrosine", ]
# Replicate previous MR analysis for gamma-glutamyltyrosine using the TwoSampleMR package
mr_ivw_fe(t2dm_gamma_glutamyltyrosine_ivs$Beta, t2dm_gamma_glutamyltyrosine_ivs$t2dm_Beta, t2dm_gamma_glutamyltyrosine_ivs$SE, t2dm_gamma_glutamyltyrosine_ivs$t2dm_SE) # Replicated
# Plot the MR results for gamma-glutamyltyrosine
dat <- data.frame(beta.exposure = t2dm_gamma_glutamyltyrosine_ivs$Beta, se.exposure = t2dm_gamma_glutamyltyrosine_ivs$SE,
                  beta.outcome = t2dm_gamma_glutamyltyrosine_ivs$t2dm_Beta, se.outcome = t2dm_gamma_glutamyltyrosine_ivs$t2dm_SE,
                  id.exposure = "gamma-glutamyltyrosine", id.outcome = "T2DM", mr_keep = TRUE, exposure = "gamma-glutamyltyrosine", outcome = "T2DM",
                  SNP = t2dm_gamma_glutamyltyrosine_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_ivw_mre"))
mr_scatter_plot(res, dat)

# Load the T2DM IVs for the metabolite: 3-methoxytyrosine
t2dm_3_methoxytyrosine_ivs <- read.table("Harmonised_T2DM_IVs/3-methoxytyrosineT2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Collect the previous MR results for 3-methoxytyrosine from the T2DM_MR_Results.tsv file
t2dm_mr_results[t2dm_mr_results$Metabolite == "3-methoxytyrosine", ]
# Replicate previous MR analysis for 3-methoxytyrosine using the TwoSampleMR package
mr_ivw_fe(t2dm_3_methoxytyrosine_ivs$Beta, t2dm_3_methoxytyrosine_ivs$t2dm_Beta, t2dm_3_methoxytyrosine_ivs$SE, t2dm_3_methoxytyrosine_ivs$t2dm_SE) # Replicated
mr_ivw_mre(t2dm_3_methoxytyrosine_ivs$Beta, t2dm_3_methoxytyrosine_ivs$t2dm_Beta, t2dm_3_methoxytyrosine_ivs$SE, t2dm_3_methoxytyrosine_ivs$t2dm_SE) # Replicated
mr_penalised_weighted_median(t2dm_3_methoxytyrosine_ivs$Beta, t2dm_3_methoxytyrosine_ivs$t2dm_Beta, t2dm_3_methoxytyrosine_ivs$SE, t2dm_3_methoxytyrosine_ivs$t2dm_SE) # Replicated
mr_weighted_mode(t2dm_3_methoxytyrosine_ivs$Beta, t2dm_3_methoxytyrosine_ivs$t2dm_Beta, t2dm_3_methoxytyrosine_ivs$SE, t2dm_3_methoxytyrosine_ivs$t2dm_SE) # Replicated
mr_egger_regression(t2dm_3_methoxytyrosine_ivs$Beta, t2dm_3_methoxytyrosine_ivs$t2dm_Beta, t2dm_3_methoxytyrosine_ivs$SE, t2dm_3_methoxytyrosine_ivs$t2dm_SE) # Replicated
# Plot the MR results for 3-methoxytyrosine
dat <- data.frame(beta.exposure = t2dm_3_methoxytyrosine_ivs$Beta, se.exposure = t2dm_3_methoxytyrosine_ivs$SE,
                  beta.outcome = t2dm_3_methoxytyrosine_ivs$t2dm_Beta, se.outcome = t2dm_3_methoxytyrosine_ivs$t2dm_SE,
                  id.exposure = "3-methoxytyrosine", id.outcome = "T2DM", mr_keep = TRUE, exposure = "3-methoxytyrosine", outcome = "T2DM",
                  SNP = t2dm_3_methoxytyrosine_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_fe", "mr_ivw_mre", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_forest_plot(res_single)

# Load the T2DM IVs for the metabolite: 1-palmitoyl-2-oleoyl-GPE (16_0_18_1)
t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs <- read.table("Harmonised_T2DM_IVs/1-palmitoyl-2-oleoyl-GPE (16_0_18_1)T2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Collect the previous MR results for 1-palmitoyl-2-oleoyl-GPE (16_0_18_1) from the T2DM_MR_Results.tsv file
t2dm_mr_results[t2dm_mr_results$Metabolite == "1-palmitoyl-2-oleoyl-GPE (16_0_18_1)", ]
# Replicate previous MR analysis for 1-palmitoyl-2-oleoyl-GPE (16_0_18_1) using the TwoSampleMR package
mr_ivw_fe(t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$SE, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_SE) # Replicated
mr_ivw_mre(t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$SE, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_SE) # Replicated
mr_penalised_weighted_median(t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$SE, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_SE) # Replicated
mr_weighted_mode(t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$SE, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_SE) # Replicated
mr_egger_regression(t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_Beta, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$SE, t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_SE) # Replicated
# Plot the MR results for 1-palmitoyl-2-oleoyl-GPE (16_0_18_1)
dat <- data.frame(beta.exposure = t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$Beta, se.exposure = t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$SE, 
                  beta.outcome = t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_Beta, se.outcome = t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$t2dm_SE, 
                  id.exposure = "1-palmitoyl-2-oleoyl-GPE (16_0_18_1)", id.outcome = "T2DM", mr_keep = TRUE, exposure = "1-palmitoyl-2-oleoyl-GPE (16_0_18_1)", outcome = "T2DM",
                  SNP = t2dm_1_palmitoyl_2_oleoyl_gpe_16_0_18_1_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_fe", "mr_ivw_mre", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_forest_plot(res_single)

# Load the T2DM IVs for the metabolite: valine
t2dm_valine_ivs <- read.table("Harmonised_T2DM_IVs/valineT2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Collect the previous MR results for valine from the T2DM_MR_Results.tsv file
t2dm_mr_results[t2dm_mr_results$Metabolite == "valine", ]
# Replicate previous MR analysis for valine using the TwoSampleMR package
mr_ivw_fe(t2dm_valine_ivs$Beta, t2dm_valine_ivs$t2dm_Beta, t2dm_valine_ivs$SE, t2dm_valine_ivs$t2dm_SE) # Replicated
# Plot the MR results for valine
dat <- data.frame(beta.exposure = t2dm_valine_ivs$Beta, se.exposure = t2dm_valine_ivs$SE, 
                  beta.outcome = t2dm_valine_ivs$t2dm_Beta, se.outcome = t2dm_valine_ivs$t2dm_SE, 
                  id.exposure = "Valine", id.outcome = "T2DM", mr_keep = TRUE, exposure = "Valine", outcome = "T2DM",
                  SNP = t2dm_valine_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_ivw_mre", "mr_penalised_weighted_median", "mr_weighted_mode", "mr_egger_regression"))
mr_scatter_plot(res, dat)

# Load the T2DM IVs for the metabolite: sphinganine
t2dm_sphinganine_ivs <- read.table("Harmonised_T2DM_IVs/sphinganineT2DMT2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Collect the previous MR results for sphinganine from the T2DM_MR_Results.tsv file
t2dm_mr_results[t2dm_mr_results$Metabolite == "sphinganine", ]
# Replicate previous MR analysis for sphinganine using the TwoSampleMR package
mr_wald_ratio(t2dm_sphinganine_ivs$Beta, t2dm_sphinganine_ivs$t2dm_Beta, t2dm_sphinganine_ivs$SE, t2dm_sphinganine_ivs$t2dm_SE) # Replicated

# Load the HbA1c IVs for the metabolite: xanthine
hba1c_xanthine_ivs <- read.table("Harmonised_HbA1c_IVs/xanthineHBA1CHBA1C_harmonised_IVs.tsv", header = TRUE, sep = "\t")
# Collect the previous MR results for xanthine from the HbA1c_MR_Results.tsv file
hba1c_mr_results[hba1c_mr_results$Metabolite == "xanthine", ]
# Replicate previous MR analysis for xanthine using the TwoSampleMR package
mr_wald_ratio(hba1c_xanthine_ivs$Beta, hba1c_xanthine_ivs$hba1c_Beta, hba1c_xanthine_ivs$SE, hba1c_xanthine_ivs$hba1c_SE) # Replicated


# Load in the significant T2DM metabolite results
sig_T2DM_metabolites <- read.table("Filtered_T2DM_Results.tsv", header = TRUE, sep = "\t")
# In a new dataframe import the "Metabolite", "Number_of_IVs", "Fixed_IVW_Esitmate" and "Fixed_IVW_SE" columns
plotting_sig_T2DM_metabolites <- sig_T2DM_metabolites[, c("Metabolite", "Number_of_IVs", "Fixed_IVW_Estimate", "Fixed_IVW_SE")]
# If the number of IVs is greater than 3, take the "Random_IVW_SE" value from the "sig_T2DM_metabolites" dataframe as the "Fixed_IVW_SE" value
plotting_sig_T2DM_metabolites$Fixed_IVW_SE[plotting_sig_T2DM_metabolites$Number_of_IVs > 3] <- sig_T2DM_metabolites$Random_IVW_SE[plotting_sig_T2DM_metabolites$Number_of_IVs > 3]
# Rename the columns to "Metabolite", "Number_of_IVs", "Beta" and "SE"
colnames(plotting_sig_T2DM_metabolites) <- c("Metabolite", "Number_of_IVs", "Beta", "SE")
# Order the dataframe by the "Beta" column in descending order
plotting_sig_T2DM_metabolites <- plotting_sig_T2DM_metabolites[order(plotting_sig_T2DM_metabolites$Beta, decreasing = TRUE), ]
# Add columns to the dataframe for the upper and lower confidence intervals
plotting_sig_T2DM_metabolites$Lower_CI <- plotting_sig_T2DM_metabolites$Beta - (1.96 * plotting_sig_T2DM_metabolites$SE)
plotting_sig_T2DM_metabolites$Upper_CI <- plotting_sig_T2DM_metabolites$Beta + (1.96 * plotting_sig_T2DM_metabolites$SE)
# Add a column called "outcome" and set it to "T2DM"
plotting_sig_T2DM_metabolites$outcome <- "T2DM"
# Add a column called "method" and set it to "IVW"
plotting_sig_T2DM_metabolites$method <- "IVW"
# Identify the minimum and maximum of the confidence intervals
lower <- min(plotting_sig_T2DM_metabolites$Lower_CI)
upper <- max(plotting_sig_T2DM_metabolites$Upper_CI)

# Read in the responses_combined.tsv file
responses_combined <- read.table("responses_combined.tsv", header = TRUE, sep = "\t")

# Make a forest plot of the significant T2DM metabolites
library(metafor)
forest(
  plotting_sig_T2DM_metabolites$Beta, 
  sei=plotting_sig_T2DM_metabolites$SE, 
  slab=plotting_sig_T2DM_metabolites$Metabolite, 
       alim=c(-1, 1), # Adjust the axis limits to CI
       main="Standard Deviation Change in Metabolite Levels\nSignificantly Associated with Log Odds Ratio for T2DM",
       xlab="Log(OR)", 
       at=log(c(0.3, 0.4705, 0.605, 0.78, 1, 1.28, 1.65, 2.12, 5)), # Change scale if needed
       xlim=c(-1.5, 1), # Set limits for x-axis based on CIs
       cex=0.6, # Adjust text size as necessary
       psize=0.8, # Point size for the plot
       refline=0, # Reference line at effect size 0
       col="blue") # Color for points and lines

# Load the significant FG metabolite results
sig_FG_metabolites <- read.table("Filtered_FG_Results.tsv", header = TRUE, sep = "\t")
# In a new dataframe import the "Metabolite", "Number_of_IVs", "Fixed_IVW_Esitmate" and "Fixed_IVW_SE" columns
plotting_sig_FG_metabolites <- sig_FG_metabolites[, c("Metabolite", "Number_of_IVs", "Fixed_IVW_Estimate", "Fixed_IVW_SE")]
# If the number of IVs is greater than 3, take the "Random_IVW_SE" value from the "sig_FG_metabolites" dataframe as the "Fixed_IVW_SE" value
plotting_sig_FG_metabolites$Fixed_IVW_SE[plotting_sig_FG_metabolites$Number_of_IVs > 3] <- sig_FG_metabolites$Random_IVW_SE[plotting_sig_FG_metabolites$Number_of_IVs > 3]
# Rename the columns to "Metabolite", "Number_of_IVs", "Beta" and "SE"
colnames(plotting_sig_FG_metabolites) <- c("Metabolite", "Number_of_IVs", "Beta", "SE")
# Order the dataframe by the "Beta" column in descending order
plotting_sig_FG_metabolites <- plotting_sig_FG_metabolites[order(plotting_sig_FG_metabolites$Beta, decreasing = TRUE), ]
# Add columns to the dataframe for the upper and lower confidence intervals
plotting_sig_FG_metabolites$Lower_CI <- plotting_sig_FG_metabolites$Beta - (1.96 * plotting_sig_FG_metabolites$SE)
plotting_sig_FG_metabolites$Upper_CI <- plotting_sig_FG_metabolites$Beta + (1.96 * plotting_sig_FG_metabolites$SE)
# Add a column called "outcome" and set it to "FG"
plotting_sig_FG_metabolites$outcome <- "FG"
# Add a column called "method" and set it to "IVW"
plotting_sig_FG_metabolites$method <- "IVW"
# Identify the minimum and maximum of the confidence intervals
lower <- min(plotting_sig_FG_metabolites$Lower_CI)
upper <- max(plotting_sig_FG_metabolites$Upper_CI)

# Make a forest plot of the significant FG metabolites
forest(
  plotting_sig_FG_metabolites$Beta, 
  sei=plotting_sig_FG_metabolites$SE, 
  slab=plotting_sig_FG_metabolites$Metabolite, 
       alim=c(-1, 1), # Adjust the axis limits to CI
       main="Standard Deviation Change in Metabolite Levels\nSignificantly Associated with a Change in FG",
       xlab="Standard Deviation Change in FG", 
       at=log(c(0.3, 0.4705, 0.605, 0.78, 1, 1.28, 1.65, 2.12, 5)), # Change scale if needed
       xlim=c(-1.5, 1), # Set limits for x-axis based on CIs
       cex=0.6, # Adjust text size as necessary
       psize=0.8, # Point size for the plot
       refline=0, # Reference line at effect size 0
       col="blue") # Color for points and lines

# Load the significant HBA1C metabolite results
sig_HBA1C_metabolites <- read.table("Filtered_HBA1C_Results.tsv", header = TRUE, sep = "\t")
# In a new dataframe import the "Metabolite", "Number_of_IVs", "Fixed_IVW_Esitmate" and "Fixed_IVW_SE" columns
plotting_sig_HBA1C_metabolites <- sig_HBA1C_metabolites[, c("Metabolite", "Number_of_IVs", "Fixed_IVW_Estimate", "Fixed_IVW_SE")]
# Rename the columns to "Metabolite", "Number_of_IVs", "Beta" and "SE"
colnames(plotting_sig_HBA1C_metabolites) <- c("Metabolite", "Number_of_IVs", "Beta", "SE")
# Order the dataframe by the "Beta" column in descending order
plotting_sig_HBA1C_metabolites <- plotting_sig_HBA1C_metabolites[order(plotting_sig_HBA1C_metabolites$Beta, decreasing = TRUE), ]
# Add columns to the dataframe for the upper and lower confidence intervals
plotting_sig_HBA1C_metabolites$Lower_CI <- plotting_sig_HBA1C_metabolites$Beta - (1.96 * plotting_sig_HBA1C_metabolites$SE)
plotting_sig_HBA1C_metabolites$Upper_CI <- plotting_sig_HBA1C_metabolites$Beta + (1.96 * plotting_sig_HBA1C_metabolites$SE)
# Add a column called "outcome" and set it to "HBA1C"
plotting_sig_HBA1C_metabolites$outcome <- "HBA1C"
# Add a column called "method" and set it to "IVW"
plotting_sig_HBA1C_metabolites$method <- "IVW"
# Identify the minimum and maximum of the confidence intervals
lower <- min(plotting_sig_HBA1C_metabolites$Lower_CI)
upper <- max(plotting_sig_HBA1C_metabolites$Upper_CI)

# Plot the significant HBA1C metabolites
forest(
  plotting_sig_HBA1C_metabolites$Beta, 
  sei=plotting_sig_HBA1C_metabolites$SE, 
  slab=plotting_sig_HBA1C_metabolites$Metabolite, 
  alim=c(-1, 1), # Adjust the axis limits to CI
  main="Standard Deviation Change in Metabolite Levels\nSignificantly Associated with a Change in HbA1c",
  xlab="Standard Deviation Change HbA1c", 
  at=log(c(0.3, 0.4705, 0.605, 0.78, 1, 1.28, 1.65, 2.12, 5)), # Change scale if needed
  xlim=c(-1.5, 1.5), # Set limits for x-axis based on CIs
  cex=0.6, # Adjust text size as necessary
  psize=0.8, # Point size for the plot
  refline=0, # Reference line at effect size 0
  col="blue") # Color for points and lines

