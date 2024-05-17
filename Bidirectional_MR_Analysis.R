library(tidyverse)

# Load the data
reverse_t2dm <- readr::read_tsv("Reverse_MR/T2DM_Reverse_MR_Results.tsv")
reverse_fg <- readr::read_tsv("Reverse_MR/FG_Reverse_MR_Results.tsv")
reverse_hba1c <- readr::read_tsv("Reverse_MR/HBA1C_Reverse_MR_Results.tsv")
filtered_t2dm <- readr::read_tsv("Filtered_T2DM_Results.tsv")
filtered_fg <- readr::read_tsv("Filtered_FG_Results.tsv")
filtered_hba1c <- readr::read_tsv("Filtered_HbA1C_Results.tsv")
metab_compids <- readr::read_tsv("Sig_Metabolites_Compid.tsv")

# Convert the ":" and "/" in metab_compids$name to "_"
metab_compids$name <- gsub(":", "_", metab_compids$name)
metab_compids$name <- gsub("/", "_", metab_compids$name)

# Select the metabolites from reverse_t2dm that have a value of less than 0.05 in any of the following columns: 
# Weighted_Mode_Pval, Weighted_Median_Pval, Random_IVW_Pval, Egger_Pval
reverse_sig_t2dm <- reverse_t2dm %>% 
  filter(Random_IVW_Pval < 0.05 | Egger_Pval < 0.05)
# Do the same for reverse_fg and reverse_hba1c
reverse_sig_fg <- reverse_fg %>% 
  filter(Random_IVW_Pval < 0.05 | Egger_Pval < 0.05)
reverse_sig_hba1c <- reverse_hba1c %>%
  filter(Random_IVW_Pval < 0.05 | Egger_Pval < 0.05)
# Save a list of the unique metabolites from the combination of the three reverse dataframes
reverse_sig_metabs <- unique(c(reverse_sig_t2dm$Metabolite, reverse_sig_fg$Metabolite, reverse_sig_hba1c$Metabolite))

# Remove all but the first column from the filtered_t2dm dataframe
filtered_t2dm <- filtered_t2dm %>% 
  select(Metabolite)
# Remove all but the first column from the filtered_fg dataframe
filtered_fg <- filtered_fg %>% 
  select(Metabolite)
# Remove all but the first column from the filtered_hba1c dataframe
filtered_hba1c <- filtered_hba1c %>% 
  select(Metabolite)
# Add the compids from metab_compids to the filtered_t2dm dataframe
filtered_t2dm <- merge(filtered_t2dm, metab_compids, by.x = "Metabolite", by.y = "name", all.x = TRUE)
# Add the compids from metab_compids to the filtered_fg dataframe
filtered_fg <- merge(filtered_fg, metab_compids, by.x = "Metabolite", by.y = "name", all.x = TRUE)
# Add the compids from metab_compids to the filtered_hba1c dataframe
filtered_hba1c <- merge(filtered_hba1c, metab_compids, by.x = "Metabolite", by.y = "name", all.x = TRUE)
# Remove the Metabolites from filtered_t2dm that are in reverse_sig_metabs list
filtered_t2dm <- filtered_t2dm[!filtered_t2dm$compid %in% reverse_sig_metabs,]
# Remove the Metabolites from filtered_fg that are in reverse_sig_metabs list
filtered_fg <- filtered_fg[!filtered_fg$compid %in% reverse_sig_metabs,]
# Remove the Metabolites from filtered_hba1c that are in reverse_sig_metabs list
filtered_hba1c <- filtered_hba1c[!filtered_hba1c$compid %in% reverse_sig_metabs,]


# Save the reverse_sig_t2dm dataframe to a tsv file
write_tsv(reverse_sig_t2dm, "Reverse_MR/Reverse_Sig_T2DM_Results.tsv")
# Save the RV_filtered_t2dm dataframe to a tsv file
write_tsv(filtered_t2dm, "RV_Filtered_T2DM_Results.tsv")
# Save the reverse_sig_fg dataframe to a tsv file
write_tsv(reverse_sig_fg, "Reverse_MR/Reverse_Sig_FG_Results.tsv")
# Save the RV_filtered_fg dataframe to a tsv file
write_tsv(filtered_fg, "RV_Filtered_FG_Results.tsv")
# Save the reverse_sig_hba1c dataframe to a tsv file
write_tsv(reverse_sig_hba1c, "Reverse_MR/Reverse_Sig_HbA1C_Results.tsv")
# Save the RV_filtered_hba1c dataframe to a tsv file
write_tsv(filtered_hba1c, "RV_Filtered_HbA1C_Results.tsv")

# Check which metabolites are overlapping between the three filtered dataframes
overlap_metabs <- Reduce(intersect, list(filtered_t2dm$Metabolite, filtered_fg$Metabolite, filtered_hba1c$Metabolite))
overlap_metabs


### Plot some examples

# Load the T2DM IVs for the metabolite: M52466
t2dm_M52466_ivs <- read.table("Reverse_MR/Harmonised_Reverse_T2DM_IVs/M52466_T2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
dat <- data.frame(beta.exposure = t2dm_M52466_ivs$Beta, se.exposure = t2dm_M52466_ivs$SE,
                  beta.outcome = t2dm_M52466_ivs$out_Beta, se.outcome = t2dm_M52466_ivs$out_SE,
                  id.exposure = "T2DM", id.outcome = "1-stearoyl-2-docosahexaenoyl-GPE (18:0/22:6)*", mr_keep = TRUE, exposure = "T2DM", outcome = "1-stearoyl-2-docosahexaenoyl-GPE (18:0/22:6)*",
                  SNP = t2dm_M52466_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
mr_funnel_plot(res_single)

# Load the T2DM IVs for the metabolite: M52710
t2dm_M52710_ivs <- read.table("Reverse_MR/Harmonised_Reverse_T2DM_IVs/M52710_T2DM_harmonised_IVs.tsv", header = TRUE, sep = "\t")
dat <- data.frame(beta.exposure = t2dm_M52710_ivs$Beta, se.exposure = t2dm_M52710_ivs$SE,
                  beta.outcome = t2dm_M52710_ivs$out_Beta, se.outcome = t2dm_M52710_ivs$out_SE,
                  id.exposure = "T2DM", id.outcome = "1-linoleoyl-2-arachidonoyl-GPC (18:2/20:4n6)*", mr_keep = TRUE, exposure = "T2DM", outcome = "1-stearoyl-2-docosahexaenoyl-GPE (18:0/22:6)*",
                  SNP = t2dm_M52710_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
mr_funnel_plot(res_single)

# Load the FG IVs for the metabolite: M48407
fg_M48407_ivs <- read.table("Reverse_MR/Harmonised_Reverse_FG_IVs/M48407_FG_harmonised_IVs.tsv", header = TRUE, sep = "\t")
dat <- data.frame(beta.exposure = fg_M48407_ivs$Beta, se.exposure = fg_M48407_ivs$SE,
                  beta.outcome = fg_M48407_ivs$out_Beta, se.outcome = fg_M48407_ivs$out_SE,
                  id.exposure = "FG", id.outcome = "dopamine sulfate (2)", mr_keep = TRUE, exposure = "FG", outcome = "dopamine sulfate (2)",
                  SNP = fg_M48407_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
mr_funnel_plot(res_single)

# Load the FG IVs for the metabolite: M48047
fg_M48047_ivs <- read.table("Reverse_MR/Harmonised_Reverse_FG_IVs/M48047_FG_harmonised_IVs.tsv", header = TRUE, sep = "\t")
dat <- data.frame(beta.exposure = fg_M48047_ivs$Beta, se.exposure = fg_M48047_ivs$SE,
                  beta.outcome = fg_M48047_ivs$out_Beta, se.outcome = fg_M48047_ivs$out_SE,
                  id.exposure = "FG", id.outcome = "X - 18886", mr_keep = TRUE, exposure = "FG", outcome = "X - 18886",
                  SNP = fg_M48047_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
mr_funnel_plot(res_single)

# Load the HBA1C IVs for the metabolite: M03147
hba1c_M03147_ivs <- read.table("Reverse_MR/Harmonised_Reverse_HbA1C_IVs/M03147_HBA1C_harmonised_IVs.tsv", header = TRUE, sep = "\t")
dat <- data.frame(beta.exposure = hba1c_M03147_ivs$Beta, se.exposure = hba1c_M03147_ivs$SE,
                  beta.outcome = hba1c_M03147_ivs$out_Beta, se.outcome = hba1c_M03147_ivs$out_SE,
                  id.exposure = "HbA1C", id.outcome = "xanthine", mr_keep = TRUE, exposure = "HbA1C", outcome = "xanthine",
                  SNP = hba1c_M03147_ivs$SNP)
res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression"))
mr_scatter_plot(res, dat)
res_single <- mr_singlesnp(dat, all_method = c("mr_ivw_mre", "mr_egger_regression"))
mr_funnel_plot(res_single)














