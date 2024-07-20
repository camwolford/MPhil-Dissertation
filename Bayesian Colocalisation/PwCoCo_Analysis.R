### Analyse the PwCoCo Results ###
# This script was used to analyse the results and filter the PwCoCo analysis.
library(tidyverse)

# Load in T2DM_Coloc_IVs.tsv
T2DM_Coloc_IVs <- read_tsv("Colocalisation/T2DM_Coloc_IVs.tsv")

# Remove rows with num_IVs < 1
T2DM_Coloc_IVs <- T2DM_Coloc_IVs %>% filter(num_IVs >= 1)

# Read in the pwococo_out_t2dm.coloc file
t2dm_coloc_results <- read_delim("Colocalisation/pwcoco_out_t2dm.coloc")
# Select the unconditioned results
t2dm_coloc_results_unconditioned <- t2dm_coloc_results %>% filter(SNP1 == "unconditioned" & SNP2 == "unconditioned") 
# Check the statistics of the nsnps column
summary(t2dm_coloc_results_unconditioned$nsnps)
# Select the metabolites from the unconditioned results that have an H4 probability > 0.8
t2dm_coloc_results_unconditioned_high_h4 <- t2dm_coloc_results_unconditioned %>% filter(H4 >= 0.8)

# Make a new dataframe called metabolite_t2dm_coloc_pass that contains the metabolites that are in the t2dm_coloc_results_unconditioned_high_h4$Dataset1 column
# the same number of times as the number of IVs in the T2DM_Coloc_IVs$num_IVs column
metabolite_t2dm_coloc_pass <- T2DM_Coloc_IVs 
# In a new column, count the number of times each metabolite is contained in t2dm_coloc_results_unconditioned_high_h4$Dataset1
metabolite_t2dm_coloc_pass$num_unconditioned_IVs <- map_int(
  metabolite_t2dm_coloc_pass$Metabolite,
  ~ sum(str_detect(t2dm_coloc_results_unconditioned_high_h4$Dataset1, fixed(.x)))
)

# In a new column, define a TRUE/FALSE value for whether the number of times each metabolite is contained in t2dm_coloc_results_unconditioned_high_h4$Dataset1
# is equal to the number of IVs in the T2DM_Coloc_IVs$num_IVs column
metabolite_t2dm_coloc_pass$unconditional_pass <- metabolite_t2dm_coloc_pass$num_unconditioned_IVs == metabolite_t2dm_coloc_pass$num_IVs

# Select the rows of Dataset1 that have low H4 probabilities
t2dm_coloc_results_unconditioned_low_h4 <- t2dm_coloc_results_unconditioned %>% filter(H4 < 0.8)

# Add a new column as TRUE/FALSE 
t2dm_coloc_results_unconditioned_low_h4$conditioned_pass <- FALSE

# Select the conditioned results
t2dm_coloc_results_conditioned <- t2dm_coloc_results %>% filter(SNP1 != "unconditioned" | SNP2 != "unconditioned")

# Remove rows with nsnps < 100
t2dm_coloc_results_conditioned <- t2dm_coloc_results_conditioned %>% filter(nsnps >= 100)

# For each row in t2dm_coloc_results_unconditioned_low_h4
for (i in 1:nrow(t2dm_coloc_results_unconditioned_low_h4)) {
  # Select the rows in t2dm_coloc_results_conditioned that have the same Dataset1 values as the current row in t2dm_coloc_results_unconditioned_low_h4
  conditioned_rows <- t2dm_coloc_results_conditioned %>% filter(Dataset1 == t2dm_coloc_results_unconditioned_low_h4$Dataset1[i])
  # If any of the rows in conditioned_rows have H4 > 0.8, then set conditioned_pass to TRUE
  t2dm_coloc_results_unconditioned_low_h4$conditioned_pass[i] <- any(conditioned_rows$H4 > 0.8)
}

# Select the metabolites from the unconditioned results that have a conditioned_pass value of TRUE
true_conditioned_pass <- t2dm_coloc_results_unconditioned_low_h4 %>% filter(conditioned_pass == TRUE)
true_conditioned_pass_025 <- t2dm_coloc_results_unconditioned_low_h4_025 %>% filter(conditioned_pass == TRUE)

# Add a conditional_pass column to the metabolite_t2dm_coloc_pass dataframe
metabolite_t2dm_coloc_pass$conditional_pass <- FALSE
metabolite_t2dm_coloc_pass$conditional_pass_025 <- FALSE
# If the difference between the num_IVs and num_unconditioned_IVs is equal to the number of times the metabolite is in true_conditioned_pass$Dataset1, then set conditional_pass to TRUE
metabolite_t2dm_coloc_pass <- metabolite_t2dm_coloc_pass %>%
  rowwise() %>%
  mutate(conditional_pass = sum(str_detect(true_conditioned_pass$Dataset1, fixed(Metabolite))) == num_IVs - num_unconditioned_IVs) %>%
  ungroup()
metabolite_t2dm_coloc_pass <- metabolite_t2dm_coloc_pass %>%
  rowwise() %>%
  mutate(conditional_pass_025 = sum(str_detect(true_conditioned_pass_025$Dataset1, fixed(Metabolite))) == num_IVs - num_unconditioned_IVs_025) %>%
  ungroup()
# If unconditioned_pass is TRUE or H4_greater is TRUE, then set conditional_pass to TRUE
metabolite_t2dm_coloc_pass$conditional_pass[metabolite_t2dm_coloc_pass$unconditional_pass] <- TRUE
metabolite_t2dm_coloc_pass$conditional_pass[metabolite_t2dm_coloc_pass$H4_greater] <- TRUE

# Count the number of metabolites that pass the unconditional and conditional tests
sum(metabolite_t2dm_coloc_pass$unconditional_pass | metabolite_t2dm_coloc_pass$conditional_pass)

# Save the metabolite_t2dm_coloc_pass dataframe to a file
write_tsv(metabolite_t2dm_coloc_pass, "Colocalisation/metabolite_t2dm_coloc_pass.tsv")


# Do all of the same for FG
# Load in FG_Coloc_IVs.tsv
FG_Coloc_IVs <- read_tsv("Colocalisation/FG_Coloc_IVs.tsv")

# Remove rows with num_IVs < 1
FG_Coloc_IVs <- FG_Coloc_IVs %>% filter(num_IVs >= 1)

# Read in the pwococo_out_fg.coloc file
fg_coloc_results <- read_delim("Colocalisation/pwcoco_out_fg.coloc")
# Select the unconditioned results
fg_coloc_results_unconditioned <- fg_coloc_results %>% filter(SNP1 == "unconditioned" & SNP2 == "unconditioned") 
# Check the statistics of the nsnps column
summary(fg_coloc_results_unconditioned$nsnps)
# Select the metabolites from the unconditioned results that have an H4 probability > 0.8
fg_coloc_results_unconditioned_high_h4 <- fg_coloc_results_unconditioned %>% filter(H4 >= 0.8)

# Make a new dataframe called metabolite_fg_coloc_pass that contains the metabolites that are in the fg_coloc_results_unconditioned_high_h4$Dataset1 column
# the same number of times as the number of IVs in the FG_Coloc_IVs$num_IVs column
metabolite_fg_coloc_pass <- FG_Coloc_IVs 
# In a new column, count the number of times each metabolite is contained in fg_coloc_results_unconditioned_high_h4$Dataset1
metabolite_fg_coloc_pass$num_unconditioned_IVs <- map_int(
  metabolite_fg_coloc_pass$Metabolite,
  ~ sum(str_detect(fg_coloc_results_unconditioned_high_h4$Dataset1, fixed(.x)))
)
# In a new column, define a TRUE/FALSE value for whether the number of times each metabolite is contained in fg_coloc_results_unconditioned_high_h4$Dataset1
# is equal to the number of IVs in the FG_Coloc_IVs$num_IVs column
metabolite_fg_coloc_pass$unconditional_pass <- metabolite_fg_coloc_pass$num_unconditioned_IVs == metabolite_fg_coloc_pass$num_IVs

# Select the rows of Dataset1 that have low H4 probabilities
fg_coloc_results_unconditioned_low_h4 <- fg_coloc_results_unconditioned %>% filter(H4 < 0.8)

# Add a new column as TRUE/FALSE 
fg_coloc_results_unconditioned_low_h4$conditioned_pass <- FALSE

# Select the conditioned results
fg_coloc_results_conditioned <- fg_coloc_results %>% filter(SNP1 != "unconditioned" | SNP2 != "unconditioned")
# Remove rows with nsnps < 100
fg_coloc_results_conditioned <- fg_coloc_results_conditioned %>% filter(nsnps >= 100)

# For each row in fg_coloc_results_unconditioned_low_h4
for (i in 1:nrow(fg_coloc_results_unconditioned_low_h4)) {
  # Select the rows in fg_coloc_results_conditioned that have the same Dataset1 values as the current row in fg_coloc_results_unconditioned_low_h4
  conditioned_rows <- fg_coloc_results_conditioned %>% filter(Dataset1 == fg_coloc_results_unconditioned_low_h4$Dataset1[i])
  # If any of the rows in conditioned_rows have H4 > 0.8, then set conditioned_pass to TRUE
  fg_coloc_results_unconditioned_low_h4$conditioned_pass[i] <- any(conditioned_rows$H4 > 0.8)
}
# Select the metabolites from the unconditioned results that have a conditioned_pass value of TRUE
true_conditioned_pass <- fg_coloc_results_unconditioned_low_h4 %>% filter(conditioned_pass == TRUE)

# Add a conditional_pass column to the metabolite_fg_coloc_pass dataframe
metabolite_fg_coloc_pass$conditional_pass <- FALSE
# If the difference between the num_IVs and num_unconditioned_IVs is equal to the number of times the metabolite is in true_conditioned_pass$Dataset1, then set conditional_pass to TRUE
metabolite_fg_coloc_pass <- metabolite_fg_coloc_pass %>%
  rowwise() %>%
  mutate(conditional_pass = sum(str_detect(true_conditioned_pass$Dataset1, fixed(Metabolite))) == num_IVs - num_unconditioned_IVs) %>%
  ungroup()
# If unconditioned_pass is TRUE or H4_greater is TRUE, then set conditional_pass to TRUE
metabolite_fg_coloc_pass$conditional_pass[metabolite_fg_coloc_pass$unconditional_pass] <- TRUE

# Count the number of metabolites that pass the unconditional and conditional tests
sum(metabolite_fg_coloc_pass$unconditional_pass | metabolite_fg_coloc_pass$conditional_pass)

# Save the metabolite_fg_coloc_pass dataframe to a file
write_tsv(metabolite_fg_coloc_pass, "Colocalisation/metabolite_fg_coloc_pass.tsv")


# Do the same for HBA1C
# Load in HBA1C_Coloc_IVs.tsv
HBA1C_Coloc_IVs <- read_tsv("Colocalisation/HBA1C_Coloc_IVs.tsv")

# Remove rows with num_IVs < 1
HBA1C_Coloc_IVs <- HBA1C_Coloc_IVs %>% filter(num_IVs >= 1)

# Read in the pwococo_out_hba1c.coloc file
hba1c_coloc_results <- read_delim("Colocalisation/pwcoco_out_hba1c.coloc")
# Select the unconditioned results
hba1c_coloc_results_unconditioned <- hba1c_coloc_results %>% filter(SNP1 == "unconditioned" & SNP2 == "unconditioned") 
# Check the statistics of the nsnps column
summary(hba1c_coloc_results_unconditioned$nsnps)
# Select the metabolites from the unconditioned results that have an H4 probability > 0.8
hba1c_coloc_results_unconditioned_high_h4 <- hba1c_coloc_results_unconditioned %>% filter(H4 >= 0.8)

# Make a new dataframe called metabolite_hba1c_coloc_pass that contains the metabolites that are in the hba1c_coloc_results_unconditioned_high_h4$Dataset1 column
# the same number of times as the number of IVs in the HBA1C_Coloc_IVs$num_IVs column
metabolite_hba1c_coloc_pass <- HBA1C_Coloc_IVs 
# In a new column, count the number of times each metabolite is contained in hba1c_coloc_results_unconditioned_high_h4$Dataset1
metabolite_hba1c_coloc_pass$num_unconditioned_IVs <- map_int(
  metabolite_hba1c_coloc_pass$Metabolite,
  ~ sum(str_detect(hba1c_coloc_results_unconditioned_high_h4$Dataset1, fixed(.x)))
)
# In a new column, define a TRUE/FALSE value for whether the number of times each metabolite is contained in hba1c_coloc_results_unconditioned_high_h4$Dataset1
# is equal to the number of IVs in the HBA1C_Coloc_IVs$num_IVs column
metabolite_hba1c_coloc_pass$unconditional_pass <- metabolite_hba1c_coloc_pass$num_unconditioned_IVs == metabolite_hba1c_coloc_pass$num_IVs

# Select the rows of Dataset1 that have low H4 probabilities
hba1c_coloc_results_unconditioned_low_h4 <- hba1c_coloc_results_unconditioned %>% filter(H4 < 0.8)

# Add a new column as TRUE/FALSE 
hba1c_coloc_results_unconditioned_low_h4$conditioned_pass <- FALSE

# Select the conditioned results
hba1c_coloc_results_conditioned <- hba1c_coloc_results %>% filter(SNP1 != "unconditioned" | SNP2 != "unconditioned")
# Remove rows with nsnps < 100
hba1c_coloc_results_conditioned <- hba1c_coloc_results_conditioned %>% filter(nsnps >= 100)

# For each row in hba1c_coloc_results_unconditioned_low_h4
for (i in 1:nrow(hba1c_coloc_results_unconditioned_low_h4)) {
  # Select the rows in hba1c_coloc_results_conditioned that have the same Dataset1 values as the current row in hba1c_coloc_results_unconditioned_low_h4
  conditioned_rows <- hba1c_coloc_results_conditioned %>% filter(Dataset1 == hba1c_coloc_results_unconditioned_low_h4$Dataset1[i])
  # If any of the rows in conditioned_rows have H4 > 0.8, then set conditioned_pass to TRUE
  hba1c_coloc_results_unconditioned_low_h4$conditioned_pass[i] <- any(conditioned_rows$H4 > 0.8)
}
# Select the metabolites from the unconditioned results that have a conditioned_pass value of TRUE
true_conditioned_pass <- hba1c_coloc_results_unconditioned_low_h4 %>% filter(conditioned_pass == TRUE)

# Add a conditional_pass column to the metabolite_hba1c_coloc_pass dataframe
metabolite_hba1c_coloc_pass$conditional_pass <- FALSE
# If the difference between the num_IVs and num_unconditioned_IVs is equal to the number of times the metabolite is in true_conditioned_pass$Dataset1, then set conditional_pass to TRUE
metabolite_hba1c_coloc_pass <- metabolite_hba1c_coloc_pass %>%
  rowwise() %>%
  mutate(conditional_pass = sum(str_detect(true_conditioned_pass$Dataset1, fixed(Metabolite))) == num_IVs - num_unconditioned_IVs) %>%
  ungroup()
# If unconditioned_pass is TRUE or H4_greater is TRUE, then set conditional_pass to TRUE
metabolite_hba1c_coloc_pass$conditional_pass[metabolite_hba1c_coloc_pass$unconditional_pass] <- TRUE

# Count the number of metabolites that pass the unconditional and conditional tests
sum(metabolite_hba1c_coloc_pass$unconditional_pass | metabolite_hba1c_coloc_pass$conditional_pass)

# Save the metabolite_hba1c_coloc_pass dataframe to a file
write_tsv(metabolite_hba1c_coloc_pass, "Colocalisation/metabolite_hba1c_coloc_pass.tsv")