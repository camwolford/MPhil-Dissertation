library(readr)
library(dplyr)

### Convert files back to tsv files ###

# Read in each of the GWAS summary statistics csv files
# Read in metabolites_cleaned.csv
metabolites <- read.csv("metabolites_cleaned.csv")
# Save this file as a tsv file
write.table(metabolites, "metabolites_cleaned.tsv", sep="\t", row.names=FALSE)
# Read in the metabolites_cleaned.tsv file
metabolites <- readr::read_tsv("metabolites_cleaned.tsv")

# Read in metabolite_gwas_associations_cleaned.csv
metabolite_gwas <- read.csv("metabolite_gwas_associations_cleaned.csv")
# Save this file as a tsv file
write.table(metabolite_gwas, "metabolite_gwas_associations_cleaned.tsv", sep="\t", row.names=FALSE)
# Read in the metabolite_gwas_associations_cleaned.tsv file
metabolite_gwas <- readr::read_tsv("metabolite_gwas_associations_cleaned.tsv")

# Read in t2dm_gwas_cleaned.csv
t2dm_gwas <- read.csv("t2dm_gwas_cleaned.csv")
# Save this file as a tsv file
write.table(t2dm_gwas, "t2dm_gwas_cleaned.tsv", sep="\t", row.names=FALSE)
# Read in the t2dm_gwas_cleaned.tsv file
t2dm_gwas <- readr::read_tsv("t2dm_gwas_cleaned.tsv")

# Read in fasting_glucose_gwas_cleaned.csv
fasting_glucose_gwas <- read.csv("fasting_glucose_gwas_cleaned.csv")
# Save this file as a tsv file
write.table(fasting_glucose_gwas, "fasting_glucose_gwas_cleaned.tsv", sep="\t", row.names=FALSE)
# Read in the fasting_glucose_gwas_cleaned.tsv file
fasting_glucose_gwas <- readr::read_tsv("fasting_glucose_gwas_cleaned.tsv")

# Read in hbA1c_gwas_cleaned.csv
hbA1c_gwas <- read.csv("hbA1c_gwas_cleaned.csv")
# Save this file as a tsv file
write.table(hbA1c_gwas, "hbA1c_gwas_cleaned.tsv", sep="\t", row.names=FALSE)
# Read in the hbA1c_gwas_cleaned.tsv file
hbA1c_gwas <- readr::read_tsv("hbA1c_gwas_cleaned.tsv")

# Read in responses_combined.csv
responses <- read.csv("responses_combined.csv")
# Save this file as a tsv file
write.table(responses, "responses_combined.tsv", sep="\t", row.names=FALSE)
# Read in the responses_combined.tsv file
responses <- readr::read_tsv("responses_combined.tsv")


### Harmonise the GWAS summary statistics ###

# Rename the fasting_glucose_gwas columns to match the metabolite_gwas columns
colnames(fasting_glucose_gwas) <- c("Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval", "SampleSize")
# Save the fasting_glucose_gwas as a tsv file
write.table(fasting_glucose_gwas, "fasting_glucose_gwas_cleaned.tsv", sep="\t", row.names=FALSE)
# Rename the hbA1c_gwas columns to match the metabolite_gwas columns
colnames(hbA1c_gwas) <- c("Chromosome", "Position", "EffectAllele", "NonEffectAllele", "EAF", "Beta", "SE", "Pval", "SampleSize")
# Save the hbA1c_gwas as a tsv file
write.table(hbA1c_gwas, "hbA1c_gwas_cleaned.tsv", sep="\t", row.names=FALSE)

# For each unique metabolite in metabolite_gwas make a new dataframe
# Test on the first metabolite
metabolite <- metabolite_gwas$Metabolite[1]
# Make a new dataframe with the rows that contain the metabolite
metabolite_gwas_metabolite <- metabolite_gwas[metabolite_gwas$Metabolite == metabolite,]

unique_metabolites <- unique(metabolite_gwas$Metabolite)

# Initialize a counter for progress reporting
counter <- 1

for(i in seq_along(unique_metabolites)) {
  metabolite <- unique_metabolites[i]

  message(sprintf("Processing %s (%d/%d)...", metabolite, counter, length(unique_metabolites)))
  
  metabolite_gwas_metabolite <- metabolite_gwas %>%
    dplyr::filter(Metabolite == metabolite)
  
  # Construct file path using metabolite
  file_path <- paste0("Individual_Metabolite_GWAS/", metabolite, "_GWAS.tsv")
  
  # Save to TSV, ensure the directory exists or add code to create it
  write.table(metabolite_gwas_metabolite, file = file_path,
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  message(sprintf("Finished processing %s. File saved to %s", metabolite, file_path))
  
  # Update the counter
  counter <- counter + 1
}


# Load the data from the saved files
test_load <- readr::read_tsv("Individual_Metabolite_GWAS/1-(1-enyl-oleoyl)-GPC (P-18_1)*_GWAS.tsv")

# what are the types of the columns
str(test_load)




# how many P.VALUE are less than 0.05x10^-8
sum(metabolite_gwas$Pval < 5e-8)

# what are the top 10 unique counts of the SNPS column
head(sort(table(metabolite_gwas$SNP), decreasing = TRUE), 20)

# Make a dataframe of the rsids with the most counts
repeated_snps <- data.frame(SNP = names(sort(table(metabolite_gwas$SNP), decreasing = TRUE)),
                            Count = as.numeric(sort(table(metabolite_gwas$SNP), decreasing = TRUE)))
# Save the dataframe to a TSV file
write.table(repeated_snps, file = "Repeated_SNPS.tsv", sep = "\t", row.names = FALSE, quote = FALSE)






