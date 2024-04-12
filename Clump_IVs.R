### Clump IVs ###
# This script is used to clump IVs based on the LD structure of the data.

# Test this on the first metabolite
# Load the data
metabolite_3 <- readr::read_tsv("Filtered_IVs/1-(1-enyl-oleoyl)-GPC (P-18_1)*_IVs.tsv")

# Clump the IVs from the third metabolite
library(ieugwasr)
library(devtools)

devtools::install_github("explodecomputer/genetics.binaRies")
genetics.binaRies::get_plink_binary()

colnames(metabolite_3)[4] = "rsid"
colnames(metabolite_3)[11] = "pval"

clumped_IVs <- ld_clump_local(metabolite_3,  clump_kb = 10000,
                        clump_r2 = 0.001,
                        clump_p = 1,
                        bfile = "1kg.v3/EUR",
                        plink_bin = genetics.binaRies::get_plink_binary())

# Rename the rsid column to SNP
colnames(clumped_IVs)[4] <- "SNP"
# Rename the pval column to Pval
colnames(clumped_IVs)[11] <- "Pval"


### Perform the clumping on all metabolites ###
# Load in each metabolite data in /Filtered_IVs
metabolite_files <- list.files("Filtered_IVs", full.names = TRUE)
# Initialise a counter for the number of metabolites
metabolite_counter <- 0
metabolites_with_more_than_1_IVs <- 0
metabolites_with_1_IVs <- 0
metabolites_with_0_IVs <- 0


for (metabolite_file in metabolite_files) {
  # Load the metabolite data
  file_path <- metabolite_file
  metabolite_data <- readr::read_tsv(file_path)
  # Check if any EffectAllele is "D" or "I"
  if (any(metabolite_data$EffectAllele == "D" | metabolite_data$EffectAllele == "I")) {
    # Remove that row
    metabolite_data <- metabolite_data[metabolite_data$EffectAllele != "D" & metabolite_data$EffectAllele != "I", ]
  }
  # Check that the metabolite data has at least 1 row
  if (nrow(metabolite_data) == 0) {
    metabolites_with_0_IVs <- metabolites_with_0_IVs + 1
    next
  }
  # Check the number of IVs
  if (nrow(metabolite_data) == 1) {
    # No need to clump
    metabolites_with_1_IVs <- metabolites_with_1_IVs + 1
    # Extract the part of the file name before the "_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_IVs")[[1]][1]
    # Add "_Clumped_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "_Clumped_IVs")
    # Save the data
    write.table(metabolite_data, paste0("Clumped_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  } else {
    # Clump the IVs
    # Convert metabolite_data to a data frame called dat with the following columns: Metabolite, Chromosome, Position, EffectAllele, rsid, EAF, pval, Beta, LowCI, UpCI, SE
    colnames(metabolite_data)[4] = "rsid"
    colnames(metabolite_data)[11] = "pval"
    clumped_IVs <- ld_clump_local(metabolite_data,  clump_kb = 10000,
                                  clump_r2 = 0.001,
                                  clump_p = 1,
                                  bfile = "1kg.v3/EUR",
                                  plink_bin = genetics.binaRies::get_plink_binary())
    # Rename the rsid column to SNP
    colnames(clumped_IVs)[4] <- "SNP"
    # Rename the pval column to Pval
    colnames(clumped_IVs)[11] <- "Pval"
    # Check the number of IVs
    if (nrow(clumped_IVs) == 1) {
      metabolites_with_1_IVs <- metabolites_with_1_IVs + 1
    } else {
      metabolites_with_more_than_1_IVs <- metabolites_with_more_than_1_IVs + 1
    }
    # Extract the part of the file name before the "_IVs"
    metabolite_name <- strsplit(basename(metabolite_file), "_IVs")[[1]][1]
    # Add "_Clumped_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "_Clumped_IVs")
    # Save the data
    write.table(clumped_IVs, paste0("Clumped_IVs/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  }
  print(paste("Processed", metabolite_counter, "metabolites"))
  print(paste("Metabolites with more than 1 IVs:", metabolites_with_more_than_1_IVs))
  print(paste("Metabolites with 1 IV:", metabolites_with_1_IVs))
  print(paste("Metabolites with 0 IVs:", metabolites_with_0_IVs))
}







