### Clump IVs Liberally ###
# This script is used to clump IVs based on the LD structure of the data.

# Clump the IVs from the third metabolite
library(ieugwasr)
library(devtools)

devtools::install_github("explodecomputer/genetics.binaRies")
genetics.binaRies::get_plink_binary()

### Perform the clumping on all metabolites ###
# Load in each metabolite data in /Filtered_IVs
metabolite_files <- list.files("Liberal_Analysis/Filtered_IVs_Liberal", full.names = TRUE)
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
    metabolite_name <- strsplit(basename(metabolite_file), "_IVs_Liberal")[[1]][1]
    # Add "_Clumped_IVs" to the file name
    metabolite_name <- paste0(metabolite_name, "_Clumped_IVs_Liberal")
    # Save the data
    write.table(metabolite_data, paste0("Liberal_Analysis/Clumped_IVs_Liberal/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
    metabolite_counter <- metabolite_counter + 1
  } else {
    # Attempt to clump the IVs
    tryCatch({
      colnames(metabolite_data)[4] = "rsid"
      colnames(metabolite_data)[11] = "pval"
      clumped_IVs <- ld_clump_local(metabolite_data, clump_kb = 10000,
                                    clump_r2 = 0.01, clump_p = 1,
                                    bfile = "1kg.v3/EUR",
                                    plink_bin = genetics.binaRies::get_plink_binary())
      # Check if clumped_IVs has rows
      if (nrow(clumped_IVs) > 0) {
        colnames(clumped_IVs)[4] <- "SNP"
        colnames(clumped_IVs)[11] <- "Pval"
        # Further processing or saving the clumped_IVs as needed
        # Example: save or analyze clumped_IVs
      }
      # Check the number of IVs
      if (nrow(clumped_IVs) == 1) {
        metabolites_with_1_IVs <- metabolites_with_1_IVs + 1
      } else {
        metabolites_with_more_than_1_IVs <- metabolites_with_more_than_1_IVs + 1
      }
      # Extract the part of the file name before the "_IVs"
      metabolite_name <- strsplit(basename(metabolite_file), "_IVs_Liberal")[[1]][1]
      # Add "_Clumped_IVs" to the file name
      metabolite_name <- paste0(metabolite_name, "_Clumped_IVs_Liberal")
      # Save the data
      write.table(clumped_IVs, paste0("Liberal_Analysis/Clumped_IVs_Liberal/", metabolite_name, ".tsv"), sep = "\t", row.names = FALSE)
      metabolite_counter <- metabolite_counter + 1
    }, error = function(e) {
      cat("Error with clumping in file", metabolite_file, ":", e$message, "\n")
    })
  }
  print(paste("Processed", metabolite_counter, "metabolites"))
  print(paste("Metabolites with more than 1 IVs:", metabolites_with_more_than_1_IVs))
  print(paste("Metabolites with 1 IV:", metabolites_with_1_IVs))
  print(paste("Metabolites with 0 IVs:", metabolites_with_0_IVs))
}







