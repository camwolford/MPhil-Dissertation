library(MendelianRandomization)
library(tidyverse)
library(R.utils)

BiocManager::install("VariantAnnotation")
remotes::install_github("mrcieu/gwasvcf")
gwasvcf::set_plink(genetics.binaRies::get_plink_binary())

# Load the coloc_filtered_metabolites_data.tsv file
metabolites_data_frame <- readr::read_tsv("coloc_filtered_metabolites_data.tsv")
# Group the metabolites by subpathway
subpathways <- unique(metabolites_data_frame$subpathway)
# Add an "M" to each metabolite compid if there isnt one
metabolites_data_frame$compid <- ifelse(grepl("^M", metabolites_data_frame$compid), metabolites_data_frame$compid, paste0("M", metabolites_data_frame$compid))

# Load the ukbb_variants_data.tsv file
ukbb_variants_data_frame <- readr::read_tsv("Mediation/ukbb_variants_data.tsv")
# Does rsid... exist in ukbb_variants_data_frame$rsid?
variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == "rs6857"]
# Does rsid... exist in haem_data_frame$rsid?
haem_data_frame$variant[haem_data_frame$variant == variant]
haem_data_frame[haem_data_frame$variant == variant, ]

all_possible_variants <- readr::read_tsv("Mediation/all_possible_variants.tsv")

# Change the columns names to variant and rsid
colnames(all_possible_variants) <- c("variant", "rsid")
# Add a row called manual_proxy
all_possible_variants$manual_proxy <- NA
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs483082"] <- "rs438811"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs73004967"] <- "rs17217098"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs174546"] <- "rs174553"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs35853021"] <- "rs261291"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs2414578"] <- "rs2414577"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs28456"] <- "rs174560"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs498936"] <- "rs1219550"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs2304128"] <- "rs12608729"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs2070895"] <- "rs1077834"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs422137"] <- "rs5812935"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs3859862"] <- "rs2006227"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs59771628"] <- "rs62276527"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs174545"] <- "rs174553"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs4144185"] <- "rs9370162"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs6968554"] <- "rs2106727"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs35246381"] <- "rs1495741"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs62492368"] <- "rs7794796"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs1868361"] <- "rs2029056"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs4921913"] <- "rs1495741"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs676457"] <- "rs545971"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs673335"] <- "rs498936"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs174574"] <- "rs174578"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs102275"] <- "rs99780"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs174577"] <- "rs174581"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs174547"] <- "rs174551"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs174544"] <- "rs174549"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs4633"] <- "rs165722"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs2727271"] <- "rs2727270"
all_possible_variants$manual_proxy[all_possible_variants$rsid == "rs6857"] <- NA

# Save all_possible_variants
write.table(all_possible_variants, "Mediation/all_possible_variants.tsv", quote = FALSE, row.names = FALSE, sep = "\t")



overall_betaX <- data.frame(all_possible_variants)
overall_betaX$heam <- NA
overall_betaX$bmi <- NA
overall_betaX$waist <- NA
overall_betaX$body_fat <- NA
overall_betaX$alchohol <- NA
overall_betaX$smoke <- NA
overall_betaX$mpv <- NA
overall_betaXse <- data.frame(all_possible_variants)
overall_betaXse$heam <- NA
overall_betaXse$bmi <- NA
overall_betaXse$waist <- NA
overall_betaXse$body_fat <- NA
overall_betaXse$alchohol <- NA
overall_betaXse$smoke <- NA
overall_betaXse$mpv <- NA

# Get each risk factor gwas summary statistics for each of the unique variants for this subpathway
# Load the haem_data.tsv file
haem_data_frame <- readr::read_tsv("Mediation/haem_data.tsv")
# Add the beta values from the haem_data_frame to the betaX$heam column
for (variant in all_possible_variants$variant) {
  if (variant %in% haem_data_frame$variant) {
    overall_betaX$heam[overall_betaX$variant == variant] <- haem_data_frame$beta[haem_data_frame$variant == variant]
    overall_betaXse$heam[overall_betaXse$variant == variant] <- haem_data_frame$se[haem_data_frame$variant == variant]
  } else {
    overall_betaX$heam[overall_betaX$variant == variant] <- NA
    overall_betaXse$heam[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the haem_data_frame
rm(haem_data_frame)

# Load the bmi_data.tsv file
bmi_data_frame <- readr::read_tsv("Mediation/bmi_data.tsv")
# Add the beta values from the bmi_data_frame to the betaX$bmi column
for (variant in all_possible_variants$variant) {
  if (variant %in% bmi_data_frame$variant) {
    overall_betaX$bmi[overall_betaX$variant == variant] <- bmi_data_frame$beta[bmi_data_frame$variant == variant]
    overall_betaXse$bmi[overall_betaXse$variant == variant] <- bmi_data_frame$se[bmi_data_frame$variant == variant]
  } else {
    overall_betaX$bmi[overall_betaX$variant == variant] <- NA
    overall_betaXse$bmi[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the bmi_data_frame
rm(bmi_data_frame)

# Load the waist_data.tsv file
waist_data_frame <- readr::read_tsv("Mediation/waist_data.tsv")
# Add the beta values from the waist_data_frame to the betaX$waist column
for (variant in all_possible_variants$variant) {
  if (variant %in% waist_data_frame$variant) {
    overall_betaX$waist[overall_betaX$variant == variant] <- waist_data_frame$beta[waist_data_frame$variant == variant]
    overall_betaXse$waist[overall_betaXse$variant == variant] <- waist_data_frame$se[waist_data_frame$variant == variant]
  } else {
    overall_betaX$waist[overall_betaX$variant == variant] <- NA
    overall_betaXse$waist[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the waist_data_frame
rm(waist_data_frame)

# Load the body_fat_data.tsv file
body_fat_data_frame <- readr::read_tsv("Mediation/body_fat_data.tsv")
# Add the beta values from the body_fat_data_frame to the betaX$body_fat column
for (variant in all_possible_variants$variant) {
  if (variant %in% body_fat_data_frame$variant) {
    overall_betaX$body_fat[overall_betaX$variant == variant] <- body_fat_data_frame$beta[body_fat_data_frame$variant == variant]
    overall_betaXse$body_fat[overall_betaXse$variant == variant] <- body_fat_data_frame$se[body_fat_data_frame$variant == variant]
  } else {
    overall_betaX$body_fat[overall_betaX$variant == variant] <- NA
    overall_betaXse$body_fat[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the body_fat_data_frame
rm(body_fat_data_frame)

# Load the alchohol_data.tsv file
alchohol_data_frame <- readr::read_tsv("Mediation/alcohol_data.tsv")
# Add the beta values from the alchohol_data_frame to the betaX$alchohol column
for (variant in all_possible_variants$variant) {
  if (variant %in% alchohol_data_frame$variant) {
    overall_betaX$alchohol[overall_betaX$variant == variant] <- alchohol_data_frame$beta[alchohol_data_frame$variant == variant]
    overall_betaXse$alchohol[overall_betaXse$variant == variant] <- alchohol_data_frame$se[alchohol_data_frame$variant == variant]
  } else {
    overall_betaX$alchohol[overall_betaX$variant == variant] <- NA
    overall_betaXse$alchohol[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the alchohol_data_frame
rm(alchohol_data_frame)

# Load the smoke_data.tsv file
smoke_data_frame <- readr::read_tsv("Mediation/smoke_data.tsv")
# Add the beta values from the smoke_data_frame to the betaX$smoke column
for (variant in all_possible_variants$variant) {
  if (variant %in% smoke_data_frame$variant) {
    overall_betaX$smoke[overall_betaX$variant == variant] <- smoke_data_frame$beta[smoke_data_frame$variant == variant]
    overall_betaXse$smoke[overall_betaXse$variant == variant] <- smoke_data_frame$se[smoke_data_frame$variant == variant]
  } else {
    overall_betaX$smoke[overall_betaX$variant == variant] <- NA
    overall_betaXse$smoke[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the smoke_data_frame
rm(smoke_data_frame)

# Load the mpv_data.tsv file
mpv_data_frame <- readr::read_tsv("Mediation/mpv_data.tsv")
# Add the beta values from the mpv_data_frame to the betaX$mpv column
for (variant in all_possible_variants$variant) {
  if (variant %in% mpv_data_frame$variant) {
    overall_betaX$mpv[overall_betaX$variant == variant] <- mpv_data_frame$beta[mpv_data_frame$variant == variant]
    overall_betaXse$mpv[overall_betaXse$variant == variant] <- mpv_data_frame$se[mpv_data_frame$variant == variant]
  } else {
    overall_betaX$mpv[overall_betaX$variant == variant] <- NA
    overall_betaXse$mpv[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the mpv_data_frame
rm(mpv_data_frame)

# Add a ApoB column to the overall_betaX and overall_betaXse data frames
overall_betaX$ApoB <- NA
overall_betaXse$ApoB <- NA

# Load the ApoB_data.tsv file
ApoB_data_frame <- readr::read_tsv("Mediation/ApoB_data.tsv")
# Add the beta values from the ApoB_data_frame to the betaX$ApoB column
for (variant in all_possible_variants$variant) {
  if (variant %in% ApoB_data_frame$variant) {
    overall_betaX$ApoB[overall_betaX$variant == variant] <- ApoB_data_frame$beta[ApoB_data_frame$variant == variant]
    overall_betaXse$ApoB[overall_betaXse$variant == variant] <- ApoB_data_frame$se[ApoB_data_frame$variant == variant]
  } else {
    overall_betaX$ApoB[overall_betaX$variant == variant] <- NA
    overall_betaXse$ApoB[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the ApoB_data_frame
rm(ApoB_data_frame)

# Add a crp column to the overall_betaX and overall_betaXse data frames
overall_betaX$crp <- NA
overall_betaXse$crp <- NA

# Load the crp_data.tsv file
crp_data_frame <- readr::read_tsv("Mediation/crp_data.tsv")
# Add the beta values from the crp_data_frame to the betaX$crp column
for (variant in all_possible_variants$variant) {
  if (variant %in% crp_data_frame$variant) {
    overall_betaX$crp[overall_betaX$variant == variant] <- crp_data_frame$beta[crp_data_frame$variant == variant]
    overall_betaXse$crp[overall_betaXse$variant == variant] <- crp_data_frame$se[crp_data_frame$variant == variant]
  } else {
    overall_betaX$crp[overall_betaX$variant == variant] <- NA
    overall_betaXse$crp[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the crp_data_frame
rm(crp_data_frame)

# Add a rbc column to the overall_betaX and overall_betaXse data frames
overall_betaX$rbc <- NA
overall_betaXse$rbc <- NA

# Load the rbc_data.tsv file
rbc_data_frame <- readr::read_tsv("Mediation/rbc_data.tsv")
# Add the beta values from the rbc_data_frame to the betaX$rbc column
for (variant in all_possible_variants$variant) {
  if (variant %in% rbc_data_frame$variant) {
    overall_betaX$rbc[overall_betaX$variant == variant] <- rbc_data_frame$beta[rbc_data_frame$variant == variant]
    overall_betaXse$rbc[overall_betaXse$variant == variant] <- rbc_data_frame$se[rbc_data_frame$variant == variant]
  } else {
    overall_betaX$rbc[overall_betaX$variant == variant] <- NA
    overall_betaXse$rbc[overall_betaXse$variant == variant] <- NA
  }
}
# Remove the rbc_data_frame
rm(rbc_data_frame)

# Save the overall_betaX and overall_betaXse dataframes
write_tsv(overall_betaX, "Mediation/overall_betaX.tsv")
write_tsv(overall_betaXse, "Mediation/overall_betaXse.tsv")


all_possible_variants <- readr::read_tsv("Mediation/all_possible_variants.tsv")

overall_betaY <- data.frame(all_possible_variants)
overall_betaY$T2DM <- NA
overall_betaYse <- data.frame(all_possible_variants)
overall_betaYse$T2DM <- NA

t2dm_data_frame <- readr::read_tsv("t2dm_gwas_cleaned.tsv")
# Make a variant column by pasting the Chromosome:Position:Effect_Allele:Noneffect_Allele
t2dm_data_frame$variant <- paste0(t2dm_data_frame$Chromosome, ":", t2dm_data_frame$Position, ":", t2dm_data_frame$EffectAllele, ":", t2dm_data_frame$NonEffectAllele)
match_flag <- FALSE
proxy_flag <- FALSE
for (variant in all_possible_variants$variant) {
  # If the variant is in the t2dm_data_frame, add the beta and se values to the overall_betaY and overall_betaYse dataframes
  if (variant %in% t2dm_data_frame$variant) {
    overall_betaY$T2DM[overall_betaY$variant == variant] <- t2dm_data_frame$Beta[t2dm_data_frame$variant == variant]
    overall_betaYse$T2DM[overall_betaYse$variant == variant] <- t2dm_data_frame$SE[t2dm_data_frame$variant == variant]
    match_flag <- TRUE
  }
  if (!match_flag) {
    rsid <- all_possible_variants$rsid[all_possible_variants$variant == variant]
    # Find a proxy for the IV
    ld_proxies <- gwasvcf::get_ld_proxies(
      rsid = rsid,
      bfile = "1kg.v3/EUR",
      searchspace = NULL,
      tag_kb = 5000,
      tag_nsnp = 5000,
      tag_r2 = 0.8,
      threads = 1,
      out = tempfile())
    # Remove proxies with ld_proxies$R < 0.8
    ld_proxies <- ld_proxies[ld_proxies$R >= 0.8, ]
    if (nrow(ld_proxies) > 0 ) {
      # Order the proxies by R in descending order
      ld_proxies <- ld_proxies[order(ld_proxies$R, decreasing = TRUE), ]
      # Make a MarkerName and reverseMarkerName column
      ld_proxies$MarkerName <- paste0(ld_proxies$CHR_B, ":", ld_proxies$BP_B, ":", ld_proxies$B1, ":", ld_proxies$B2)
      for (j in 1:nrow(ld_proxies)) {
        if (ld_proxies$MarkerName[j] %in% t2dm_data_frame$variant) {
          overall_betaY$T2DM[overall_betaY$variant == variant] <- t2dm_data_frame$Beta[t2dm_data_frame$variant == ld_proxies$MarkerName[j]]
          overall_betaYse$T2DM[overall_betaY$variant == variant] <- t2dm_data_frame$SE[t2dm_data_frame$variant == ld_proxies$MarkerName[j]]
          proxy_flag <- TRUE
          break
        }
      }
    }
  }
  if (!proxy_flag && !match_flag) {
    # If the variant is not in the haem_data_frame, set the beta value to NA
    overall_betaY$T2DM[overall_betaY$variant == variant] <- NA
    overall_betaYse$T2DM[overall_betaY$variant == variant] <- NA
  }
  match_flag <- FALSE
  proxy_flag <- FALSE
}
# Remove the t2dm_data_frame
rm(t2dm_data_frame)

# Save the overall_betaY and overall_betaYse dataframes
write_tsv(overall_betaY, "Mediation/overall_betaY.tsv")
write_tsv(overall_betaYse, "Mediation/overall_betaYse.tsv")




# Load the coloc_filtered_metabolites_data.tsv file
metabolites_data_frame <- readr::read_tsv("coloc_filtered_metabolites_data.tsv")
# Group the metabolites by subpathway
subpathways <- unique(metabolites_data_frame$subpathway)
# Add an "M" to each metabolite compid if there isnt one
metabolites_data_frame$compid <- ifelse(grepl("^M", metabolites_data_frame$compid), metabolites_data_frame$compid, paste0("M", metabolites_data_frame$compid))


# Load the all_possible_variants dataframe
all_possible_variants <- readr::read_tsv("Mediation/all_possible_variants.tsv")
# Load the overall_betaX and overall_betaXse dataframes
overall_betaX <- readr::read_tsv("Mediation/overall_betaX.tsv")
overall_betaXse <- readr::read_tsv("Mediation/overall_betaXse.tsv")


# Load the ukbb_variants_data.tsv file
ukbb_variants_data_frame <- readr::read_tsv("Mediation/ukbb_variants_data.tsv")

# Load the haem_data.tsv file
haem_data_frame <- readr::read_tsv("Mediation/haem_data.tsv")
# Grab the rows that have NA in the overall_betaX$heam and overall_betaXse$heam
na_rows <- which(is.na(overall_betaX$heam))
# For each row with NA in the overall_betaX$heam and overall_betaXse$heam
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$heam[overall_betaX$variant == variant] <- haem_data_frame$beta[haem_data_frame$variant == proxy_variant]
  overall_betaXse$heam[overall_betaXse$variant == variant] <- haem_data_frame$se[haem_data_frame$variant == proxy_variant]
}
# Remove the haem_data_frame
rm(haem_data_frame)

# Load the bmi_data.tsv file
bmi_data_frame <- readr::read_tsv("Mediation/bmi_data.tsv")
# Grab the rows that have NA in the overall_betaX$bmi and overall_betaXse$bmi
na_rows <- which(is.na(overall_betaX$bmi))
# For each row with NA in the overall_betaX$bmi and overall_betaXse$bmi
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$bmi[overall_betaX$variant == variant] <- bmi_data_frame$beta[bmi_data_frame$variant == proxy_variant]
  overall_betaXse$bmi[overall_betaXse$variant == variant] <- bmi_data_frame$se[bmi_data_frame$variant == proxy_variant]
}
# Remove the bmi_data_frame
rm(bmi_data_frame)

# Load the waist_data.tsv file
waist_data_frame <- readr::read_tsv("Mediation/waist_data.tsv")
# Grab the rows that have NA in the overall_betaX$waist and overall_betaXse$waist
na_rows <- which(is.na(overall_betaX$waist))
# For each row with NA in the overall_betaX$waist and overall_betaXse$waist
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$waist[overall_betaX$variant == variant] <- waist_data_frame$beta[waist_data_frame$variant == proxy_variant]
  overall_betaXse$waist[overall_betaXse$variant == variant] <- waist_data_frame$se[waist_data_frame$variant == proxy_variant]
}
# Remove the waist_data_frame
rm(waist_data_frame)

# Load the body_fat_data.tsv file
body_fat_data_frame <- readr::read_tsv("Mediation/body_fat_data.tsv")
# Grab the rows that have NA in the overall_betaX$body_fat and overall_betaXse$body_fat
na_rows <- which(is.na(overall_betaX$body_fat))
# For each row with NA in the overall_betaX$body_fat and overall_betaXse$body_fat
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$body_fat[overall_betaX$variant == variant] <- body_fat_data_frame$beta[body_fat_data_frame$variant == proxy_variant]
  overall_betaXse$body_fat[overall_betaXse$variant == variant] <- body_fat_data_frame$se[body_fat_data_frame$variant == proxy_variant]
}
# Remove the body_fat_data_frame
rm(body_fat_data_frame)

# Load the alchohol_data.tsv file
alchohol_data_frame <- readr::read_tsv("Mediation/alcohol_data.tsv")
# Grab the rows that have NA in the overall_betaX$alchohol and overall_betaXse$alchohol
na_rows <- which(is.na(overall_betaX$alchohol))
# For each row with NA in the overall_betaX$alchohol and overall_betaXse$alchohol
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$alchohol[overall_betaX$variant == variant] <- alchohol_data_frame$beta[alchohol_data_frame$variant == proxy_variant]
  overall_betaXse$alchohol[overall_betaXse$variant == variant] <- alchohol_data_frame$se[alchohol_data_frame$variant == proxy_variant]
}
# Remove the alchohol_data_frame
rm(alchohol_data_frame)

# Load the smoke_data.tsv file
smoke_data_frame <- readr::read_tsv("Mediation/smoke_data.tsv")
# Grab the rows that have NA in the overall_betaX$smoke and overall_betaXse$smoke
na_rows <- which(is.na(overall_betaX$smoke))
# For each row with NA in the overall_betaX$smoke and overall_betaXse$smoke
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$smoke[overall_betaX$variant == variant] <- smoke_data_frame$beta[smoke_data_frame$variant == proxy_variant]
  overall_betaXse$smoke[overall_betaXse$variant == variant] <- smoke_data_frame$se[smoke_data_frame$variant == proxy_variant]
}
# Remove the smoke_data_frame
rm(smoke_data_frame)

# Load the mpv_data.tsv file
mpv_data_frame <- readr::read_tsv("Mediation/mpv_data.tsv")
# Grab the rows that have NA in the overall_betaX$mpv and overall_betaXse$mpv
na_rows <- which(is.na(overall_betaX$mpv))
# For each row with NA in the overall_betaX$mpv and overall_betaXse$mpv
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$mpv[overall_betaX$variant == variant] <- mpv_data_frame$beta[mpv_data_frame$variant == proxy_variant]
  overall_betaXse$mpv[overall_betaXse$variant == variant] <- mpv_data_frame$se[mpv_data_frame$variant == proxy_variant]
}
# Remove the mpv_data_frame
rm(mpv_data_frame)

# Load the ApoB_data.tsv file
ApoB_data_frame <- readr::read_tsv("Mediation/ApoB_data.tsv")
# Grab the rows that have NA in the overall_betaX$ApoB and overall_betaXse$ApoB
na_rows <- which(is.na(overall_betaX$ApoB))
# For each row with NA in the overall_betaX$ApoB and overall_betaXse$ApoB
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$ApoB[overall_betaX$variant == variant] <- ApoB_data_frame$beta[ApoB_data_frame$variant == proxy_variant]
  overall_betaXse$ApoB[overall_betaXse$variant == variant] <- ApoB_data_frame$se[ApoB_data_frame$variant == proxy_variant]
}
# Remove the ApoB_data_frame
rm(ApoB_data_frame)

# Load the crp_data.tsv file
crp_data_frame <- readr::read_tsv("Mediation/crp_data.tsv")
# Grab the rows that have NA in the overall_betaX$crp and overall_betaXse$crp
na_rows <- which(is.na(overall_betaX$crp))
# For each row with NA in the overall_betaX$crp and overall_betaXse$crp
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$crp[overall_betaX$variant == variant] <- crp_data_frame$beta[crp_data_frame$variant == proxy_variant]
  overall_betaXse$crp[overall_betaXse$variant == variant] <- crp_data_frame$se[crp_data_frame$variant == proxy_variant]
}
# Remove the crp_data_frame
rm(crp_data_frame)

# Load the rbc_data.tsv file
rbc_data_frame <- readr::read_tsv("Mediation/rbc_data.tsv")
# Grab the rows that have NA in the overall_betaX$rbc and overall_betaXse$rbc
na_rows <- which(is.na(overall_betaX$rbc))
# For each row with NA in the overall_betaX$rbc and overall_betaXse$rbc
for (i in na_rows) {
  variant <- overall_betaX$variant[i]
  # Use the manual proxy
  manual_proxy <- overall_betaX$manual_proxy[overall_betaX$variant == variant]
  proxy_variant <- ukbb_variants_data_frame$variant[ukbb_variants_data_frame$rsid == manual_proxy]
  overall_betaX$rbc[overall_betaX$variant == variant] <- rbc_data_frame$beta[rbc_data_frame$variant == proxy_variant]
  overall_betaXse$rbc[overall_betaXse$variant == variant] <- rbc_data_frame$se[rbc_data_frame$variant == proxy_variant]
}
# Remove the rbc_data_frame
rm(rbc_data_frame)

# Save overall_betaX and overall_betaXse
readr::write_tsv(overall_betaX, "Mediation/overall_betaX.tsv")
readr::write_tsv(overall_betaXse, "Mediation/overall_betaXse.tsv")


# Load the overall_betaY and overall_betaYse dataframes
overall_betaY <- readr::read_tsv("Mediation/overall_betaY.tsv")
overall_betaYse <- readr::read_tsv("Mediation/overall_betaYse.tsv")

t2dm_data_frame <- readr::read_tsv("t2dm_gwas_cleaned.tsv")
# Make a variant column by pasting the Chromosome:Position:Effect_Allele:Noneffect_Allele
t2dm_data_frame$variant <- paste0(t2dm_data_frame$Chromosome, ":", t2dm_data_frame$Position, ":", t2dm_data_frame$EffectAllele, ":", t2dm_data_frame$NonEffectAllele)

# Load the unqiue colocalised metabolite IVs
coloc_metab_variants <- readr::read_tsv("coloc_metab_variants.tsv")

# The only variant data missing from overall_betaY and overall_betaYse is for rs6857.
# This SNP does not have any good proxies but should have a match in the dataset.
# The variant is 19:45392254:T:C, check if it is in the dataset
t2dm_data_frame$Beta[t2dm_data_frame$variant == "19:45392254:T:C"]
# The alleles are switched in the dataset, so I need to manually add it and multiply the beta by -1.
# Manually add the data to the overall_betaY and overall_betaYse dataframes
overall_betaY$T2DM[overall_betaY$variant == "19:45392254:C:T"] <- (t2dm_data_frame$Beta[t2dm_data_frame$variant == "19:45392254:T:C"]*-1)
overall_betaYse$T2DM[overall_betaYse$variant == "19:45392254:C:T"] <- t2dm_data_frame$SE[t2dm_data_frame$variant == "19:45392254:T:C"]

# Save the overall_betaY and overall_betaYse dataframes
readr::write_tsv(overall_betaY, "Mediation/overall_betaY.tsv")
readr::write_tsv(overall_betaYse, "Mediation/overall_betaYse.tsv")





# Load the coloc_filtered_metabolites_data.tsv file
metabolites_data_frame <- readr::read_tsv("coloc_filtered_metabolites_data.tsv")
# Group the metabolites by subpathway
subpathways <- unique(metabolites_data_frame$subpathway)
# Add an "M" to each metabolite compid if there isnt one
metabolites_data_frame$compid <- ifelse(grepl("^M", metabolites_data_frame$compid), metabolites_data_frame$compid, paste0("M", metabolites_data_frame$compid))


# Load the all_possible_variants dataframe
all_possible_variants <- readr::read_tsv("Mediation/all_possible_variants.tsv")
# Load the overall_betaX and overall_betaXse dataframes
overall_betaX <- readr::read_tsv("Mediation/overall_betaX.tsv")
overall_betaXse <- readr::read_tsv("Mediation/overall_betaXse.tsv")
# Load the overall_betaY and overall_betaYse dataframes
overall_betaY <- readr::read_tsv("Mediation/overall_betaY.tsv")
overall_betaYse <- readr::read_tsv("Mediation/overall_betaYse.tsv")

# Make a dataframe with the unique variants as the rows
all_possible_variants$NAs <- NA
# If there is an NA in overall_betaX, overall_betaXse, overall_betaY, and overall_betaYse dataframes, set the NAs column to TRUE
for (variant in all_possible_variants$variant) {
  if (any(is.na(overall_betaX$heam[overall_betaX$variant == variant])) | any(is.na(overall_betaXse$heam[overall_betaXse$variant == variant])) | any(is.na(overall_betaY$T2DM[overall_betaY$variant == variant])) | any(is.na(overall_betaYse$T2DM[overall_betaYse$variant == variant]))) {
    all_possible_variants$NAs[all_possible_variants$variant == variant] <- TRUE
  }
}
# If the value is not TRUE, set it to FALSE
all_possible_variants$NAs[is.na(all_possible_variants$NAs)] <- FALSE

# Remove rows from the overall_betaX, overall_betaXse, overall_betaY, and overall_betaYse dataframes where all_possible_variants$NAs is TRUE
overall_betaX <- overall_betaX[!all_possible_variants$NAs, ]
overall_betaXse <- overall_betaXse[!all_possible_variants$NAs, ]
overall_betaY <- overall_betaY[!all_possible_variants$NAs, ]
overall_betaYse <- overall_betaYse[!all_possible_variants$NAs, ]

# Load in the T2DM_files
T2DM_files <- list.files("Harmonised_T2DM_IVs", full.names = TRUE)

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

# Add a column to all_possible_variants to store proxy information
all_possible_variants$proxy <- NA

iv_chromosome_old <- 0

# For each subpathway
for (subpathway in subpathways){
  # Load the unique variants for this subpathway
  unique_variants <- readr::read_tsv(paste0("Mediation/", "unique_variants_", subpathway, ".tsv"))
  # Get the metabolites in the subpathway
  metabolites <- metabolites_data_frame$name[metabolites_data_frame$subpathway == subpathway]
  # Initialise metabolites_removed
  metabolites_removed <- c()
  
  # Make a dataframe with the unique variants as the rows
  betaX <- data.frame(unique_variants)
  betaXse <- data.frame(unique_variants)
  for (metabolite in metabolites) {
    betaX[metabolite] <- NA
  }
  for (metabolite in metabolites) {
    betaXse[metabolite] <- NA
  }
  # Add the columns from overall_betaX to betaX by matching the variant
  betaX <- merge(betaX, overall_betaX, by = "variant", all.x = TRUE)
  betaXse <- merge(betaXse, overall_betaXse, by = "variant", all.x = TRUE)
  
  # Remove rows from betaX and betaXse that are not in overall_betaX$variant and overall_betaXse$variant
  betaX <- betaX[!is.na(betaX$heam), ]
  betaXse <- betaXse[!is.na(betaXse$heam), ]
  
  # Order the betaX and betaXse dataframes by the variant column values before the first ":"
  # For each row in betaX and betaXse, get the chromosome by splitting the variant column by ":" and taking the first element
  for (variant in betaX$variant) {
    betaX$chromosome[betaX$variant == variant] <- as.numeric(strsplit(variant, ":")[[1]][1])
    betaXse$chromosome[betaXse$variant == variant] <- as.numeric(strsplit(variant, ":")[[1]][1])
  }
  betaX <- betaX[order(betaX$chromosome), ]
  betaXse <- betaXse[order(betaXse$chromosome), ]
  
  for (metabolite in metabolites) {
    # Get the IVs of each metabolite
    metabolite_file <- metabolite_file <- T2DM_files[str_detect(T2DM_files, fixed(metabolite, ignore_case = TRUE))]
    metabolite_data <- readr::read_tsv(metabolite_file)
    # Change the MarkerName column _ to : and remove chr
    metabolite_data$MarkerName <- gsub("_", ":", metabolite_data$MarkerName)
    metabolite_data$MarkerName <- gsub("chr", "", metabolite_data$MarkerName)
    
    # Check if one of the metabolite_data$MarkerName is in the betaX$variant
    if (any(metabolite_data$MarkerName %in% betaX$variant)) {
      # Get the compid of the metabolite
      compid <- metabolites_data_frame$compid[metabolites_data_frame$name == metabolite]
      
      for (variant in betaX$variant) {
        match_flag <- FALSE
        proxy_flag <- FALSE
        # Get the chromosome of the variant, it is the part before the first ":" in the variant
        iv_chromosome <- strsplit(variant, ":")[[1]][1]
        # If the iv_chromosome is the different than iv_chromosome, re-read the necessary files
        if (iv_chromosome != iv_chromosome_old) {
          # Re-read the necessary files
          raw_gz_file <- paste0("Reverse_MR/Raw_Metab_GWAS_Data/", compid, "/", "INTERVAL_MRC-Epi_",  compid, "_sorted_chr_", iv_chromosome, ".tbl.gz")
          # If the file is not already unzipped, unzip it
          tryCatch({
            R.utils::gunzip(raw_gz_file, overwrite = TRUE)
          }, error = function(e) {
            message("An error occurred while unzipping the file: ", e)
          })
          raw_tbl_file <- paste0("Reverse_MR/Raw_Metab_GWAS_Data/", compid, "/","INTERVAL_MRC-Epi_", compid, "_sorted_chr_", iv_chromosome, ".tbl")
          # If the file is already unzipped, read it and re-zip it
          combined_data <- readr::read_tsv(raw_tbl_file)
          # Re-zip the file
          R.utils::gzip(raw_tbl_file, overwrite = TRUE)
        }
        
        # Add chr to the variant
        modified_variant <- paste0("chr", variant)
        
        # If the modified_variant is in combined_data$MarkerName
        if (modified_variant %in% combined_data$MarkerName){
          # Get the beta value for the metabolite
          betaX[betaX$variant == variant, metabolite] <- combined_data$Effect[combined_data$MarkerName == modified_variant]
          # Get the standard error value for the metabolite
          betaXse[betaXse$variant == variant, metabolite] <- combined_data$StdErr[combined_data$MarkerName == modified_variant]
          match_flag <- TRUE
        } 
        if (!match_flag) {
          if (!is.na(all_possible_variants$proxy[all_possible_variants$variant == variant])) {
            proxy_variant <- all_possible_variants$proxy[all_possible_variants$variant == variant]
            
            betaX[betaX$variant == variant, metabolite] <- combined_data$Effect[combined_data$MarkerName == proxy_variant]
            betaXse[betaXse$variant == variant, metabolite] <- combined_data$StdErr[combined_data$MarkerName == proxy_variant]
          } else {
            rsid <- all_possible_variants$rsid[all_possible_variants$variant == variant]
            # Find a proxy for the IV
            ld_proxies <- gwasvcf::get_ld_proxies(
              rsid = rsid,
              bfile = "1kg.v3/EUR",
              searchspace = NULL,
              tag_kb = 5000,
              tag_nsnp = 5000,
              tag_r2 = 0.9,
              threads = 1,
              out = tempfile())
            # Remove proxies with ld_proxies$R < 0.9
            ld_proxies <- ld_proxies[ld_proxies$R >= 0.9, ]
            if (nrow(ld_proxies) > 0 ) {
              # Order the proxies by R in descending order
              ld_proxies <- ld_proxies[order(ld_proxies$R, decreasing = TRUE), ]
              # Make a MarkerName and reverseMarkerName column
              ld_proxies$MarkerName <- paste0("chr", ld_proxies$CHR_B, ":", ld_proxies$BP_B, ":", ld_proxies$B1, ":", ld_proxies$B2)
              for (j in 1:nrow(ld_proxies)) {
                if (ld_proxies$MarkerName[j] %in% combined_data$MarkerName) {
                  betaX[betaX$variant == variant, metabolite] <- combined_data$Effect[combined_data$MarkerName == ld_proxies$MarkerName[j]]
                  betaXse[betaXse$variant == variant, metabolite] <- combined_data$StdErr[combined_data$MarkerName == ld_proxies$MarkerName[j]]
                  proxy_flag <- TRUE
                  all_possible_variants$proxy[all_possible_variants$variant == variant] <- ld_proxies$MarkerName[j]
                  break
                }
              }
            }
          }
        }
        if (!proxy_flag && !match_flag) {
          betaX[betaX$variant == variant, metabolite] <- NA
          betaXse[betaXse$variant == variant, metabolite] <- NA
        }
        match_flag <- FALSE
        proxy_flag <- FALSE
        iv_chromosome_old <- iv_chromosome
      }
    }
  }
  
  # Remove the chromosome column
  betaX <- betaX[, -ncol(betaX)]
  betaXse <- betaXse[, -ncol(betaXse)]
  
  # Make a dataframe with the unique variants as the rows
  betaY <- data.frame(unique_variants)
  betaYse <- data.frame(unique_variants)
  # Add the columns from overall_betaY to betaY by matching the variant
  betaY <- merge(betaY, overall_betaY, by = "variant", all.x = TRUE)
  betaYse <- merge(betaYse, overall_betaYse, by = "variant", all.x = TRUE)
  
  # Save the subpathway specific betaX betaXse betaY and betaYse
  write_tsv(betaX, file = paste0("Mediation/betaX", subpathway, ".tsv"))
  write_tsv(betaXse, file = paste0("Mediation/betaXse", subpathway, ".tsv"))
  write_tsv(betaY, file = paste0("Mediation/betaY", subpathway, ".tsv"))
  write_tsv(betaYse, file = paste0("Mediation/betaYse", subpathway, ".tsv"))
  
}