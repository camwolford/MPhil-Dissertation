library(readr)

# Reading the contents of TSV file using read_tsv() method
metabolite_gwas_associations<-readr::read_tsv("metabolite-gwas-associations.tsv")
metabolites<-readr::read_tsv("metabolites.tsv")

# remove the firstAuthor, publicationDate, journal, title, ssApiFlag, agreeToCc0, pubmedId, initialSampleDescription and replicateSampleDescription columns
metabolites<-metabolites[,-c(1,3,4,5,13,14,16,17,18)]
# what are the unique values of the "bgTraits" column?
unique(metabolites$bgTraits)
# remove the "bgTraits" column
metabolites<-metabolites[,-c(4)]
# what are the unique values of the "summaryStatistics" column?
unique(metabolites$summaryStatistics)
# remove the "summaryStatistics" column
metabolites<-metabolites[,-c(7)]
# How many of the reportedTraits are have a value of 0 in the associationCount column?
sum(metabolites$associationCount==0)
# what are the types of the columns in the metabolites data frame?
str(metabolites)
# Save the cleaned data to a new file called "metabolites_cleaned.csv"
write.csv(metabolites, file = "metabolites_cleaned.csv", row.names = FALSE)

# remove the 1,2,3,4,5,6,7,9,10,33,35,36,37 columns from the metabolite_gwas_associations data frame
metabolite_gwas_associations<-metabolite_gwas_associations[,-c(1,2,3,4,5,6,7,9,10,33,35,36,37)]
# what are the unique values of the "REPORTED GENE(S)" column?
unique(metabolite_gwas_associations$`REPORTED GENE(S)`)
# remove the "REPORTED GENE(S)" column
metabolite_gwas_associations<-metabolite_gwas_associations[,-c(5)]
# what are the unique values of the "P-VALUE (TEXT)" column?
unique(metabolite_gwas_associations$`P-VALUE (TEXT)`)
# remove the "P-VALUE (TEXT)" column
metabolite_gwas_associations<-metabolite_gwas_associations[,-c(20)]
# what are the unique values of the "CNV" column?
unique(metabolite_gwas_associations$CNV)
# remove the "CNV" column
metabolite_gwas_associations<-metabolite_gwas_associations[,-c(22,23)]
# Change the name of the "DISEASE/TRAIT" column to "Metabolite"
colnames(metabolite_gwas_associations)[1]<-"Metabolite"

# how many Metabolites from the metabolite_gwas_associations data frame are in the reportedTraits column of the metabolites data frame?
sum(metabolites$reportedTrait %in% metabolite_gwas_associations$Metabolite)
# 713 Metabolites from the metabolite_gwas_associations data frame are in the reportedTraits column of the metabolites data frame
# what are the types of the columns in the metabolite_gwas_associations data frame?
str(metabolite_gwas_associations)
# what are the unique values of the "MERGED" column?
unique(metabolite_gwas_associations$MERGED)
# value counts of the "MERGED" column
table(metabolite_gwas_associations$MERGED)
# Look at the one row with a value of "1" in the "MERGED" column
metabolite_gwas_associations[metabolite_gwas_associations$MERGED==1,]
# Save the cleaned data to a new file called "metabolite_gwas_associations_cleaned.csv"
write.csv(metabolite_gwas_associations, file = "metabolite_gwas_associations_cleaned.csv", row.names = FALSE)

# Import the EUR_Metal_LDSC-CORR_Neff.v2.txt file
t2dm_gwas<-readr::read_tsv("EUR_Metal_LDSC-CORR_Neff.v2.txt")
# what are the types of the columns in the t2dm_gwas data frame?
str(t2dm_gwas)
# Save the cleaned data to a new file called "t2dm_gwas_cleaned.csv"
write.csv(t2dm_gwas, file = "t2dm_gwas_cleaned.csv", row.names = FALSE)


library(R.utils)
# gunzip the fasting_glucose_gwas.tsv.gz file
gunzip("fasting_glucose_gwas.tsv.gz", "new_fasting_glucose_gwas.tsv")
# Read in the new_fasting_glucose_gwas.tsv file
fasting_glucose_gwas<-readr::read_tsv("fasting_glucose_gwas.tsv")
# gunzip the HbA1c_gwas.tsv.gz file
gunzip("HbA1c_gwas.tsv.gz", "new_HbA1c_gwas.tsv")
# Read in the new_HbA1c_gwas.tsv file
hbA1c_gwas<-readr::read_tsv("new_HbA1c_gwas.tsv")
# Save the cleaned data to a new file called "hbA1c_gwas_cleaned.csv"
write.csv(hbA1c_gwas, file = "hbA1c_gwas_cleaned.csv", row.names = FALSE)
# Save the cleaned data to a new file called "fasting_glucose_gwas_cleaned.csv"
write.csv(fasting_glucose_gwas, file = "fasting_glucose_gwas_cleaned.csv", row.names = FALSE)


# Read in the responses.tsv file
responses<-readr::read_tsv("responses.tsv")
# Read in the responsesns.tsv file
responsesns<-readr::read_tsv("responsesns.tsv")

# make a new column in responses and called source and set it to "responses"
responses$source<-"responses"
# make a new column in responsesns and called source and set it to "responsesns"
responsesns$source<-"responsesns"
# add a column called key to resonsesns and set it to NA
responsesns$key<-NA
# combine the responses and responsesns data frames
responses_combined<-rbind(responses, responsesns)
# what are the types of the columns in the responses_combined data frame?
str(responses_combined)
# Save the cleaned data to a new file called "responses_combined.csv"
write.csv(responses_combined, file = "responses_combined.csv", row.names = FALSE)



### New Metabolite Association Data ###
# Read in the metabolite_sup_associations.txt file
metabolite_sup_associations<-read.delim("metabolite_sup_associations.txt")
# Rename the columns to the values of the first row
colnames(metabolite_sup_associations)<-metabolite_sup_associations[1,]
# Remove the first row
metabolite_sup_associations<-metabolite_sup_associations[-1,]
# Remove columns 1, 2, 3, 11, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 27,28,29,30,31,33,34,35,36,37,38,39,40,41
metabolite_sup_associations<-metabolite_sup_associations[,-c(1,2,3,11,12,13,14,15,17,18,19,20,22,23,24,25,27,28,29,30,31,33,34,35,36,37,38,39,40,41)]
# Rename the biochemical column to "Metabolite"
colnames(metabolite_sup_associations)[4]<-"SNP"
colnames(metabolite_sup_associations)[5]<-"EffectAllele"
colnames(metabolite_sup_associations)[6]<-"NonEffectAllele"
colnames(metabolite_sup_associations)[7]<-"Metabolite"
colnames(metabolite_sup_associations)[8]<-"EAF"
colnames(metabolite_sup_associations)[9]<-"Beta"
colnames(metabolite_sup_associations)[10]<-"SE"
colnames(metabolite_sup_associations)[11]<-"Pval"
str(metabolite_sup_associations)
# Convert columns 1,2,8,9,10,11 to numeric
metabolite_sup_associations$Chromosome<-as.numeric(metabolite_sup_associations$Chromosome)
metabolite_sup_associations$Position<-as.numeric(metabolite_sup_associations$Position)
metabolite_sup_associations$EAF<-as.numeric(metabolite_sup_associations$EAF)
metabolite_sup_associations$Beta<-as.numeric(metabolite_sup_associations$Beta)
metabolite_sup_associations$SE<-as.numeric(metabolite_sup_associations$SE)
metabolite_sup_associations$Pval<-as.numeric(metabolite_sup_associations$Pval)
# Change the Pval values to 10^-Pval
metabolite_sup_associations$Pval_new<-(10^(-1*metabolite_sup_associations$Pval))
metabolite_sup_associations$Pval<-metabolite_sup_associations$Pval_new
metabolite_sup_associations<-metabolite_sup_associations[,-12]
# unique values of the Metabolite column
length(unique(metabolite_sup_associations$Metabolite))
# how many NA values in the Metabolite column
sum(is.na(metabolite_sup_associations$Metabolite))

# Count how many NAs are in each column
colSums(is.na(metabolite_sup_associations))
# Remove any rows with NAs
metabolite_sup_associations<-metabolite_sup_associations[complete.cases(metabolite_sup_associations),]
# Number of unique values in the Metabolite column
length(unique(metabolite_sup_associations$Metabolite))

# Replace any / with _ in the Metabolite column
metabolite_sup_associations$Metabolite<-gsub(":", "_", metabolite_sup_associations$Metabolite)

# Save the data to "metabolite_gwas_associations_cleaned.tsv"
write.table(metabolite_sup_associations, "metabolite_gwas_associations_cleaned.tsv", sep="\t", row.names=FALSE)
# Load the same file to check if it was saved correctly
metabolite_sup_associations<-readr::read_tsv("metabolite_gwas_associations_cleaned.tsv")









