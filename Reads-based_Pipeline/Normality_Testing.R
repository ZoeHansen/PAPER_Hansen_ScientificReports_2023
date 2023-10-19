########################################################

# Aim 2: Normality Testing (Case, Control, Follow-Up and Case vs. Follow)

########################################################

# Load libraries and data

library(tidyverse)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(vegan)
library(calibrate)

pheno <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/ERIN_Metagenomes_Metadata_60_CaseFollowPairs_wControls.csv',
                  header = TRUE) 

args = read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_ALL_TYPES_gene_level.csv', 
                header = TRUE)

microbiome = read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/GEnorm_ERIN_kaiju_CaseControlFollow_SPECIES.csv',
                      header=TRUE)

### If data needs to be transposed (species):

microbiome$Species <- microbiome$Species %>%
  replace_na("Unknown")

micro <- ddply(microbiome, "Species", numcolwise(sum))

micro<- micro[, colSums(micro[,-1] != 0) > 0] 
micro <- micro[rowSums(micro[,-1] != 0) > 0,]

data <- micro %>%
  gather(key = key, value = value, 2:ncol(micro)) %>%
  spread(key=names(micro)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

data <- args

data<- data[, colSums(data[,-1] != 0) > 0] 
data <- data[rowSums(data[,-1] != 0) > 0,] 

data[is.na(data)] <- 0

########################################################
# Selecting relevant metadata for analysis
########################################################

# Select the relevant metadata fields to include in analysis
health <- pheno %>%
  dplyr::select(ER_ID, Case.status, Pathogen, Antibiotics, Reassigned.Pair.ID) #%>%#, Case.Follow_ID)
#  filter(!grepl('Control', Case.status))

health$Antibiotics[health$Antibiotics == ""] <- "No"

IDs_use <- health %>%
  dplyr::select(ER_ID)


# Merge metadata with AGS-normalized abundances
#Full Case_Control Campylobacter dataset with metadata of interest
data_2 <- left_join(IDs_use, data, by = 'ER_ID')

# Remove any lingering zeroes to avoid downstream errors
data_2 <- data_2[, colSums(data_2 != 0) > 0]
data_2 <- data_2[rowSums(data_2 != 0) > 0,] 


data_cc <- left_join(health, data_2, by = "ER_ID")

# If needed, remove any samples from analysis here
data_cc <- data_cc %>%
  filter(!grepl('Control', Case.status)) #%>%
#  filter(!grepl('FollowUp',Case.status))%>%
#  drop_na(.,Case.Follow_ID)

data_cc$ER_ID <- factor(data_cc$ER_ID)
data_cc$Case.status <- factor(data_cc$Case.status)
data_cc$Reassigned.Pair.ID <- factor(data_cc$Reassigned.Pair.ID)
data_cc$Pathogen <- factor(data_cc$Pathogen)
data_cc$Antibiotics <- factor(data_cc$Antibiotics)

data_cc[is.na(data_cc)] <- 0
data_cc <- data_cc[, colSums(data_cc != 0) > 0]
data_cc <- data_cc[rowSums(data_cc != 0) > 0,] 


########################################################
# Calculate and plot within-sample (Alpha) Diversity
########################################################

# Richness (of genes) across our samples
r.c = specnumber(data_cc[,-c(1:5)])

# Shannon diversity
h.c = diversity(data_cc[,-c(1:5)], index = 'shannon')

# Calculate Pielou's Evenness index
pielou.c = h.c/log(r.c)

### Case Control Follow
# Combine alpha diversity data and Case status information
div.c=tibble(data_cc$ER_ID, data_cc$Case.status, data_cc$Pathogen, data_cc$Reassigned.Pair.ID, r.c, h.c, pielou.c)
colnames(div.c)=c("ER_ID", "Case.status", "Pathogen","Pair","Richness", "Shannon", "Pielou")

div.c.cf <- div.c %>%
  filter(!grepl('Control', Case.status))


########################################################
# Perform Normality Test with Shapiro-Wilks & Visualization
########################################################

# Density Plots

ggdensity(div.c$Richness,
          main='Density Plot - Richness',
          xlab='Richness')

ggdensity(div.c$Shannon,
          main='Density Plot - Shannon Diversity',
          xlab='Shannon Diversity')

ggdensity(div.c$Pielou,
          main='Density Plot - Evenness',
          xlab='Evenness')


# Q-Q Plots

ggqqplot(div.c$Richness,
         main='Q-Q Plot - Richness')

ggqqplot(div.c$Shannon,
         main='Q-Q Plot - Shannon Diversity')

ggqqplot(div.c$Pielou,
         main='Q-Q Plot - Evenness')


# Shapiro-Wilks Test

shapiro.test(div.c$Richness)

shapiro.test(div.c$Shannon)

shapiro.test(div.c$Pielou)

