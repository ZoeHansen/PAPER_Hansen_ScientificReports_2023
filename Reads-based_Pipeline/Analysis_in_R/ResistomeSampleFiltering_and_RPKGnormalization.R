########## R Code for Data-wrangling with Tidyverse ##########

# Author : Zoe Hansen

# This code is designed to wrangle our large datasets from the ERIN metagenomics study. I've included a number
# of sections in this code which wrangles the data in a certain way for a certain step in our analysis

# In a number of these sections, there is code for plotting our data. In others, it is just data organization

############## Loading the Data #########################

library(tidyverse)
library(plyr)

# We will use the 'SamDedupResistome' files since these have been deduplicated. If needed, the original
# resistome files are available. 

resistome <- read.csv('D://Resistome/Reads_based/SamDedup_AMR_analytic_matrix.csv', header=TRUE)

# Transpose our gene dataframe to make a column of sample IDs (helps with merging later)
all_seq <- resistome %>%
  gather(key = key, value = value, 2:ncol(resistome)) %>%
  spread(key=names(resistome)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(.,ER_ID=key)

############## Import Phenotypic Data ################

# Here, we will import a data file with all of the phenotypic information/metadata
meta <- read.csv('D://ERIN_Metagenomes_Metadata.csv', header=TRUE)

############## Filtering data to exclude certain samples ##########################

# The 'ERIN_Metagenomes_Metadata.csv' has been modified to include only samples which should be included
# in the downstream analysis
      # i.e. samples with <50,000 reads, those that did not sequence/align well, highly contaminated, etc. 
      # have been removed from the sample list

# In order to filter these out from our resistome data, we need to use the ER_ID variable from our metadata
# file to isolate these samples

ID_use <- meta %>%
  select(ER_ID)

# Now, we perform the filtering of the SamDedupResistome file(s):
for (i in ID_use){
  resistome_filt <- all_seq %>%
    dplyr::filter(ER_ID %in% i)
}

# Now, we have a dataframe containing all of the samples of interest (n = 259), which exclude the following:
# - duplicates
# - read counts <50,000
# - samples with the case status 'Missing'
# - samples which were excluded in Aim One (Campylobacter) due to high contamination, low quality sequencing/alignment
# - any samples not included in our final spreadsheet 

############## Normalization by Genome Equivalents ############

# Load in the .csv file with the ouptut from MicrobeCensus
ags <- read.csv('D://MicrobeCensus_Output/mc_output_amrplusplus.csv', header=TRUE)

ags$ER_ID = as.character(ags$ER_ID)

# The suggested metric for normalization from MicrobeCensus is RPKG - reads per kilobase per genome equivalent
# To accomplish this, we must take (# reads mapped to gene)/ (gene length) / (# genome equivalents)

# Unfortunately, the RPKG metric is only calculable for the 'Gene' level data; for all other levels (e.g. 
# Group, Class, Type), read counts will be normalized by the number of genome equivalents in each sample

######## Gene Length ########
# To find the gene length, I have used Tablet with the annotations.csv file from MEGARes 2.0. 

length <- read.csv('D://megares_gene_lengths.csv', header=TRUE)

# Now, we will divide the gene count values by the length of the gene. To do this, we need to use the 
# non-transposed version that contains all genes in the first column:

resistome.wlength <- left_join(resistome, length, by = 'Gene')

resistome.length.norm <- resistome.wlength %>%
  mutate_at(vars(c(-Gene, -Length, -Length_kb)), list(~(. / Length_kb))) %>%
  select(-Length, -Length_kb)

# To match our other file, we need to transpose:
resistome.length.norm <- resistome.length.norm %>%
  gather(key = key, value = value, 2:ncol(resistome.length.norm)) %>%
  spread(key=names(resistome.length.norm)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  rename(., ER_ID=key)

### We need to filter out unneeded samples from this dataframe
ID_use <- meta %>%
  select(ER_ID)

# Now, we perform the filtering of the SamDedupResistome file(s):
for (i in ID_use){
  resistome.length.norm.filt <- resistome.length.norm %>%
    dplyr::filter(ER_ID %in% i)
}

######## Normalization w/ Genome Equivalents ########
# Now, we will use our genome equivalents metric to complete the normalization
# We must then use our MicrobeCensus data to isolate the number of genome equivalents:
genome_equiv <- ags %>%
  select(ER_ID, GenomeEquivalents) 

#resistome.wGE <- left_join(resistome.length.norm.filt, genome_equiv, by = 'ER_ID')
resistome.wGE <- left_join(resistome_filt, genome_equiv, by = 'ER_ID') # For Group, Class, and Type

rpkg_normalized <- resistome.wGE %>%
  mutate_at(vars(c(-GenomeEquivalents, -ER_ID)), list(~(. / GenomeEquivalents))) %>%
  select(-GenomeEquivalents) 


write.csv(rpkg_normalized, 'D://Resistome/Reads_based/ERIN_GEnorm_resistome_fullpath.csv', row.names = FALSE)

########## Other Levels ##########
# Now that we have our normalized reads, we can use the group_by() function to get counts relevant to the 
# group, mechanism, class, and type levels. 

rpkg_normalized <- read.csv('I://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_resistome_gene.csv',
                            heade=TRUE)

rpkg_t <- rpkg_normalized %>%
  gather(key = key, value = value, 2:ncol(rpkg_normalized)) %>%
  spread(key=names(rpkg_normalized)[1], value = 'value') %>%
  dplyr::rename(., ARG_name=key)

full_path <- separate(rpkg_t, col= ARG_name, into=c('Gene','Type','Class','Mechanism','Group','ConfirmationNeeded'), sep="\\|")

snps2 <- na.omit(full_path)
write.csv(snps2, 'I://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_RequiresSNPconfirmation_GEnorm_fullpath.csv',
          row.names = FALSE)

snps2_genes <- snps2$Gene
ex_snps2 <- full_path %>%
  filter(!(Gene %in% snps2_genes))
write.csv(ex_snps, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_RequiresSNPconfirmation_EXCLUDED_GEnorm_fullpath.csv',
          row.names=FALSE)

### Finding gene, group, mechanism, class, and type level counts for all samples

full_path2 <- full_path %>%
  select(., -ConfirmationNeeded)
full_path2<- full_path2[, colSums(full_path2[,-1] != 0) > 0] 
full_path2 <- full_path2[rowSums(full_path2[,-1] != 0) > 0,]

genes <- ddply(full_path2, "Gene", numcolwise(sum))
genes<- genes[, colSums(genes[,-1] != 0) > 0] 
genes <- genes[rowSums(genes[,-1] != 0) > 0,]

operons <- ddply(full_path2,"Group",numcolwise(sum))
operons<- operons[, colSums(operons[,-1] != 0) > 0] 
operons <- operons[rowSums(operons[,-1] != 0) > 0,] 

mechs <- ddply(full_path2,"Mechanism",numcolwise(sum))
mechs<- mechs[, colSums(mechs[,-1] != 0) > 0] 
mechs <- mechs[rowSums(mechs[,-1] != 0) > 0,] 

classes <- ddply(full_path2,"Class",numcolwise(sum))
classes<- classes[, colSums(classes[,-1] != 0) > 0] 
classes <- classes[rowSums(classes[,-1] != 0) > 0,] 

types <- ddply(full_path2, "Type", numcolwise(sum))
types <- types[, colSums(types[,-1] != 0) >0]
types <- types[rowSums(types[,-1] != 0) >0, ]

# Transpose each level and save to a CSV file for storage/later use
full_path_t <- types %>%
  gather(key = key, value = value, 2:ncol(types)) %>%
  spread(key=names(types)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

write.csv(full_path_t, 'D://Resistome/Reads_based/ERIN_GEnorm_ALL_TYPES_type_level.csv', row.names = FALSE)
