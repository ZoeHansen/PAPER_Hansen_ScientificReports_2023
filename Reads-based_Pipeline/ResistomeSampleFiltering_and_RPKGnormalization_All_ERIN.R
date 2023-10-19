########## R Code for Data-wrangling with Tidyverse ##########

# Author : Zoe Hansen
# Last modified : 2021.10.04

# This code is designed to wrangle our large datasets from the ERIN metagenomics study. I've included a number
# of sections in this code which wrangles the data in a certain way for a certain step in our analysis

# In a number of these sections, there is code for plotting our data. In others, it is just data organization


############## Loading the Data #########################

# First, load libraries

library(tidyverse)
library(plyr)

# Since I've already merged the resistome/microbiome data for all runs in Python, we just need to load
# this data:

# We will use the 'SamDedupResistome' files since these have been deduplicated. If needed, the original
# resistome files are available. 

resistome <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/SamDedup_AMR_analytic_matrix.csv',
                      header=TRUE)

# Since this file contains "|" as a delimiter between identification levels, we can separate them as needed
# to get our count information

broken <- separate(resistome, col= GENE, into=c('Gene','Type','Class','Mechanism','Group','ConfirmationNeeded'), sep="\\|")

# Isolate all genes that "Require SNP Confirmation"
snps <- na.omit(broken)
write.csv(snps, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_RequiresSNPconfirmation_fullpath.csv',
          row.names = FALSE)

# Creating a matrix that excludes genes that "Require SNP Confirmation
snp_genes <- snps$Gene

ex_snps <- broken %>%
  filter(!(Gene %in% snp_genes))
write.csv(ex_snps, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_RequiresSNPconfirmation_EXCLUDED_fullpath.csv',
          row.names=FALSE)


# Transpose our gene dataframe to make a column of sample IDs (helps with merging later)
all_seq <- resistome %>%
  gather(key = key, value = value, 2:ncol(resistome)) %>%
  spread(key=names(resistome)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(.,ER_ID=key)


############## Import Phenotypic Data ################

# Here, we will import a data file with all of the phenotypic information/metadata

meta <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/SecondAnalysis_MEGARes_v1/ERIN_Metagenomes_Metadata.csv',
                 header=TRUE)


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
ags <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/SecondAnalysis_MEGARes_v1//MicrobeCensus_Output/mc_output_amrplusplus.csv',
                   header=TRUE)

ags$ER_ID = as.character(ags$ER_ID)


# The suggested metric for normalization from MicrobeCensus is RPKG - reads per kilobase per genome equivalent
# To accomplish this, we must take (# reads mapped to gene)/ (gene length) / (# genome equivalents)

# Unfortunately, the RPKG metric is only calculable for the 'Gene' level data; for all other levels (e.g. 
# Group, Class, Type), read counts will be normalized by the number of genome equivalents in each sample


######## Gene Length (2nd Analysis Only) ########
# To find the gene length, I have used Tablet with the annotations.csv file from MEGARes 2.0. -- this was only 
# completed in the "Second Analysis" which actually has MEGARes v1.0 data. 

length <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/Second_Analysis/megares_gene_lengths.csv',
                   header=TRUE)

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


write.csv(rpkg_normalized, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_resistome_fullpath.csv',
          row.names = FALSE)

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

full_path_t <- types %>%
  gather(key = key, value = value, 2:ncol(types)) %>%
  spread(key=names(types)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

write.csv(full_path_t, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_ALL_TYPES_type_level.csv',
          row.names = FALSE)

##### Selecting for DRUGS and Multi-Compounds relevant to DRUGS ######
drugs_multi <- full_path2 %>%
  filter(., Type == 'Drugs'| Type == 'Multi-compound')%>%
  filter(., Class != 'Biocide_and_metal_resistance')

dm_gene <- ddply(drugs_multi, "Gene", numcolwise(sum))
dm_gene<- dm_gene[, colSums(dm_gene[,-1] != 0) > 0] 
dm_gene <- dm_gene[rowSums(dm_gene[,-1] != 0) > 0,]

dm_group <- ddply(drugs_multi, "Group", numcolwise(sum))
dm_group<- dm_group[, colSums(dm_group[,-1] != 0) > 0] 
dm_group <- dm_group[rowSums(dm_group[,-1] != 0) > 0,]

dm_mech <- ddply(drugs_multi, "Mechanism", numcolwise(sum))
dm_mech<- dm_mech[, colSums(dm_mech[,-1] != 0) > 0] 
dm_mech <- dm_mech[rowSums(dm_mech[,-1] != 0) > 0,]

dm_class <- ddply(drugs_multi, "Class", numcolwise(sum))
dm_class<- dm_class[, colSums(dm_class[,-1] != 0) > 0] 
dm_class <- dm_class[rowSums(dm_class[,-1] != 0) > 0,]

drugs_multi_t <- dm_class %>%
  gather(key = key, value = value, 2:ncol(dm_class)) %>%
  spread(key=names(dm_class)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

write.csv(drugs_multi_t, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_DRUG_and_MULTI_class_level.csv',
          row.names = FALSE)

###### Separating by MEGARESv2 Types ######

drugs <- full_path2 %>%
  filter(., Type == 'Drugs')
drugs<- drugs[, colSums(drugs[,-1] != 0) > 0] 
drugs <- drugs[rowSums(drugs[,-1] != 0) > 0,] 

drugs1 <- ddply(drugs,"Gene",numcolwise(sum))
drugs1<- drugs1[, colSums(drugs1[,-1] != 0) > 0] 
drugs1 <- drugs1[rowSums(drugs1[,-1] != 0) > 0,] 

drugs2 <- ddply(drugs_sub,"Group",numcolwise(sum))
drugs2<- drugs2[, colSums(drugs2[,-1] != 0) > 0] 
drugs2 <- drugs2[rowSums(drugs2[,-1] != 0) > 0,] 

drugs3 <- ddply(drugs_sub,"Mechanism",numcolwise(sum))
drugs3<- drugs3[, colSums(drugs3[,-1] != 0) > 0] 
drugs3 <- drugs3[rowSums(drugs3[,-1] != 0) > 0,] 

drugs4 <- ddply(drugs_sub,"Class",numcolwise(sum))
drugs4<- drugs4[, colSums(drugs4[,-1] != 0) > 0] 
drugs4 <- drugs4[rowSums(drugs4[,-1] != 0) > 0,] 

drugs_t <- drugs3 %>%
  gather(key = key, value = value, 2:ncol(drugs3)) %>%
  spread(key=names(drugs3)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

write.csv(drugs_t, 'I://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_DRUGS_mechlevel.csv',
          row.names = FALSE)


metals <- full_path2 %>% 
  filter(., Type == 'Metals')
metals<- metals[, colSums(metals[,-1] != 0) > 0] 
metals <- metals[rowSums(metals[,-1] != 0) > 0,] 

metals1 <- ddply(metals,"Gene",numcolwise(sum))
metals1<- metals1[, colSums(metals1[,-1] != 0) > 0] 
metals1 <- metals1[rowSums(metals1[,-1] != 0) > 0,] 

metals2 <- ddply(metals,"Group",numcolwise(sum))
metals2<- metals2[, colSums(metals2[,-1] != 0) > 0] 
metals2 <- metals2[rowSums(metals2[,-1] != 0) > 0,] 

metals3 <- ddply(metals,"Mechanism",numcolwise(sum))
metals3<- metals3[, colSums(metals3[,-1] != 0) > 0] 
metals3 <- metals3[rowSums(metals3[,-1] != 0) > 0,] 

metals4 <- ddply(metals,"Class",numcolwise(sum))
metals4<- metals4[, colSums(metals4[,-1] != 0) > 0] 
metals4 <- metals4[rowSums(metals4[,-1] != 0) > 0,] 

metals_t <- metals4 %>%
  gather(key = key, value = value, 2:ncol(metals4)) %>%
  spread(key=names(metals4)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

write.csv(metals_t, 'I://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_METALS_classlevel.csv',
          row.names = FALSE)



biocides <- full_path2 %>%
  filter(., Type == 'Biocides')
biocides<- biocides[, colSums(biocides[,-1] != 0) > 0] 
biocides <- biocides[rowSums(biocides[,-1] != 0) > 0,] 

biocides1 <- ddply(biocides,"Gene",numcolwise(sum))
biocides1<- biocides1[, colSums(biocides1[,-1] != 0) > 0] 
biocides1 <- biocides1[rowSums(biocides1[,-1] != 0) > 0,] 

biocides2 <- ddply(biocides,"Group",numcolwise(sum))
biocides2<- biocides2[, colSums(biocides2[,-1] != 0) > 0] 
biocides2 <- biocides2[rowSums(biocides2[,-1] != 0) > 0,] 

biocides3 <- ddply(biocides,"Mechanism",numcolwise(sum))
biocides3<- biocides3[, colSums(biocides3[,-1] != 0) > 0] 
biocides3 <- biocides3[rowSums(biocides3[,-1] != 0) > 0,] 

biocides4 <- ddply(biocides,"Class",numcolwise(sum))
biocides4<- biocides4[, colSums(biocides4[,-1] != 0) > 0] 
biocides4 <- biocides4[rowSums(biocides4[,-1] != 0) > 0,] 

biocides_t <- biocides4 %>%
  gather(key = key, value = value, 2:ncol(biocides4)) %>%
  spread(key=names(biocides4)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

write.csv(biocides_t, 'I://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_BIOCIDES_classlevel.csv',
          row.names = FALSE)

multi_comp <- full_path2 %>%
  filter(., Type == 'Multi-compound')
multi_comp<- multi_comp[, colSums(multi_comp[,-1] != 0) > 0] 
multi_comp <- multi_comp[rowSums(multi_comp[,-1] != 0) > 0,] 

multi_comp1 <- ddply(multi_comp,"Gene",numcolwise(sum))
multi_comp1<- multi_comp1[, colSums(multi_comp1[,-1] != 0) > 0] 
multi_comp1 <- multi_comp1[rowSums(multi_comp1[,-1] != 0) > 0,] 

multi_comp2 <- ddply(multi_comp,"Group",numcolwise(sum))
multi_comp2<- multi_comp2[, colSums(multi_comp2[,-1] != 0) > 0] 
multi_comp2 <- multi_comp2[rowSums(multi_comp2[,-1] != 0) > 0,] 

multi_comp3 <- ddply(multi_comp,"Mechanism",numcolwise(sum))
multi_comp3<- multi_comp3[, colSums(multi_comp3[,-1] != 0) > 0] 
multi_comp3 <- multi_comp3[rowSums(multi_comp3[,-1] != 0) > 0,] 

multi_comp4 <- ddply(multi_comp,"Class",numcolwise(sum))
multi_comp4<- multi_comp4[, colSums(multi_comp4[,-1] != 0) > 0] 
multi_comp4 <- multi_comp4[rowSums(multi_comp4[,-1] != 0) > 0,] 

multi_comp_t <- multi_comp4 %>%
  gather(key = key, value = value, 2:ncol(multi_comp4)) %>%
  spread(key=names(multi_comp4)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

write.csv(multi_comp_t, 'I://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_MULTICOMP_classlevel.csv',
          row.names = FALSE)

