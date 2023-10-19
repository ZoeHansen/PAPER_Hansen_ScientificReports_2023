#####################################
# R Code for Data-wrangling with Tidyverse - Microbiome
# Author : Zoe Hansen
#####################################
# This code is designed to wrangle our large datasets from the ERIN
# metagenomics study. Here, I use the concatenated output from the 
# taxonomic classifier, Kaiju, and perform basic data wrangling and
# normalization (by # genome equivalents). 

# This code is used exclusively for organization and generating new
# CSV files for later use; there is not code for plotting/visualization here.

############## Loading the Data #########################

library(tidyverse)
library(plyr)

# We will use the 'SamDedupResistome' files since these have been deduplicated. If needed, the original
# resistome files are available. 

microbiome <- read.csv('D://Microbiome/Reads_based/ERIN_CaseControlFollow_kaiju_SPECIES.csv', header=TRUE)
microbiome <- microbiome[, colSums(microbiome != 0) > 0] 
microbiome <- microbiome[rowSums(microbiome != 0) > 0,] 

# Transpose our gene dataframe to make a column of sample IDs (helps with merging later)
all_seq <- setNames(data.frame(t(microbiome[,-1])), microbiome[,1])
all_seq <- all_seq %>%
  rownames_to_column(., "ER_ID")

### Import phenotypic data 

meta <- read.csv('D://ERIN_Metagenomes_Metadata.csv', header = TRUE)

############## Filtering data to exclude certain samples ##########################

# The 'ERIN_Metagenomes_Metadata.csv' includes only samples relevant to the downstream analysis
# i.e. samples with <50,000 reads, those that did not sequence/align well, highly contaminated, etc. have been removed from the sample list

# In order to filter these out from our resistome data, we need to use the ER_ID variable from our metadata
# file to isolate these samples

ID_use <- meta %>%
  dplyr::select(ER_ID)

# Now, we perform the filtering of the microbiome file(s):
all_seq_use <- left_join(ID_use, all_seq, by = 'ER_ID')
all_seq_use <- all_seq_use[, colSums(all_seq_use != 0) > 0] 
all_seq_use <- all_seq_use[rowSums(all_seq_use != 0) > 0,] 

# Now, we have a dataframe containing all of the samples of interest (n = 259), which exclude the following:
# - duplicates
# - read counts <50,000
# - samples with the case status 'Missing'
# - samples which were excluded in Aim One (Campylobacter) due to high contamination, low quality sequencing/alignment
# - any samples not included in our final spreadsheet 

write.csv(all_seq_use, 'D://Microbiome/Reads_based/ERIN_kaiju_reads_fullpath_casefollow_pairs.csv', row.names=FALSE)

############## Normalization by Genome Equivalents ############

# Load in the .csv file with the ouptut from MicrobeCensus
ags<- read.csv('D://MicrobeCensus_Output/mc_output_amrplusplus.csv', header=TRUE)
ags$ER_ID = as.character(ags$ER_ID)

# The suggested metric for normalization from MicrobeCensus is RPKG - reads per kilobase per genome equivalent
# To accomplish this, we must take (# reads mapped to gene)/ (gene length) / (# genome equivalents)

# Unfortunately, the RPKG metric is only calculable for the 'Gene' level resistome data; for all 
# other data (e.g. microbiome taxonomic levels) we will normalize by number of genome equivalents

# Now, we will use our genome equivalents metric to complete the normalization
# We must then use our MicrobeCensus data to isolate the number of genome equivalents:

genome_equiv <- ags %>%
  dplyr::select(ER_ID, GenomeEquivalents) 
genome_equiv_use <- left_join(ID_use, genome_equiv, by='ER_ID')

microbiome.wGE <- left_join(genome_equiv_use, all_seq_use, by = 'ER_ID')


# Perform the normalization by # genome equivalents in a sample
ge_normalized <- microbiome.wGE %>%
  mutate_at(vars(c(-GenomeEquivalents, -ER_ID)), list(~(. / GenomeEquivalents)))%>%
  dplyr::select(-GenomeEquivalents) 

########### Taxonomic Levels ###########

# Now that we have our normalized reads, we can isolate individual taxonomic levels of interest from
# our dataframes. 

ge_t <- ge_normalized %>%
  gather(key = key, value = value, 2:ncol(ge_normalized)) %>%
  spread(key=names(ge_normalized)[1], value = 'value') %>%
  dplyr::rename(., taxa_name=key)

# KINGDOM (isolate from PHYLUM file)
king <- ddply(full_path, "Kingdom", numcolwise(sum))
king<- king[, colSums(king[,-1] != 0) > 0] 
king <- king[rowSums(king[,-1] != 0) > 0,]
king[1, 1] <- "Cannot_be_assigned"

# PHYLUM file
full_path <- separate(ge_t, col= taxa_name, into=c('Kingdom','Phylum','extra'), sep="\\;")
full_path[179,2] <- 'cannot be assigned a (non-viral) phylum'
full_path[208,2] <- 'unclassified'
phy <- full_path %>%
  select(., -c(Kingdom, extra))

# CLASS file
full_path <- separate(ge_t, col= taxa_name, into=c('Kingdom','Phylum','Class','extra'), sep="\\;")
full_path[122,3] <- 'cannot be assigned a (non-viral) class'
full_path[207,3] <- 'unclassified'
class <- full_path %>%
  select(., -c(Kingdom, Phylum, extra))

# ORDER file
full_path <- separate(ge_t, col= taxa_name, into=c('Kingdom','Phylum','Class','Order','extra'), sep="\\;")
full_path[263,4] <- 'cannot be assigned a (non-viral) order'
full_path[495,4] <- 'unclassified'
ord <- full_path %>%
  select(., -c(Kingdom, Phylum, Class, extra))

# FAMILY file
full_path <- separate(ge_t, col= taxa_name, into=c('Kingdom','Phylum','Class','Order','Family','extra'), sep="\\;")
full_path[600,5] <- 'cannot be assigned a (non-viral) family'
full_path[1066,5] <- 'unclassified'
fam <- full_path %>%
  select(., -c(Kingdom, Phylum, Class, Order, extra))

# GENUS file
full_path <- separate(ge_t, col= taxa_name, into=c('Kingdom','Phylum','Class','Order','Family','Genus','extra'), sep="\\;")
full_path[3429,6] <- 'cannot be assigned a (non-viral) genus'
full_path[4297,6] <- 'unclassified'
genus <- full_path %>%
  select(., -c(Kingdom, Phylum, Class, Order, Family, extra))

# SPECIES file
full_path <- separate(ge_t, col= taxa_name, into=c('Kingdom','Phylum','Class','Order','Family','Genus','Species','extra'), sep="\\;")
full_path[35646,7] <- 'cannot be assigned a (non-viral) species'
full_path[37373,7] <- 'unclassified'
species <- full_path %>%
  select(., -c(Kingdom, Phylum, Class, Order, Family, Genus, extra))

# Transpose (if needed -- not ideal due to larg number of columns)
full_path_t <- phy %>%
  gather(key = key, value = value, 2:ncol(phy)) %>%
  spread(key=names(phy)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

write.csv(species, 'D://Microbiome/Reads_based/GEnorm_ERIN_kaiju_CaseControlFollow_SPECIES.csv', row.names = FALSE)


##### Grouping by a taxonomic rank in a "full_path" file #######
#If trying to 'group_by' when using a full_path file (not the method used above)
# the following code is equipped to do so:

class <- ddply(full_path, "Class", numcolwise(sum))
class<- class[, colSums(class[,-1] != 0) > 0] 
class <- class[rowSums(class[,-1] != 0) > 0,]
class[1, 1] <- "Cannot_be_assigned"

# For example, we needed to do this for the "KINGDOM" rank above, as Kaiju does not allow
# direct output of kingdom information :( 

