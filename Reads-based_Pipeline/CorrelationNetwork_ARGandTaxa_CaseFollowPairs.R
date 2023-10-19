#################################################

# Correlation Networks - Constructing CSV files for Gephi
# AIM 2 - Microbiome and Resistome

# Author: Zoe Hansen (via Karla Vasco)
# Last Modified: 2021.10.27

#################################################

### Load libraries and data ###

library(tidyverse)
library(Hmisc)
library(plyr)


##################### Case - Follow Global Network #########################

# Data should contain samples in the rows and metadata fields, args, and taxa in the columns

# Metadata
meta <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/ERIN_Metagenomes_Metadata_60_CaseFollowPairs.csv', 
                 header = TRUE)

meta_ids <- meta %>%
  select(ER_ID, Case.status)

# Normalized Actual abundance ARG data
args <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_ALL_TYPES_group_level_CaseFollowPairs.csv',
                 header = TRUE)

args_t <- args %>%
  gather(key = key, value = value, 2:ncol(args)) %>%
  spread(key=names(args)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)


# Isolate ARGs that occur in at least 10% samples (~50% of total ARGs in Case/Follow)
args_sum <- args_t %>%
  mutate(Abundance = rowSums(args_t[,-1] !=0 ))

args_keep <- args_sum[which(args_sum$Abundance > 12),] %>%
  select(-Abundance)

args_sub <- args_keep %>%
  gather(key = key, value = value, 2:ncol(args_keep)) %>%
  spread(key=names(args_keep)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

# Subset ARGs file to contain only 60 Case-Follow Pairs
args.cf <- left_join(meta_ids, args_sub, by='ER_ID')
args.cf <- args.cf[, colSums(args.cf[,-c(1:2)] != 0) > 0] 
args.cf <- args.cf[rowSums(args.cf[,-c(1:2)] != 0) > 0,]



# Normalized Actual abundance Taxonomic data
taxonomy <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/GEnorm_ERIN_kaiju_60CaseFollowPairs_GENUS.csv',
                 header=TRUE)

# Fix and Transpose our taxa file

taxonomy$Genus <- taxonomy$Genus %>%
  replace_na("Unknown")

taxa <- ddply(taxonomy, "Genus", numcolwise(sum))

## For some reason, my 'taxa_cutoff' code is removing the most abundant genera (e.g. Escherichia)....
# Maybe subset based on total abundance? To do this, we should split case and follow up though probably :(
# So, maybe split cases and follow-ups then perform a cutoff based on the total abundances in each respective group? 
# Do this for ARGs too. Have fun. 

# Isolate taxa that occur in 20% number of samples (results in ~50% of total genera in Case/Follow)
### THIS WORKS #####

taxa_sum <- taxa %>%
  mutate(Abundance = rowSums(taxa[,-1] !=0 ))

taxa_keep <- taxa_sum[which(taxa_sum$Abundance > 20),] %>%
  select(-Abundance)


taxa_t <- taxa_keep %>%
  gather(key = key, value = value, 2:ncol(taxa_keep)) %>%
  spread(key=names(taxa_keep)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

taxa.cf <- left_join(meta_ids, taxa_t, by = 'ER_ID')
taxa.cf <- taxa.cf[, colSums(taxa.cf[,-c(1:2)] != 0) > 0] 
taxa.cf <- taxa.cf[rowSums(taxa.cf[,-c(1:2)] != 0) > 0,]



### Merging all of the data together ###

taxa_args <- left_join(taxa.cf, args.cf, by=c('ER_ID', 'Case.status')) 

meta_taxa_args <- left_join(meta_ids, taxa_args, by=c('ER_ID', 'Case.status'))


##### Create the "nodes" dataframe #####

### CASES ### (Repeat for FollowUps)

# ARGs
args_meta <- left_join(meta_ids, args.cf, by = c('ER_ID','Case.status'))

ARGs_case <- args_meta %>%
  filter(Case.status == 'FollowUp')%>%
  select(-Case.status)

ARGs_case.t <- ARGs_case %>%
  gather(key = key, value = value, 2:ncol(ARGs_case)) %>%
  spread(key=names(ARGs_case)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., Group=key)

ARGs_case.group <- ARGs_case.t %>%
  mutate(Abundance = rowSums(.[,-1]))%>%
  mutate(Type='ARG')%>%
  select(Type, Group, Abundance)


# Taxa
taxa_meta <- left_join(meta_ids, taxa.cf, by = c('ER_ID','Case.status'))

taxa_case <- taxa_meta %>%
  filter(Case.status == 'FollowUp')%>%
  select(-Case.status)

taxa_case.t <- taxa_case %>%
  gather(key = key, value = value, 2:ncol(taxa_case)) %>%
  spread(key=names(taxa_case)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., Group=key)

taxa_case.group <- taxa_case.t %>%
  mutate(Abundance = rowSums(.[,-1]))%>%
  mutate(Type='Taxa')%>%
  select(Type, Group, Abundance)


### Join the taxa and ARG data
nodes_ARGs_taxa <- dplyr::union(ARGs_case.group, taxa_case.group)

nrow(nodes_ARGs_taxa)


#Adding id numbers required for Gephi (corresponding to node id)
nodes_gephi <- nodes_ARGs_taxa %>%
  mutate(ID = 1:nrow(.))%>% # Create an ID associated with each ARG or Taxa
  mutate(Label = Group)%>%
  select(Label, ID, Abundance, everything()) # Reordering columns

#Saving final file as csv
write_csv(nodes_gephi, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_Global/FollowUp/nodes_GENUS_ARGgroup_FollowUp_cutoff12samples_20220126.csv')


##### Create the "Edges" dataframe #####

### CASE (repeat for FollowUp)

#Abundance matrix for a group (Here is samples from the location "Punta Espinoza", but it could be "cases" or "controls")
abundances_group <- meta_taxa_args %>% 
  filter(Case.status == 'FollowUp')%>%
  select(3:ncol(.))%>% # removing metadata from table to isolate abundances
  select(order(colnames(.))) # order the column names


abundances_group[is.na(abundances_group)] <- 0 #adding 0 if there is no value

abd.mat <- as.matrix(abundances_group[sapply(abundances_group, is.numeric)]) #ensure data are numeric

### Calculating correlations and formatting tables for Gephi

correlation_group<- rcorr(abd.mat, type = 'spearman')

#Changing table format
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat) # takes the upper triangle of the correlation matrix
  data.frame( # creates a new dataframe
    row = rownames(cormat)[row(cormat)[ut]], # takes rownames of the upper triangle for the correlation matrix
    column = rownames(cormat)[col(cormat)[ut]], # compiles column names of the upper triangle of correlation matrix
    cor  =(cormat)[ut], # records correlation values of the row x column values comparison
    p = pmat[ut] # records p-value of the correlation between row x column values
  )
}

edges_group = flattenCorrMatrix(correlation_group$r, correlation_group$P)


# Perform p-value adjustment for multiple comparisons:

q_val <- p.adjust(edges_group$p, method='BH')

edges_group_adj<- cbind(edges_group, q_val)

#Renaming columns, adding correlation type, and filtering correlations higher than 0.75
edges_group_clean <- edges_group_adj %>% 
  dplyr::rename(Label = row, Target = column, Correlation = cor, p_value=p) %>% 
  mutate(Type = "undirected") %>% 
  filter(Correlation >= 0.80)#%>%  # set the Spearman correlation cutoff as 0.80
#  filter(q_value < 0.01)

#Changing labels for id numbers (corresponding to node id)
edges_group_gephi <- left_join(edges_group_clean,nodes_gephi, by = "Label") %>%
  dplyr::rename(Source = ID) %>% 
  select(Source, Target, Correlation, Type.x) %>% 
  dplyr::rename(Label = Target)

edges_gephi <- left_join(edges_group_gephi,nodes_gephi, by = "Label") %>% 
  select(Source, ID, Correlation, Type.x) %>% 
  dplyr::rename(Target = ID, Type = Type.x)


##Exporting edge's file for Gephi

#Saving edges file as csv
write_csv(edges_group_clean, "D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_Global/FollowUp/edges_GENUS_ARGgroup_FollowUp_labels_80cutoff_50percTaxaARG_20220126.csv")

write_csv(edges_gephi, "D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_Global/FollowUp/edges_GENUS_ARGgroup_FollowUp_GEPHI_nolabels_80cutoff_50percTaxaARG_20220126.csv")





############ Beta-lactam Resistance Genes ###################

megares <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/megares_full_annotations_v2.csv',
                    header=TRUE)

# Isolate the betalactam group information
betalactams <- megares %>%
  filter(class == 'betalactams') %>%
  select(group) %>%
  group_by(group)

# combine duplicate groups to form a dataframe
betalactams.groups <- as.data.frame(unique(betalactams$group))
colnames(betalactams.groups) <- 'Group'

############################# Beta-lactam Method #1 #######################################

# This method strictly subsets beta-lactam ARGs from the overall ARG set. 
# This limits our analysis to only beta-lactam ARGs and not necessarily other co-occuring ARGs
# of interest (for this, see Method #2)
######### Subset Cases with highest abundance ############

# Subset our cases based on containment of beta-lactam ARGs to use in paired analysis

args.cases <- args.cf %>%
#  filter(Case.status == 'Case') %>%
  select(-Case.status)

args.case.t <- args.cases %>% 
  gather(key=key, value=value, 2:ncol(args.cases)) %>%
  spread(key=names(args.cases)[1], value='value') %>%
  dplyr::rename(., ER_ID=key)

# Isolate beta-lactam ARGs in cases only
case.betalactam <- left_join(betalactams.groups, args.case.t, by=c('Group'='ER_ID'))

# Remove ARGs that weren't in our dataset (NAs)
case.betalactam <- na.omit(case.betalactam)


# Try to determine which cases have highest abundance of beta-lactams
# Take subset containing these cases

betalactam.sums <- as.data.frame(colSums(case.betalactam[,-1]))
colnames(betalactam.sums) <- 'Sum'

top.beta.abd <- betalactam.sums %>%
  rownames_to_column(., 'ER_ID') %>%
  arrange(.,desc(Sum))%>%
  slice(1:10)

top.beta.ids <- top.beta.abd %>%
  select(ER_ID)


########## Find associated Follow-up to for the Pair #############

# Make the case subset to include in analysis and find Pair ID:

meta.pair <- meta %>%
  select(ER_ID, Case.status, Case.Follow_ID, Pathogen)

meta.case.sub <- left_join(top.beta.ids, meta.pair, by='ER_ID')


meta.pairs.sub <- left_join(meta.case.sub, meta.pair, by='Case.Follow_ID')

meta.pairs.sub <- meta.pairs.sub %>%
  select(-ER_ID.x, -Case.status.x, -Pathogen.x)%>%
  dplyr::rename(., ER_ID=ER_ID.y, Case.status=Case.status.y, Pathogen=Pathogen.y)


col_order <- c("ER_ID", "Case.status", "Case.Follow_ID","Pathogen")

meta.pairs.sub <- meta.pairs.sub[, col_order]


# Subset our ARGs list based on the pairs we are interested in 
args.pair <- subset(args.cf, ER_ID %in% meta.pairs.sub$ER_ID) %>%
  select(-Case.status)

args.pair.t <- args.pair %>% 
  gather(key=key, value=value, 2:ncol(args.pair)) %>%
  spread(key=names(args.pair)[1], value='value') %>%
  dplyr::rename(., ER_ID=key)

# Isolate beta-lactam ARGs in our 10 pairs
pair.betalactam <- left_join(betalactams.groups, args.pair.t, by=c('Group'='ER_ID'))

# Remove ARGs that show 'NA'
pair.betalactam <- na.omit(pair.betalactam)

pair.betalactam.t <- pair.betalactam %>% 
  gather(key=key, value=value, 2:ncol(pair.betalactam)) %>%
  spread(key=names(pair.betalactam)[1], value='value') %>%
  dplyr::rename(., ER_ID=key)


# Want to figure out how to get rid of columns that sum to 0....
# colSums is not working...ugh. 

pairs.final <- left_join(meta.pairs.sub, pair.betalactam.t, by = "ER_ID")

write.csv(pairs.final, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_CaseFollowPairs_Betalactams/Betalactam_ARGs_Top10Pairs_viaCases_20220122.csv',
          row.names=FALSE)

################# Create Correlation Matrices for Pairs ###################

# Metadata
meta <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/ERIN_Metagenomes_Metadata_60_CaseFollowPairs.csv', 
                 header = TRUE)

meta_ids <- meta %>%
  select(ER_ID, Case.status)


# If reading in Top-10 pair information:

pair.beta <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_CaseFollowPairs_Betalactams/Betalactam_ARGs_Top10Pairs_viaCases_20220122.csv',
                      header=TRUE)

# If reading in all pair information: 

pair.beta <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_CaseFollowPairs_Betalactams/Betalactam_ARGs_AllCaseFollowPairs_20220123.csv',
                      header=TRUE)

# Same taxa code from above: 
# Normalized Actual abundance Taxonomic data
taxonomy <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/GEnorm_ERIN_kaiju_CaseControlFollow_GENUS.csv',
                     header=TRUE)

# Fix and Transpose our taxa file

taxonomy$Genus <- taxonomy$Genus %>%
  replace_na("Unknown")

taxa <- ddply(taxonomy, "Genus", numcolwise(sum))

taxa_t <- taxa %>%
  gather(key = key, value = value, 2:ncol(taxa)) %>%
  spread(key=names(taxa)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

taxa.cf <- left_join(meta_ids, taxa_t, by = 'ER_ID')
taxa.cf <- taxa.cf[, colSums(taxa.cf[,-c(1:2)] != 0) > 0] 
taxa.cf <- taxa.cf[rowSums(taxa.cf[,-c(1:2)] != 0) > 0,]

# Isolate taxa that occur in 20% number of samples (results in ~50% of total genera in Case/Follow)
taxa_cutoff <- taxa.cf
taxa_cutoff <- taxa_cutoff[, colSums(taxa_cutoff[,-c(1:2)] != 0) > 20] 

colSums(taxa.cf[,-c(1:2)])



# For now, can create correlations for these betalactam ARGs and taxa from above

beta_taxa_args <- left_join(pair.beta, taxa_cutoff, by=c('ER_ID'))
  

##### Create the "nodes" dataframe for each pair #####

### For Subset of Pairs only ###
# Retrieve ER_IDs from pair file:

pair.ids <- pairs.final %>%
  select(ER_ID, Case.Follow_ID)



### CASES ### (Repeat for FollowUps)

# Betalactam ARGs

beta_pair <- pair.beta %>%
#  filter(Case.Follow_ID == '3') %>%
  filter(Case.status == 'FollowUp')%>%
  select(-Case.status, -Case.Follow_ID)

beta_pair.t <- beta_pair %>%
  gather(key = key, value = value, 2:ncol(beta_pair)) %>%
  spread(key=names(beta_pair)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., Group=key)

beta_pair.group <- beta_pair.t %>%
  mutate(Abundance = rowSums(.[,-1]))%>%
  mutate(Type='ARG')%>%
  select(Type, Group, Abundance)


# Taxa
taxa_meta <- left_join(meta_ids, taxa_cutoff, by = c('ER_ID'))

taxa_pair <- taxa_meta %>%
#  filter(Case.Follow_ID == '3') %>%
  filter(Case.status == 'FollowUp')%>%
  select(-Case.status)

taxa_pair.t <- taxa_pair %>%
  gather(key = key, value = value, 2:ncol(taxa_pair)) %>%
  spread(key=names(taxa_pair)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., Group=key)

taxa_pair.group <- taxa_pair.t %>%
  mutate(Abundance = rowSums(.[,-1]))%>%
  mutate(Type='Taxa')%>%
  select(Type, Group, Abundance)


### Join the taxa and ARG data
nodes_beta_taxa <- dplyr::union(beta_pair.group, taxa_pair.group)

nrow(nodes_beta_taxa)


#Adding id numbers required for Gephi (corresponding to node id)
nodes_pair_gephi <- nodes_beta_taxa %>%
  mutate(ID = 1:nrow(.))%>% # Create an ID associated with each ARG or Taxa
  mutate(Label = Group)%>%
  select(Label, ID, Abundance, everything()) # Reordering columns

#Saving final file as csv
write_csv(nodes_pair_gephi, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_CaseFollowPairs_Betalactams/All_CaseFollowPairs/nodes_GENUS_BETALACTAM_AllCaseFollowPairs_FOLLOWUPS_actualabundance_20220123.csv')


##### Create the "Edges" dataframe for each pair #####

### CASE (repeat for FollowUp)

#Abundance matrix for a group (Here is samples from the location "Punta Espinoza", but it could be "cases" or "controls")
abundances_group.pair <- beta_taxa_args %>%
#  filter(Case.Follow_ID == '3')%>%
  filter(Case.status == 'FollowUp')%>%
  select(4:ncol(.))%>% # removing metadata from table to isolate abundances
  select(order(colnames(.))) # order the column names


abundances_group.pair[is.na(abundances_group.pair)] <- 0 #adding 0 if there is no value

abd.mat.pair <- as.matrix(abundances_group.pair[sapply(abundances_group.pair, is.numeric)]) #ensure data are numeric

### Calculating correlations and formatting tables for Gephi

correlation_group.pair<- rcorr(abd.mat.pair, type = 'spearman')

#Changing table format
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat) # takes the upper triangle of the correlation matrix
  data.frame( # creates a new dataframe
    row = rownames(cormat)[row(cormat)[ut]], # takes rownames of the upper triangle for the correlation matrix
    column = rownames(cormat)[col(cormat)[ut]], # compiles column names of the upper triangle of correlation matrix
    cor  =(cormat)[ut], # records correlation values of the row x column values comparison
    p = pmat[ut] # records p-value of the correlation between row x column values
  )
}

edges_group.pair = flattenCorrMatrix(correlation_group.pair$r, correlation_group.pair$P)


# Perform p-value adjustment for multiple comparisons:

q_val <- p.adjust(edges_group.pair$p, method='BH')

edges_group_pair_adj<- cbind(edges_group.pair, q_val)

#Renaming columns, adding correlation type, and filtering correlations higher than 0.75
edges_group_pair_clean <- edges_group_pair_adj %>% 
  dplyr::rename(Label = row, Target = column, Correlation = cor, p_value=p) %>% 
  mutate(Type = "undirected") %>% 
  filter(Correlation >= 0.80)#%>%  # set the Spearman correlation cutoff as 0.80
#  filter(q_value < 0.01)

#Changing labels for id numbers (corresponding to node id)
edges_group_pair_gephi <- left_join(edges_group_pair_clean, nodes_pair_gephi, by = "Label") %>%
  dplyr::rename(Source = ID) %>% 
  select(Source, Target, Correlation, Type.x) %>% 
  dplyr::rename(Label = Target)

edges_pair_gephi <- left_join(edges_group_pair_gephi,nodes_pair_gephi, by = "Label") %>% 
  select(Source, ID, Correlation, Type.x) %>% 
  dplyr::rename(Target = ID, Type = Type.x)


##Exporting edge's file for Gephi

#Saving edges file as csv
write_csv(edges_group_pair_clean, "D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_CaseFollowPairs_Betalactams/All_CaseFollowPairs/edges_GENUS_BETALACTAM_AllCaseFollowPairs_FOLLOWUPS_labels_80cutoff_actualabundance_20220123.csv")

write_csv(edges_pair_gephi, "D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_CaseFollowPairs_Betalactams/All_CaseFollowPairs/edges_GENUS_GENUS_BETALACTAM_AllCaseFollowPairs_FOLLOWUPS_GEPHI_nolabels_80cutoff_actualabundance_20220123.csv")









############################## Beta-lactam Method #2 ######################################

# This method attempts to extract any connection involving a beta-lactam-relevant ARG from the 
# global network dataset generated above. 

# Can use the same nodes file for input to GEPHI -- create new column designating which ARGs are beta-lactams
# so that we can color-code in Gephi

beta_nodes <- nodes_gephi 
  
beta_nodes$Type[beta_nodes$Label %in% betalactams.groups$Group] <- 'Betalactam'

write.csv(beta_nodes, "D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_Betalactam/FollowUp/nodes_GENUS_BetalactamARGgroup_FollowUp_50percTaxaARG_20220126.csv",
          row.names=FALSE)

# Use the edges_group_clean variable to extract all relevant connections with beta-lactam ARGs
# Use the betalactams.group variable for extracting these


beta.pair.edges <- edges_group_clean %>%
  filter(Label %in% betalactams.groups$Group | Target %in% betalactams.groups$Group)

#Changing labels for id numbers (corresponding to node id)
edges_beta_group_gephi <- left_join(beta.pair.edges,nodes_gephi, by = "Label") %>%
  dplyr::rename(Source = ID) %>% 
  select(Source, Target, Correlation, Type.x) %>% 
  dplyr::rename(Label = Target)

edges_beta_group_gephi <- left_join(edges_beta_group_gephi,nodes_gephi, by = "Label") %>% 
  select(Source, ID, Correlation, Type.x) %>% 
  dplyr::rename(Target = ID, Type = Type.x)

write.csv(beta.pair.edges, "D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_Betalactam/FollowUp/edges_GENUS_BETALACTAM_CONNECTIONS_FollowUp_labels_80cutoff_20220126.csv",
          row.names=FALSE)

write.csv(edges_beta_group_gephi, "D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_Betalactam/FollowUp/edges_GENUS_BETALACTAM_CONNECTIONS_FollowUp_GEPHI_nolabels_80cutoff_20220126.csv",
          row.names=FALSE)



######################### Networks by Pathogen ##############################

### Note: we can only perform the network construction for Salmonella and Campylobacter cases/follow-ups
# This is because the Spearman Rank correlation test requires >4 observations, and we have fewer than that
# for our shigella and STEC numbers (4 and 3, respectively).


# Data should contain samples in the rows and metadata fields, args, and taxa in the columns

# Metadata
meta <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/ERIN_Metagenomes_Metadata_60_CaseFollowPairs.csv', 
                 header = TRUE)

meta_ids <- meta %>%
  select(ER_ID, Case.status, Pathogen)

# Normalized Actual abundance ARG data
args <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_ALL_TYPES_group_level_CaseFollowPairs.csv',
                 header = TRUE)

args_t <- args %>%
  gather(key = key, value = value, 2:ncol(args)) %>%
  spread(key=names(args)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)


# Isolate ARGs that occur in at least 10% samples (~50% of total ARGs in Case/Follow)
args_sum <- args_t %>%
  mutate(Abundance = rowSums(args_t[,-1] !=0 ))

args_keep <- args_sum[which(args_sum$Abundance > 12),] %>%
  select(-Abundance)

args_sub <- args_keep %>%
  gather(key = key, value = value, 2:ncol(args_keep)) %>%
  spread(key=names(args_keep)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

# Subset ARGs file to contain only 60 Case-Follow Pairs
args.cf <- left_join(meta_ids, args_sub, by='ER_ID')
args.cf <- args.cf[, colSums(args.cf[,-c(1:2)] != 0) > 0] 
args.cf <- args.cf[rowSums(args.cf[,-c(1:2)] != 0) > 0,]



# Normalized Actual abundance Taxonomic data
taxonomy <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/GEnorm_ERIN_kaiju_60CaseFollowPairs_GENUS.csv',
                     header=TRUE)

# Fix and Transpose our taxa file

taxonomy$Genus <- taxonomy$Genus %>%
  replace_na("Unknown")

taxa <- ddply(taxonomy, "Genus", numcolwise(sum))

# Isolate taxa that occur in 20% number of samples (results in ~50% of total genera in Case/Follow)
### THIS WORKS #####

taxa_sum <- taxa %>%
  mutate(Abundance = rowSums(taxa[,-1] !=0 ))

taxa_keep <- taxa_sum[which(taxa_sum$Abundance > 20),] %>%
  select(-Abundance)


taxa_t <- taxa_keep %>%
  gather(key = key, value = value, 2:ncol(taxa_keep)) %>%
  spread(key=names(taxa_keep)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

taxa.cf <- left_join(meta_ids, taxa_t, by = 'ER_ID')
taxa.cf <- taxa.cf[, colSums(taxa.cf[,-c(1:2)] != 0) > 0] 
taxa.cf <- taxa.cf[rowSums(taxa.cf[,-c(1:2)] != 0) > 0,]



### Merging all of the data together ###

taxa_args <- left_join(taxa.cf, args.cf, by=c('ER_ID', 'Case.status','Pathogen')) 

meta_taxa_args <- left_join(meta_ids, taxa_args, by=c('ER_ID', 'Case.status','Pathogen'))


##### Create the "nodes" dataframe #####

### CASES ### (Repeat for FollowUps)

# ARGs
args_meta <- left_join(meta_ids, args.cf, by = c('ER_ID','Case.status','Pathogen'))

ARGs_case <- args_meta %>%
  filter(Pathogen == 'Campylobacter (CA)')%>%
  filter(Case.status == 'Case')%>%
  select(-Case.status, -Pathogen)

ARGs_case.t <- ARGs_case %>%
  gather(key = key, value = value, 2:ncol(ARGs_case)) %>%
  spread(key=names(ARGs_case)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., Group=key)

ARGs_case.group <- ARGs_case.t %>%
  mutate(Abundance = rowSums(.[,-1]))%>%
  mutate(Type='ARG')%>%
  select(Type, Group, Abundance)


# Taxa
taxa_meta <- left_join(meta_ids, taxa.cf, by = c('ER_ID','Case.status','Pathogen'))

taxa_case <- taxa_meta %>%
  filter(Pathogen == 'Campylobacter (CA)')%>%
  filter(Case.status == 'Case')%>%
  select(-Case.status, -Pathogen)

taxa_case.t <- taxa_case %>%
  gather(key = key, value = value, 2:ncol(taxa_case)) %>%
  spread(key=names(taxa_case)[1], value = 'value') %>%
  replace(is.na(.), 0)%>%
  dplyr::rename(., Group=key)

taxa_case.group <- taxa_case.t %>%
  mutate(Abundance = rowSums(.[,-1]))%>%
  mutate(Type='Taxa')%>%
  select(Type, Group, Abundance)


### Join the taxa and ARG data
nodes_ARGs_taxa <- dplyr::union(ARGs_case.group, taxa_case.group)

nrow(nodes_ARGs_taxa)


#Adding id numbers required for Gephi (corresponding to node id)
nodes_gephi <- nodes_ARGs_taxa %>%
  mutate(ID = 1:nrow(.))%>% # Create an ID associated with each ARG or Taxa
  mutate(Label = Group)%>%
  select(Label, ID, Abundance, everything()) # Reordering columns

#Saving final file as csv
write.csv(nodes_gephi, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_Global/Pathogen/Shigella/CASE_nodes_GENUS_ARGgroup_Shigella_20220127.csv',
          row.names=FALSE)


##### Create the "Edges" dataframe #####

### CASE (repeat for FollowUp)

#Abundance matrix for a group (Here is samples from the location "Punta Espinoza", but it could be "cases" or "controls")
abundances_group <- meta_taxa_args %>% 
  filter(Pathogen == 'Campylobacter (CA)')%>%
  filter(Case.status == 'Case')%>%
  select(4:ncol(.))%>% # removing metadata from table to isolate abundances
  select(order(colnames(.))) # order the column names


abundances_group[is.na(abundances_group)] <- 0 #adding 0 if there is no value

abd.mat <- as.matrix(abundances_group[sapply(abundances_group, is.numeric)]) #ensure data are numeric

### Calculating correlations and formatting tables for Gephi

correlation_group<- rcorr(abd.mat, type = 'spearman')

#Changing table format
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat) # takes the upper triangle of the correlation matrix
  data.frame( # creates a new dataframe
    row = rownames(cormat)[row(cormat)[ut]], # takes rownames of the upper triangle for the correlation matrix
    column = rownames(cormat)[col(cormat)[ut]], # compiles column names of the upper triangle of correlation matrix
    cor  =(cormat)[ut], # records correlation values of the row x column values comparison
    p = pmat[ut] # records p-value of the correlation between row x column values
  )
}

edges_group = flattenCorrMatrix(correlation_group$r, correlation_group$P)


# Perform p-value adjustment for multiple comparisons:

q_val <- p.adjust(edges_group$p, method='BH')

edges_group_adj<- cbind(edges_group, q_val)

#Renaming columns, adding correlation type, and filtering correlations higher than 0.75
edges_group_clean <- edges_group_adj %>% 
  dplyr::rename(Label = row, Target = column, Correlation = cor, p_value=p) %>% 
  mutate(Type = "undirected") %>% 
  filter(Correlation >= 0.80)#%>%  # set the Spearman correlation cutoff as 0.80
#  filter(q_value < 0.01)

#Changing labels for id numbers (corresponding to node id)
edges_group_gephi <- left_join(edges_group_clean,nodes_gephi, by = "Label") %>%
  dplyr::rename(Source = ID) %>% 
  select(Source, Target, Correlation, Type.x) %>% 
  dplyr::rename(Label = Target)

edges_gephi <- left_join(edges_group_gephi,nodes_gephi, by = "Label") %>% 
  select(Source, ID, Correlation, Type.x) %>% 
  dplyr::rename(Target = ID, Type = Type.x)


##Exporting edge's file for Gephi

#Saving edges file as csv
write.csv(edges_group_clean, "D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_Global/Pathogen/Campylobacter/CASE_edges_GENUS_ARGgroup_Campylobacter_labels_80cutoff_50percTaxaARG_20220127.csv",
          row.names=FALSE)

write.csv(edges_gephi, "D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/CorrelationNetwork/Jan_2022_Global/Pathogen/Campylobacter/CASE_edges_GENUS_ARGgroup_Campylobacter_GEPHI_nolabels_80cutoff_50percTaxaARG_20220127.csv",
          row.names=FALSE)
