#################################################

# Aim 2: Relative Abundance - Microbiome - Case, Control, Follow

# Author: Zoe Hansen
# Last Modified: 2022.01.14

#################################################
# Load libraries and data

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(scales)
library(viridis)

microbiome <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/GEnorm_ERIN_kaiju_CaseControlFollow_SPECIES.csv',
                       header = TRUE)

# Since some fields have NA, we need to redefine and combine them
microbiome$Class <- microbiome$Species %>%
  replace_na("Unknown")

micro <- ddply(microbiome, "Species", numcolwise(sum))
micro<- micro[, colSums(micro[,-1] != 0) > 0] 
micro <- micro[rowSums(micro[,-1] != 0) > 0,]

# We then need to transpose to obtain samples in the rows and taxa in the columns
micro_data <- micro %>%
  gather(key = key, value = value, 2:ncol(micro)) %>%
  spread(key=names(micro)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

micro_data <- micro_data[, colSums(micro_data != 0) > 0]
micro_data <- micro_data[rowSums(micro_data !=0)>0,]


# Import metadata

meta <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/ERIN_Metagenomes_Metadata_60_CaseFollowPairs.csv',
                 header = TRUE)

meta <- meta %>%
  dplyr::select(ER_ID, Case.status, Pathogen, Case.Follow_ID)%>%
  filter(!grepl('Control',Case.status)) %>%
  #  filter(!grepl('FollowUp', Case.status))%>%
  drop_na()


#################################################
# Calculate relative abundance of ARGs in samples
#################################################

# This was used to transpose the hierarchical labels of Type, Class, Group
#classt = setNames(data.frame(t(class_data[,-1])), class_data[,1])
#classt <- classt %>%
#  rownames_to_column('ER_ID')

# Add a column with the total reads/genome equivalents for each sample
micro_data1 <- micro_data %>%
  mutate(., SampTotal = rowSums(micro_data[,-1]))


# Create a new dataframe with the relative abundance information
ra <- cbind(micro_data1$ER_ID, micro_data1[, -c(1,40024)] / micro_data1$SampTotal)
colnames(ra)[1] <- 'ER_ID'


# Merge our relative abundance dataframe with metadata and order by Case.status
ra_df <- left_join(meta, ra, by = "ER_ID") %>%
  arrange(., Case.status)


####### Rename Select Columns ######
ra_df <-dplyr::rename(ra_df, 'Multi-compound' = 'Multi.compound')

ra_df <- ra_df[, colSums(ra_df != 0)>0]

ra_df_t <- ra_df %>%
  gather(key = key, value = value, 2:ncol(ra_df)) %>%
  spread(key=names(ra_df)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

write.csv(ra_df_t, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/Abundance/RelativeAbundance/RelativeAbundance_GENUS_CaseFollowPairs.csv',
          row.names = FALSE)

#################################################
# Prepare data for plotting
#################################################

ra_df <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/Abundance/RelativeAbundance/RelativeAbundance_GENUS_CaseFollowPairs.csv',
                  header=TRUE, check.names = FALSE)

# Extract Case status data to append later (the melt function requires only one variable to collapse upon)
ra_df_ordered <- ra_df %>%
  arrange(.,Case.status, Case.Follow_ID)#%>%
#  filter(.,Case.status == 'FollowUp')%>%
#  arrange(.,Pathogen)

Health <- ra_df_ordered$Case.status
Bact <- ra_df_ordered$Pathogen

# Remove the case status data from the dataframe to prepare for melting 
ra_df.cc <- ra_df_ordered%>%
  dplyr::select(., -c(Case.status, Pathogen, Case.Follow_ID))

# Create an arbitrary sequential numbering system to avoid spacing in the x-axis of the plot
id.num <- seq(1,120,1)
id.num.case <- seq(1,60,1)

# Melt our ra_df.cc variable, and attach the Health and Num variables to the dataframe
ra_df.long <- melt(ra_df.cc, id.vars = 'ER_ID', variable.name = 'Kingdom')
ra_df.long$Case.status = rep(Health, times = (ncol(ra_df.cc)-1))
ra_df.long$Num = rep(id.num, times = (ncol(ra_df.cc)-1))

ra_df.long <- melt(ra_df.cc, id.vars = 'ER_ID', variable.name = 'Resistance_Class')
ra_df.long$Bacteria = rep(Bact, times = (ncol(ra_df.cc)-1))
ra_df.long$Num = rep(id.num.case, times = (ncol(ra_df.cc)-1))


#################################################
# Plot relative abundance 
#################################################

# Designate a color palette (used later in our ggplot function)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

### Create the plot

# Facet by Health Status
ggplot(data = ra_df.long, aes(x = Num, y = value, fill = Kingdom))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'D', 
                     guide=guide_legend(nrow=1))+
  scale_x_discrete('Num', name = 'Health Status')+
  scale_y_continuous(expand = c(0.01,0))+
  facet_wrap( ~ Case.status, strip.position = 'bottom', scales = 'free_x')+ 
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Kingdom')+
  xlab('\nHealth Status\n')+
  ylab('Relative Abundance per Sample\n')


# Facet by Pathogen
ggplot(data = ra_df.long, aes(x = Num, y = value, fill = Resistance_Class))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'D', 
                     guide=guide_legend(nrow=6))+
  scale_x_discrete('Num', name = 'Pathogen')+
  scale_y_continuous(expand = c(0.01,0))+
  facet_wrap( ~ Bacteria, strip.position = 'bottom', scales = 'free_x')+ 
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Resistance Class')+
  xlab('\nPathogen\n')+
  ylab('Relative Abundance per Sample\n')

#############################################
# Plot Top 10 Phyla
#############################################

### PHYLUM ###

# Case

case.phylum <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/Abundance/RelativeAbundance/AverageRelativeAbundance_PHYLUM_Top10_CasesOnly.csv',
                        header=TRUE)
case.phylum <- case.phylum %>%
  select(-AvgRA) %>%
  filter(ER_ID != 'Case.status')


case.p.t <- case.phylum %>%
  gather(key = key, value = value, 2:ncol(case.phylum)) %>%
  spread(key=names(case.phylum)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

case.p.t <-dplyr::rename(case.p.t, 'Cannot be assigned' = 'cannot be assigned a (non-viral) phylum')
case.p.t <-dplyr::rename(case.p.t, 'Unclassified' = 'unclassified')

case.p.wID <- left_join(meta, case.p.t, by='ER_ID')

case.p.wID <- case.p.wID %>%
  filter(Case.status != 'FollowUp')%>%
  arrange(., Case.Follow_ID)

case.p.long <- melt(case.p.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'),
                    variable.name = 'Phylum')

case.p.long$value <- as.numeric(case.p.long$value)


# Follow-Ups

follow.phylum <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/Abundance/RelativeAbundance/AverageRelativeAbundance_PHYLUM_Top10_FollowUpsOnly.csv',
                           header=TRUE)
follow.phylum <- follow.phylum %>%
  select(-AvgRA) %>%
  filter(ER_ID != 'Case.status')

follow.p.t <- follow.phylum %>%
  gather(key = key, value = value, 2:ncol(follow.phylum)) %>%
  spread(key=names(follow.phylum)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

follow.p.t <-dplyr::rename(follow.p.t, 'Cannot be assigned' = 'cannot be assigned a (non-viral) phylum')
follow.p.t <-dplyr::rename(follow.p.t, 'Unclassified' = 'unclassified')

follow.p.wID <- left_join(meta, follow.p.t, by='ER_ID')

follow.p.wID <- follow.p.wID %>%
  filter(Case.status != 'Case')%>%
  arrange(., Case.Follow_ID)

follow.p.long <- melt(follow.p.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'), 
                      variable.name = 'Phylum')


follow.p.long$value <- as.numeric(follow.p.long$value)

combined.p.long <- rbind(case.p.long, follow.p.long)


ggplot(data = combined.p.long, aes(x = Case.Follow_ID, y = value, fill = Phylum))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE, option = 'H')+
  scale_fill_viridis(discrete = TRUE, option = 'H', 
                     guide=guide_legend(nrow=4))+
  scale_x_discrete('Num', name = NULL)+
  scale_y_continuous(expand = c(0.01,0))+
#  facet_wrap( ~ Case.status, strip.position = 'bottom', scales = 'free_x')+ 
  facet_wrap( ~Case.status, strip.position='left', ncol=1, scales = 'free_y')+
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.spacing = unit(0.5, 'cm'),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Phylum')+
#  xlab('\nHealth Status\n')+
  ylab('Relative Abundance per Sample\n')


###### Plotting Case/Follow Top 10 Phylum Separately #######

case.p.plot <- ggplot(data = case.p.long, aes(x = Num, y = value, fill = Phylum))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'D', 
                     guide=guide_legend(nrow=5))+
  scale_x_discrete('Num', name = 'Case')+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Phylum')+
  ylab('Relative Abundance per Sample\n')


follow.p.plot <- ggplot(data = follow.p.long, aes(x = Num, y = value, fill = Phylum))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'C', 
                     guide=guide_legend(nrow=5))+
  scale_x_discrete('Num', name = 'Follow-Ups')+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Phylum')+
  ylab('Relative Abundance per Sample\n')


#Combine Plots #

ggarrange(case.p.plot + rremove("ylab"),
          follow.p.plot + rremove("ylab"),
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 2, nrow = 1)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))


#############################################
# Plot Top 25 Genera
#############################################

# Cases

case.genus <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/Abundance/RelativeAbundance/AverageRelativeAbundance_GENUS_Top25_CasesOnly.csv',
                        header=TRUE)

case.genus <- case.genus %>%
  select(-AvgRA) %>%
  filter(ER_ID != 'Case.status')

case.g.t <- case.genus %>%
  gather(key = key, value = value, 2:ncol(case.genus)) %>%
  spread(key=names(case.genus)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

case.g.t <-dplyr::rename(case.g.t, 'Cannot be assigned' = 'cannot be assigned a (non-viral) genus')
case.g.t <-dplyr::rename(case.g.t, 'Unclassified' = 'unclassified')

case.g.wID <- left_join(meta, case.g.t, by='ER_ID')

case.g.wID <- case.g.wID %>%
  filter(Case.status != 'FollowUp')%>%
  arrange(., Case.Follow_ID)

case.g.long <- melt(case.g.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'), 
                      variable.name = 'Genus')

case.g.long$value <- as.numeric(case.g.long$value)


# Follow-Ups

follow.genus <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/Abundance/RelativeAbundance/AverageRelativeAbundance_GENUS_Top25_FollowUpsOnly.csv',
                       header=TRUE)

follow.genus <- follow.genus %>%
  select(-AvgRA) %>%
  filter(ER_ID != 'Case.status')

follow.g.t <- follow.genus %>%
  gather(key = key, value = value, 2:ncol(follow.genus)) %>%
  spread(key=names(follow.genus)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

follow.g.t <-dplyr::rename(follow.g.t, 'Cannot be assigned' = 'cannot be assigned a (non-viral) genus')
follow.g.t <-dplyr::rename(follow.g.t, 'Unclassified' = 'unclassified')

follow.g.wID <- left_join(meta, follow.g.t, by='ER_ID')

follow.g.wID <- follow.g.wID %>%
  filter(Case.status != 'Case')%>%
  arrange(., Case.Follow_ID)

follow.g.long <- melt(follow.g.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'), 
                      variable.name = 'Genus')


follow.g.long$value <- as.numeric(follow.g.long$value)

combined.g.long <- rbind(case.g.long, follow.g.long)


ggplot(data = combined.g.long, aes(x = Case.Follow_ID, y = value, fill = Genus))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE, option = 'H')+
  scale_fill_viridis(discrete = TRUE, option = 'H', 
                     guide=guide_legend(nrow=7))+
  scale_x_discrete('Num', name = 'Health Status')+
  scale_y_continuous(expand = c(0.01,0))+
  facet_wrap( ~ Case.status, strip.position = 'bottom', scales = 'free_x')+ 
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Genus')+
  xlab('\nHealth Status\n')+
  ylab('Relative Abundance per Sample\n')



###### Plotting Case/Follow Top 25 Genera Separately ######

case.g.plot<- ggplot(data = case.g.long, aes(x = Num, y = value, fill = Genus))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'D', 
                     guide=guide_legend(nrow=9))+
  scale_x_discrete('Num', name = 'Case')+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Genus')+
  ylab('Relative Abundance per Sample\n')

follow.g.plot <- ggplot(data = follow.g.long, aes(x = Num, y = value, fill = Genus))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'C', 
                     guide=guide_legend(nrow=9))+
  scale_x_discrete('Num', name = 'Follow-Ups')+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Genus')+
  ylab('Relative Abundance per Sample\n')


#Combine Plots #

ggarrange(case.g.plot + rremove("ylab"),
          follow.g.plot + rremove("ylab"),
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 2, nrow = 1)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))


#############################################
# Plot Top 10 Genera
#############################################

# Cases

case.genus10 <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/Abundance/RelativeAbundance/AverageRelativeAbundance_GENUS_Top10_CasesOnly.csv',
                       header=TRUE)

case.genus10 <- case.genus10 %>%
  select(-AvgRA) %>%
  filter(ER_ID != 'Case.status')

case.g.t10 <- case.genus10 %>%
  gather(key = key, value = value, 2:ncol(case.genus10)) %>%
  spread(key=names(case.genus10)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

case.g.t10 <-dplyr::rename(case.g.t10, 'Cannot be assigned' = 'cannot be assigned a (non-viral) genus')

case.g.wID10 <- left_join(meta, case.g.t10, by='ER_ID')

case.g.wID10 <- case.g.wID10 %>%
  filter(Case.status != 'FollowUp')%>%
  arrange(., Case.Follow_ID)
#  arrange(., Pathogen, Case.Follow_ID)

case.g.long10 <- melt(case.g.wID10, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'), 
                      variable.name = 'Genus')

case.g.long10$value <- as.numeric(case.g.long10$value)

# Add "Num" if reordering for Pathogens
#id.num.case <- seq(1,60,1)
#case.g.long10$Num = rep(id.num.case, times = (ncol(case.g.long10)-1))


# Follow-Ups

follow.genus10 <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Microbiome/Reads_based/Abundance/RelativeAbundance/AverageRelativeAbundance_GENUS_Top10_FollowUpsOnly.csv',
                         header=TRUE)

follow.genus10 <- follow.genus10 %>%
  select(-AvgRA) %>%
  filter(ER_ID != 'Case.status')

follow.g.t10 <- follow.genus10 %>%
  gather(key = key, value = value, 2:ncol(follow.genus10)) %>%
  spread(key=names(follow.genus10)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

follow.g.t10 <-dplyr::rename(follow.g.t10, 'Cannot be assigned' = 'cannot be assigned a (non-viral) genus')
follow.g.t10 <-dplyr::rename(follow.g.t10, 'Unclassified' = 'unclassified')

follow.g.wID10 <- left_join(meta, follow.g.t10, by='ER_ID')

follow.g.wID10 <- follow.g.wID10 %>%
  filter(Case.status != 'Case')%>%
  arrange(.,Case.Follow_ID)
#arrange(., Pathogen, Case.Follow_ID)

follow.g.long10 <- melt(follow.g.wID10, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'), 
                      variable.name = 'Genus')


follow.g.long10$value <- as.numeric(follow.g.long10$value)

# Add "Num" if reordering for Pathogens
#id.num.follow <- seq(1,60,1)
#follow.g.long10$Num = rep(id.num.follow, times = (ncol(follow.g.long10)-1))

combined.g.long10 <- rbind(case.g.long10, follow.g.long10)


ggplot(data = combined.g.long10, aes(x = Case.Follow_ID, y = value, fill = Genus))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE, option = 'H')+
  scale_fill_viridis(discrete = TRUE, option = 'H', 
                     guide=guide_legend(nrow=4))+
  scale_x_discrete('Case.Follow_ID', name = NULL)+
  scale_y_continuous(expand = c(0.01,0))+
#  facet_wrap( ~ Case.status, strip.position = 'bottom', scales = 'free_x')+ 
  facet_wrap( ~Case.status, strip.position='left', ncol=1, scales = 'free_y')+
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        strip.text = element_text(size=12),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.spacing = unit(0.5, "cm"),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 11), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Genus')+
#  xlab('\nHealth Status\n')+
  ylab('Relative Abundance per Sample\n')


###### Plotting Case/Follow Top 10 Genera Separately ###########

case.g.plot<- ggplot(data = case.g.long, aes(x = Num, y = value, fill = Genus))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'D', 
                     guide=guide_legend(nrow=5))+
  scale_x_discrete('Num', name = 'Cases')+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Genus')+
  ylab('Relative Abundance per Sample\n')




follow.g.plot<-ggplot(data = follow.g.long, aes(x = Num, y = value, fill = Genus))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'C', 
                     guide=guide_legend(nrow=5))+
  scale_x_discrete('Num', name = 'Follow-Ups')+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'bottom', 
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size = 11),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        #        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'Genus')+
  ylab('Relative Abundance per Sample\n')


#Combine Plots #

ggarrange(case.g.plot + rremove("ylab"),
          follow.g.plot + rremove("ylab"),
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 2, nrow = 1)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))




