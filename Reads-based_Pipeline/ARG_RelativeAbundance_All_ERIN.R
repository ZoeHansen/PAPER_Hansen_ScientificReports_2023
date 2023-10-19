#################################################

# Aim 2: Relative Abundance of ARGs - Case, Control, Follow

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

class_data <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/ERIN_GEnorm_DRUG_and_MULTI_class_level.csv',
                       header = TRUE)
class_data <- class_data[, colSums(class_data != 0) > 0]


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
class_data1 <- class_data %>%
  mutate(., SampTotal = rowSums(class_data[,-1]))


# Create a new dataframe with the relative abundance information
ra <- cbind(class_data1$ER_ID, class_data1[, -c(1,81)] / class_data1$SampTotal)
colnames(ra)[1] <- 'ER_ID'


# Merge our relative abundance dataframe with metadata and order by Case.status
ra_df <- left_join(meta, ra, by = "ER_ID") %>%
  arrange(., Case.status)


####### Rename Select Columns ######
ra_df <-dplyr::rename(ra_df, 'Multi-compound' = 'Multi.compound')

# Rename a few of the columns for formatting (Class level - Drugs and Multi-compounds)
ra_df <- dplyr::rename(ra_df, 'Beta-lactams' = 'betalactams')
ra_df <- dplyr::rename(ra_df, 'CAP' = 'Cationic_antimicrobial_peptides')
ra_df <- dplyr::rename(ra_df, "Drug, Biocide & Metal resistance" = 'Drug_and_biocide_and_metal_resistance')
ra_df <- dplyr::rename(ra_df, "Drug & Biocide resistance" = 'Drug_and_biocide_resistance')
ra_df <- dplyr::rename(ra_df, "Drug & Metal resistance" = 'Drug_and_metal_resistance')
ra_df <- dplyr::rename(ra_df, "MDR" = 'Multi.drug_resistance')

# If Metals/Biocides are included:
ra_df <- dplyr::rename(ra_df, "Acetate resistance" = 'Acetate_resistance')
ra_df <- dplyr::rename(ra_df, "Acid resistance" = 'Acid_resistance')
ra_df <- dplyr::rename(ra_df, "Arsenic resistance" = 'Arsenic_resistance')
ra_df <- dplyr::rename(ra_df, "Biguanide resistance" = 'Biguanide_resistance')
ra_df <- dplyr::rename(ra_df, "Biocide & Metal resistance" = 'Biocide_and_metal_resistance')
ra_df <- dplyr::rename(ra_df, "Chromium resistance" = 'Chromium_resistance')
ra_df <- dplyr::rename(ra_df, "Copper resistance" = 'Copper_resistance')
ra_df <- dplyr::rename(ra_df, "Gold resistance" = 'Gold_resistance')
ra_df <- dplyr::rename(ra_df, "Iron resistance" = 'Iron_resistance')
ra_df <- dplyr::rename(ra_df, "Mercury resistance" = 'Mercury_resistance')
ra_df <- dplyr::rename(ra_df, "Multi-biocide resistance" = 'Multi.biocide_resistance')
ra_df <- dplyr::rename(ra_df, "Multi-metal resistance" = 'Multi.metal_resistance')
ra_df <- dplyr::rename(ra_df, "Nickel resistance" = 'Nickel_resistance')
ra_df <- dplyr::rename(ra_df, "Paraquat resistance" = 'Paraquat_resistance')
ra_df <- dplyr::rename(ra_df, "Peroxide resistance" = 'Peroxide_resistance')
ra_df <- dplyr::rename(ra_df, "Phenicol compound resistance" = 'Phenolic_compound_resistance')
ra_df <- dplyr::rename(ra_df, "QACs resistance" = 'Quaternary_Ammonium_Compounds_.QACs._resistance')
ra_df <- dplyr::rename(ra_df, 'Sodium resistance' = 'Sodium_resistance')
ra_df <- dplyr::rename(ra_df, 'Tellurium resistance' = 'Tellurium_resistance')
ra_df <- dplyr::rename(ra_df, 'Zinc resistance' = 'Zinc_resistance')

ra_df <- ra_df[, colSums(ra_df != 0)>0]

write.csv(ra_df, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/RelativeAbundance_DRUGS_AND_MULTI_MechanismLevel_CaseFollowPaired.csv',
          row.names = FALSE)

#################################################
# Prepare data for plotting
#################################################

ra_df <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/DRUGS_and_MULTI/RelativeAbundance_DRUG_and_MULTI_ClassLevel_CaseFollowPaired.csv',
                  header=TRUE, check.names = FALSE)

# Extract Case status data to append later (the melt function requires only one variable to collapse upon)
ra_df_ordered <- ra_df %>%
#  arrange(.,Case.status, Case.Follow_ID)%>%
  filter(.,Case.status == 'FollowUp')%>%
  arrange(.,Bacteria)

Health <- ra_df_ordered$Case.status
Bact <- ra_df_ordered$Bacteria

# Remove the case status data from the dataframe to prepare for melting 
ra_df.cc <- ra_df_ordered%>%
  dplyr::select(., -c(Case.status, Bacteria, Case.Follow_ID))

# Create an arbitrary sequential numbering system to avoid spacing in the x-axis of the plot
id.num <- seq(1,120,1)
id.num.case <- seq(1,60,1)

# Melt our ra_df.cc variable, and attach the Health and Num variables to the dataframe
ra_df.long <- melt(ra_df.cc, id.vars = 'ER_ID', variable.name = 'Resistance_Class')
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
ggplot(data = ra_df.long, aes(x = Num, y = value, fill = Resistance_Class))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
#  scale_fill_manual(values = getPalette(17), 
#                    guide = guide_legend(nrow=4))+
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
  labs(fill = 'Resistance Class')+
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


####################################################
# Plotting Top 10 ARG Classes 
####################################################
case.class <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/ALL_TYPES/AverageRelativeAbundance_ALL_TYPES_ClassLevel_Top10_Cases.csv',
                       header=TRUE)

case.class <- case.class %>%
  filter(ER_ID != 'Case.status')

case.c.t <- case.class %>%
  gather(key = key, value = value, 2:ncol(case.class)) %>%
  spread(key=names(case.class)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

case.c.wID <- left_join(meta, case.c.t, by='ER_ID')

case.c.wID <- case.c.wID %>%
  filter(Case.status != 'FollowUp')%>%
  arrange(., Case.Follow_ID)

# Use Case.Follow_ID as the number to order the columns in the plot
case.c.long <- melt(case.c.wID, id.vars = c('ER_ID', 'Case.status','Pathogen','Case.Follow_ID'), 
                    variable.name = 'CLASS')

case.c.long$value <- as.numeric(case.c.long$value)


# Follow-Ups

follow.class <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/ALL_TYPES/AverageRelativeAbundance_ALL_TYPES_ClassLevel_Top10_FollowUps.csv',
                         header=TRUE)

follow.class <- follow.class %>%
  filter(ER_ID != 'Case.status')

follow.c.t <- follow.class %>%
  gather(key = key, value = value, 2:ncol(follow.class)) %>%
  spread(key=names(follow.class)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

follow.c.wID <- left_join(meta, follow.c.t, by='ER_ID')

follow.c.wID <- follow.c.wID %>%
  filter(Case.status != 'Case')%>%
  arrange(., Case.Follow_ID)

# Use Case.Follow_ID to order the columns in the plot
follow.c.long <- melt(follow.c.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'),
                      variable.name = 'CLASS')


follow.c.long$value <- as.numeric(follow.c.long$value)


combined.c.long <- rbind(case.c.long, follow.c.long)


ggplot(data = combined.c.long, aes(x = Case.Follow_ID, y = value, fill = CLASS))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE, option = 'C')+
  scale_fill_viridis(discrete = TRUE, option = 'C', 
                     guide=guide_legend(nrow=5))+
  scale_x_discrete('Num', name = 'Health Status')+
  scale_y_continuous(expand = c(0.01,0))+
  facet_wrap( ~ Case.status, strip.position = 'bottom', scales = 'free_x')+ 
#  facet_wrap( ~Case.status, strip.position='left', ncol=1, scales = 'free_y')+
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
#        panel.spacing = unit(0.5,'cm'),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 12), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'ARG Class')+
  xlab('\nHealth Status\n')+
  ylab('Relative Abundance per Sample\n')


####################################################
# Plotting Top 25 ARG Groups 
####################################################

### Cases

# Cases

case.group <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/ALL_TYPES/AverageRelativeAbundance_ALL_TYPES_GroupLevel_Top25_Cases.csv',
                       header=TRUE)

case.group <- case.group %>%
  filter(ER_ID != 'Case.status')

case.g.t <- case.group %>%
  gather(key = key, value = value, 2:ncol(case.group)) %>%
  spread(key=names(case.group)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

case.g.wID <- left_join(meta, case.g.t, by='ER_ID')

case.g.wID <- case.g.wID %>%
  filter(Case.status != 'FollowUp')%>%
  arrange(., Case.Follow_ID)

# Use Case.Follow_ID as the number to order the columns in the plot
case.g.long <- melt(case.g.wID, id.vars = c('ER_ID', 'Case.status','Pathogen','Case.Follow_ID'), 
                    variable.name = 'GROUP')

case.g.long$value <- as.numeric(case.g.long$value)


# Follow-Ups

follow.group <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/ALL_TYPES/AverageRelativeAbundance_ALL_TYPES_GroupLevel_Top25_FollowUps.csv',
                         header=TRUE)

follow.group <- follow.group %>%
  filter(ER_ID != 'Case.status')

follow.g.t <- follow.group %>%
  gather(key = key, value = value, 2:ncol(follow.group)) %>%
  spread(key=names(follow.group)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

follow.g.t <-dplyr::rename(follow.g.t, 'APH3-DPRIME' = 'APH3.DPRIME')

follow.g.wID <- left_join(meta, follow.g.t, by='ER_ID')

follow.g.wID <- follow.g.wID %>%
  filter(Case.status != 'Case')%>%
  arrange(., Case.Follow_ID)

# Use Case.Follow_ID to order the columns in the plot
follow.g.long <- melt(follow.g.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'),
                      variable.name = 'GROUP')


follow.g.long$value <- as.numeric(follow.g.long$value)


combined.long <- rbind(case.g.long, follow.g.long)


ggplot(data = combined.long, aes(x = Case.Follow_ID, y = value, fill = GROUP))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'D', 
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
  labs(fill = 'ARG Group')+
  xlab('\nHealth Status\n')+
  ylab('Relative Abundance per Sample\n')


###### Plotting Case/Follow Top 25 GroupsSeparately ######

#To plot these separately (with different color schemes):
ggplot(data = case.g.long, aes(x = Case.Follow_ID, y = value, fill = GROUP))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_fill_manual(values = colors, 
                    guide = guide_legend(nrow=5))+
  #  scale_color_viridis(discrete = TRUE)+
  #  scale_fill_viridis(discrete = TRUE, option = 'D', 
  #                     guide=guide_legend(nrow=5))+
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
  labs(fill = 'ARG Group')+
  ylab('Relative Abundance per Sample\n')



ggplot(data = follow.g.long, aes(x = Case.Follow_ID, y = value, fill = GROUP))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_fill_manual(values = colors$Color, 
                    guide = guide_legend(nrow=5))+
#  scale_color_viridis(discrete = TRUE)+
#  scale_fill_viridis(discrete = TRUE, option = 'C', 
#                     guide=guide_legend(nrow=5))+
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
  labs(fill = 'ARG Group')+
  ylab('Relative Abundance per Sample\n')


#Combine Plots #

ggarrange(case.g.plot + rremove("ylab"),
          follow.g.plot + rremove("ylab"),
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 2, nrow = 1)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))

####################################################
# Plotting Top 10 ARG Groups 
####################################################

### Cases

# Cases

case.group <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/ALL_TYPES/AverageRelativeAbundance_ALL_TYPES_GroupLevel_Top10_Cases.csv',
                       header=TRUE)

case.group <- case.group %>%
  filter(ER_ID != 'Case.status')

case.g.t <- case.group %>%
  gather(key = key, value = value, 2:ncol(case.group)) %>%
  spread(key=names(case.group)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

case.g.wID <- left_join(meta, case.g.t, by='ER_ID')

case.g.wID <- case.g.wID %>%
  filter(Case.status != 'FollowUp')%>%
  arrange(., Case.Follow_ID)


case.g.long <- melt(case.g.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'),
                    variable.name = 'GROUP')

case.g.long$value <- as.numeric(case.g.long$value)


# Follow-Ups

follow.group <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/ALL_TYPES/AverageRelativeAbundance_ALL_TYPES_GroupLevel_Top10_FollowUps.csv',
                         header=TRUE)

follow.group <- follow.group %>%
  filter(ER_ID != 'Case.status')

follow.g.t <- follow.group %>%
  gather(key = key, value = value, 2:ncol(follow.group)) %>%
  spread(key=names(follow.group)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

follow.g.wID <- left_join(meta, follow.g.t, by='ER_ID')

follow.g.wID <- follow.g.wID %>%
  filter(Case.status != 'Case')%>%
  arrange(., Case.Follow_ID)


follow.g.long <- melt(follow.g.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'),
                      variable.name = 'GROUP')


follow.g.long$value <- as.numeric(follow.g.long$value)


combined.long <- rbind(case.g.long, follow.g.long)


ggplot(data = combined.long, aes(x = Case.Follow_ID, y = value, fill = GROUP))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'D', 
                     guide=guide_legend(nrow=2))+
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
        panel.spacing = unit(0.5,'cm'),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  labs(fill = 'ARG Group')+
  xlab('\nHealth Status\n')+
  ylab('Relative Abundance per Sample\n')



###### Plotting Case/Follow Top 10 Groups Separately ########

case.g.plot<-ggplot(data = case.g.long, aes(x = Num, y = value, fill = GROUP))+
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
  labs(fill = 'ARG Group')+
  ylab('Relative Abundance per Sample\n')


follow.g.plot<-ggplot(data = follow.g.long, aes(x = Num, y = value, fill = GROUP))+
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
  labs(fill = 'ARG Group')+
  ylab('Relative Abundance per Sample\n')


#Combine Plots #

ggarrange(case.g.plot + rremove("ylab"),
          follow.g.plot + rremove("ylab"),
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 2, nrow = 1)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))




####################################################
# Plotting Top 10 ARG Mechanisms
####################################################

### Cases

# Cases

case.mech <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/ALL_TYPES/AverageRelativeAbundance_ALL_TYPES_MechanismLevel_Top10_Cases.csv',
                       header=TRUE)

case.mech <- case.mech %>%
  select(-AvgRA)%>%
  filter(ER_ID != 'Case.status')

case.m.t <- case.mech %>%
  gather(key = key, value = value, 2:ncol(case.mech)) %>%
  spread(key=names(case.mech)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

case.m.t <- dplyr::rename(case.m.t, 'Aminoglycoside-resistant 16S ribosomal subunit protein' = 'Aminoglycoside.resistant_16S_ribosomal_subunit_protein')
case.m.t <- dplyr::rename(case.m.t, 'Copper resistance protein' = 'Copper_resistance_protein')
case.m.t <- dplyr::rename(case.m.t, 'Drug, Biocide & Metal RND efflux pumps' = 'Drug_and_biocide_and_metal_RND_efflux_pumps')
case.m.t <- dplyr::rename(case.m.t, 'Drug & Biocide MFS efflux pumps' = 'Drug_and_biocide_MFS_efflux_pumps')
case.m.t <- dplyr::rename(case.m.t, 'Drug & Biocide RND efflux pumps' = 'Drug_and_biocide_RND_efflux_pumps')
case.m.t <- dplyr::rename(case.m.t, 'Fluoroquinolone-resistant DNA topoisomerases' = 'Fluoroquinolone.resistant_DNA_topoisomerases')
case.m.t <- dplyr::rename(case.m.t, 'Fosfomycin target mutation' = 'Fosfomycin_target_mutation')
case.m.t <- dplyr::rename(case.m.t, 'Macrolide-resistant 23S rRNA mutation' = 'Macrolide.resistant_23S_rRNA_mutation')
case.m.t <- dplyr::rename(case.m.t, 'Multi-metal resistance protein' = 'Multi.metal_resistance_protein')
case.m.t <- dplyr::rename(case.m.t, 'Rifampin-resistant beta-subunit (RNA polymerase RpoB)' = 'Rifampin.resistant_beta.subunit_of_RNA_polymerase_RpoB')


case.m.wID <- left_join(meta, case.m.t, by='ER_ID')

case.m.wID <- case.m.wID %>%
  filter(Case.status != 'FollowUp')%>%
  arrange(., Case.Follow_ID)


case.m.long <- melt(case.m.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'), 
                    variable.name = 'MECHANISM')


case.m.long$value <- as.numeric(case.m.long$value)


# Follow-Ups

follow.mech <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Resistome/Reads_based/Abundance/RelativeAbundance/ALL_TYPES/AverageRelativeAbundance_ALL_TYPES_MechanismLevel_Top10_FollowUps.csv',
                         header=TRUE)

follow.mech<- follow.mech %>%
  select(-AvgRA)%>%
  filter(ER_ID != 'Case.status')

follow.m.t <- follow.mech %>%
  gather(key = key, value = value, 2:ncol(follow.mech)) %>%
  spread(key=names(follow.mech)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

follow.m.t <- dplyr::rename(follow.m.t, '23S rRNA methyltransferases' = 'X23S_rRNA_methyltransferases')
follow.m.t <- dplyr::rename(follow.m.t, 'Aminoglycoside-resistant 16S ribosomal subunit protein' = 'Aminoglycoside.resistant_16S_ribosomal_subunit_protein')
follow.m.t <- dplyr::rename(follow.m.t, 'Aminoglycoside O-nucleotidyltransferases' = 'Aminoglycoside_O.nucleotidyltransferases')
follow.m.t <- dplyr::rename(follow.m.t, 'Class A Beta-lactamases' = 'Class_A_betalactamases')
follow.m.t <- dplyr::rename(follow.m.t, 'Drug & Biocide RND efflux pumps' = 'Drug_and_biocide_RND_efflux_pumps')
follow.m.t <- dplyr::rename(follow.m.t, 'Macrolide-resistant 23S rRNA mutation' = 'Macrolide.resistant_23S_rRNA_mutation')
follow.m.t <- dplyr::rename(follow.m.t, 'MLS efflux pumps' = 'MLS_resistance_MFS_efflux_pumps')
follow.m.t <- dplyr::rename(follow.m.t, 'Multi-metal resistance protein' = 'Multi.metal_resistance_protein')
follow.m.t <- dplyr::rename(follow.m.t, 'Tetracycline inactivation enzymes' = 'Tetracycline_inactivation_enzymes')
follow.m.t <- dplyr::rename(follow.m.t, 'Tetracycline-resistant Ribosomal protection proteins' = 'Tetracycline_resistance_ribosomal_protection_proteins')


follow.m.wID <- left_join(meta, follow.m.t, by='ER_ID')

follow.m.wID <- follow.m.wID %>%
  filter(Case.status != 'Case')%>%
  arrange(., Case.Follow_ID)


follow.m.long <- melt(follow.m.wID, id.vars = c('ER_ID','Case.status','Pathogen','Case.Follow_ID'), 
                      variable.name = 'MECHANISM')

follow.m.long$value <- as.numeric(follow.m.long$value)

combined.m.long <- rbind(case.m.long, follow.m.long)


ggplot(data = combined.m.long, aes(x = Case.Follow_ID, y = value, fill = MECHANISM))+
  geom_bar(stat = 'identity', width = 1, position = 'stack')+ 
  scale_color_viridis(discrete = TRUE)+
  scale_fill_viridis(discrete = TRUE, option = 'D', 
                     guide=guide_legend(nrow=8))+
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
  labs(fill = 'ARG Mechanism')+
  xlab('\nHealth Status\n')+
  ylab('Relative Abundance per Sample\n')



###### Plotting Case/Follow Top 10 Mechanism Separately #######

case.m.plot<-ggplot(data = case.m.long, aes(x = Num, y = value, fill = MECHANISM))+
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
  labs(fill = 'ARG \nMechanism')+
  ylab('Relative Abundance per Sample\n')

follow.m.plot<-ggplot(data = follow.m.long, aes(x = Num, y = value, fill = MECHANISM))+
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
  labs(fill = 'ARG \nMechanism')+
  ylab('Relative Abundance per Sample\n')


#Combine Plots #

ggarrange(case.m.plot,
          follow.m.plot + rremove("ylab"),
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 2, nrow = 1)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))
