##############################################

# Aim 2 - Actual Abundance of ARGs

##############################################
library(tidyverse)
library(ggplot2)
library(ggpubr)
#################################################
# Top 10 ARG Classes
#################################################

# Cases

case.class <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/AverageAbundance_ALL_TYPES_ClassLevel_Top10_Cases.csv',
                  header = TRUE)

case.c.t <- case.class %>%
  gather(key = key, value = value, 2:ncol(case.class)) %>%
  spread(key=names(case.class)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

case.c.t <-dplyr::rename(case.c.t, 'Drug & Biocide Resistance' = 'Drug_and_biocide_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Multi-metal Resistance' = 'Multi-metal_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Drug, Biocide & Metal Resistance' = 'Drug_and_biocide_and_metal_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Copper resistance' = 'Copper_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Multi-biocide Resistance' = 'Multi-biocide_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Beta-lactams' = 'betalactams')

case.c <- case.c.t %>%
  gather(key = key, value = value, 2:ncol(case.c.t)) %>%
  spread(key=names(case.c.t)[1], value = 'value') %>%
  dplyr::rename(., Class=key)

case.c$CASE <- as.numeric(case.c$CASE)

case.c.plot <- ggplot(data = case.c, aes(x =reorder(Class, -CASE), y = CASE, fill = Class))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18, face='bold', hjust=0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Cases', fill = 'ARG Group')+
  xlab('ARG Class')+
  ylab('Average Abundance per Sample\n')


# Follow-Ups

follow.class <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/AverageAbundance_ALL_TYPES_ClassLevel_Top10_FollowUps.csv',
                       header = TRUE)

follow.c.t <- follow.class %>%
  gather(key = key, value = value, 2:ncol(follow.class)) %>%
  spread(key=names(follow.class)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

follow.c.t <-dplyr::rename(follow.c.t, 'Drug & Biocide Resistance' = 'Drug_and_biocide_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Multi-metal Resistance' = 'Multi-metal_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Drug, Biocide & Metal Resistance' = 'Drug_and_biocide_and_metal_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Copper resistance' = 'Copper_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'CAP' = 'Cationic_antimicrobial_peptides')
follow.c.t <-dplyr::rename(follow.c.t, 'Beta-lactams' = 'betalactams')

follow.c <- follow.c.t %>%
  gather(key = key, value = value, 2:ncol(follow.c.t)) %>%
  spread(key=names(follow.c.t)[1], value = 'value') %>%
  dplyr::rename(., Class=key)

follow.c$FOLLOWUP <- as.numeric(follow.c$FOLLOWUP)

follow.c.plot <- ggplot(data = follow.c, aes(x =reorder(Class, -FOLLOWUP), y = FOLLOWUP, fill = Class))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18, face='bold', hjust=0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title='Follow-Ups',fill = 'ARG Group')+
  xlab('ARG Class')+
  ylab('Average Abundance per Sample\n')


# Combine plots #
ggarrange(case.c.plot,
          follow.c.plot + rremove("ylab"),
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 2, nrow = 1)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))


#################################################
# Top 25 ARG Classes
#################################################

# Cases

case.class <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/AverageAbundance_ALL_TYPES_ClassLevel_Top25_Cases.csv',
                       header = TRUE)

case.c.t <- case.class %>%
  gather(key = key, value = value, 2:ncol(case.class)) %>%
  spread(key=names(case.class)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

case.c.t <-dplyr::rename(case.c.t, 'Drug & Biocide Resistance' = 'Drug_and_biocide_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Multi-metal Resistance' = 'Multi-metal_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Drug, Biocide & Metal Resistance' = 'Drug_and_biocide_and_metal_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Copper resistance' = 'Copper_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Multi-biocide Resistance' = 'Multi-biocide_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Beta-lactams' = 'betalactams')
case.c.t <-dplyr::rename(case.c.t, 'Biocide & Metal Resistance' = 'Biocide_and_metal_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Acid Resistance' = 'Acid_resistance')
case.c.t <-dplyr::rename(case.c.t, 'CAP' = 'Cationic_antimicrobial_peptides')
case.c.t <-dplyr::rename(case.c.t, 'Arsenic Resistance' = 'Arsenic_resistance')
case.c.t <-dplyr::rename(case.c.t, 'MDR' = 'Multi-drug_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Nickel Resistance' = 'Nickel_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Zinc Resistance' = 'Zinc_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Acetate Resistance' = 'Acetate_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Sodium Resistance' = 'Sodium_resistance')
case.c.t <-dplyr::rename(case.c.t, 'Peroxide Resistance' = 'Peroxide_resistance')

case.c <- case.c.t %>%
  gather(key = key, value = value, 2:ncol(case.c.t)) %>%
  spread(key=names(case.c.t)[1], value = 'value') %>%
  dplyr::rename(., Class=key)

case.c$CASE <- as.numeric(case.c$CASE)

case.c.plot <- ggplot(data = case.c, aes(x =reorder(Class, -CASE), y = CASE, fill = Class))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        plot.title = element_text(size=18, face='bold', hjust=0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Cases', fill = 'ARG Group')+
  xlab('ARG Class')+
  ylab('Average Abundance per Sample\n')


# Follow-Ups

follow.class <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/AverageAbundance_ALL_TYPES_ClassLevel_Top25_FollowUps.csv',
                         header = TRUE)

follow.c.t <- follow.class %>%
  gather(key = key, value = value, 2:ncol(follow.class)) %>%
  spread(key=names(follow.class)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

follow.c.t <-dplyr::rename(follow.c.t, 'Drug & Biocide Resistance' = 'Drug_and_biocide_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Multi-metal Resistance' = 'Multi-metal_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Drug, Biocide & Metal Resistance' = 'Drug_and_biocide_and_metal_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Copper resistance' = 'Copper_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'CAP' = 'Cationic_antimicrobial_peptides')
follow.c.t <-dplyr::rename(follow.c.t, 'Beta-lactams' = 'betalactams')
follow.c.t <-dplyr::rename(follow.c.t, 'Multi-biocide Resistance' = 'Multi-biocide_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Biocide & Metal Resistance' = 'Biocide_and_metal_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Acid Resistance' = 'Acid_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Arsenic Resistance' = 'Arsenic_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'MDR' = 'Multi-drug_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Nickel Resistance' = 'Nickel_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Zinc Resistance' = 'Zinc_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Acetate Resistance' = 'Acetate_resistance')
follow.c.t <-dplyr::rename(follow.c.t, 'Sodium Resistance' = 'Sodium_resistance')

follow.c <- follow.c.t %>%
  gather(key = key, value = value, 2:ncol(follow.c.t)) %>%
  spread(key=names(follow.c.t)[1], value = 'value') %>%
  dplyr::rename(., Class=key)

follow.c$FOLLOWUP <- as.numeric(follow.c$FOLLOWUP)

follow.c.plot <- ggplot(data = follow.c, aes(x =reorder(Class, -FOLLOWUP), y = FOLLOWUP, fill = Class))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        plot.title = element_text(size=18, face='bold', hjust=0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title='Follow-Ups',fill = 'ARG Group')+
  xlab('ARG Class')+
  ylab('Average Abundance per Sample\n')


# Combine plots #
ggarrange(case.c.plot + rremove('xlab'),
          follow.c.plot,
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 1, nrow = 2)+
  theme(plot.margin = margin(0.5,0.1,0.1,1.5, "cm"))

#################################################
# Top 10 ARG Groups
#################################################

# Cases

case.group <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/AverageAbundance_ALL_TYPES_GroupLevel_Top10_Cases.csv',
                       header = TRUE)

case.group$Group <- factor(case.group$Group, levels = c('MLS23S',
                                                            'RPOB',
                                                            'A16S',
                                                            'TUFAB',
                                                            'MDTB',
                                                            'ACRD',
                                                            'MGTA',
                                                            'MDTC',
                                                            'ACRB',
                                                            'PARC'))

case.group$CASE <- as.numeric(case.group$CASE)

case.g.plot <- ggplot(data = case.group, aes(x =Group, y = CASE, fill = Group))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18, face='bold', hjust=0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Cases', fill = 'ARG Group')+
  xlab('ARG Group')+
  ylab('Average Abundance per Sample\n')


# Follow-Ups

follow.group <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/AverageAbundance_ALL_TYPES_GroupLevel_Top10_FollowUps.csv',
                       header = TRUE)

follow.group$Group <- factor(follow.group$Group, levels = c('MLS23S',
                                                        'TETQ',
                                                        'A16S',
                                                        'CFX',
                                                        'RPOB',
                                                        'ERMF',
                                                        'MEFE',
                                                        'TETW',
                                                        'RRSH',
                                                        'CAP16S'))

follow.group$FOLLOWUP <- as.numeric(follow.group$FOLLOWUP)

follow.g.plot <- ggplot(data = follow.group, aes(x =Group, y = FOLLOWUP, fill = Group))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18, face='bold', hjust=0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(size = 12),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Follow-Ups', fill = 'ARG Group')+
  xlab('ARG Group')+
  ylab('Average Abundance per Sample\n')


# Combine plots #
ggarrange(case.g.plot + rremove('xlab'),
          follow.g.plot,
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 1, nrow = 2)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))


#################################################
# Top 25 ARG Groups
#################################################

# Cases

case.group <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/AverageAbundance_ALL_TYPES_GroupLevel_Top25_Cases.csv',
                       header = TRUE)

case.group$Group <- factor(case.group$Group, levels = c('MLS23S',
                                                        'RPOB',
                                                        'A16S',
                                                        'TUFAB',
                                                        'MDTB',
                                                        'ACRD',
                                                        'MGTA',
                                                        'MDTC',
                                                        'ACRB',
                                                        'PARC',
                                                        'CAP16S', 'ACRF',
                                                        'GYRA','RRSC','RRSH',
                                                        'LPDT','MDTF','CPXAR',
                                                        'GYRBA','COPA','CUSA',
                                                        'ZNTA','YJCG','PARE','PTSL'))

case.group$CASE <- as.numeric(case.group$CASE)

case.g.plot <- ggplot(data = case.group, aes(x =Group, y = CASE, fill = Group))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18, face='bold', hjust=0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle=45, hjust=1, size = 12),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Cases', fill = 'ARG Group')+
  xlab('ARG Group')+
  ylab('Average Abundance per Sample\n')


# Follow-Ups

follow.group <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/AverageAbundance_ALL_TYPES_GroupLevel_Top25_FollowUps.csv',
                         header = TRUE)

follow.group$Group <- factor(follow.group$Group, levels = c('MLS23S',
                                                            'TETQ',
                                                            'A16S',
                                                            'CFX',
                                                            'RPOB',
                                                            'ERMF',
                                                            'MEFE',
                                                            'TETW',
                                                            'RRSH',
                                                            'CAP16S',
                                                            'ANT6','RRSC','TETX',
                                                            'TUFAB','BEXA','ACRB',
                                                            'TET16S','LNUA','O23S',
                                                            'GYRBA','CBLA','CPXAR',
                                                            'MDTB','MDTF','MEFA'))

follow.group$FOLLOWUP <- as.numeric(follow.group$FOLLOWUP)

follow.g.plot <- ggplot(data = follow.group, aes(x =Group, y = FOLLOWUP, fill = Group))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18, face='bold', hjust=0.5),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle=45, hjust=1,size = 12),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Follow-Ups', fill = 'ARG Group')+
  xlab('ARG Group')+
  ylab('Average Abundance per Sample\n')


# Combine plots #
ggarrange(case.g.plot + rremove('xlab'),
          follow.g.plot,
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 1, nrow = 2)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))



#################################################
# Top 10 ARG Classes Combined
#################################################

class.combined <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/2022_01_14/AverageAbundance_ALL_TYPES_ClassLevel_CaseFollowPairs_Top10combined.csv',
                         header=TRUE)

class.combined$Class <- factor(class.combined$Class, levels = c('Other',
                                                                  'Drug and biocide resistance',
                                                                  'Multi-metal resistance',
                                                                  'MLS',
                                                                  'Drug, biocide, and metal resistance',
                                                                  'Aminoglycosides',
                                                                  'Copper resistance',
                                                                  'Fluoroquinolones',
                                                                  'Multi-biocide resistance',
                                                                  'Beta-lactams',
                                                                  'Tetracyclines',
                                                                  'CAP'))

cf_colors <- c('cyan4', 'darkorchid3')

ggplot(data = class.combined, aes(x = Class, y = value, fill = Case.status))+
  geom_bar(stat = 'identity', width = 0.8, position = position_dodge(0.9))+ 
  scale_color_manual(values = cf_colors)+
  scale_fill_manual(values=cf_colors)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'),
        plot.margin = margin(1,1,1,1,'cm'))+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  guides(fill=guide_legend(nrow=2))+
  labs(fill = 'Health Status')+
  ylab('Average Abundance per Genome Equivalent\n')+
  xlab('ARG Class')

#################################################
# Top 10 ARG Groups Combined
#################################################

group.combined <- read.csv('D://Resistome/Reads_based/Abundance/ActualAbundance/ALL_TYPES/2022_01_14/AverageAbundance_ALL_TYPES_GroupLevel_CaseFollowPairs_Top10combined.csv',
                           header=TRUE)

group.combined$Group <- factor(group.combined$Group, 
                               levels = c('MLS23S',
                                          'RPOB',
                                          'A16S',
                                          'TUFAB',
                                          'MDTB',
                                          'ACRD',
                                          'TETQ',
                                          'MGTA',
                                          'MDTC',
                                          'ACRB',
                                          'PARC',
                                          'CAP16S',
                                          'RRSH',
                                          'CFX',
                                          'ERMF',
                                          'MEFE',
                                          'TETW'))

cf_colors <- c('cyan4', 'darkorchid3')

ggplot(data = group.combined, aes(x = Group, y = value, fill = Case.status))+
  geom_bar(stat = 'identity', width = 0.8, position = position_dodge(0.9))+ 
  scale_color_manual(values = cf_colors)+
  scale_fill_manual(values=cf_colors)+
  scale_y_continuous(expand = c(0.01,0))+
  theme(panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=12,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        strip.text.x = element_text(size =10), 
        axis.line = element_line(colour = 'black'),
        plot.margin = margin(1,1,1,1,'cm'))+
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))+
  guides(fill=guide_legend(nrow=2))+
  labs(fill = 'Health Status')+
  ylab('Average Abundance per Genome Equivalent\n')+
  xlab('ARG Group')
