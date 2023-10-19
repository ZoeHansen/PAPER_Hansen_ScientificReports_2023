##############################################

# Actual Abundance of Taxa

##############################################
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(viridis)
#################################################
# Top 10 PHYLA
#################################################

### PHYLUM ###

# Cases

case.phylum <- read.csv('D://Microbiome/Reads_based/Abundance/ActualAbundance/AverageAbundance_Top10_PHYLUM_CaseOnly.csv',
                        header=TRUE, stringsAsFactors = FALSE)

case.phylum$Phylum <- factor(case.phylum$Phylum, levels = case.phylum$Phylum[order(-case.phylum$value)])

case.p.plot<- ggplot(data = case.phylum, aes(x = Phylum, y = value, fill = Phylum))+
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
        plot.title = element_text(size=18, hjust=0.5, face='bold'),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Cases',fill = 'Phylum')+
  ylab('Average Abundance\n')


# Follow-Ups

follow.phylum <- read.csv('D://Microbiome/Reads_based/Abundance/ActualAbundance/AverageAbundance_Top10_PHYLUM_FollowUpOnly.csv',
                        header=TRUE, stringsAsFactors = FALSE)

follow.phylum$Phylum <- factor(follow.phylum$Phylum, levels = follow.phylum$Phylum[order(-follow.phylum$value)])

follow.p.plot<- ggplot(data = follow.phylum, aes(x = Phylum, y = value, fill = Phylum))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+ 
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE, option='C')+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18, hjust=0.5, face='bold'),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Follow-Ups',fill = 'Phylum')+
  ylab('Average Abundance\n')

# Combine plots #
ggarrange(case.p.plot + rremove('xlab'),
          follow.p.plot,
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 1, nrow = 2)+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))

#################################################
# Top 10 GENUS
#################################################

# Case

case.genus <- read.csv('D://Microbiome/Reads_based/Abundance/ActualAbundance/AverageAbundance_Top10_GENUS_CasesOnly.csv',
                       header=TRUE)

case.genus$Genus<- factor(case.genus$Genus, levels = case.genus$Genus[order(-case.genus$CASE)])

case.g.plot<- ggplot(data = case.genus, aes(x = Genus, y = CASE, fill = Genus))+
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
        plot.title = element_text(size=18, hjust=0.5, face='bold'),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Cases',fill = 'Genus')+
  ylab('Average Abundance\n')

# Follow-Ups

follow.genus <- read.csv('D://Microbiome/Reads_based/Abundance/ActualAbundance/AverageAbundance_Top10_GENUS_FollowUpsOnly.csv',
                          header=TRUE, stringsAsFactors = FALSE)

follow.genus$Genus <- factor(follow.genus$Genus, levels = follow.genus$Genus[order(-follow.genus$FOLLOW)])

follow.g.plot<- ggplot(data = follow.genus, aes(x = Genus, y = FOLLOW, fill = Genus))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+ 
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE, option='C')+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18, hjust=0.5, face='bold'),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Follow-Ups',fill = 'Genus')+
  ylab('Average Abundance\n')

# Combine plots #
ggarrange(case.g.plot + rremove('xlab'),
          follow.g.plot,
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 1, nrow = 2)+
  theme(plot.margin = margin(0.5,0.1,0.1,2.5, "cm"))


#################################################
# Top 25 GENUS
#################################################

# Case

case.genus <- read.csv('D://Microbiome/Reads_based/Abundance/ActualAbundance/AverageAbundance_Top25_GENUS_CasesOnly.csv',
                       header=TRUE)

case.genus$Genus<- factor(case.genus$Genus, levels = case.genus$Genus[order(-case.genus$CASE)])

case.g.plot<- ggplot(data = case.genus, aes(x = Genus, y = CASE, fill = Genus))+
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
        plot.title = element_text(size=18, hjust=0.5, face='bold'),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 50, hjust = 1, size = 12),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Cases',fill = 'Genus')+
  ylab('Average Abundance\n')

# Follow-Ups

follow.genus <- read.csv('D://Microbiome/Reads_based/Abundance/ActualAbundance/AverageAbundance_Top25_GENUS_FollowUpsOnly.csv',
                         header=TRUE, stringsAsFactors = FALSE)

follow.genus$Genus <- factor(follow.genus$Genus, levels = follow.genus$Genus[order(-follow.genus$FOLLOW)])

follow.g.plot<- ggplot(data = follow.genus, aes(x = Genus, y = FOLLOW, fill = Genus))+
  geom_bar(stat = 'identity', width = 0.9, position = 'dodge')+ 
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE, option = 'C')+
  scale_y_continuous(expand = c(0.01,0))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, 'lines'),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18, hjust=0.5, face='bold'),
        axis.text.y = element_text(size = 10), 
        axis.text.x = element_text(angle = 50, hjust = 1, size = 12),
        strip.text.x = element_text(size =14), 
        axis.line = element_line(colour = 'black'))+
  guides(fill=guide_legend(nrow=2))+
  labs(title = 'Follow-Ups',fill = 'Genus')+
  ylab('Average Abundance\n')

# Combine plots #
ggarrange(case.g.plot + rremove('xlab'),
          follow.g.plot,
          labels = c("A", "B"),
          vjust = 0.10,
          ncol = 1, nrow = 2)+
  theme(plot.margin = margin(0.5,0.1,0.1,2.5, "cm"))


#################################################
# Top 10 PHYLA Combined
#################################################

phyla.combined <- read.csv('D://Microbiome/Reads_based/Abundance/ActualAbundance/AverageAbundance_Top10_PHYLUM_CaseFollowCombined.csv',
                           header=TRUE)

phyla.combined$Phylum[5]<- 'Cannot be assigned'
phyla.combined$Phylum[6]<- 'Cannot be assigned'
phyla.combined$Phylum[11]<- 'Unclassified'
phyla.combined$Phylum[12]<- 'Unclassified'

phyla.combined$Phylum <- factor(phyla.combined$Phylum, 
                               levels = c('Bacteroidetes',
                                          'Proteobacteria',
                                          'Cannot be assigned',
                                          'Firmicutes',
                                          'Preplasmiviricota',
                                          'Unclassified',
                                          'Verrucomicrobia',
                                          'Actinobacteria',
                                          'Uroviricota',
                                          'Other'))

cf_colors <- c('cyan4', 'darkorchid3')

ggplot(data = phyla.combined, aes(x = Phylum, y = value, fill = Case.status))+
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
  ylab('Average Abundance per Sample\n')+
  xlab('Phylum')

#################################################
# Top 10 GENUS Combined
#################################################

genus.combined <- read.csv('D://Microbiome/Reads_based/Abundance/ActualAbundance/AverageAbundance_Top10_GENUS_CaseFollowCombined.csv',
                           header=TRUE)

# To exclude "cannot be assigned"
#genus.combined <- genus.combined %>%
#  filter(!grepl('Cannot be assigned', Genus))

genus.combined$Genus <- factor(genus.combined$Genus, 
                                levels = c('Cannot be assigned',
                                           'Bacteroides',
                                           'Other',
                                           'Salmonella',
                                           'Escherichia',
                                           'Alistipes',
                                           'Pseudomonas',
                                           'Parabacteroides',
                                           'Prevotella',
                                           'Faecalibacterium',
                                           'Akkermansia',
                                           'Phocaeicola',
                                           'Klebsiella',
                                           'Unclassified'))

cf_colors <- c('cyan4', 'darkorchid3')

ggplot(data = genus.combined, aes(x = Genus, y = value, fill = Case.status))+
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
  ylab('Average Abundance per Sample\n')+
  xlab('Genus')
