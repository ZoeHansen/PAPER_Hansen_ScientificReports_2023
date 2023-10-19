#################################################
# Sequencing Metrics 
# Author: Zoe Hansen
#################################################
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(outliers)

#################################################
# Exploring Average Genome Size (AGS) and Genome Equivalents (GE)
#################################################

meta <- read.csv('D://ERIN_Metagenomes_Metadata_60_CaseFollowPairs_wControls.csv', header=TRUE)

# Use the metadata dataframe from above to extract the Average Genome Size (AGS) and Genome Equivalents (GE)
ags <- meta %>%
  select(ER_ID, Case.status, Avg_GS, GenomeEquivalents, Pathogen)#%>%
#  filter(!grepl('Control',Case.status))

#Determine means
c_means <- ags %>%
  group_by(Case.status)%>%
  summarise(MeanAGS = mean(Avg_GS), MeanGE = mean(GenomeEquivalents))

# Convert to "long" format
ags_long <- melt(ags, id.vars=c('ER_ID','Case.status', 'Pathogen'))

#Generate the Plot
cf.pal = c('cyan4','darkorchid4')
labels <- c(Avg_GS = 'Average Genome Size', GenomeEquivalents = "Genome Equivalents")
cf_comparisons <- list(c('Case', 'FollowUp'))

ags_plot<-ggplot(data=ags_long, aes(x=Case.status, y=value))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Case.status, color=Case.status))+
  facet_wrap(~variable, scales="free_y",labeller=labeller(variable=labels)) +
  scale_color_manual(values=cf.pal)+
  theme_bw(base_size = 12) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text.x = element_text(size=12, hjust=0.5,vjust=0.5),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size=12),
        panel.spacing.x = unit(2, 'lines'),
        plot.margin=unit(c(1,1,0.5,0.5), 'cm'))+
  labs(
    x = 'Health Status\n',
    y = 'Value\n'
  )+
  stat_pvalue_manual(means.to.plot,   label = "p.adj", tip.length = 0.01)


# Compare Means
wilcox.test(Avg_GS ~ Case.status, data = ags, paired = TRUE)
wilcox.test(GenomeEquivalents ~ Case.status, data = ags, paired = TRUE)

means.to.plot <- tibble::tribble(
  ~group1,    ~group2,         ~p.adj,  ~y.position, ~variable,
  "Case",    "FollowUp",    0.086,       8.0e+06,    "Avg_GS",
  "Case",    "FollowUp",    0.029,       950,    "GenomeEquivalents",
)

########################################################
# Exploring Quality/Filtering Statistics from trimming and host removal
########################################################

##### Trimmomatic Statistics 

trim <- read.csv('D://Sequencing_Metrics/AmrPlusPlus/TrimmomaticStats.csv',
                 header=TRUE)

trim_meta <- meta %>%
  select(ER_ID, Case.status, Pathogen)


trim.filt <- left_join(trim_meta, trim, by='ER_ID')
trim.filt <- trim.filt %>%
  dplyr::select(-c(ForwardOnlySurviving, ReverseOnlySurviving,BothSurviving, PercentSurviving)) %>%
  filter(!grepl('Control', Case.status))

#Determine means
trim_means <- trim.filt %>%
  group_by(Case.status)%>%
  summarise(MeanInput = mean(NumberOfInputReads), MeanSurvive = mean(TotalSurviving),
            MeanDrop = mean(Dropped))

trim_minmax <- trim.filt %>%
  group_by(Case.status)%>%
  summarise(MinInput=min(NumberOfInputReads), MaxInput=max(NumberOfInputReads), MinSurvive=min(TotalSurviving),
            MaxSurvive=max(TotalSurviving), MinDrop=min(Dropped), MaxDrop=max(Dropped))

# Convert to "long" format
trim_long <- melt(trim.filt, id.vars=c('ER_ID','Case.status', 'Pathogen'))
                 
#Generate the Plot
labels <- c(NumberOfInputReads = 'Input Reads', TotalSurviving='Surviving Reads', PercentSurviving='Percent Surviving')

survive.plot<-ggplot(data=trim.filt, aes(x=Case.status, y=TotalSurviving))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Case.status, color=Case.status))+
  #facet_wrap(~Case.status, scales="free_y",labeller=labeller(variable=labels)) +
  scale_color_manual(values=cf.pal)+
  theme_bw(base_size = 12) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text.x=element_text(size=12, hjust=0.5, vjust=0.5),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size=12),
        panel.spacing.x = unit(2, 'lines'),
        plot.margin=unit(c(1,1,0.5,0.5), 'cm'))+
  labs(
    x = 'Health Status\n',
    y = 'Number of Surviving Reads\n'
  )+
  stat_pvalue_manual(means.to.plot,   label = "p.adj", tip.length = 0.01)

# Wilcoxon Test - Case-Follow Pairs
wilcox.test(NumberOfInputReads ~ Case.status, data = trim.filt, paired = TRUE)
wilcox.test(TotalSurviving ~ Case.status, data = trim.filt, paired = TRUE)

means.to.plot <- tibble::tribble(
  ~group1,    ~group2,         ~p.adj,  ~y.position, ~variable,
  "Case",    "FollowUp",    0.5147,     8.0e+06,    "NumberOfInputReads",
  "Case",    "FollowUp",    0.5147,     8.0e+06,    "TotalSurviving"
)

###### Host Removal Statistics #####

host <- read.csv('D://Sequencing_Metrics/AmrPlusPlus/HostRemovalStats.csv',
                 header=TRUE)
host.filt <- left_join(ags, host, by='ER_ID')
host.filt <- host.filt %>%
  dplyr::select(-c(Avg_GS, GenomeEquivalents, Pathogen, NumberOfInputReads))%>%
  filter(!grepl('Control',Case.status))

#Determine means
host_means <- host.filt %>%
  group_by(Case.status)%>%
  summarise(MeanMapped = mean(Mapped), MeanUnmapped = mean(Unmapped))

host_minmax <- host.filt %>%
  group_by(Case.status)%>%
  summarise(MinMapped=min(Mapped), MaxMapped=max(Mapped), MinUnmapped=min(Unmapped), MaxUnmapped=max(Unmapped))

# Convert to "long" format
host_long <- melt(host.filt, id.vars=c('ER_ID','Case.status'))

#Generate the Plot
labels <- c(Mapped = 'Mapped Reads (Host)', Unmapped='Unmapped Reads (Non-host)')

ggplot(data=host.filt, aes(x=Case.status, y=Unmapped))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Case.status, color=Case.status))+
#  facet_wrap(~variable, scales="free_y",labeller=labeller(variable=labels)) +
#  scale_color_brewer(palette = 'Dark2') +
  scale_color_manual(values = cf.pal)+
  theme_bw(base_size = 12) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text.x=element_text(size=12),#, angle=30, hjust=1),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 14),
        strip.text.x = element_text(size=12),
        panel.spacing.x = unit(2, 'lines'),
        plot.margin=unit(c(1,1,0.5,0.5), 'cm'))+
  labs(
    x = 'Health Status\n',
    y = 'Number of Non-Host Reads\n'
  )+
  stat_pvalue_manual(means.to.plot,   label = "p.adj", tip.length = 0.01)

# Wilcoxon Test - Case-Follow Pairs
wilcox.test(Unmapped ~ Case.status, data = host.filt, paired = TRUE)

means.to.plot <- tibble::tribble(
  ~group1,    ~group2,         ~p.adj,  ~y.position, ~variable,
  "Case",    "FollowUp",    0.04,       1.80e+07,    "Unmapped",
)

