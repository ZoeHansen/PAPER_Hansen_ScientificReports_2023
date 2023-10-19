########################################################
# ACC Quantification Code

# Author: Zoe Hansen
# Last modified: 2022.01.12
########################################################

# This code uses the "ACCs_per_sample" information retrieved from a Python script
# Each sample has an associated # of ACCs and this code is designed to visualize
# these differences. 

# Load libraries and data

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(vegan)
library(calibrate)

accs <- read.csv('D://ACC_analysis/ERIN_total_ACCs_per_sample_CaseFollowPairs_wControls.csv',
                 header = TRUE)

meta = read.csv('D://ERIN_Metagenomes_Metadata_60_CaseFollowPairs_wControls.csv', 
                 header = TRUE)

meta_sub <- meta %>%
  dplyr::select(ER_ID, Case.status, Pathogen)

acc_data <- left_join(meta_sub, accs, by=c('ER_ID', 'Case.status'))

# Remove samples that did not align properly/contain ACCs; subset cases and follow-ups only
acc_sub <- acc_data %>%
  filter(!grepl('ER0535', ER_ID))%>%
  filter(!grepl('ER0611', ER_ID))%>%
  filter(!grepl('ER0612', ER_ID)) %>%
  filter(!grepl('Control', Case.status))

acc.means <- acc_sub %>%
  group_by(Case.status, Pathogen) %>%
  summarise(MeanACC=mean(num_accs))

acc.minmax <- acc_sub %>%
  group_by(Case.status, Pathogen)%>%
  summarise(MinACC=min(num_accs), MaxACC=max(num_accs))

case.follow.comp <- list(c('Case','FollowUp'))
colors.cf <- c('cyan4','darkorchid3')

ggplot(data=acc_sub, aes(x=Case.status, y=num_accs, color=Case.status))+
  geom_boxplot() + 
  geom_jitter(aes(shape=Case.status, color=Case.status))+
  scale_color_manual(values=colors.cf) +
  scale_fill_manual(values=colors.cf)+
  facet_wrap(~Case.status, scales="free") +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text.x=element_text(size=12, angle=30, hjust=0.5, vjust=0.55),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        panel.spacing.x = unit(1, 'lines'))+
  labs(
    x = '\nHealth Status',
    y = 'Number of ACCs\n'
  )+
#  stat_pvalue_manual(means.to.plot,   label = "p.adj", tip.length = 0.01)
  stat_compare_means(comparisons = path.comp,
                     label = 'p.format')


# Wilcoxon Test - Case-Follow Pairs
wilcox.test(num_accs ~ Case.status, data = acc_sub, paired = TRUE)

means.to.plot <- tibble::tribble(
  ~group1,    ~group2,         ~p.adj,  ~y.position, ~variable,
  "Case",    "FollowUp",    2.744e-07,     850,   "num_accs"
)





