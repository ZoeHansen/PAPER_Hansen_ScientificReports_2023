########################################################

# Aim 2: Analyzing Assembly Quality (QUAST Output)

########################################################
# Load libraries and data

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(vegan)
library(calibrate)

data <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Sequencing_Metrics/QUAST_MEGAHIT/ERIN_QUAST_assembly_statistics.csv',
                  header = TRUE)
data_t <- data %>%
  gather(key = key, value = value, 2:ncol(data)) %>%
  spread(key=names(data)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

pheno = read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/ERIN_Metagenomes_Metadata_60_CaseFollowPairs_wControls.csv', 
                header = TRUE)

pheno_sub <- pheno %>%
  dplyr::select(ER_ID, Case.status, Pathogen)

quast <- left_join(pheno_sub, data_t, by='ER_ID')

########################################################
# Selecting relevant metadata for analysis
########################################################

# If needed, remove any samples from analysis here
quast <- quast %>%
  filter(!grepl('Control', Case.status)) 

# Remove any lingering zeroes to avoid downstream errors

quast$Case.status <- factor(quast$Case.status)
quast$Pathogen <- factor(quast$Pathogen)

########################################################
# Calculate and plot within-sample (Alpha) Diversity
########################################################

# Isolate variables of interest

quast_sub <- quast %>%
  dplyr::select(ER_ID, Case.status, Pathogen, `# contigs`, `Total length`,
                N50,`Avg. coverage depth`)%>%
#  dplyr::select(ER_ID, Case.status, Pathogen, `# contigs`, `Total length`,
#                `GC (%)`, N50,`Mapped (%)`, `Avg. coverage depth`)%>%
  filter(!grepl('ER0535', ER_ID))%>%
  filter(!grepl('ER0611', ER_ID))%>%
  filter(!grepl('ER0612', ER_ID))

names(quast_sub)[4] <- 'Num.contigs'
names(quast_sub)[5] <- 'Total_length'
names(quast_sub)[7] <- 'Avg.coverage_depth'

# Contig lengths

### Plot alpha diversity ###

# Determine the Mean, Min, and Max for these measures
# Health Status
div.means <- quast_sub %>%
  group_by(Case.status) %>%
  summarise(MeanContigs = mean(Num.contigs), MeanTotalLength = mean(Total_length),
            MeanGC = mean(GC_percent), MeanN50=mean(N50), MeanPercentMapped=mean(MappedReads_percent),
            MeanAvgCoverageDepth=mean(Avg.coverage_depth))


# Reshape the data to "long" format
div.quast.long=melt(quast_sub, id.vars=c("ER_ID", "Case.status", 'Pathogen'))


# Comparisons
health_comparisons <- list(c("Case", "Control"), c("FollowUp", "Control"), c("Case", "FollowUp"))
colors.health = c('cyan4', 'darkorange3', 'darkorchid3')

case.follow.comp <- list(c('Case','FollowUp'))
colors.cf <- c('cyan4','darkorchid3')


labels <- c(Num.contigs='Number of contigs', Total_length='Total Contig Length', 
            N50='N50', Avg.coverage_depth='Avg. Coverage Depth')
#labels <- c(Num.contigs='Number of contigs', Total_length='Total Contig Length', 
#            GC_percent='GC content (%)', N50='N50', MappedReads_percent='Reads Mapped (%)', 
#            Avg.coverage_depth='Avg. Coverage Depth')


div.quast.long$variable <- factor(div.quast.long$variable, 
                                  levels=c("Num.contigs", "Total_length", "N50", "Avg.coverage_depth"),
                                  ordered = TRUE)


new_order <- c("Num.contigs", "Total_length", "N50", "Avg.coverage_depth")
div.quast.long <- div.quast.long %>% 
  mutate(variable=variable %>% fct_relevel(new_order))

# Plot alpha diversity values for Case Status
quast.plot<- ggplot(data=div.quast.long, aes(x=Case.status, y=value, color=Case.status))+
  geom_boxplot() + 
#  geom_point(aes(group=Pair, fill=Case.status),size=2,
#             position = position_dodge(0.2), shape=21)+
  geom_jitter(aes(shape=Case.status, color=Case.status))+
#  geom_line(aes(group=Pair), position = position_dodge(0.2), color='gray44')+ # For connecting paired data
  facet_wrap(~variable, scales="free", 
             labeller=labeller(variable=labels)) +
  scale_color_manual(values=colors.health) +
  scale_fill_manual(values=colors.health)+
  theme_bw(base_size = 10) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text.x=element_text(size=12),#, angle=40, hjust=1),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=18),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        panel.spacing.x = unit(1, 'lines'))+
  labs(
    x = '\nHealth Status',
    y = 'Assembly Statistic\n'
  )+
  stat_pvalue_manual(means.to.plot,   label = "p.adj", tip.length = 0.01)
#  stat_compare_means(comparisons = health_comparisons,
#                     label = 'p.format')


# Get the stats (Case, Control, Follow)
compare_means(Num.contigs ~ Case.status, quast_sub, method = 'wilcox.test')
compare_means(Total_length ~ Case.status, quast_sub, method = 'wilcox.test')
compare_means(GC_percent ~ Case.status, quast_sub, method = 'wilcox.test')
compare_means(MappedReads_percent ~ Case.status, quast_sub, method = 'wilcox.test')
compare_means(N50 ~ Case.status, quast_sub, method = 'wilcox.test')
compare_means(Avg.coverage_depth ~ Case.status, quast_sub, method = 'wilcox.test')


# Case-Follow Pairs
wilcox.test(Num.contigs ~ Case.status, data = quast_sub, paired = TRUE)
wilcox.test(Total_length ~ Case.status, data = quast_sub, paired = TRUE)
wilcox.test(GC_percent ~ Case.status, data = quast_sub, paired = TRUE)
wilcox.test(MappedReads_percent ~ Case.status, data = quast_sub, paired = TRUE)
wilcox.test(N50 ~ Case.status, data = quast_sub, paired = TRUE)
wilcox.test(Avg.coverage_depth ~ Case.status, data = quast_sub, paired = TRUE)

means.to.plot <- tibble::tribble(
  ~group1,    ~group2,         ~p.adj,  ~y.position, ~variable,
  "Case",    "FollowUp",      4.9e-5,     2.0e+05,   "Num.contigs",
  "Case",    "Control",       1.0e-7,     1.8e+05,   "Num.contigs",
  "FollowUp", "Control",      0.069,      1.65e+05,   "Num.contigs",
  "Case",    "FollowUp",      2.2e-4,     3.0e+08,  "Total_length",
  "Case",    "Control",       1.5e-5,     2.5e+08,  "Total_length",
  "FollowUp", "Control",      0.28,       2.88e+08,  "Total_length",
  "Case",    "FollowUp",      0.036,      7.8e+04,   "N50",
  "Case",    "Control",       1.5e-5,     6.2e+04,   "N50",
  "FollowUp", "Control",      0.0042,     7.0e+04,   "N50",
  "Case",    "FollowUp",      0.011,      71,        "Avg.coverage_depth",
  "Case",    "Control",       0.00063,    60,        "Avg.coverage_depth",
  "FollowUp", "Control",      0.38,       65,        "Avg.coverage_depth"
)


ggsave('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/Sequencing_Metrics/QUAST_MEGAHIT/QUAST_assembly_stats_CaseControlFollow_updated_20230818.eps',
       plot=quast.plot,
       width=20, height=30, 
       units='cm',dpi=300, device='eps')
