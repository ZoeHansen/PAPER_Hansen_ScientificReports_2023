#################################################

# MMUPHin - Differentially Abundant Taxa & ARGs

#################################################

library(MMUPHin)
library(magrittr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(vegan)

# Load data (for this, we will use relative abundance/proportions since MMUPHin doesn't like our normalized abundances :( ))

# Note: metadata should contain samples as rows and feature data should contain samples as columns

# Microbiome
abd.data <- read.csv('D://Microbiome/Reads_based/Abundance/RelativeAbundance/RelativeAbundance_GENUS_CaseFollowPairs.csv', header=TRUE)

# Resistome
abd.data <- read.csv('D://Resistome/Reads_based/Abundance/RelativeAbundance/ALL_TYPES/RelativeAbundance_ALL_TYPES_GroupLevel_CaseFollowPaired.csv')

# Need to remove a couple of columns and transpose resistome RA (SKIP for taxonomy!!)
abd.data <- abd.data %>%
  dplyr::select(.,-c(Case.status, Case.Follow_ID, Bacteria))

abd.data <- abd.data %>%
  gather(key = key, value = value, 2:ncol(abd.data)) %>%
  spread(key=names(abd.data)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

# Convert taxonomy/ARGs to rownames 
rownames(abd.data) <- abd.data$ER_ID
abd.data <- abd.data %>%
  dplyr::select(-ER_ID)

# Make sure to sort ER_ID by descending before importing (*insert eye roll here*)
meta.data <- read.csv('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/ERIN_Metagenomes_Metadata_60_CaseFollowPairs.csv',
                      header=TRUE)

meta.data <- meta.data %>%
  dplyr::select('ER_ID','Case.status','Pathogen','Run','Case.Follow_ID',
                'Followup.days','Residence.type','Age.years','Hospital',
                'Avg_GS','GenomeEquivalents','Antibiotics', 'Gender','Year')

meta.data$Antibiotics[meta.data$Antibiotics == ""] <- "No"
meta.data$Residence.type[meta.data$Residence.type == ''] <- 'Unknown'
meta.data$Hospital[meta.data$Hospital == ''] <- "Unknown"

meta.data$Run <- factor(meta.data$Run, levels=c('1','2','3','4'))
meta.data$Year <- as.numeric(meta.data$Year)
meta.data$Followup.days <- as.integer(meta.data$Followup.days)
meta.data$Case.status <- factor(meta.data$Case.status, levels=c('Case','FollowUp'))
meta.data$Antibiotics <- factor(meta.data$Antibiotics, levels=c('No','Yes'))

rownames(meta.data) <- meta.data$ER_ID

#### Batch Correction ####

# First, we will attempt to correct for batch effects due to Sequencing Run (only necessary for microbiome data)

fit_adjust_batch <- adjust_batch(feature_abd=abd.data,    # abundance matrix
                                 batch='Run',             # the batch variable requiring correction
                                 covariates='Case.status',# covariates to adjust for (i.e. Case vs. Follow-up status)
                                 data=meta.data,          # metadata file
                                 control=list(verbose=FALSE)) # other specifiers, if desired

# Isolate adjusted abundance table
micro_abd_adj <- fit_adjust_batch$feature_abd_adj

# Perform PERMANOVA to examine variability due to Sequencing Run before and after adjustment

dist_before <- vegdist(t(abd.data), method='bray')
dist_after <- vegdist(t(micro_abd_adj), method='bray')

fit_adonis_before <- adonis(dist_before ~ Run, meta.data)
fit_adonis_after <- adonis(dist_after ~ Run, meta.data)

print(fit_adonis_before)
print(fit_adonis_after)

# Note: compare amount of variability explained using the R2 measure; check for significance
# This may depend on the resolution of taxa (e.g. phylum vs. genus vs. species)

#### Meta-analytical Differential Abundance ####

# 'lm_meta()' combines the functionality of MaAsLin2 and vegan to retrieve differentially abundant
# features between groups

# General setup of model
fit_lm_meta <- lm_meta(feature_abd=abd.data,       # designate abundance data
                       batch='Run',                # list the variable to correct for (batch)
                       exposure='Case.status',     # group variable
                       covariates=c(''),           # Other environmental variables to account for in the model
                       data=meta.data,             # metadata file
                       control=list(verbose=FALSE))# other specifications (if desired)

#meta.data$Run <- as.numeric(meta.data$Run)

# Run the model
fit_lm_meta <- lm_meta(feature_abd=abd.data,       
                       batch='Run',                
                       exposure='Case.status',     
                       covariates=c('Age.years','Avg_GS','GenomeEquivalents','Year',
                                    'Antibiotics'),           
                       data=meta.data,             
                       control=list(verbose=FALSE))

# Microbiome: envfit() pulled out:
  # age in years
  # sequencing run (batch)
  # Average genome size
  # genome equivalents
  # year of sampling
  # use of Abx
  # Follow-up days and hospital were nearly significant -- could include

# Resistome: envfit() pulled out: 
  # age in years
  # year of sampling
  # genome equivalents
  # residence type, follow-up days and hospitals were trending towards significance

# Unfortunately, MaAsLin2 is not allowing categorical covariates with more than two levels
# because of a reference-call issue :(
# I can't figure out how to include these with numeric identifiers as a work-around...it's 
# definitely limiting. 

# View the 'meta_fits' output
meta_fits <- fit_lm_meta$meta_fits

meta_fits_signif <- meta_fits %>%
  filter(qval.fdr < 0.05)         # pull out significant p-adjusted values

write.csv(meta_fits_signif, 'D://MMUPHin/MMUPHin_ARG_GENUS_CaseStatus_CaseFollowPairs.csv', row.names=FALSE)

# Plot the significant findings (differential abundance)
meta_fits_signif %>%
  filter(abs(coef) > 0.065) %>%
  arrange(coef) %>% 
  mutate(feature = factor(feature, levels = feature)) %>% 
  ggplot(aes(y = coef, x = feature, fill=coef>0)) +
  geom_bar(stat = "identity") +
  coord_flip()+
  scale_fill_manual(values=c("cyan4","darkorchid3"),
                    labels=c("Case","FollowUp"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=10))+
  theme(legend.title = element_text(face='bold',size=12),
        legend.text= element_text(size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(colour = 'black',linetype = 'dashed'),
        panel.spacing = unit(0,'lines'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.25,0.25,0,0.25), 'cm'))+
  labs(
    x = '\nGenus\n',
    y = '\nCoefficient\n',
    fill = 'Health Status')


##### Investigating Continuous Population Structure #####

# Identify variables driving continuous structure

###### Cases & Follow-Ups Separate ######
### Cases ### 

case.fit_continuous <- continuous_discover(feature_abd = case.abd.adj,
                                      batch = "Run",
                                      data = case.meta,
                                      control = list(var_perc_cutoff = 0.75,
                                                     verbose = FALSE))

case.loading <- data.frame(feature = rownames(case.fit_continuous$consensus_loadings),
                      loading1 = case.fit_continuous$consensus_loadings[, 1])

case.loading %>%
  arrange(-abs(loading1)) %>%
  slice(1:20) %>%
  arrange(loading1) %>%
  mutate(feature = factor(feature, levels = feature)) %>%
  ggplot(aes(x = feature, y = loading1)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle("Features with top loadings")

case.mds <- cmdscale(d = D_case)
colnames(case.mds) <- c("Axis1", "Axis2")
as.data.frame(case.mds) %>% 
  mutate(score1 = case.fit_continuous$consensus_scores[, 1]) %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = score1)) +
  geom_point() +
  coord_fixed()


### Follow-Ups ###
follow.fit_continuous <- continuous_discover(feature_abd = follow.abd.adj,
                                           batch = "Run",
                                           data = follow.meta,
                                           control = list(var_perc_cutoff = 0.75,
                                                          verbose = FALSE))

follow.loading <- data.frame(feature = rownames(follow.fit_continuous$consensus_loadings),
                           loading1 = follow.fit_continuous$consensus_loadings[, 1])

follow.loading %>%
  arrange(-abs(loading1)) %>%
  slice(1:20) %>%
  arrange(loading1) %>%
  mutate(feature = factor(feature, levels = feature)) %>%
  ggplot(aes(x = feature, y = loading1)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  ggtitle("Features with top loadings")

follow.mds <- cmdscale(d = D_follow)
colnames(follow.mds) <- c("Axis1", "Axis2")
as.data.frame(follow.mds) %>% 
  mutate(score1 = follow.fit_continuous$consensus_scores[, 1]) %>% 
  ggplot(aes(x = Axis1, y = Axis2, color = score1)) +
  geom_point() +
  coord_fixed()


##### Cases + Follow-Ups Together ######
cf.fit_continuous <- continuous_discover(feature_abd = micro_abd_adj,
                                             batch = "Run",
                                             data = meta.data,
                                             control = list(var_perc_cutoff = 0.4,
                                                            verbose = FALSE))

cf.loading <- data.frame(feature = rownames(cf.fit_continuous$consensus_loadings),
                             loading1 = cf.fit_continuous$consensus_loadings[, 1])

shape=rep(1, nrow(meta.data))
shape[meta.data$Case.status=='Case']='Case'
shape[meta.data$Case.status=='FollowUp']='FollowUp'
shape[meta.data$Antibiotics=='Yes']='Received Antibiotics'

cf.loading.data <- cf.loading %>%
  arrange(-abs(loading1)) %>%
  slice(1:30) %>%
  arrange(loading1) %>%
  mutate(feature = factor(feature, levels = feature))

write.csv(cf.loading.data, 'D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/MMUPHin/MMUPHin_ARG_CLASS_LoadingData_CaseStatus_CaseFollowPairs.csv',
          row.names=FALSE)

arg_loading<-ggplot(data=cf.loading.data, aes(x = feature, y = loading1, fill=loading1>0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values=c("blue4","firebrick4"),
                    labels=c("Score < 0","Score > 0"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=10))+
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(colour = 'black',linetype = 'dashed'),
        panel.spacing = unit(0,'lines'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.5,0.25,0,0.25), 'cm'),
        plot.title=element_text(face='bold',hjust=0.5, vjust=2))+
  ggtitle("Features (ARG Tradeoffs)")+
  labs(
    x = '\nARG Group\n',
    y = '\nLoading 1\n')

cf.mds <- cmdscale(d = dist_after)
colnames(cf.mds) <- c("Axis1", "Axis2")
cf.mds.data <- as.data.frame(cf.mds) %>% 
  mutate(score1 = cf.fit_continuous$consensus_scores[, 1])

write.csv(cf.mds.data, 'D://MMUPHin/MMUPHin_ARG_CLASS_MDS_ScoringData_CaseStatus_CaseFollowPairs.csv', row.names=FALSE)

arg_mds<-ggplot(data=cf.mds.data, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color=score1, fill=score1, shape=shape),
             size=3) +
  scale_shape_manual(values=c(16, 15, 17))+
  guides(fill=guide_colorbar(title = "Score",
                             label.position='right',
                             title.position='top',
                             title.vjust=1,
                             frame.color='black',
                             barwidth=1.5,
                             barheight=15),
         shape=guide_legend(override.aes=list(size=5)),
         color='none',
         size='none')+
  scale_color_viridis(discrete = FALSE, option = 'H')+
  scale_fill_viridis(discrete = FALSE, option = 'H')+
  coord_fixed()+
  theme(legend.title = element_text(face='bold',size=12),
        legend.text= element_text(size=12),
        legend.key= element_rect(fill=NA),
        legend.key.size = unit(0.75, "cm"),
        legend.spacing.y = unit(0.5, 'cm'),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0,'lines'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.25,0.25,0,0.25), 'cm'))+
  labs(fill = "Score",
       shape = "Health Status",
       x='\nMDS1',
       y='MDS2\n')

##### Create Combination Plot #####

comb_plot<-ggarrange(genus_loading,
          genus_mds,
          arg_loading,
          arg_mds,
          labels = c("A", "B","C","D"),
          font.label=list(color="black",size=24),
          ncol = 2, nrow = 2,
          align='hv')+
  theme(plot.margin = margin(0.5,0.1,0.1,0.5, "cm"))

ggsave('MMUPHin_LoadingsOrdination_CombinedPlot_GENUS_GROUP_v2.eps',
       width=50, height=45, units='cm', dpi=300,
       device='eps')

