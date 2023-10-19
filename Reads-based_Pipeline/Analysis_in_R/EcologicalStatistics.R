########################################################

# Performing Ecological Statistics (Case, Control, Follow-Up and Case vs. Follow)

########################################################
# Load libraries and data

library(tidyverse)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(vegan)
library(calibrate)

pheno <- read.csv('D://ERIN_Metagenomes_Metadata_60_CaseFollowPairs_wControls.csv', header = TRUE) 
args = read.csv('D://Resistome/Reads_based/ERIN_GEnorm_ALL_TYPES_gene_level.csv', header = TRUE)
microbiome = read.csv('D://Microbiome/Reads_based/GEnorm_ERIN_kaiju_CaseControlFollow_SPECIES.csv', header=TRUE)

### If data needs to be transposed (species):

microbiome$Species <- microbiome$Species %>%
  replace_na("Unknown")

micro <- ddply(microbiome, "Species", numcolwise(sum))

### If data needs to be transposed (genus):
microbiome$Genus <- microbiome$Genus %>%
  replace_na('Unknown')

micro <- ddply(microbiome, 'Genus', numcolwise(sum))
micro<- micro[, colSums(micro[,-1] != 0) > 0] 
micro <- micro[rowSums(micro[,-1] != 0) > 0,]

# Assign 'data' variable to microbiome set
data <- micro %>%
  gather(key = key, value = value, 2:ncol(micro)) %>%
  spread(key=names(micro)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

# Overwrite 'data' variable for ARG data (change axis names, labels, limits, etc. accordingly)
# data <- args

data<- data[, colSums(data[,-1] != 0) > 0] 
data <- data[rowSums(data[,-1] != 0) > 0,] 
data[is.na(data)] <- 0

########################################################
# Selecting relevant metadata for analysis
########################################################

# Select the relevant metadata fields to include in analysis
health <- pheno %>%
  dplyr::select(ER_ID, Case.status, Pathogen, Antibiotics, Reassigned.Pair.ID) %>%
  filter(!grepl('Control', Case.status))

health$Antibiotics[health$Antibiotics == ""] <- "No"

IDs_use <- health %>%
  dplyr::select(ER_ID)

enviro <- pheno %>%
  dplyr::select(ER_ID, Case.status, Followup.days, Age.years,Gender, Residence.type, Reassigned.Pair.ID,
                Run, Avg_GS, GenomeEquivalents, Year, Hospital, County, Stool.type) %>%
  filter(!grepl('Control', Case.status))

enviro$Stool.type[enviro$Stool.type == ''] <- 'Unknown'
enviro$Hospital[enviro$Hospital == ''] <- "Unknown"

meta <- left_join(health, enviro, by = c('ER_ID','Reassigned.Pair.ID'))
meta <- meta %>%
  filter(!grepl('Control', Case.status))
meta[is.na(meta)] <- 0

# Merge metadata with AGS-normalized abundances
#Full Case_Control Campylobacter dataset with metadata of interest
data_2 <- left_join(IDs_use, data, by = 'ER_ID')

data_cc <- left_join(health, data_2, by = "ER_ID")

data_cc$ER_ID <- factor(data_cc$ER_ID)
data_cc$Case.status <- factor(data_cc$Case.status)
data_cc$Reassigned.Pair.ID <- factor(data_cc$Reassigned.Pair.ID)
data_cc$Pathogen <- factor(data_cc$Pathogen)
data_cc$Antibiotics <- factor(data_cc$Antibiotics)
data_cc[is.na(data_cc)] <- 0
data_cc <- data_cc[, colSums(data_cc != 0) > 0]
data_cc <- data_cc[rowSums(data_cc != 0) > 0,] 

########################################################
# Calculate and plot within-sample (Alpha) Diversity
########################################################

# Richness (of genes) across our samples
r.c = specnumber(data_cc[,-c(1:5)])

# Shannon diversity
h.c = diversity(data_cc[,-c(1:5)], index = 'shannon')

# Calculate Pielou's Evenness index
pielou.c = h.c/log(r.c)

# Combine alpha diversity data and Case status information
div.c=tibble(data_cc$ER_ID, data_cc$Case.status, data_cc$Pathogen, data_cc$Reassigned.Pair.ID, r.c, h.c, pielou.c)
colnames(div.c)=c("ER_ID", "Case.status", "Pathogen","Pair","Richness", "Shannon", "Pielou")

### Plot alpha diversity ###

# Reshape the data to "long" format
div.c.long=melt(div.c, id.vars=c("ER_ID", "Case.status", 'Pathogen','Pair'))

# Determine the Mean, Min, and Max for these measures
# Health Status
div.means <- div.c %>%
  dplyr::group_by(Case.status) %>%
  dplyr::summarise(MeanRich = mean(Richness), MeanShannon = mean(Shannon), MeanPielou = mean(Pielou)) 

div.minmax <- div.c %>%
  dplyr::group_by(Case.status) %>%
  dplyr::summarise(MinRich = min(Richness), MaxRich = max(Richness), MinS = min(Shannon), 
            MaxS = max(Shannon), MinP = min(Pielou), MaxP = max(Pielou))

# Comparisons
case.follow.comp <- list(c('Case','FollowUp'))
colors.cf <- c('cyan4','darkorchid3')
labels <- c(Richness='Richness', Shannon='Shannon Index', Pielou='Pielou Evenness')

# Plot alpha diversity values for Case Status
alpha.div.plot<- ggplot(data=div.c.long, aes(x=Case.status, y=value))+
  geom_boxplot(outlier.shape = NA) + 
#  geom_point(aes(group=Pair, fill=Case.status),size=2,
#             position = position_dodge(0.2), shape=21)+
  geom_jitter(aes(color=Case.status, shape=Case.status))+
  facet_wrap(~factor(variable, levels=c('Richness','Shannon','Pielou')),
                     scales="free", labeller=labeller(variable=labels)) +
  scale_color_manual(values=colors.health) +
  scale_fill_manual(values=colors.health)+
  theme_bw(base_size = 10) +
  theme(legend.position = 'none', 
        axis.line = element_line(colour = 'black'),
        axis.text.x=element_text(size=12, hjust=0.5),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14, face='bold', hjust=0.5),
        axis.text.y.right = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        strip.text.y = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        panel.spacing.x = unit(1, 'lines'))+
  labs(
    x = '\nHealth Status',
    y = 'Alpha Diversity Value\n',
  )+
  stat_pvalue_manual(means.to.plot,   label = "p.adj", tip.length = 0.01)

### Wilcoxon Signed-Rank Test (for Case-Follow Pairs Only) 
cf.r <- wilcox.test(Richness ~ Case.status, data = div.c, paired = TRUE)
cf.s <- wilcox.test(Shannon ~ Case.status, data = div.c, paired = TRUE)
cf.e <- wilcox.test(Pielou ~ Case.status, data = div.c, paired = TRUE)

### ARGS
means.to.plot <- tibble::tribble(
  ~group1,    ~group2,         ~p.adj,  ~y.position, ~variable,
  "Case",    "FollowUp",     1.9e-12,  560,       "Richness",
  "Case",    "FollowUp",     1.9e-13,  6.4,       "Shannon",
  "Case",    "FollowUp",     3.5e-10,  1.2,     "Pielou",
)

### Microbiome
means.to.plot <- tibble::tribble(
  ~group1,    ~group2,         ~p.adj,  ~y.position, ~variable,
  "Case",    "FollowUp",     3.5e-8,     16500,       "Richness",
  "Case",    "FollowUp",     3.4e-6,     4,       "Shannon",
  "Case",    "FollowUp",     2.6e-5,     0.5,     "Pielou",
)

########################################################
# Calculate and plot across-sample (Beta) diversity
########################################################

data_cc <- data_cc[, colSums(data_cc != 0) > 0]
data_cc <- data_cc[rowSums(data_cc != 0) > 0,] 

#Calculate Bray-curtis dissimilarity. 
bc.d.c=vegdist(data_cc[,-c(1:5)], method="bray")

#Map desired colors to the Case statuses to make our future legend for the PCoA.  
Class=rep('black' ,nrow(data_cc))
Class[data_cc$Case.status=="Case"]= 'cyan4'
Class[data_cc$Case.status=='FollowUp']= 'darkorchid3'

shape=rep(1, nrow(data_cc))
shape[data_cc$Case.status=='Case']=21
shape[data_cc$Case.status=='FollowUp']=22
shape[data_cc$Antibiotics=='Yes']=24

# Bacteria
bact=rep('black',nrow(data_cc))
bact[data_cc$Pathogen=='Campylobacter (CA)']='dodgerblue2'
bact[data_cc$Pathogen=='Salmonella (SA)']='firebrick2'
bact[data_cc$Pathogen=='Shigella (SH)']='goldenrod2'
bact[data_cc$Pathogen=='STEC (EC)']='mediumpurple2'

# Sequencing Run 
runid=rep('black',nrow(data_cc))
runid[env.sub$Run=='1']='springgreen'
runid[env.sub$Run=='2']='purple'
runid[env.sub$Run=='3']='gold'
runid[env.sub$Run=='4']='tomato'

######## PCoA #########
### Perform PCoA 
bc.c.pcoa=cmdscale(bc.d.c, eig=TRUE, x.ret=TRUE, k=3) 

# Calculate percent variance explained by each axes 1 and 2
ax1.v.bc.c=bc.c.pcoa$eig[1]/sum(bc.c.pcoa$eig)
ax2.v.bc.c=bc.c.pcoa$eig[2]/sum(bc.c.pcoa$eig)
ax3.v.bc.c=bc.c.pcoa$eig[3]/sum(bc.c.pcoa$eig)

# Add environmental factors if desired
env.sub=meta
envEF.bc=envfit(bc.c.pcoa, env.sub, permutations = 999, na.rm = TRUE)
envEF.bc

### PLOTTING
#### Plot the ordination for Case Status for Axes 1&2 ####
plot(bc.c.pcoa$points[,1],bc.c.pcoa$points[,2],
     cex=1.5,
     pch=shape,
     bg=Class,
     xlab= paste("PCoA1: ",100*round(ax1.v.bc.c,3),"% var. explained",sep=""), 
     ylab= paste("PCoA2: ",100*round(ax2.v.bc.c,3),"%var. explained",sep=""));#,
beta.div.plot <- recordPlot()

# Add environmental variables
plot(envEF.bc, p.max=0.01, col='black', cex = 1)#,
     labels=list(vectors=paste(c('Follow-up Days','Age (years)'))))


### Ellipses (run code following base plot generation)
# Case Follow
ordiellipse(bc.c.pcoa, groups=data_cc$Case.status, 
            col=c('cyan4','darkorchid3'),
            draw = 'lines', kind = 'sd', conf = 0.95, alpha = 0.05, label = TRUE, lwd = 2)

# Bacterial Pathogen 
ordiellipse(bc.c.pcoa, groups=data_cc$Pathogen, 
            col=c('dodgerblue2','firebrick2','goldenrod2','mediumpurple2'),
            draw = 'lines', kind = 'sd', conf = 0.95, alpha = 0.05, label = TRUE, lwd = 2)

# Sequencing run
ordiellipse(bc.c.pcoa, groups=env.sub$Run, draw = 'lines', 
            col=c('springgreen','purple','gold','tomato'),
            kind = 'ehull', conf = 0.95, alpha = 0.05, label = TRUE, lwd = 2)


### Legend (run code following base plot generation)
# Case Follow
legend('topleft', c('Case','FollowUp','Antibiotics'), pch = c(21,22,24), 
       pt.bg=c('cyan4', 'darkorchid3','black'), lty=0)

# Pathogen
legend('topleft',c('Campylobacter (CA)','Salmonella (SA)','Shigella (SH)','STEC (EC)','Case','FollowUp','Antibiotics'), 
       pch=c(21,21,21,21,21,22,24), pt.bg=c('dodgerblue2','firebrick2','goldenrod2','mediumpurple2','black','black','black'))

# Sequencing run
legend('topleft', c('Run 1','Run 2','Run 3','Run 4','Case','FollowUp','Antibiotics'), pch = c(21,21,21,21,21,22,24), 
       pt.bg=c('springgreen','purple','gold','tomato','black', 'black','black'), lty=0)

### Text Labels
#Add ID labels to points (if desired)
textxy(X=bc.c.pcoa$points[,1], Y=bc.c.pcoa$points[,2],labs=data_cc$Reassigned.Pair.ID, 
       cex=0.6, pos=2)


############## COMBINE BASE AND GGPLOT (alpha and beta diversity) #####################

alpha.beta.comb<- plot_grid(alpha.div.plot, beta.div.plot, 
          rel_heights = c(1, 0.7),
          ncol=1, nrow=2,
          labels = c("A", "B"))

ggsave('D://Manning_ERIN/ERIN_FullDataset_AIM_TWO/ThirdAnalysis_MEGARes_v2/microbiome_alpha_beta_div_CaseControlFollow_updated_20230819.eps',
       plot=alpha.beta.comb,
       width=25, height=38, 
       units='cm',dpi=300, device='eps')


##### PCoA EnvFit - Investigate species or genus ####
 
# Intrinsic variables (species or genus)
pcoa.spp.fit <- envfit(bc.c.pcoa, data_cc[,c(6:ncol(data_cc))], permutations=999)
head(pcoa.spp.fit)

pcoa.spp.df <- data.frame((pcoa.spp.fit$vectors)$arrows, (pcoa.spp.fit$vectors)$r,
                          (pcoa.spp.fit$vectors)$pvals)

write.csv(pcoa.spp.df, 'D://Microbiome/Reads_based/EcologicalStatistics/ERIN_PCoA_EnvFit_intrinsic_VECTORS_GENUS_CaseFollowPairs.csv',
          row.names=TRUE)
write.csv(pcoa.spp.df, 'D://Resistome/Reads_based/Ecological_Statistics/ALL_TYPES/ERIN_PCoA_EnvFit_intrinsic_VECTORS_GROUP_CaseFollowPairs.csv',
          row.names=TRUE)

# Plot Intrinsic Variables

# First, subset based on R2 value (otherwise, way too many arrows)
#__FUNCTION: select.envfit__# (found from https://www.researchgate.net/post/How_do_I_set_cutoff_r_values_for_plotting_arrows_from_function_envfit_in_R)
# function (select.envit) filters the resulting list of function (envfit) based on their p values. This allows to display only significant values in the final plot.

select.envfit<-function(fit, r.select){ #needs two sorts of input: fit= result of envfit, r.select= numeric, correlation minimum threshold
  for (i in 1:length(fit$vectors$r)) { #run for-loop through the entire length of the column r in object fit$vectors$r starting at i=1
    if (fit$vectors$r[i]<r.select) { #Check wether r<r.select, i.e. if the correlation is weaker than the threshold value. Change this Parameter for r-based selection
      fit$vectors$arrows[i,]=NA #If the above statement is TRUE, i.e. r is smaller than r.select, then the coordinates of the vectors are set to NA, so they cannot be displayed
      i=i+1 #increase the running parameter i from 1 to 2, i.e. check the next value in the column until every value has been checked
    } #close if-loop
  } #close for-loop
  return(fit) #return fit as the result of the function
} #close the function

pcoa.spp.fit2 <-select.envfit(pcoa.spp.fit, r.select=0.82) # Perform select.envfit on dataset
plot(pcoa.spp.fit2, p.max=0.001, col = "black", cex = 0.7)#, arrow.mul=0.8)


# Extrinsic variables (Environmental factors)
pcoa.env.fit <- envfit(bc.c.pcoa, env.sub, permutations = 999, na.rm = TRUE)
pcoa.vect.df <- data.frame((pcoa.env.fit$vectors)$arrows, (pcoa.env.fit$vectors)$r,
                          (pcoa.env.fit$vectors)$pvals)
pcoa.fact.df <- data.frame((pcoa.env.fit$factors)$r,
                           (pcoa.env.fit$factors)$pvals)

write.csv(pcoa.fact.df,'D://Microbiome/Reads_based/EcologicalStatistics/ERIN_PCoa_EnvFit_extrinsic_FACTORS_genus_CaseFollowPairs.csv',
          row.names = TRUE)

########################################################
# PERMANOVA
########################################################
# Centroid
a.c = adonis(bc.d.c~data_cc$Case.status, distance=TRUE, permutations=1000)

# Disperson
b.c=betadisper(bc.d.c, group=data_cc$Case.status)
permutest(b.c)
plot(b.c)
boxplot(b.c)

# Post-hoc Tukey's Honestly Significant Difference 
TukeyHSD(b.c, which = "group", ordered = FALSE,conf.level = 0.95)

# Performing a pairwise PERMANOVA:

### For some reason, this code isn't working for me anymore :( Troubleshoot later
# (from P. Martinez Arbizu: https://github.com/pmartinezarbizu/pairwiseAdonis)

pairwise.adonis <- function(x,factors, sim.method, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] 
                ~ factors[factors %in% c(as.character(co[1,elem]),(co[2,elem]))] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

PW.Adonis= pairwise.adonis(x, factors, sim.method = 'bray', p.adjust.m='bonferroni')
