#################################################

# ANCOM-BC for Differentially Abundant Taxa & ARGs

#################################################

# ANCOM-BC requires a phyloseq object in order to run. The first section of code produces this
# phyloseq object using metagenomic data. 

# This code was found in a tutorial for formatting metagenomic data to be used in 
# Phyloseq. The original markdown can be found here: https://mvuko.github.io/meta_phyloseq/

# Load libraries and data

library(tidyverse)
library(phyloseq)
library(ANCOMBC)

#### Load and create your 'OTU' table #####

# Microbiome
otu <- read.csv('D://Microbiome/Reads_based/GEnorm_ERIN_kaiju_60CaseFollowPairs_PHYLUM.csv', header=TRUE)

# Resistome (requires transposition)
otu <- read.csv('D://Resistome/Reads_based/ERIN_GEnorm_ALL_TYPES_class_level_CaseFollowPairs.csv', header=TRUE)

otu <- otu %>%
  gather(key = key, value = value, 2:ncol(otu)) %>%
  spread(key=names(otu)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

otu<- otu[, colSums(otu[,-1] != 0) > 0] 
otu <- otu[rowSums(otu[,-1] != 0) > 0,]

# Double check the structure
str(otu)

# OTU IDs must be the rownames
rownames(otu)<- otu$ER_ID

# Double check to see that your read counts match that expected
count <- colSums(otu[, c(2:ncol(otu))], na.rm = TRUE)
count

###### Load and create your Taxonomy table #####
taxa<- otu %>% 
  dplyr::select(ER_ID)  
str(taxa)

# Taxa need to be recognized as factors 
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
str(taxa)

# Add the first column of OTU tables so that Taxa and OTU tables match
taxa<- cbind(otu$ER_ID, taxa)

#colnames(taxa)[1]<- "Genus"
colnames(taxa)[1] <- 'ARG_class'
str(taxa)

# Rownames must be consistent across tables
rownames(taxa)<- taxa$ARG_class

#now you can delete the first columns of both OTU and taxa tables
# because they are now the rownames, and therefore redundant in the first column
otu<- otu %>% 
  dplyr::select(-ER_ID)

taxa<- taxa %>% 
  dplyr::select(-ER_ID)

###### Load your metadata #####

meta <- read.csv('D://ERIN_Metagenomes_Metadata_60_CaseFollowPairs.csv', header = TRUE)

meta <- meta %>%
  dplyr::select(ER_ID, Case.status, Pathogen, Case.Follow_ID, Run)

colnames(meta) <- c('ER_ID','Case.status','Pathogen', 'PairID','Seq_Run')
str(meta)

# Make sample names into the rownames
rownames(meta)<- meta$ER_ID

#and delete the first column because it is now redundant
meta<- meta %>% 
  dplyr::select(-ER_ID)

# Clarify your factor levels 
meta$Case.status<- factor(meta$Case.status, levels = c("Case","FollowUp"))
meta$Pathogen <- factor(meta$Pathogen, levels = c('Campylobacter (CA)','Salmonella (SA)',
                                                  'Shigella (SH)', 'STEC (EC)'))

##### Running Phyloseq #####

# OTU and Taxonomy tables must be matrices
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)

#transform data to phyloseq objects
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples<- sample_data(meta)

#and put them in one object
phylo_object<- phyloseq(phylo_OTU, phylo_TAX, phylo_samples)

#check if everything looks good
sample_sums(phylo_object)       #should sum up the number of all reads per sample
sample_names(phylo_object)      #sample names
rank_names(phylo_object)        #taxa levels  
sample_variables(phylo_object)  #factors

otu_table(phylo_object)[1:3, 1:2]
taxa_names(phylo_object)[1:5]

##### Running ANCOM-BC #####

ancom_out <- ancombc(phyloseq = phylo_object,
                     formula = 'Case.status', # + Seq_Run',
                     p_adj_method = "BH", # Benjamini-Hochberg correction (FDR)
                     zero_cut = 0.9,      # removes taxa that have a greater proportion of zeroes than 0.9
                     lib_cut = 0,         # samples with librar sizes < lib_cut are excluded
                     group = 'Case.status', # group variable in metadata
                     struc_zero = TRUE,  # whether to detect structural zeroes
                     neg_lb = FALSE,      # whether to classify taxon as structural zero using lower bound
                     tol = 1e-05,         # iteration convergence tolerance for E-M algorithm
                     max_iter = 100,      # max # iterations for E-M algorithm
                     conserve = TRUE,     # whether to use conservative variance estimate of test statistic (recommended for small 'n' or large # diff. abd. taxa)
                     alpha = 0.05,        # level of significance
                     global = TRUE)       # whether to perform global test (across groups)


# Note: for our microbiome data, I included sequencing run as a confounder/covariate in the formula
# This does not need to be included for the resistome since envfit() didn't pull it out for now

results <- ancom_out$res  # isolate the results output from ANCOM-BC

ancom_results <- results %>%
  as_tibble(results)

ancom_signif <- ancom_results %>%  # Subset the significant adjusted p-values (q-values)
  filter(q_val < 0.05)

ancom_signif_beta <- as.data.frame(ancom_signif$beta)
ancom_signif_se <- as.data.frame(ancom_signif$se)
ancom_signif_w <- as.data.frame(ancom_signif$W)
ancom_signif_pval <- as.data.frame(ancom_signif$p_val)
ancom_signif_qval <- as.data.frame(ancom_signif$q_val)
ancom_signif_da <- as.data.frame(ancom_signif$diff_abn)

ancom_signif_names <- as.data.frame(rownames(ancom_signif$beta)) 
#colnames(ancom_signif_names) <- 'Genus'
colnames(ancom_signif_names) <- 'ARG_class'

ancom_signif_final <- as.data.frame(cbind(ancom_signif_names, 
                                          ancom_signif_beta,
                                          ancom_signif_se,
                                          ancom_signif_w,
                                          ancom_signif_pval,
                                          ancom_signif_qval,
                                          ancom_signif_da))

write.csv(ancom_signif_final, "D://ANCOM-BC/ANCOMBC_ARG_CLASS_CaseStatusRun_CaseFollowPairs.csv",
          row.names = FALSE)

# Note: the coefficients in the 'beta' column refer to fold-change values relative to the reference
# In this circumstance, "Cases" are our reference level, and so any positive FC would indicate that
# that particular taxon is more highly represented in Follow-ups than Cases
# Conversely, a negative 'beta' value would indicate that the FC in Follow-ups is exp(-beta) 
# compared to cases (meaning that taxon is more abundant in cases). 


###############################################################################

########### Plotting ANCOM-BC Output #################

ancom_signif_final <- read.csv('D://ANCOM-BC/ANCOMBC_ARG_GROUP_CaseStatusRun_CaseFollowPairs.csv', header=TRUE)

ancom.data <- ancom_signif_final %>%
  filter(Case.statusFollowUp.qval<0.05)%>%
  arrange(-abs(Case.statusFollowUp.beta)) %>%
  slice(1:40) %>%
  arrange(Case.statusFollowUp.beta) %>%
  mutate(ARG_group = factor(ARG_group, levels = ARG_group))

ggplot(data=ancom.data, aes(x = ARG_group, y = Case.statusFollowUp.beta, fill=Case.statusFollowUp.beta>0)) +
  geom_bar(stat = "identity") +
  coord_flip() +
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
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.5,0.25,0,0.25), 'cm'),
        plot.title=element_text(face='bold',hjust=0.5, vjust=2))+
  labs(
    x = '\nARG Group\n',
    y = '\nCoefficient\n',
    fill = 'Health Status')
