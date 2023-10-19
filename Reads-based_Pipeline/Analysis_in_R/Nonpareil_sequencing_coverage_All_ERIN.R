#################################################

# Nonpareil - Sequencing Coverage (Reads-based pipeline)

#################################################
# Load libraries

library(Nonpareil)
library(tidyverse)
library(ggplot2)

#################################################
# Assessing Coverage with Nonpareil
#################################################
# Dataset Preparation

# Nonpareil output
non_out <- read.csv('D://Sequencing_Metrics/Nonpareil/ERIN_Nonpareil_npo_output_files.csv', header = TRUE, as.is = TRUE)
non_out$File <- as.character(non_out$File)
non_out$ER_ID <- as.character(non_out$ER_ID)
non_out$color <- as.character(non_out$color)

# Metadata file used to extract relevant sample information
meta <- read.csv('D://ERIN_Metagenomes_Metadata_60_CaseFollowPairs_wControls.csv', header = TRUE)

coverage <- left_join(meta, non_out, by = 'ER_ID')

# Plot Nonpareil coverage curves
# (Figure S1)

nps<- Nonpareil.curve.batch(files=coverage$File[1:91], col=coverage$colors, labels='', plot=TRUE, 
                            plot.opts=list(plot.observed=FALSE), star = 95)


# Show the estimated values of Nonpareil
non_output <- as.data.frame(print(nps))
rownames(non_output) <- coverage$ER_ID

write.csv(non_output, 'D://Sequencing_Metrics/Nonpareil/nonpareil_statistics_CaseFollowPairs_20220128.csv', row.names=FALSE)

# Show summary for 'kappa'
summary(non_output$kappa)

# Show current coverage of the dataset (as %)
summary(non_output$C)

# Extract Nonpareil sequence diversity (Nd) index
summary(non_output$diversity)

# Extract actual sequencing effort (in Gbp)
summary(non_output$LR)

# Extract sequencing effort for nearly complete coverage (in Gbp)
summary(non_output$LRstar)

# Predict coverage for a sequencing effort of 10Gbp
summary(sapply(nps$np.curves, predict, 10e9))
