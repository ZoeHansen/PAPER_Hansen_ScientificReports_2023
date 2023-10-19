########################################################

# Mantel Test

########################################################

library(tidyr)
library(vegan)
library(plyr)

pair.mat <- read.csv('D://PairwiseDissimilarity_Tests/ERIN_CaseFollowPairs_binary_matrix_MantelTest.csv', header=TRUE, row.names = 1)

# Change 'pair.mat' variable to a matrix
pair.mat <- as.matrix(pair.mat)

#### Microbiome ####
micro.abd.data <- read.csv('D://Microbiome/Reads_based/GEnorm_ERIN_kaiju_60CaseFollowPairs_SPECIES.csv', header=TRUE)

# If microbiome data: 
micro.abd.data$ER_ID <- micro.abd.data$ER_ID %>%
  replace_na("Unknown")

micro <- ddply(micro.abd.data, "ER_ID", numcolwise(sum))

micro<- micro[, colSums(micro[,-1] != 0) > 0] 
micro <- micro[rowSums(micro[,-1] != 0) > 0,]

data <- micro %>%
  gather(key = key, value = value, 2:ncol(micro)) %>%
  spread(key=names(micro)[1], value = 'value') %>%
  dplyr::rename(., ER_ID=key)

data <- data[, colSums(data != 0) > 0]
data <- data[rowSums(data != 0) > 0,] 


# Calculate Bray-Curtis Dissimilarity

micro.bc.dist=vegdist(data[,-1], method="bray")
micro.bc.dist.matrix <- as.matrix(micro.bc.dist)
colnames(micro.bc.dist.matrix) <- rownames(micro.bc.dist.matrix) <- data[['ER_ID']]

# Perform the Mantel Test
micro_mantel <- mantel(micro.bc.dist.matrix, pair.mat, 
                       method='spearman', 
                       permutations=9999,
                       na.rm=TRUE)

##### Resistome #####

arg.abd.data <- read.csv('D://Resistome/Reads_based/ERIN_GEnorm_ALL_TYPES_gene_level_CaseFollowPairs.csv', header=TRUE)

arg.abd.data <- arg.abd.data[, colSums(arg.abd.data != 0) > 0]
arg.abd.data <- arg.abd.data[rowSums(arg.abd.data != 0) > 0,]
arg.bc.dist = vegdist(arg.abd.data[,-1], method='bray')
arg.bc.dist.matrix <- as.matrix(arg.bc.dist)
colnames(arg.bc.dist.matrix) <- rownames(arg.bc.dist.matrix) <- arg.abd.data[['ER_ID']]

arg_mantel <- mantel(arg.bc.dist.matrix, pair.mat,
                     method = 'spearman',
                     permutations=9999,
                     na.rm=TRUE)

##### Microbiome vs. Resistome #####

arg_micro_mantel <- mantel(micro.bc.dist.matrix, arg.bc.dist.matrix,
                           method='spearman',
                           permutations=9999,
                           na.rm=TRUE)

#############################################################################
################ Plot Histogram of Dissimilarity Measures ###################

# Retrieve dissimilarity values for each pairwise comparison

### Binary Pairs Matrix 

pairs.pairwise.dist <- t(combn(colnames(pair.mat),2))
pairs.pairwise.df <- data.frame(pairs.pairwise.dist, dist=pair.mat[pairs.pairwise.dist])
colnames(pairs.pairwise.df) <- c('sample1','sample2','dist')

# Extract case-follow pair comparisons

case.follow.pairs <- pairs.pairwise.df %>%
  dplyr::filter(dist == 1) %>%
  dplyr::select(-dist)

### Microbiome

micro.pairwise.dist <- t(combn(colnames(micro.bc.dist.matrix), 2))
micro.pairwise.df <- data.frame(micro.pairwise.dist, dist=micro.bc.dist.matrix[micro.pairwise.dist])
colnames(micro.pairwise.df) <- c('sample1','sample2','dist')

# Calculate the mean for pairwise distances across samples
all.mean <- mean(micro.pairwise.df$dist)


# Extract actual paired Case-Follow samples and calculate mean
cf.pairs.dist <- dplyr::left_join(case.follow.pairs, micro.pairwise.df, by=c('sample1','sample2'))
cf.pairs.mean <- mean(cf.pairs.dist$dist)

# Perform t-test to statistically compare means
t.test(micro.pairwise.df$dist, cf.pairs.dist$dist, 
       var.equal = FALSE, # specify that the variance will differ between groups
       conf.level = 0.95,
       alternative = 'greater') # the is that 'x' has a greater mean than 'y'

# Set the breaks for the histogram
hist.breaks <- seq(0,1, by=0.05)

hist(micro.pairwise.df$dist, 
     breaks=hist.breaks,
     ylim=c(0,1300),
     col='lightgray',
     main='Pairwise Bray-Curtis Dissimilarity',
     xlab= 'Pairwise Dissimilarity')

# Add lines to show the different means
lines(x=c(cf.pairs.mean,cf.pairs.mean), 
      y=c(cf.pairs.mean,1200),
      col='purple4',
      lwd=3)

lines(x=c(all.mean, all.mean),
      y=c(all.mean, 1100),
       col='firebrick2',
       lwd=3)

# Add values of mean above lines
text(x=cf.pairs.mean, y=1250,
     labels=round(cf.pairs.mean,3))
text(x=all.mean, y=1150,
     labels=round(all.mean,3))

# Add legend
legend(0.7, 800,
       c('Paired Mean','Overall Mean'),
       lty=1, 
       col=c('purple4','firebrick2'),
       lwd=3)


### Resistome

arg.pairwise.dist <- t(combn(colnames(arg.bc.dist.matrix), 2))
arg.pairwise.df <- data.frame(arg.pairwise.dist, dist=arg.bc.dist.matrix[arg.pairwise.dist])
colnames(arg.pairwise.df) <- c('sample1','sample2','dist')

# Calculate the mean for pairwise distances across samples
all.mean <- mean(arg.pairwise.df$dist)

# Extract acutal paired Case-Follow samples and calculate mean
cf.pairs.dist <- dplyr::left_join(case.follow.pairs, arg.pairwise.df, by=c('sample1','sample2'))
cf.pairs.mean <- mean(cf.pairs.dist$dist)

# Perform t-test to statistically compare means
t.test(arg.pairwise.df$dist, cf.pairs.dist$dist, 
       var.equal = FALSE, # specify that the variance will differ between groups
       conf.level = 0.95,
       alternative = 'greater') # the is that 'x' has a greater mean than 'y'

# Set the breaks for the histogram
hist.breaks <- seq(0,1, by=0.025)

hist(arg.pairwise.df$dist, 
     breaks=hist.breaks,
     ylim=c(0,1000),
     col='lightgray',
     main='Pairwise Bray-Curtis Dissimilarity',
     xlab= 'Pairwise Dissimilarity')

# Add lines to show the different means
lines(x=c(cf.pairs.mean,cf.pairs.mean), 
      y=c(cf.pairs.mean,1200),
      col='purple4',
      lwd=3)

lines(x=c(all.mean, all.mean),
      y=c(all.mean, 1100),
      col='firebrick2',
      lwd=3)

# Add values of mean above lines
text(x=cf.pairs.mean, y=1250,
     labels=round(cf.pairs.mean,3))
text(x=all.mean, y=1150,
     labels=round(all.mean,3))

# Add legend
legend(0.2, 1000,
       c('Paired Mean','Overall Mean'),
       lty=1, 
       col=c('purple4','firebrick2'),
       lwd=3)

