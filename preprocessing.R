# The goal of this project is to experiment with classifiers to predict cancer type from the expression data

# The dataset consists of gene expression values extracted for 801 tumor samples
# across 5 cancer types (labels correspond to cancer types). The data
# was extracted for 20,532 genes using RNA-seq (Illumina HiSeq platform).
# and then it was normalized to FPKM (fragments per kilobase of transcript per million mapped reads)

# https://gdc.cancer.gov/about-data/data-harmonization-and-generation/genomic-data-harmonization-0#Pipelines

# 1: Data cleaning and exploration
# 2: Pre-processing: Reduce the dimension of the dataset using feature selection
# 3: Modeling: Use SVM to build a classifier

library(tidyverse)
library(ggplot2)
library(caret)
library(ggfortify)
setwd('/Users/ericaspada/Desktop/projects/PANCAN')

# read in raw data files
exp_values = read_csv('data.csv')
labels = read_csv('labels.csv')

# examine the files
print(str(exp_values))
print(str(labels))
head(exp_values)
head(labels)

# merge the labels file into the expression values file
data = merge(labels, exp_values, by.x='X1', by.y='X1')
data[1:25,1:25] # examine

# many genes have 0 values, get the percent of 0s per gene
pct_zeros = as.data.frame(sapply(data, function(x) {round(length(which(x==0))/length(x)*100,digits=2)}))
pct_zeros = cbind(row.names(pct_zeros), pct_zeros)
colnames(pct_zeros) = c('gene', 'pct')
pct_zeros = as_tibble(pct_zeros)

# remove all genes that have 100% zeros since these won't contribute to the classification
to_remove = pct_zeros %>%
  filter(pct==100.00) %>%
  select(gene)
to_remove = pull(to_remove,gene) #pull creates a vector
data_cleaned = data %>%
  select(-all_of(to_remove))



# look at the data for a few random samples on a qq-norm plot and make histograms
sample_1 = as.numeric(data_cleaned[2,3:20266])
qqnorm_plot_1 = qqnorm(sample_1)
hist_1 = hist(sample_1, breaks=30)
sample_2 = as.numeric(data_cleaned[56,3:20266])
qqnorm_plot_2 = qqnorm(sample_2)
hist_2 = hist(sample_2,breaks=30)
sample_3 = as.numeric(data_cleaned[24,3:20266])
qqnorm_plot_3 = qqnorm(sample_3)
hist_3 = hist(sample_3, breaks=30) 
# looks like a negatively skewed distribution

# Use PCA to visualize the samples in two dimensions
# TO DO: does it need to be log transformed?
prin_comp = prcomp(data_cleaned[,3:20266], scale=FALSE)
summary(prin_comp)$importance[1:2,1:20] # examine
pc1_pc2_plot = autoplot(prin_comp, data=data_cleaned,colour= 'Class') # note: ggfortify package is needed to use autoplot with prcomp objects
print(pc1_pc2_plot)

# Implement the FCBF filter approach (fast correlation based filter)
# FCBF uses symmetrical uncertainty as the correlation measure
# it uses the concept of 'predominant correlation' to select features
library(FCBF)

# transpose the df (this package expects the rows to be features)
data_t = t(data_cleaned)
data_t[1:25,1:25] # examine
print(dim(data_t))

# First step is to discretize the expression values (required for symmetrical uncertainty)
# the discretize() function breaks the values up into thirds (low,high,high)
discr_exp = as.data.frame(discretize_exprs(data_t[3:20266,])) # make a new df of discretized values
discr_exp[1:25,1:25]
class_labels = data_t[2,] # extract the class labels
class_labels = as.factor(class_labels)
features_fcbf_1 = fcbf(discr_exp, class_labels, thresh=.25,verbose=TRUE, samples_in_rows=FALSE)
# this results in 36 features
su_plot(discr_exp,class_labels_all) # histogram of feature counts and correlation
# .15 looks like a better threshold for this dataset
features_fcbf_2 = fcbf(discr_exp, class_labels, thresh=.15,verbose=TRUE, samples_in_rows=FALSE)
# this results in 98 features
data_fcbf = data_cleaned[,c(2,features_fcbf_2$index)] # subset the data to only include the selected features
data_fcbf[1:25,1:25] # examine
