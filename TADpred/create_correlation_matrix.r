# DESCRIPTION: creates a file-backed probe-probe correlation matrix 
# using the TAD clusters as blocks. Also creates a vector called ld,
# wich is used in further steps.

# INPUT: TAD clusters, test data OR methylation data from separate set

# OUTPUT: TAD block correlation matrix, ld data

# EXECUTE: using the slurm script since this step takes a while and 
# needs a lot of RAM. So dont execute anything in this script,
# just fill in paths etc.

##############
##### SETTINGS

# Only use file names. The paths will be set in the slurm script.

# Methylation data.
# Dataset used to calculate correlations. Can be your test set or a separate
# set. Must be the following format:
#
#                 sample1       sample2      sample3
# cg21870274 -0.002014552 -0.0497791585  0.019817845
# cg23917638  0.039038936 -0.0839184964 -0.010086821
# cg11422233  0.002169797  0.0173729978 -0.003790025

# Can be beta or M values, or residuals.
methylation_data_file<- ("residuals_betas_allMarkosCovariables-withSex_wholeSample.Robj")

# TAD clusters
# Can be any kind of clusters as long as they have the following columns:
#
#        cpg  cluster
# cg05638471       1
# cg22696995       1
# cg12612065       2
#
# and the cluster numbering needs to be sequential without any gaps.
# If you have tens of thousands of clusters, you might want to comment out
# the logging info in the loop: "if (i %% 100 == 0)" etc to avoid extremely
# long logging files.
tad_clusters_file<- ("tad_clusters.RData")

# You also need to make sure that the variable names correspond to the ones 
# used in the script.

# So once the data is loaded:
load(methylation_data_file)
load(tad_clusters_file)

# You need to assign the correct variable names to the ones that will be
# used in the script. So dont change the names methylation_data and cluster_data. 
methylation_data<- residuals_all_chr
cluster_data<- tad_clusters

# Remove your old variables for cleanliness
rm(residuals_all_chr, tad_clusters)

# Save paths

# The absolute path for the file-backed matrix must be hardcoded here.
# Should be your analysis path plus /corr-data
corr_data_path<- paste("/tsd/p33/data/durable/characters/kristofs/",
                       "TADpred_TOP_SCZ/corr-data", sep="")

# File names for correlation matrix and ld (RData)
matrix_path<- "correlation_matrix.RData"
ld_path<- "ld.RData"

###############################################################################
# The rest should work as it is.

# install packages if not installed
if (!requireNamespace("bigsnpr", quietly = TRUE)) {
  install.packages("bigsnpr")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

library(Matrix)
library(bigsnpr)
library(dplyr)

# version check 
print("bigsnpr version:")
packageVersion("bigsnpr")

# logging
print("methylation data:")
str(methylation_data)

# transpose residuals
data_t<- t(methylation_data)
print("transposed:")
str(data_t)

# not needed anymore
rm(methylation_data)

print("clusters:")
str(cluster_data)

# create the matrix

print("number of clusters:")
(max_cluster <- max(cluster_data$cluster))

tmp <- tempfile(tmpdir = "corr-data")

singletons<- vector()

print("starting loop...")
for (i in 1:max_cluster) {
  if (i %% 100 == 0) {
    cat("Number of clusters processed:", i, "out of", max_cluster, "\n")
  }
  
  # Extract cpg entries where cluster == i
  cluster_cpgs<- cluster_data$cpg[cluster_data$cluster == i]
  
  if (length(cluster_cpgs) == 1) {
    singletons <- c(singletons, cluster_cpgs)
  } else {
    # Create a subset DataFrame using the probes in cluster_cpgs
    subset_df<- data_t[, cluster_cpgs, drop = FALSE]
    
    correlations<- cor(subset_df)
    correlations <- as(correlations, "sparseMatrix")
    
    if (i == 1) {
      print("cluster number:")
      print(i)
      print("cluster cpgs:")
      print(str(cluster_cpgs))
      print("subset_df:")
      print(str(subset_df))
      print("correlations:")
      print(str(correlations))
      ld <- Matrix::colSums(correlations^2)
      print("attempting to initiate SFBM...")
      corr <- as_SFBM(correlations, tmp, compact = TRUE)
      print("first cluster processed. SFBM created.")
      print("corr:")
      print(str(corr))
    } else {
      ld <- c(ld, Matrix::colSums(correlations^2))
      corr$add_columns(correlations, nrow(corr))
    }
  } 
}
print("number of singletons:")
print(length(singletons))

singleton_matrix <- sparseMatrix(i = integer(0), j = integer(0), x = double(0),
                                 dims = c(length(singletons), length(singletons)))
colnames(singleton_matrix)<- singletons
rownames(singleton_matrix)<- singletons
diag(singleton_matrix) <- 1

ones_vector <- rep(1, length(singletons))
names(ones_vector) <- singletons

ld <- c(ld, ones_vector)
corr$add_columns(singleton_matrix, nrow(corr))

print("length of final ld:")
length(ld)
print("final corr:")
str(corr)

print("changing backing file name")

print("current file name:")
(file_name <- basename(corr$backingfile))

print("new tmp path:")
print(corr_data_path)

new_backingfile <- file.path(corr_data_path, file_name)

corr$backingfile <- new_backingfile

print("saving matrix and ld...")
save(corr, file = matrix_path)
save(ld, file = ld_path)

print("thank_you_and_goodnight")

