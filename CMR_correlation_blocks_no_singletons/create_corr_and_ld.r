# creates a co-methyl clustered sparse file-backed correlation matrix without singletons
# should be run in slurm: slurm_create_corr_l.sh

# install package if not installed
if (!requireNamespace("bigsnpr", quietly = TRUE)) {
  install.packages("bigsnpr")
}

# install package if not installed
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

library(Matrix)
library(bigsnpr)
library(dplyr)

print("Matrix version:")
packageVersion("Matrix")
print("bigsnpr version:")
packageVersion("bigsnpr")
print("dplyr:")
packageVersion("dplyr")

# load subset of test data
residuals_path<- "non_singleton_residuals.RData"
load(residuals_path)

print("residuals:")
str(non_singleton_residuals)

# transpose data
residuals_t<- t(non_singleton_residuals)

print("transposed:")
str(residuals_t)

# not needed anymore
rm(non_singleton_residuals)

# load non-singleton cmrs from CoMeBack
co_meth_regions_path<- "non_singletons.RData"
load(co_meth_regions_path)

co_meth_regions<- non_singletons

print("co meth regions:")
str(co_meth_regions)

# create the matrix

print("number of clusters:")
(max_cluster <- max(co_meth_regions$CoMethylCluster))

tmp <- tempfile(tmpdir = "corr-data")

singletons<- vector()

print("starting loop...")
for (i in 1:max_cluster) {
  if (i %% 1000 == 0) {
    cat("Number of clusters processed:", i, "out of", max_cluster, "\n")
  }
  
  # Extract clustercpg entries where CoMethylCluster == i
  cluster_cpgs<- co_meth_regions$clustercpg[co_meth_regions$CoMethylCluster == i]
  
  if (length(cluster_cpgs) == 1) {
    singletons <- c(singletons, cluster_cpgs)
  } else {
    # Create a subset DataFrame using the probes in cluster_cpgs
    subset_df<- residuals_t[, cluster_cpgs, drop = FALSE]
    
    correlations<- cor(subset_df)
    correlations <- as(correlations, "sparseMatrix")
    
    if (i == 1) {
      print("cluster number:")
      print(i)
      print("cluster cpgs:")
      print(cluster_cpgs)
      print("subset_df:")
      dim(subset_df)
      str(subset_df)
      print("correlations:")
      str(correlations)
      correlations
      ld <- Matrix::colSums(correlations^2)
      print("attempting to initiate SFBM")
      corr <- as_SFBM(correlations, tmp, compact = TRUE)
      print("first cluster processed. SFBM created.")
      print("corr:")
      str(corr)
      corr
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

# change backing file name
print("changing backing file name")

print("file name:")
(file_name <- basename(corr$backingfile))

print("new tmp path:")
(new_tmp_path<- paste("/tsd/p33/data/durable/characters/kristofs/",
                      "LDpred_pipeline_1/Step3_create_correaltion_matrix/",
                      "corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/",
                      "no_cometh_singletons/corr-data", sep=""
                      ))

new_backingfile <- file.path(new_tmp_path, file_name)

corr$backingfile <- new_backingfile

# save the corr and ld
save(ld, corr, file = "corr_and_ld.RData")

