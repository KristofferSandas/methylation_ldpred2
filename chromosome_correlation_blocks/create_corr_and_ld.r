# creates a  sparse file-backed correlation matrix 
# should be run in slurm - takes a long time

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

# load test data
residuals_path<- "residuals_subset_ordered.RData"
load(residuals_path)

print("residuals:")
str(residuals_subset_ordered)

# transpose data
residuals_t<- t(residuals_subset_ordered)

print("transposed:")
str(residuals_t)

# not needed anymore
rm(residuals_subset_ordered)

# load chromosome clusters
# variable name from CMR run kept to minimize changes in code
co_meth_regions_path<- "chr_clusters_EPIC.RData"
load(co_meth_regions_path)

co_meth_regions<- chr_clusters_EPIC

print("chromosome clusters:")
str(co_meth_regions)

# create the matrix

print("number of chromosomes (should be 22):")
(max_cluster <- max(co_meth_regions$CoMethylCluster))

tmp <- tempfile(tmpdir = "corr-data")

singletons<- vector()

print("starting loop...")
for (i in 1:max_cluster) {
  if (i %% 1 == 0) {
    cat("Number of chromosomes processed:", i, "out of", max_cluster, "\n")
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
      print(str(cluster_cpgs))
      print("subset_df:")
      print(str(subset_df))
      print("correlations:")
      str(correlations)
      ld <- Matrix::colSums(correlations^2)
      print("attempting to initiate SFBM")
      corr <- as_SFBM(correlations, tmp, compact = TRUE)
      print("first cluster processed. SFBM created.")
      print("corr:")
      str(corr)
    } else {
      ld <- c(ld, Matrix::colSums(correlations^2))
      corr$add_columns(correlations, nrow(corr))
    }
  } 
}
print("number of singletons (should be 0:")
print(length(singletons))

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
                      "chromosome_clusters/",
                      "corr-data", sep=""
                      ))

new_backingfile <- file.path(new_tmp_path, file_name)

corr$backingfile <- new_backingfile

# save the corr and ld
save(ld, corr, file = "corr_and_ld.RData")

