# this should be run in slurm, requires lots of RAM
# creates CMR clusters from the test data

# install comeback package if not installed
if (!requireNamespace("comeback", quietly = TRUE)) {
  print("comeback not installed. installing..")
  install.packages("comeback_1.0.1.tar", repos = NULL, type="source")
}

library(comeback)
print("comeback version check:")
packageVersion("comeback")

# Additional functions from CoMeBack
source("MRS_func.R")

# load test data
load("residuals_subset_ordered.RData")
print("Residuals:")
str(residuals_subset_ordered)

# Transpose data
betas_t<- t(residuals_subset_ordered)
print("Transposed residuals:")
dim(betas_t)

rm(residuals_subset_ordered)

# CoMeBack settings
print("running cmr():")
print('corlo = 0.2, maxprbdst = 100000, corlodst=800, Iarray = "EPIC"')
cmrs = cmr(Mdata = betas_t, corlo = 0.2, maxprbdst = 100000, corlodst=800, Iarray = "EPIC")
    
print("total nr of probes allocated to clusters:")
(alloc_probes<- length(unlist(cmrs)))
print("fraction of probes allocated to clusters:")
alloc_probes/ncol(betas_t)

print("saving cmrs...")
save(cmrs, file = "cmrs.RData")
print("cmrs saved.")

print("running GenCoMeRegion()")
co_meth_regions = GenCoMeRegion(cmrs = cmrs, beta = betas_t, Overlap = F)

print("dim co_meth_regions")
dim(co_meth_regions)

print("saving co_meth_regions...")
save(co_meth_regions, file = "co_meth_regions.RData")
print("cmrs co_meth_regions.")
print("thank_you_and_good_night")
