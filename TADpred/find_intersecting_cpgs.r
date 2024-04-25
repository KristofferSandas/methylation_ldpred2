# DESCRIPTION: finds the CpGs that are present in all datasets: 
# training data, test data and TAD map. If you are using a separate set 
# to create the correlation matrix, this set needs to be included here
# as well. You might have to improvise in this script depending on formats.

# INPUT: training data (for example meta-EWAS), test data (for example residuals
# from TOP sample), CpG TAD map 

# OUTPUT: a vector of CpG IDs present in all datasets.

# EXECUTE: with this script.

library(yaml)

# Import paths from config.
config <- yaml::read_yaml("config.yaml")

# Load test data
load(config$test_data_path)
# Inspect
str(residuals_all_chr)

# Load TAD map
load("CpG_TAD_maps/cpg_tad_map.RData")
# Inspect
str(cpg_tad_map)

# Load training data
meta_ewas<- read.csv(config$training_data_path)

# create a subset of columns we need from the meta-ewas
# (LDpred authors recommend using fixed effects.)
cols_to_keep <- c("ProbeID", 
                  "All_Effect_Random", 
                  "All_Effect_SE_Random",
                  "All_P_Random")
meta_ewas<- meta_ewas[, cols_to_keep, drop = FALSE]

# Inspect
str(meta_ewas)

# Remove NA values

# check for NA values in meta_ewas
any(is.na(meta_ewas))

# You can remove with na.omit()
# but if you want some more control of what is happening 
# you can do the following:

# Find the locations of the NA values:
na_matrix <- is.na(meta_ewas)
(num_na <- sum(is.na(meta_ewas))) # number of NA fields
(na_indices <- which(na_matrix, arr.ind = TRUE)) # where they are

# remove the rows with NA values. Might not work for removing 
# multiple rows.
meta_ewas<- meta_ewas[-na_indices[1], , drop = FALSE]

# inspect
str(meta_ewas)

# check for NA values in residuals
any(is.na(residuals_all_chr))

# check for NA values in TAD map
any(is.na(cpg_tad_map))

# NA values can be removed as above.

# Find the intersecting CpGs by ID.
# If you have an additional set for the correlation 
# matrix, include it as well.
intersect_cpgs<- intersect(meta_ewas$ProbeID, rownames(residuals_all_chr))
intersect_cpgs<- intersect(intersect_cpgs, cpg_tad_map$cpg)

# Inspect
str(intersect_cpgs)

# Save to your analysis directory
save_path<- file.path(config$analysis_path, "intersect_cpgs.RData")
# Create dir if it doesnt exist.
if (!dir.exists(dirname(save_path))) {
  dir.create(dirname(save_path), recursive = TRUE)
}
save(intersect_cpgs, file = save_path)

# Clean variables
rm(config, meta_ewas, residuals_all_chr, cpg_tad_map, 
   intersect_cpgs, save_path, na_matrix, na_indices, num_na,
   cols_to_keep)
