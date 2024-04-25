# DESCRIPTION: This script clusters the CpGs in your analysis into TADs.
# The resulting cluster table is used to create the correlation matrix in 
# the next step.

# INPUT: TAD map, intersecting CpGs vector

# OUTOUT: TAD clusters

# EXECUTE: with this script

library(yaml)
library(dplyr)

###########
##### FILES

# Import paths from config.
config<- yaml::read_yaml("config.yaml")

# Load TAD map
load("CpG_TAD_maps/cpg_tad_map.RData")
# Inspect
str(cpg_tad_map)

# Load intersecting CpGs
intersect_cpgs_path<- file.path(config$analysis_path, "intersect_cpgs.RData")
load(intersect_cpgs_path)
# Inspect
str(intersect_cpgs)

# Path for saving the result file to later
save_path<- file.path(config$analysis_path, "tad_clusters.RData")

############################################################
##### Some QC of the TAD map before clustering analysis data.
##### Can be skipped.

# Check the number of TAD regions.
(num_unique <- length(unique(cpg_tad_map$TADid)))

# Check the sizes of the TAD regions
frequencies<- as.vector(table(cpg_tad_map$TADid))
summary(frequencies) # Stats for number of CpGs in TADs

#################
##### FILTER CpGs

# Keep only the CpGs from your intersection vector.
tad_clusters<- cpg_tad_map %>% filter(cpg %in% intersect_cpgs)

# Check if there are CpG duplicates in the TAD map. There is usually one.
nrow(tad_clusters) == length(intersect_cpgs)

# If there are duplicates, remove them. This method is for removing 
# a single duplicate and may not work for multiple.
duplicates <- subset(tad_clusters, duplicated(tad_clusters$cpg))
duplicates$cpg
(rows_with_duplicates <- subset(tad_clusters, cpg == duplicates$cpg))
tad_clusters<- tad_clusters[-as.integer(rownames(rows_with_duplicates)[1]),]
# Double check that it worked
nrow(tad_clusters) == length(intersect_cpgs)

# Check the number of TAD regions remaining.
(num_unique2 <- length(unique(tad_clusters$TADid)))

# Check the sizes of the remaining TAD regions
# Usually some of the biggest TADs are removed,
# which is good for creating the matrix later.
frequencies2<- as.vector(table(tad_clusters$TADid))
summary(frequencies2) # Stats for number of CpGs in TADs

############################
##### ASSIGNING CLUSTER IDS. 
# For matrix creation to work we need to assign sequential
# cluster IDs. The TAD IDs skip numbers.

# Make sure the TAD IDs are ascending
tad_clusters<- tad_clusters[order(tad_clusters$TADid), ]
(is_ascending <- all(diff(tad_clusters$TADid) >= 0))

# Add sequential cluster IDs
tad_clusters$cluster<- with(rle(tad_clusters$TADid),rep(seq_along(lengths), lengths))

# Basic QC of the cluster IDs. 

# should be the same 
length(tad_clusters$cluster) == length(intersect_cpgs)
max(tad_clusters$cluster) == num_unique2

# Special check to see that cluster IDs are sequential.
# We're getting paranoid here but if the cluster IDs are not
# sequential the matrix creation will stop midway, and it
# takes a bit of time so we want to be sure.
check_df <- data.frame(number = numeric(),freq = numeric())
prev_number <- tad_clusters$cluster[1]
prev_freq <- 1
for (i in 2:nrow(tad_clusters)) {
  current_number <- tad_clusters$cluster[i]
  
  if (current_number != prev_number) {
    check_df <- rbind(check_df, data.frame(number = prev_number, freq = prev_freq))
    prev_freq <- 1
    prev_number <- current_number
  } else {
    prev_freq <- prev_freq + 1
  }
}
check_df <- rbind(check_df, data.frame(number = prev_number, freq = prev_freq))

# If this is true it means all cluster IDs are sequential
(is_sequential <- all(diff(check_df$number) == 1))

# Check how many singleton clusters there are
rows_with_1 <- check_df[check_df$freq == 1, ]
(nr_singletons<- length(rows_with_1$number))

#######################
##### Save the clusters
save(tad_clusters, file = save_path)

# Clean up
rm(config, cpg_tad_map, intersect_cpgs, save_path, num_unique, num_unique2,
   frequencies, frequencies2, tad_clusters, duplicates, rows_with_duplicates,
   is_ascending, check_df, prev_freq, prev_number, is_sequential, rows_with_1,
   nr_singletons, current_number, i, intersect_cpgs_path)



