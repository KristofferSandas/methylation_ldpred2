
# load test data:
load("residuals_subset_ordered.RData")

print("residuals:")
str(residuals_subset_ordered)

ids_to_map<- rownames(residuals_subset_ordered)

########################################### debug line
#ids_to_map<- sample(ids_to_map, 500, replace = FALSE)

print("ids to map:")
str(ids_to_map)

rm(residuals_subset_ordered)

# load chromosome data from i450K and EPIC
load("sysdata.rda")

print("sysdata:")
str(init_data)

df_450<- init_data$I450K_Manifest
df_450<- df_450[grepl("^(1[0-9]?|2[0-2]?|[1-9])$", df_450$CHR), ]

print("chromosomes in 450:")
print(unique(df_450$CHR))

print("450:")
str(df_450)

df_EPIC<- init_data$EPIC_Manifest
df_EPIC<- df_EPIC[grepl("^(1[0-9]?|2[0-2]?|[1-9])$", df_EPIC$CHR), ]

print("chromosomes in EPIC:")
print(unique(df_EPIC$CHR))

print("EPIC:")
str(df_EPIC)

rm(init_data)

# match CpGs to chromosomes in 450K
print("matching 450...")
not_found_count_450 <- 0

chr_clusters_450 <- data.frame(clustercpg = character(),
                               CoMethylCluster = numeric(), 
                               stringsAsFactors = FALSE)

# Iterate through each ID in ids_to_map
for (id in ids_to_map) {
  # Check if the ID is in row names of df_450
  if (id %in% rownames(df_450)) {
    # Get the corresponding CHR number from df_450
    chr_number <- df_450[id, "CHR"]
    
    # Add the ID and CHR number to chr_clusters dataframe
    chr_clusters_450 <- rbind(chr_clusters_450, 
                              data.frame(clustercpg = id, 
                                         CoMethylCluster = as.numeric(chr_number)))
  } else {
    # Increment the counter for IDs not found
    not_found_count_450 <- not_found_count_450 + 1
  }
}

print("450 clusters created:")
print(str(chr_clusters_450))

print("chromosomes:")
print(unique(chr_clusters_450$CoMethylCluster))

print("any NA values:")
any(is.na(chr_clusters_450))

print("number of residual ids not found in 450:")
print(not_found_count_450)

rm(df_450)

# match CpGs to chromosomes in EPIC
print("matching EPIC...")
not_found_count_EPIC <- 0

chr_clusters_EPIC <- data.frame(clustercpg = character(),
                                CoMethylCluster = numeric(), 
                                stringsAsFactors = FALSE)

# Iterate through each ID in ids_to_map
for (id in ids_to_map) {
  # Check if the ID is in row names of df_450
  if (id %in% rownames(df_EPIC)) {
    # Get the corresponding CHR number from df_450
    chr_number <- df_EPIC[id, "CHR"]
    
    # Add the ID and CHR number to chr_clusters dataframe
    chr_clusters_EPIC <- rbind(chr_clusters_EPIC, 
                               data.frame(clustercpg = id, 
                                          CoMethylCluster = as.numeric(chr_number)))
  } else {
    # Increment the counter for IDs not found
    not_found_count_EPIC <- not_found_count_EPIC + 1
  }
}

print("EPIC clusters created:")
print(str(chr_clusters_EPIC))

print("chromosomes:")
print(unique(chr_clusters_EPIC$CoMethylCluster))

print("any NA values:")
any(is.na(chr_clusters_EPIC))

print("number of residual ids not found in EPIC:")
print(not_found_count_EPIC)

rm(df_EPIC)

print("saving 450..")
save(chr_clusters_450, file = "chr_clusters_450.RData")
print("done")

print("saving EPIC..")
save(chr_clusters_EPIC, file = "chr_clusters_EPIC.RData")
print("done")

print("thank_you_and_good_night")
