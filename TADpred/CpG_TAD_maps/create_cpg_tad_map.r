# DESCRIPTION: This script maps CpGs from the illumina manifests
# to the TAD regions in the scaffold and saves the results as an RData obj.

# INPUT: tad scaffold, illumina manifest(s)

# OUTPUT: CpG TAD map

# EXECUTE: Use the slurm script CpG_TAD_maps/slurm_create_cpg_tad_maps.sh
# This script does not take long but might need some extra RAM (20G should
# be more than enough). The outfile from slurm works well as a log. 

##### SETTINGS:

# Use only file names, the paths are added in the SLURM script.

# Load tad scaffold data
# Select hg19 for I450K/EPIC1 or hg38 for EPIC2
tad_data <- read.table("tad_scaffold_hg19.txt", header = FALSE)

# Load illumina manifest data
load("I450K_Manifest.RData")
load("EPIC_Manifest.RData")

# If you are using other manifests make sure they have cgIDs as 
# rownames and contain columns CHR and MAPINFO

# Example:
#str(EPIC_Manifest)
# 'data.frame':	X obs. of  2 variables:
# $ CHR    : chr
# $ MAPINFO: int

#head(EPIC_Manifest)
#             CHR MAPINFO
# cg07881041  19  5236016
# cg18478105  20  61847650
# cg23229610  1   6841125

# Rename manifests for name consistency in further processing
# So dont change names df_450 / df_EPIC.
df_450<- I450K_Manifest
df_EPIC<- EPIC_Manifest

# Set to "df_EPIC" if EPIC or EPIC2 is used.
# Set to "df_450" if I450K is used.
# Set to "both" if both EPIC1 and I450K should be used.
# "both" finds the most CpGs, but EPIC2 and I450K cant 
# be combined since they use different genome builds.
use_manifest<- "both"

# Create a unique name for the RData file where the
# TAD map will be saved so you can tell different maps
# apart.
TAD_map_file_name<- "cpg_tad_map.RData"

###############################################################################
# The rest should work as it is.

# Name columns and add IDs to all TADs
colnames(tad_data) <- c("chrom", "start", "end")
tad_data <- cbind(TADid = as.numeric(rownames(tad_data)), tad_data)

# Some logging
cat("tad data:\n")
str(tad_data)
cat("Any NA values:")
print(any(is.na(tad_data)))

cat("\nManifest(s) used:", use_manifest, "\n")

if (use_manifest == "both") {
  cat("\nI450K:\n")
  str(df_450)
  cat("Any NA values:")
  print(any(is.na(df_450)))
  cat("\nEPIC:\n")
  str(df_EPIC)
  cat("Any NA values:")
  print(any(is.na(df_EPIC)))
  
} else if (use_manifest == "df_450") {
  cat("\nI450K:\n")
  str(df_450)
  cat("Any NA values:")
  print(any(is.na(df_450)))
  
} else if (use_manifest == "df_EPIC") {  
  cat("\nEPIC:\n")
  str(df_EPIC)
  cat("Any NA values:")
  print(any(is.na(df_EPIC)))
  
} else {
  print("Error: Invalid value for **use_manifest**")
}

# Create empty df for TAD map
cpg_tad_map <- data.frame(cpg = character(), 
                          TADid = numeric(), 
                          stringsAsFactors = FALSE)

# Keep track of the number of TADs that no CpGs are mapped to  
empty_TADS<- 0

# Loop through each chromosome
for (chrom in 1:22) {
  
  cat("\n########## PROCESSING CHROMOSOME", chrom, "##########\n\n")
  
  # Filter tad_data for the current chromosome
  tad_chrom <- tad_data[tad_data$chr == paste0("chr", chrom), ]
  
  # logging
  cat("number of TADs in chromosome", chrom, ":", nrow(tad_chrom), "\n")
  #cat("ranging from TADid", tad_chrom$TADid[1],"to", tad_chrom$TADid[nrow(tad_chrom)])
  print("failsafe: all TAD entries labeled:")
  print(unique(tad_chrom$chrom))
  
  if (use_manifest == "both") {
    # Combines both EPIC1 and I450K
    
    # Filter df_450 for the current chromosome
    print("parsing out chromosome from I450K...")
    df_450_chrom <- df_450[df_450$CHR == chrom, ]
    #print("df:")
    #print(str(df_450_chrom))
    
    print("failsafe: all 450 entries labeled:")
    print(unique(df_450_chrom$CHR))
    
    # Filter df_EPIC for the current chromosome
    print("parsing out chromosome from EPIC...")
    df_EPIC_chrom <- df_EPIC[df_EPIC$CHR == chrom, ]
    #print("df:")
    #print(str(df_EPIC_chrom))
    
    print("failsafe: all EPIC entries labeled:")
    print(unique(df_EPIC_chrom$CHR))
    
    print("combining 450 and EPIC CpGs...")
    common_row_names <- intersect(rownames(df_450_chrom), rownames(df_EPIC_chrom))
    #print("common row names:")
    #print(str(common_row_names))
    
    df_450_filtered <- df_450_chrom[!rownames(df_450_chrom) %in% common_row_names, ]
    
    #print("entries from 450 not appearing in EPIC:")
    #print(str(df_450_filtered))
    
    combined_df <- rbind(df_450_filtered, df_EPIC_chrom)
    #print("combined_df:")
    #print(str(combined_df))
    
    cat("\nSUMMARY FOR CHROMOSOME", chrom, ":\n")
    print("nr of cpgs in I450K:")
    print(nrow(df_450_chrom))
    print("nr of cpgs in EPIC:")
    print(nrow(df_EPIC_chrom))
    print("intersect:")
    print(length(common_row_names))
    print("cpgs unique to I450K:")
    print(nrow(df_450_filtered))
    print("unique I450K cpgs combined with EPIC:")
    print(nrow(combined_df))
    print("should be the same as:")
    print(nrow(df_EPIC_chrom) + nrow(df_450_filtered))
    
  } else if (use_manifest == "df_450") {
    # Uses only I450K
    
    # Filter df_450 for the current chromosome
    print("parsing out chromosome from I450K...")
    df_450_chrom <- df_450[df_450$CHR == chrom, ]
    #print("df:")
    #print(str(df_450_chrom))
    
    print("failsafe: all 450 entries labeled:")
    print(unique(df_450_chrom$CHR))
    
    cat("\nSUMMARY FOR CHROMOSOME", chrom, ":")
    print("nr of cpgs in I450K:")
    print(nrow(df_450_chrom))
    
  } else if (use_manifest == "df_EPIC") {
    # Uses only EPIC
    
    # Filter df_EPIC for the current chromosome
    print("parsing out chromosome from EPIC...")
    df_EPIC_chrom <- df_EPIC[df_EPIC$CHR == chrom, ]
    #print("df:")
    #print(str(df_EPIC_chrom))
    
    print("failsafe: all EPIC entries labeled:")
    print(unique(df_EPIC_chrom$CHR))
    
    cat("\nSUMMARY FOR CHROMOSOME", chrom, ":")
    print("nr of cpgs in EPIC:")
    print(nrow(df_EPIC_chrom))
    
  } else {
    print("Error: Invalid value for **use_manifest**")
  }
  
  cat("\n")
  print("mapping TADs...")
  cat("\n")
  # Loop through each row in tad_chrom
  for (i in 1:nrow(tad_chrom)) {
    
    # Create range from start to end
    range_start <- tad_chrom$start[i]
    range_end <- tad_chrom$end[i]
    
    # Filter manifest(s) for MAPINFO within the range
    cpgs_within_range <- combined_df[combined_df$MAPINFO >= range_start & combined_df$MAPINFO <= range_end, ]
    #str(cpgs_within_range)
    
    if (nrow(cpgs_within_range) == 0) {
      cat("NOTE: ", tad_chrom$chrom[i], 
          tad_chrom$start[i], tad_chrom$end[i], 
          "matched to zero CpGs.\n\n")
      empty_TADS<- empty_TADS + 1
      next  # Move to the next iteration of the for loop
    }
    
    # Add the matching rownames to cpg_tad_map
    cpg_tad_map <- rbind(cpg_tad_map, data.frame(cpg = rownames(cpgs_within_range), TADid = tad_chrom$TADid[i]))
  }
  
  cat("mapping TADs for chromosome", chrom, "done.\n")
  print("cpgs in tad map so far:")
  print(nrow(cpg_tad_map))
}

cat("\n")
print("all chromosomes processed")
print("number of empty TADs:")
print(empty_TADS)

print("saving the tad map...")

save(cpg_tad_map, file = TAD_map_file_name)

# A string I use for logging to quickly 
# check that the entire script executed.
print("thank_you_and_goodnight")


