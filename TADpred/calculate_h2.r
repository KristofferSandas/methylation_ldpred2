# DESCRIPTION: calculates the h2 value needed for further analysis.
# This is one of the parameters determining the genetic architechture 
# of the phenotype studied.
# h2 in LDpred represents heritability. With methylation data this 
# concept does not apply. The math still works, but it is not uncommon
# to see h2-values of 10 or even >100 so dont be shocked.

# INPUT: df_beta, ld

# OUTPUT: h2 

# EXECUTE: from this script

library(bigsnpr)
library(dplyr)

packageVersion("bigsnpr")
# should be >= v1.11.4

###########
##### FILES

# Import paths from config.
config<- yaml::read_yaml("config.yaml")

# Load ld
ld_path<- file.path(config$analysis_path, "ld.RData")
load(ld_path)
# Inspect
str(ld)

# Load df_beta
df_beta_path<- file.path(config$analysis_path, "df_beta.RData")
load(df_beta_path)
# Inspect
str(df_beta)

# Where to save the h2 later
h2_save_path<- file.path(config$analysis_path, "h2.RData")

###################################
##### ORDER df-beta ACCORDING TO LD

# Generate an index vector based on the order of 
# appearance of rsid values in ld
ids<- df_beta$rsid
index <- match(ids, names(ld))
any(is.na(index)) # just to be sure

# Order df_beta by the index vector
df_beta_ordered <- df_beta[order(index), ]

# double check that they are in the same order
head(df_beta_ordered$rsid)
head(names(ld))
all(df_beta_ordered$rsid == names(ld))

##################
##### CALCULATE h2

ldsc <- with(df_beta_ordered, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                      sample_size = n_eff, blocks = NULL, ncores=1))
# Inspect
ldsc

# Save h2
save(ldsc, file = h2_save_path)

# Clean up variables
rm(config, df_beta, df_beta_path, df_beta_ordered, h2_save_path,
   ids, index, ld, ld_path, ldsc)
