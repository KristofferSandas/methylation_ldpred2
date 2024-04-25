# DESCRIPTION: This script transforms the training data into the format needed
# for LDpred. 

# INPUT: training data, test data, intersecting CpGs, ld, case/control info for 
# training data

# OUTPUT: df_beta

# EXECUTE: from this script

library(bigsnpr)
library(yaml)
library(dplyr)

packageVersion("bigsnpr")
# should be >= v1.11.4

###############
##### LOAD DATA

# Import paths from config.
config<- yaml::read_yaml("config.yaml")

# Load training data
meta_ewas<- read.csv(config$training_data_path)
# Inspect
str(meta_ewas)
  
# Load test data
load(config$test_data_path)
# Inspect
str(residuals_all_chr)

# Load ld
ld_path<- file.path(config$analysis_path, "ld.RData")
load(ld_path)
# Inspect
str(ld)

# Load intersect CpGs
intersect_path<- file.path(config$analysis_path, "intersect_cpgs.RData")
load(intersect_path)
# Inspect
str(intersect_cpgs)

# Path for saving the results later
df_beta_save_path<- file.path(config$analysis_path, "df_beta.RData")

##########################
##### DATA PREP AND QC

# Examine if intersect_cpgs and ld contain the same CpGs
# Everyting needs to be ordered according to names(ld) from now on.
# We will order things just before they are used to avoid confusion
# about what order things are in.
ordered_intersect_cpgs<- intersect_cpgs[match(names(ld), intersect_cpgs)]
str(ordered_intersect_cpgs)
str(names(ld))
all(ordered_intersect_cpgs == names(ld))
# If they are not the same something may have gone wrong. But it might still
# be ok depending on the clusters used.

# create a subset of columns we need from the meta-ewas
# LDpred authors recommend using fixed effects.
cols_to_keep <- c("ProbeID", 
                  "All_Effect_Random", 
                  "All_Effect_SE_Random",
                  "All_P_Random")
meta_ewas<- meta_ewas[, cols_to_keep, drop = FALSE]

# Filter by intersecting CpGs
meta_ewas <- meta_ewas %>% filter(ProbeID %in% intersect_cpgs)

# inspect
str(meta_ewas)

# calculate effective sample size for df_beta
# According to the formula
# n_eff = 4 / (1 / n_control + 1 / n_case)
# for each study in the meta analysis.
# Then sum the n_eff from all studies in meta-ewas
# cases/controls in example meta-ewas:
# UCL:    353∕322
# ABR:    414∕433
# IoPPN:  290/203
n_eff_UCL<- 4 / ((1/322) + (1/353))
n_eff_ABR<- 4 / ((1/433) + (1/414))
n_eff_IoPPN<- 4 / ((1/203) + (1/290))
(n_eff_value = round(n_eff_UCL + n_eff_ABR + n_eff_IoPPN))

# Filter test data by intersecting CpGs
test_data<- residuals_all_chr[intersect_cpgs, ]
# Inspect
str(test_data)

#################################################
###### making sumstats compatible with snp_match()

# Placeholder value for alleles "a1" and "a0"
a1<- "B"
a0<- "D"

# the number of rows in the training data
len1<- dim(meta_ewas)[1]

# reformat the training data
sumstats<- data.frame(
  rsid = meta_ewas$ProbeID,
  chr = as.integer(rep(0, len1)),
  pos = as.integer(rep(0, len1)),
  a1 = rep(a1, len1),
  a0 = rep(a0, len1),
  beta = meta_ewas$All_Effect_Random,
  beta_se = meta_ewas$All_Effect_SE_Random,
  N = as.integer(rep(n_eff_value, len1)),
  p = meta_ewas$All_P_Random,
  n_eff = as.integer(rep(n_eff_value, len1)) 
)

# inspect
str(sumstats)

# Creating the "map" variable for snp_match()

# number of rows in test data
# (len1 == len2) should be true
len2<- dim(test_data)[1]

# Create the "map" dataframe
map <- data.frame(
  chr = as.integer(rep(0, len2)),
  rsid = rownames(test_data),
  pos = as.integer(rep(0, len2)),
  a1 = rep(a1, len2),
  a0 = rep(a0, len2)
)

# Inspect
str(map)

####################
##### create df_beta 
df_beta <- snp_match(sumstats, map, 
                     strand_flip = FALSE, 
                     join_by_pos = FALSE,
                     remove_dups = FALSE,
                     return_flip_and_rev = TRUE
)

# The message displayed should say: 
# "X variants to be matched. X variants have been matched."
# I.e. all variants should be matched. This function is 
# from LDpred and thus talks about flipped and reversed
# variants which does not apply in methylation so it can
# be ignored.

# Inspect
str(df_beta)

# check for NA values
any(is.na(df_beta))

# Save df_beta
save(df_beta, file = df_beta_save_path)

# Clean up
rm(a0, a1, cols_to_keep, config, df_beta, df_beta_save_path, 
   intersect_cpgs, intersect_path, ld, ld_path, len1, len2,
   map, meta_ewas, n_eff_ABR, n_eff_IoPPN, n_eff_UCL, n_eff_value,
   ordered_intersect_cpgs, residuals_all_chr, sumstats, test_data)
