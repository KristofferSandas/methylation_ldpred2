# This script finds the CpGs with the highest prior and posterior effect sizes.
# It maps the CpGs to the genes they are associated with.
# Finally, it identifies the genes present in the posterior data set which 
# are not present in the prior. These are the genes associated with the 
# CpGs that the posterior model found that the prior did not.

setwd("/home/monkey/bergen_project")

# get top cpgs with highest absolute prior effect sizes
prior_data <- read.csv("SCZ_notop_sex_stratified_meta_analysis.csv")
prior_data <- prior_data[, c("ProbeID", "All_Effect_Random")]
prior_data <- prior_data[order(abs(prior_data$All_Effect_Random), decreasing = TRUE), ]
str(prior_data)
prior_top_cpgs <- head(prior_data, 10000)
str(prior_top_cpgs)

# get top cpgs with highest absolute posterior effect sizes
load('ld.RData') 
str(ld)
load('beta_auto.RData')
str(beta_auto)
post_data <- data.frame(ProbeID = names(ld), beta_auto = beta_auto)
str(post_data)
post_data <- post_data[order(abs(post_data$beta_auto), decreasing = TRUE), ]
post_top_cpgs <- head(post_data, 10000)
str(post_top_cpgs)

# save as csv
write.csv(prior_top_cpgs, file = 'prior_top_cpgs.csv', row.names = FALSE)
write.csv(post_top_cpgs, file = 'post_top_cpgs.csv', row.names = FALSE)

prior_top_cpgs<- read.csv("prior_top_cpgs.csv")
post_top_cpgs<- read.csv("post_top_cpgs.csv")

# Get illumina manifest
manifest<- read.csv("humanmethylation450_15017482_v1-2(1).csv", skip = 7, header = TRUE)
str(manifest)
head(manifest,5)

manifest<- manifest[, c("IlmnID", "UCSC_RefGene_Name")]
str(manifest)

# map prior CpGs to genes:
matching_indices <- match(prior_top_cpgs$ProbeID, manifest$IlmnID, nomatch = 0)
matching_UCSC_RefGene_Names <- manifest$UCSC_RefGene_Name[matching_indices]
prior_top_cpgs$UCSC_RefGene_Names <- matching_UCSC_RefGene_Names
# remove doubles in name col
prior_top_cpgs$UCSC_RefGene_Names <- gsub(";.*$", "", prior_top_cpgs$UCSC_RefGene_Names)

str(prior_top_cpgs)

# map posterior CpGs to genes:
matching_indices <- match(post_top_cpgs$ProbeID, manifest$IlmnID, nomatch = 0)
matching_UCSC_RefGene_Names <- manifest$UCSC_RefGene_Name[matching_indices]
post_top_cpgs$UCSC_RefGene_Names <- matching_UCSC_RefGene_Names
# remove doubles in name col
post_top_cpgs$UCSC_RefGene_Names <- gsub(";.*$", "", post_top_cpgs$UCSC_RefGene_Names)

str(post_top_cpgs)

# check for na vals
(any(is.na(prior_top_cpgs)))
(any(is.na(post_top_cpgs)))

# remove rows with empty cols
prior_top_cpgs <- prior_top_cpgs[prior_top_cpgs$UCSC_RefGene_Names != "", ]
str(prior_top_cpgs)
post_top_cpgs <- post_top_cpgs[post_top_cpgs$UCSC_RefGene_Names != "", ]
str(post_top_cpgs)

# Keep the top 6500 genes from both priors and posteriors
prior_top_cpgs <- head(prior_top_cpgs, 6500)
post_top_cpgs <- head(post_top_cpgs, 6500)

# save 
write.csv(prior_top_cpgs, file = 'prior_top_cpgs.csv', row.names = FALSE)
write.csv(post_top_cpgs, file = 'post_top_cpgs.csv', row.names = FALSE)

# find the genes present in the posterior list that are not present in the prior list = "new genes"
in_post_but_not_in_prior <- post_top_cpgs$UCSC_RefGene_Names[!post_top_cpgs$UCSC_RefGene_Names %in% prior_top_cpgs$UCSC_RefGene_Names]
in_post_but_not_in_prior<- unique(in_post_but_not_in_prior)
str(in_post_but_not_in_prior)

# save
write.csv(in_post_but_not_in_prior, file = 'all_new_in_post.csv', row.names = FALSE)



