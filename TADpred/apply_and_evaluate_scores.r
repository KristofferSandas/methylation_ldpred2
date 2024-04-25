# DESCRIPTION: Apply the weights to the test data and evaluate the 
# scores using phenotype data. 

# INPUT: beta_auto (posterior weights), ld, test data, phenotype data

# OUTPUT: risk scores, results, model

# EXECUTE: from this script

library(bigsnpr)
library(fmsb)

###########
##### FILES

# Import paths from config.
config<- yaml::read_yaml("config.yaml")

# Load beta_auto
beta_auto_path<- file.path(config$analysis_path, "beta_auto.RData")
load(beta_auto_path)
# Inspect
str(beta_auto)

# Load test data 
load(config$test_data_path)
# inspect
str(residuals_all_chr)

# Load ld for ordering
ld_path<- file.path(config$analysis_path, "ld.RData")
load(ld_path)
# Inspect
str(ld)

# Load phenotype data
load(config$pheno_path)
# Inspect
str(pheno)

# Remove NA values from phenotpe data
pheno = na.omit(pheno)
str(pheno)

# Save path for results
results_save_path<- file.path(config$analysis_path, "model_results.RData")

#####################################
##### ORDER TEST DATA ACCORDING TO LD

# We first filter the test data by ld since we already checked
# that ld and inersecting_cpgs are identical in an earlier script.
intersect_ids <- intersect(rownames(residuals_all_chr), names(ld))
str(intersect_ids)

test_data <- residuals_all_chr[intersect_ids, ]
str(test_data)

# Then we order the test data according to ld
order_ids<- names(ld)
test_data_ordered <- test_data[match(order_ids, rownames(test_data)), ]

# Check that the order is the same
all(rownames(test_data_ordered) == names(ld))
head(rownames(test_data_ordered))
head(names(ld))

#################################
##### CALCULATING THE RISK SCORES

# Matrix-vector multiplication
risk_scores <- t(test_data_ordered) %*% beta_auto

# Transform to a vector and name using samples
risk_scores_vector <- as.vector(risk_scores)
names(risk_scores_vector) <- rownames(risk_scores)
str(risk_scores_vector)

# Double check that the samples are in the same order in the 
# risk scores and the phenotype data
all(names(risk_scores_vector) == pheno$Sample_Name)

##############################
##### EVALUATE THE RISK SCORES

# LDpred usually tests partial correlation
(partial_corr<- pcor(risk_scores_vector, factor(pheno$diagnosis), NULL))
# correlation lower_CI upper_CI

# Evaluate using logistic regression
logReg_results = glm(formula = factor(pheno$diagnosis) ~ risk_scores_vector,
                                       family = "binomial")

(logReg_summary<- summary(logReg_results))

(R2<- NagelkerkeR2(logReg_results))

# Save model results
save(risk_scores_vector, logReg_results, logReg_summary, R2, partial_corr,
     file = results_save_path)

# Clean up
rm(beta_auto, beta_auto_path, config, intersect_ids, ld, ld_path,
   logReg_results, logReg_summary, order_ids, partial_corr, pheno,
   R2, residuals_all_chr, results_save_path, risk_scores, risk_scores_vector,
   test_data, test_data_ordered)
