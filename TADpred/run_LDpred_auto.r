# DESCRIPTION: Runs LDpred-auto using all the data created so far.
# Here is where the posterior weights are produced in a vector called 
# beta_auto. The raw scores are also saved in multi_auto for reference.
# QC plots on the model are compiled in a pdf file.

# INPUT: df_beta, ld, h2, correlation matrix

# OUTPUT: beta_auto, multi_auto, qc pdf

# EXECUTE: from slurm script slurm_run_LDpred_auto.sh

##############
##### SETTINGS

# Use only file names, the paths need to be added in the 
# slurm script. 

# Load df_beta
load("df_beta.RData")

# Load correlation matrix
load("correlation_matrix.RData")

# Load ld
load("ld.RData")

# Load h2
load("h2.RData")

# Save paths
multi_auto_save_path<- ("multi_auto.RData")
beta_auto_save_path<- ("beta_auto.RData")
qc_pdf_save_path<- ("model_qc.pdf")

################################################################
# The rest should work as it is

library(bigsnpr)
library(ggplot2)
library(gridExtra)

# install packages if not installed
if (!requireNamespace("bigsnpr", quietly = TRUE)) {
  install.packages("bigsnpr")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}

# version check 
print("bigsnpr version:")
packageVersion("bigsnpr")

# Logging
print("correlation matrix:")
str(corr)

print("h2:")
(h2_est <- ldsc[["h2"]])

print('ld:')
str(ld)

print("df_beta:")
str(df_beta)

print("ordering df_beta according to ld...")
ids<- df_beta$rsid
index <- match(ids, names(ld))
print("any NA values in index:")
any(is.na(index)) # just to be sure

# Order df_beta by the index vector
df_beta_ordered <- df_beta[order(index), ]

print("all CpGs in df_beta and ld are in the same location:")
all(df_beta_ordered$rsid == names(ld))
print("first entries:")
print("df_beta:")
head(df_beta_ordered$rsid)
print("ld:")
head(names(ld))

print("running LDpred...")
coef_shrink <- 0.95
# run LDpred 
multi_auto <- snp_ldpred2_auto(
  corr, df_beta_ordered, h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30), ncores = 1,
  # use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
  allow_jump_sign = FALSE, shrink_corr = coef_shrink)

print("run completed.")
print("multi_auto:")
str(multi_auto, max.level = 1)

print("multi_auto[1]:")
str(multi_auto[[1]], max.level = 1)

print("creating QC plots...")

pdf(qc_pdf_save_path) # open print to pdf

# Loop through the items
for (i in 1:30) {
  auto <- multi_auto[[i]]
  
  # Create the first plot
  plot1 <- qplot(y = auto$path_p_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p")
  
  # Create the second plot
  plot2 <- qplot(y = auto$path_h2_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2")
  
  # Combine plots
  combined_plot <- plot_grid(plot1, plot2, ncol = 1, align = "hv")
  
  # Print the combined plot
  print(combined_plot)
}
dev.off() # close print to pdf
print("pdf creation complete.")

print("filtering markov chains by quality...")
range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
print("range (all should be between 0-2):")
print(range)

print("which chains to keep (out of 30):")
(keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))

beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
print("beta_auto:")
str(beta_auto)

# Saving results
print("saving results...")
save(multi_auto, file = multi_auto_save_path)
save(beta_auto, file = beta_auto_save_path)
print("results saved.")

print("thank_you_and_goodnight")

