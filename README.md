# masters_project_2024

The directory *TADpred* contains a prototype pipeline for TADpred, including a tutorial and data files needed.

All analyses were performed following the pipeline in the TADpred directory, with additional extra steps for creating the different matrices in directories *chromosome_correlation_blocks*, *CMR_correlation_blocks* and *CMR_correlation_blocks_no_singletons*.

The *pdf files* are the plots from the Gibbs sampler QC.

*prior_top_cpgs* and *posterior_top_cpgs* contain the genes associated with the CpGs with the highest absolute effect sizes for both the prior and posterior training data. *all_new_in_posterior_scores* are the genes appearing in the posterior list that do not appear in the prior. These files were created with the *bioanalysis.r* script.


