
# load test data
residuals_path<- paste("/tsd/p33/data/durable/characters/kristofs/",
                       "LDpred_pipeline_1/Step1_create_df_beta/",
                       "residuals_subset_ordered.RData", sep="")

load(residuals_path)

str(residuals_subset_ordered)

# load CMRs from CoMeBack
cometh_regions_path<- paste("/tsd/p33/data/durable/characters/kristofs/",
                            "LDpred_pipeline_1/Step2_run_comeback/",
                            "corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/",
                            "co_meth_regions.RData", sep="")

load(cometh_regions_path)

str(co_meth_regions)
# 67921

# Find where the singetons begin
# using logging info from CoMeBack run
co_meth_regions[67911:67931,]
co_meth_regions[67921,]

non_singletons<- co_meth_regions[1:67921,]
head(non_singletons)
tail(non_singletons)
str(non_singletons)

# save all non-singleton CMRs
non_singletons_save_path<- paste("/tsd/p33/data/durable/characters/kristofs/",
                                 "LDpred_pipeline_1/Step3_create_correaltion_matrix/",
                                 "corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/",
                                 "no_cometh_singletons/non_singletons.RData", sep="")

save(non_singletons, file = non_singletons_save_path)

# create the corresponding subset of the test data
intersect_ids <- intersect(rownames(residuals_subset_ordered), non_singletons$clustercpg)
str(intersect_ids)

non_singleton_residuals <- residuals_subset_ordered[intersect_ids, ]
str(non_singleton_residuals)
any(is.na(non_singleton_residuals)) # F

non_singleton_residuals_save_path<- paste("/tsd/p33/data/durable/characters/kristofs/",
                                 "LDpred_pipeline_1/Step3_create_correaltion_matrix/",
                                 "corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/",
                                 "no_cometh_singletons/non_singleton_residuals.RData", sep="")

save(non_singleton_residuals, file = non_singleton_residuals_save_path)

