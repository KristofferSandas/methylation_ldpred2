#!/bin/bash

#SBATCH --account=p33_tsd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=4:00:00

#SBATCH --output=/tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step2_run_comeback/corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/slurm-%j.out

set -o errexit

module purge 
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

## Copy input files to the work directory:
cp /tsd/p33/data/durable/characters/kristofs/For_Kristoffer/comeback_package/comeback_1.0.1.tar $SCRATCH
cp /tsd/p33/data/durable/characters/kristofs/For_Kristoffer/Rscripts/MRS_func.R $SCRATCH
cp /tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step1_create_df_beta/residuals_subset_ordered.RData $SCRATCH
cp /tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step2_run_comeback/corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/run_comeback.r $SCRATCH

## Make sure the results are copied back to the submit directory:
chkfile cmrs.RData co_meth_regions.RData 

## Run Rscript:
cd $SCRATCH

Rscript run_comeback.r
