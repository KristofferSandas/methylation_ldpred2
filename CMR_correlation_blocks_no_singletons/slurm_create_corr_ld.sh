#!/bin/bash

#SBATCH --account=p33_tsd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=5:00:00

#SBATCH --output=/tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step3_create_correaltion_matrix/corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/no_cometh_singletons/slurm-%j.out

set -o errexit

module purge 
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

## Copy input files to the work directory:
cp /tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step3_create_correaltion_matrix/corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/no_cometh_singletons/non_singleton_residuals.RData $SCRATCH
cp /tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step3_create_correaltion_matrix/corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/no_cometh_singletons/non_singletons.RData $SCRATCH
cp /tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step3_create_correaltion_matrix/corlo_02_maxprbdst_100k_corlodst=800_Iarray_EPIC/no_cometh_singletons/create_corr_and_ld.r $SCRATCH

## Make sure the results are copied back to the submit directory:
chkfile corr_and_ld.RData corr-data 

## Run Rscript:
cd $SCRATCH

Rscript create_corr_and_ld.r

