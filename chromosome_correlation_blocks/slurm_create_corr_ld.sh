#!/bin/bash

#SBATCH --account=p33_tsd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --time=48:00:00

#SBATCH --output=/tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step3_create_correaltion_matrix/chromosome_clusters/slurm-%j.out

set -o errexit

module purge 
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

## Copy input files to the work directory:
cp /tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step1_create_df_beta/residuals_subset_ordered.RData $SCRATCH
cp /tsd/p33/data/durable/characters/kristofs/locus_info_tool/chr_clusters_EPIC.RData $SCRATCH
cp /tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step3_create_correaltion_matrix/chromosome_clusters/create_corr_and_ld.r $SCRATCH

## Make sure the results are copied back to the submit directory:
chkfile corr_and_ld.RData corr-data 

## Run Rscript:
cd $SCRATCH

Rscript create_corr_and_ld.r

