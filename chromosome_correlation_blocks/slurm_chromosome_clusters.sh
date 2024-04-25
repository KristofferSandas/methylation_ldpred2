#!/bin/bash

#SBATCH --account=p33_tsd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=72:00:00

#SBATCH --output=/tsd/p33/data/durable/characters/kristofs/locus_info_tool/slurm-%j.out

set -o errexit

module purge 
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

## Copy input files to the work directory:
cp /tsd/p33/data/durable/characters/kristofs/For_Kristoffer/comeback_package/comeback/R/sysdata.rda $SCRATCH
cp /tsd/p33/data/durable/characters/kristofs/LDpred_pipeline_1/Step1_create_df_beta/residuals_subset_ordered.RData $SCRATCH
cp /tsd/p33/data/durable/characters/kristofs/locus_info_tool/create_chromosome_clusters.r $SCRATCH

## Make sure the results are copied back to the submit directory:
chkfile chr_clusters_450.RData chr_clusters_EPIC.RData 

## Run Rscript:
cd $SCRATCH

Rscript create_chromosome_clusters.r
