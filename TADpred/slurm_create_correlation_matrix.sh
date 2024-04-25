#!/bin/bash

#SBATCH --account=p33_tsd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=01:00:00

#SBATCH --output=slurm-%j.out

## IMPORTANT: Execute this script *from your analysis directory* to make
## sure all files are in the correct place: $ sbatch /path/to/Rproject/dir/TADpred/slurm_create_correlation_matrix.sh

## The main path where your TADpred.Rproj file is
TADpred_path="/tsd/p33/data/durable/characters/kristofs/TADpred"

## The main path where your analysis is
analysis_path="/tsd/p33/data/durable/characters/kristofs/TADpred_TOP_SCZ"

set -o errexit

module purge 
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

## Copies input files to the work directory:

## Test data
## In an alternate reality this path could be imported from config. No such luck.
cp /tsd/p33/data/durable/characters/kristofs/For_Kristoffer/Residuals/all_samples/residuals_betas_allMarkosCovariables-withSex_wholeSample.Robj $SCRATCH

## TAD map
cp "$analysis_path/tad_clusters.RData" $SCRATCH

## Rscript
cp "$TADpred_path/create_correlation_matrix.r" $SCRATCH

## Makes sure the results are copied back to the submit directory.
## These file names must correspond to the save file names in the rscript.
chkfile correlation_matrix.RData ld.RData corr-data 

## Executes the rscript:
cd $SCRATCH
Rscript create_correlation_matrix.r

