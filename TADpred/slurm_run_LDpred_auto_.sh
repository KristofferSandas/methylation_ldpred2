#!/bin/bash

#SBATCH --account=p33_tsd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=01:00:00

#SBATCH --output=slurm-%j.out

## IMPORTANT: Execute this script *from your analysis directory* to make
## sure all files are in the correct place: $ sbatch /path/to/Rproject/dir/TADpred/slurm_run_LDpred_auto.sh

## The main path where your TADpred.Rproj file is
TADpred_path="/tsd/p33/data/durable/characters/kristofs/TADpred"

## The main path where your analysis is
analysis_path="/tsd/p33/data/durable/characters/kristofs/TADpred_TOP_SCZ"

set -o errexit

module purge 
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

## Copies input files to the work directory:
cp "$analysis_path/df_beta.RData" $SCRATCH
cp "$analysis_path/correlation_matrix.RData" $SCRATCH
cp "$analysis_path/ld.RData" $SCRATCH
cp "$analysis_path/h2.RData" $SCRATCH
cp "$TADpred_path/run_LDpred_auto.r" $SCRATCH

## Makes sure the results are copied back to the submit directory:
chkfile multi_auto.RData beta_auto.RData model_qc.pdf

## Run the r script:
cd $SCRATCH
Rscript run_LDpred_auto.r

