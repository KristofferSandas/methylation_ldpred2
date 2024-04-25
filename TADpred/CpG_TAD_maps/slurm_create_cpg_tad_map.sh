#!/bin/bash

#SBATCH --account=p33_tsd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=00:10:00

#SBATCH --output=slurm-%j.out

## 10 min with 20G RAM should be enough.
## It is best to execute this script standing in its directory CpG_TAD_maps.
## That way the system moves the result files and saves the log file in
## the correct directory.

## The main path where your TADpred.Rproj file is
TADpred_path="/tsd/p33/data/durable/characters/kristofs/TADpred"

set -o errexit

module purge 
module load R-bundle-Bioconductor/3.15-foss-2021b-R-4.2.0

## Copies input files to the work directory (SCRATCH).
## Make sure the correct manifests and TAD scaffolds are specified.
cp "$TADpred_path/illumina_manifests/I450K_Manifest.RData" $SCRATCH
cp "$TADpred_path/illumina_manifests/EPIC_Manifest.RData" $SCRATCH
cp "$TADpred_path/TAD_scaffolds/tad_scaffold_hg19.txt" $SCRATCH
cp "$TADpred_path/CpG_TAD_maps/create_cpg_tad_map.r" $SCRATCH

## Makes sure the results are copied back to the submit directory.
## This file name must correspond to the TAD_map_file_name in the rscript create_cpg_tad_map.r
chkfile cpg_tad_map.RData

## Executes the rscript:
cd $SCRATCH
Rscript create_cpg_tad_map.r
