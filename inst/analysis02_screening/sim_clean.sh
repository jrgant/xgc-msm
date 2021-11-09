#!/bin/bash
#SBATCH -J CleanEpi
#SBATCH--time=2:00:00
#SBATCH -p batch
#SBATCH--mem=150GB
#SBATCH -n 64
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu

module load R/4.0.3

cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/analysis02_screening/sim_clean.r --vanilla
