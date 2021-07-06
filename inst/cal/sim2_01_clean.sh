#!/bin/bash
#SBATCH -J Sim2-CleanLHS
#SBATCH--time=1:00:00
#SBATCH -p batch
#SBATCH--mem=150GB
#SBATCH -n 208
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu

module load R/4.0.3

cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/cal/sim2_01_clean.r --vanilla
