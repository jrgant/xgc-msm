#!/bin/bash
#SBATCH -J Sim3-CleanLHS
#SBATCH--time=1:00:00
#SBATCH -p batch
#SBATCH--mem=150GB
#SBATCH -n 208
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu

module load R/4.0.3

cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/cal/sim3_01_clean.r --vanilla
