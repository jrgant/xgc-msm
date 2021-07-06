#!/bin/bash
#SBATCH -J Main-CleanEpi
#SBATCH--time=2:00:00
#SBATCH -p batch
#SBATCH--mem=150GB
#SBATCH -n 208
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu

module load R/4.0.3

cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/analysis01_epi/03_main_clean.r --vanilla
