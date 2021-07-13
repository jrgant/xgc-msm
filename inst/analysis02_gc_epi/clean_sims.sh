#!/bin/bash
#SBATCH -J CleanEpiSim
#SBATCH--time=1:00:00
#SBATCH -p batch
#SBATCH--mem=150GB
#SBATCH -n 208
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu

module load R/4.0.3

cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/analysis02_gc_epi/clean_sims.r --vanilla
