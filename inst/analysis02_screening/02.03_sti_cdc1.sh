#!/bin/bash
#SBATCH -J ScreenEpi-STI_CDC1
#SBATCH --time=3:00:00
#SBATCH -p batch
#SBATCH --mem=3GB
#SBATCH -n 1
#SBATCH --array=1-1000
#SBATCH -o ScreenEpi-STI_CDC1_ARRAY-%A_JOB-%J_SIMNO-%4a.log
#SBATCH --export=ALL,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,STI_SCREEN_TYPE=cdc,STI_SCREEN_KISS_EXPOSURE=FALSE,SIMDIR=~/scratch/ScreenEpi-STI_CDC1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu
module load R/4.0.3
cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/analysis02_screening/02.00_main.r --vanilla
