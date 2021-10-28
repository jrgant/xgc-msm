#!/bin/bash
#SBATCH -J ScreenEpi-STI_BURNIN
#SBATCH --time=3:00:00
#SBATCH -p batch
#SBATCH --mem=3GB
#SBATCH -n 1
#SBATCH --array=1-1000
#SBATCH -o ScreenEpi-STI_BURNIN_ARRAY-%A_JOB-%J_SIMNO-%4a.log
#SBATCH --export=ALL,NSIMS=1,NSTEPS=3120,ARRIVE_RATE_ADD_PER20K=1.285,EPI_RUN_TYPE=burnin,STI_SCREEN_TYPE=base,STI_SCREEN_KISS_EXPOSURE=FALSE,STARTTIME=1,SIMDIR=~/scratch/ScreenEpi-STI_BURNIN
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu
module load R/4.0.3
cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/analysis02_screening/02.00_epi.r --vanilla
