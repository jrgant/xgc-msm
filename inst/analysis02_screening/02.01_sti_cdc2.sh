#!/bin/bash
#SBATCH -J ScreenEpi-STI_CDC2
#SBATCH --time=3:00:00
#SBATCH -p batch
#SBATCH --mem=3GB
#SBATCH -n 1
#SBATCH --array=1-1000
#SBATCH -o ScreenEpi-STI_CDC2_ARRAY-%A_JOB-%J_SIMNO-%4a.log
#SBATCH --export=ALL,NSIMS=1,NSTEPS=260,ARRIVE_RATE_ADD_PER20K=1.285,EPI_RUN_TYPE=scenario,STI_SCREEN_TYPE=cdc,STI_SCREEN_KISS_EXPOSURE=TRUE,STARTTIME=3121,SIMDIR=~/scratch/ScreenEpi-STI_CDC2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu
module load R/4.0.3
cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/analysis02_screening/02.00_epi.r --vanilla
