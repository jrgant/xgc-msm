#!/bin/bash
#SBATCH -J Main-Episim
#SBATCH --time=3:00:00
#SBATCH -p batch
#SBATCH --mem=3GB
#SBATCH -n 1
#SBATCH --array=1-1000
#SBATCH -o MainSim_ARRAY-%A_JOB-%J_SIMNO-%4a.log
#SBATCH --export=ALL,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,SIMDIR=~/scratch/MainEpi
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu
module load R/4.0.3
cd ~/data/jgantenb/xgcmsm/
Rscript ./inst/analysis01_epi/02_main.r --vanilla
