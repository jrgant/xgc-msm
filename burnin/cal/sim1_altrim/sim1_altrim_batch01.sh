#!/bin/bash
#SBATCH -J Sim1-LHS-XGC
#SBATCH --time=2:00:00
#SBATCH -p batch
#SBATCH --mem=3GB
#SBATCH -n 1
#SBATCH --array=1-1000
#SBATCH -o LHS-Sim1_ARRAY-%A_JOB-%J_SIMNO-%4a.log
#SBATCH --export=ALL,NSIMS=1,NSTEPS=3120,ARRIVE_RATE_ADD_PER20K=1.285,SIMDIR=~/scratch/sim1_altrim
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu
module load R/4.0.3
cd ~/data/jgantenb/xgcmsm/
Rscript ./burnin/cal/sim1_02_altrim_lhs.r --vanilla
