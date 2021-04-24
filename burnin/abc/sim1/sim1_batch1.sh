#!/bin/bash
#SBATCH -J Sim1-LHS-XGC
#SBATCH--time=3:00:00
#SBATCH -p batch
#SBATCH--mem=200GB
#SBATCH -n 1
#SBATCH--array=1-1000
#SBATCH -o ~/scratch/sim1/LHS-Sim1_ARRAY-%A_JOB-%J_SIMNO-%4a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu
module load R/4.0.3
cd ~/data/jgantenb/xgcmsm/
Rscript ./burnin/abc/sim1_02_lhs.r --vanilla
