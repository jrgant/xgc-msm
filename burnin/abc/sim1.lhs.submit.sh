#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu

module load R/4.0.3
cd ~/data/jgantenb/xgcmsm/
Rscript ./burnin/abc/sim1.lhs.r --vanilla
