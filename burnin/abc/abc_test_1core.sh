#!/bin/bash
#
#SBATCH -J ABC-calibration
#SBATCH -t 24:00:00
#SBATCH -p batch
#SBATCH --mem 350GB
#SBATCH -n 1

#SBATCH -o ABC-JobID-%J.log
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu
#

echo ""
echo ""
echo ""

module load R/4.0.3
cd ~/data/jgantenb/xgc-msm

export TAIL_LENGTH=52
export NSTEPS=200
export SIMS_BELOW_TOL=100

echo ""
echo ""

echo "==============================================================================="
echo " SETTINGS "
echo "==============================================================================="
echo ""

echo "Tail length:"
echo $TAIL_LENGTH
echo ""

echo "Time steps:"
echo $NSTEPS
echo ""

echo "Simulations below tolerance:"
echo $SIMS_BELOW_TOL
echo ""

echo "=============================================================================="
echo " DEBUG SCRIPT "
echo "=============================================================================="
echo ""

cat ./burnin/abc/sim.prep.debug.r

echo "=============================================================================="
echo " ABC SCRIPT "
echo "=============================================================================="
echo ""

cat ./burnin/abc/sim.prep.r

Rscript ./burnin/abc/sim.prep.debug.r --vanilla
