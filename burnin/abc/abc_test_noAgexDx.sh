#!/bin/bash
#
#SBATCH -J ABC-calibration
#SBATCH -t 24:00:00
#SBATCH -p batch
#SBATCH --mem 200GB
#SBATCH -n 96
#SBATCH -o ABC-JobID-%J.log
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu
#

echo ""
echo ""
echo ""

module load R/4.0.3
cd ~/data/jgantenb/xgc-msm/burnin/abc

export TAIL_LENGTH=52
export NSTEPS=1040
export SIMS_BELOW_TOL=200


echo ""
echo ""

echo "==============================================================================="
echo " CONFIRM SETTINGS "
echo "==============================================================================="
echo ""

echo "The following environmental variables are currently set..."
echo ""

echo "Tail length (time steps):"
echo $TAIL_LENGTH
echo ""

echo "Number of time steps:"
echo $NSTEPS
echo ""

echo "Simulations below tolerance:"
echo $SIMS_BELOW_TOL
echo ""
echo "* For the Lenormand algorithm, the value above is multiplied by n_alpha, which"
echo "is set to the EasyABC package default 0.5."
echo ""
echo ""

echo "=============================================================================="
echo " ABC SCRIPT "
echo "=============================================================================="
echo ""

cat sim.prep.noAgeDx.r

echo ""
echo ""

Rscript sim.prep.noAgeDx.r --vanilla
