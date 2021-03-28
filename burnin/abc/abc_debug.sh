#!/bin/bash
#
#SBATCH -o ABC-JobID-%J.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu

echo ""
echo ""
echo ""

module load R/4.0.3
cd ~/data/jgantenb/xgcmsm/burnin/abc

export TAIL_LENGTH=52
export NSTEPS=3120
export SIMS_BELOW_TOL=500


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

echo " DEBUG ---------------------------------------------------------------------- "
cat echo $SIMFILE

echo " ABC PROCEDURE -------------------------------------------------------------- "
cat $SIMFILE

echo ""
echo ""

Rscript $SIMFILE --vanilla
