#!/bin/bash
sed -e 's/1-1000/1-100/' -e 's/SIMDIR=sim1/SIMDIR=~\/scratch\/sim1_ArrivePlus_1/' sim1_batch1.sh | sbatch
sed -e 's/1-1000/1-100/' -e 's/20K=1/20K=1.1/' -e 's/SIMDIR=sim1/SIMDIR=~\/scratch\/sim1_ArrivePlus_1.1/' sim1_batch1.sh | sbatch
sed -e 's/1-1000/1-100/' -e 's/20K=1/20K=1.25/' -e 's/SIMDIR=sim1/SIMDIR=~\/scratch\/sim1_ArrivePlus_1.25/' sim1_batch1.sh | sbatch
sed -e 's/1-1000/1-100/' -e 's/20K=1/20K=1.5/' -e 's/SIMDIR=sim1/SIMDIR=~\/scratch\/sim1_ArrivePlus_1.5/' sim1_batch1.sh | sbatch
sed -e 's/1-1000/1-100/' -e 's/20K=1/20K=2/' -e 's/SIMDIR=sim1/SIMDIR=~\/scratch\/sim1_ArrivePlus_2/' sim1_batch1.sh | sbatch
