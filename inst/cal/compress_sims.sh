#!/bin/bash
#SBATCH -J Compress-Sims
#SBATCH --time=3:00:00
#SBATCH -p batch
#SBATCH --mem=50GB
#SBATCH -n 64
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jrgant@brown.edu

XZ_OPT=-9T0 tar Jcf $COMPDIR.tar.xz $COMPDIR
