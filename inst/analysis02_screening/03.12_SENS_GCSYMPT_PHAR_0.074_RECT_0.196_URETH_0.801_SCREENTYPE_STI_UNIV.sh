sbatch -J SENS_03.12_GCSYMPT_PHAR_0.074_RECT_0.196_URETH_0.801_SCREENTYPE_STI_UNIV -o SENS_03.12_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.12_GCSYMPT_PHAR_0.074_RECT_0.196_URETH_0.801_SCREEN_STI_UNIV,STI_SCREEN_TYPE=universal,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.0741294650685293,RECT_GC_SYMPT_PROB_ALTPARAM=0.195718435856286,URETH_GC_SYMPT_PROB_ALTPARAM=0.801203002239751 02.05_sti_univ.sh
