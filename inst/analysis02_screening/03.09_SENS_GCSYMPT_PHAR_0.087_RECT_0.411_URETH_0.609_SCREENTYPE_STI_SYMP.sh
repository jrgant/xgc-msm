sbatch -J SENS_03.09_GCSYMPT_PHAR_0.087_RECT_0.411_URETH_0.609_SCREENTYPE_STI_SYMP -o SENS_03.09_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.09_GCSYMPT_PHAR_0.087_RECT_0.411_URETH_0.609_SCREEN_STI_SYMP,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.0866430887169464,RECT_GC_SYMPT_PROB_ALTPARAM=0.410861819430695,URETH_GC_SYMPT_PROB_ALTPARAM=0.608705628106292 02.02_sti_symp.sh
