sbatch -J SENS_02.06_GCSYMPT_PHAR_0.092_RECT_0.234_URETH_0.89_SCREENTYPE_STI_SYMP -o SENS_02.06_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00  --export=ALL,SIMDIR=~/scratch/SENS_02.06_GCSYMPT_PHAR_0.092_RECT_0.234_URETH_0.89_SCREEN_STI_SYMP,EPI_RUN_TYPE=scenario,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.0923847475482076,RECT_GC_SYMPT_PROB_ALTPARAM=0.234367460478201,URETH_GC_SYMPT_PROB_ALTPARAM=0.89040481694669 02.01_sti_symp.sh