sbatch -J SENS_03.07_GCSYMPT_PHAR_0.105_RECT_0.319_URETH_0.919_SCREENTYPE_STI_UNIV -o SENS_03.07_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.07_GCSYMPT_PHAR_0.105_RECT_0.319_URETH_0.919_SCREEN_STI_UNIV,STI_SCREEN_TYPE=universal,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.105140577072503,RECT_GC_SYMPT_PROB_ALTPARAM=0.318698484167601,URETH_GC_SYMPT_PROB_ALTPARAM=0.919198100673061 02.05_sti_univ.sh
