sbatch -J SENS_03.11_GCSYMPT_PHAR_0.057_RECT_0.295_URETH_0.865_SCREENTYPE_STI_CDC2 -o SENS_03.11_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.11_GCSYMPT_PHAR_0.057_RECT_0.295_URETH_0.865_SCREEN_STI_CDC2,STI_SCREEN_TYPE=cdc,STI_SCREEN_KISS_EXPOSURE=TRUE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.0570771782122352,RECT_GC_SYMPT_PROB_ALTPARAM=0.294725463776314,URETH_GC_SYMPT_PROB_ALTPARAM=0.865360672013233 02.04_sti_cdc2.sh
