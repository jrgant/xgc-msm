sbatch -J SENS_03.08_GCSYMPT_PHAR_0.042_RECT_0.374_URETH_0.68_SCREENTYPE_STI_CDC2 -o SENS_03.08_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.08_GCSYMPT_PHAR_0.042_RECT_0.374_URETH_0.68_SCREEN_STI_CDC2,STI_SCREEN_TYPE=cdc,STI_SCREEN_KISS_EXPOSURE=TRUE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.0422251377242339,RECT_GC_SYMPT_PROB_ALTPARAM=0.373678522718566,URETH_GC_SYMPT_PROB_ALTPARAM=0.67971773087135 02.04_sti_cdc2.sh
