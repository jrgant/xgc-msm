sbatch -J SENS_03.03_PrEP_2.5_SCREENTYPE_STI_CDC2 -o SENS_03.03_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.03_PrEP_2.5_SCREEN_STI_CDC2,STI_SCREEN_TYPE=cdc,STI_SCREEN_KISS_EXPOSURE=TRUE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=2.5 02.04_sti_cdc2.sh
