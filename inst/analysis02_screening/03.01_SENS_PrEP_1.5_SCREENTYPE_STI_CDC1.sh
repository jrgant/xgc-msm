sbatch -J SENS_03.01_PrEP_1.5_SCREENTYPE_STI_CDC1 -o SENS_03.01_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.01_PrEP_1.5_SCREEN_STI_CDC1,STI_SCREEN_TYPE=cdc,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=1.5 02.03_sti_cdc1.sh
