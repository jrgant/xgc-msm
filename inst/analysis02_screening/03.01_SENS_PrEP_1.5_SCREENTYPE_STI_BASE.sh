sbatch -J SENS_03.01_PrEP_1.5_SCREENTYPE_STI_BASE -o SENS_03.01_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.01_PrEP_1.5_SCREEN_STI_BASE,STI_SCREEN_TYPE=base,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=1.5 02.01_sti_base.sh
