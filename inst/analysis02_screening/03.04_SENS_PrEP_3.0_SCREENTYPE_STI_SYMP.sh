<<<<<<< HEAD
sbatch -J SENS_03.04_PrEP_3.0_SCREENTYPE_STI_SYMP -o SENS_03.04_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.04_PrEP_3.0_SCREEN_STI_SYMP,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=3 02.02_sti_symp.sh
=======
sbatch -J SENS_03.04_PrEP_3.0_SCREENTYPE_STI_SYMP -o SENS_03.04_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00  --export=ALL,SIMDIR=~/scratch/SENS_03.04_PrEP_3.0_SCREEN_STI_SYMP,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=3 02.02_sti_symp.sh
>>>>>>> 86312dc4f8f46245f236e655d929d0e96957145f
