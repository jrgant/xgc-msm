sbatch -J SENS_02.03_PrEP_2.0_SCREENTYPE_STI_UNIV -o SENS_02.03_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00  --export=ALL,SIMDIR=~/scratch/SENS_02.03_PrEP_2.0_SCREEN_STI_UNIV,STI_SCREEN_TYPE=universal,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PREP_SCALE_ALTPARAM=2 02.01_sti_univ.sh
