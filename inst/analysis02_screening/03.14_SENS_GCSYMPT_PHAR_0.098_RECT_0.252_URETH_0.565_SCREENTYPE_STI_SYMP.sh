sbatch -J SENS_03.14_GCSYMPT_PHAR_0.098_RECT_0.252_URETH_0.565_SCREENTYPE_STI_SYMP -o SENS_03.14_ARRAY-%A_JOB-%J_SIMNO-%4a.log --export=ALL,SIMDIR=~/scratch/SENS_03.14_GCSYMPT_PHAR_0.098_RECT_0.252_URETH_0.565_SCREEN_STI_SYMP,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.0983004049458697,RECT_GC_SYMPT_PROB_ALTPARAM=0.251940310516416,URETH_GC_SYMPT_PROB_ALTPARAM=0.564744788061711 02.02_sti_symp.sh