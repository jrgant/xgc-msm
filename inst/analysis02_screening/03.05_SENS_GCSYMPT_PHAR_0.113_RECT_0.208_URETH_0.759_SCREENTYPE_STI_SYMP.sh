sbatch -J SENS_03.05_GCSYMPT_PHAR_0.113_RECT_0.208_URETH_0.759_SCREENTYPE_STI_SYMP -o SENS_03.05_ARRAY-%A_JOB-%J_SIMNO-%4a.log -t 5:00:00 --export=ALL,SIMDIR=~/scratch/SENS_03.05_GCSYMPT_PHAR_0.113_RECT_0.208_URETH_0.759_SCREEN_STI_SYMP,STI_SCREEN_TYPE=symptomatic,STI_SCREEN_KISS_EXPOSURE=FALSE,NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,PHAR_GC_SYMPT_PROB_ALTPARAM=0.1132304593585,RECT_GC_SYMPT_PROB_ALTPARAM=0.207907975351219,URETH_GC_SYMPT_PROB_ALTPARAM=0.758776287365799 02.02_sti_symp.sh
