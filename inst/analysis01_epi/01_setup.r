################################################################################
## SETUP ##
################################################################################

library(pacman)
p_load(EpiModelHIV, rlecuyer)

selected_paramsets <- readRDS(here::here("inst", "cal", "simid_sel.rds"))
an01_path <- here::here("inst", "analysis01_epi")


################################################################################
## SAMPLE PARAMETER SETS ##
################################################################################

nsim_spec <- 1000

set.seed(1971)
paramset_sampled <- sample(selected_paramsets, nsim_spec, replace = TRUE)

saveRDS(
  paramset_sampled,
  here::here(an01_path, "paramset_sampled.rds")
)


################################################################################
## SETUP RANDOM NUMBER STREAMS ##
################################################################################

streams <- paste0("str", sprintf("%04d", 1:nsim_spec))
.lec.CreateStream(streams)

saveRDS(
  .lec.Random.seed.table,
  here::here(an01_path, "seeds.rds")
)


################################################################################
## MAKE BATCH SCRIPT ##
################################################################################

# This function writes a batch script to submit a job array.
make_batch_script <- function(jobname, walltime, partition, mem,
                              ncores, array, log_fullpath, batchid,
                              nsims, nsteps, add_arrivals, simdir) {

  sb <- "#SBATCH"

  specs <- paste(
    "#!/bin/bash",
    paste(  sb, "-J", jobname       ),
    paste0( sb, " --time=", walltime ),
    paste(  sb, "-p", partition     ),
    paste0( sb, " --mem=", mem       ),
    paste(  sb, "-n", ncores        ),
    paste0( sb, " --array=", array   ),
    paste(  sb, "-o", log_fullpath  ),
    paste0(
      sb,
      " --export=ALL,NSIMS=", nsims,
      ",NSTEPS=", nsteps,
      ",ARRIVE_RATE_ADD_PER20K=", add_arrivals,
      ",SIMDIR=", simdir
    ),
    paste(  sb, "--mail-type=ALL"   ),
    paste(  sb, "--mail-user=jrgant@brown.edu"),
    "module load R/4.0.3",
    "cd ~/data/jgantenb/xgcmsm/",
    "Rscript ./inst/analysis01_epi/02_main.r --vanilla",
    sep = "\n"
  )

  sdir <- here::here("inst", "analysis01_epi")
  writeLines(
    text = specs,
    con = file.path(sdir, paste0("main_epi", batchid, ".sh"))
  )

}

make_batch_script(
  jobname = "Main-Episim",
  walltime = "3:00:00",
  partition = "batch",
  mem = "3GB",
  ncores = 1,
  array = paste0("1-", nsim_spec),
  log_fullpath = "MainSim_ARRAY-%A_JOB-%J_SIMNO-%4a.log",
  batchid = 1,
  nsims = 1,
  nsteps = 52 * 65,
  add_arrivals = 1.285,
  simdir = "~/scratch/MainEpi"
)
