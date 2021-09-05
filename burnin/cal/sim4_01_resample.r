# SETUP ------------------------------------------------------------------------

pacman::p_load(
  xgcmsm,
  methods,
  EpiModelHIV,
  data.table,
  magrittr,
  rms,
  stringr,
  pscl,
  lhs,
  rlecuyer
)

inputs <- readRDS(here::here("inst", "cal", "main_analysis_inputs.rds"))


# RANDOMLY RESAMPLE FROM SELECTED INPUTS ---------------------------------------

set.seed(1066)
sampled_list <- sample(names(inputs), 1000, replace = TRUE)

caldir <- here::here("burnin", "cal")

if (!file.exists(file.path(caldir, "sim4"))) {
  dir.create(file.path(caldir, "sim4"))
}

saveRDS(sampled_list, file = file.path(caldir, "sim4", "sampled_sim4.rds"))


# CREATE RNG STREAMS FOR SIMULATIONS -------------------------------------------

## use the rlecuyer package to create independent random numbers streams
## for each simulation.

streams <- paste0("str", sprintf("%04d", seq_along(sampled_list)))
.lec.CreateStream(streams)

saveRDS(.lec.Random.seed.table, file.path(caldir, "sim4", "seeds_sim4.rds"))


# MAKE BATCH SCRIPTS -----------------------------------------------------------

## # make batch scripts to submit arrays
## # SLURM limit is 1200 job requests at a time
nbatches <- 1
start_index <- seq(1, length(sampled_list), length(sampled_list) / nbatches)

arrays <- sapply(
  start_index,
  function(.x) paste0(.x, "-", .x + (length(sampled_list) / nbatches) - 1)
)

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
    "Rscript ./burnin/cal/sim4_02_lhs.r --vanilla",
    sep = "\n"
  )

  sdir <- here::here("burnin", "cal", "sim4")
  writeLines(
    text = specs,
    con = file.path(sdir, paste0("sim4_batch", batchid, ".sh"))
  )

}

for (i in seq_along(arrays)) {
  make_batch_script(
    jobname = "Sim4-LHS-XGC",
    walltime = "3:00:00",
    partition = "batch",
    mem = "3GB",
    ncores = 1,
    array = arrays[i],
    log_fullpath = "LHS-Sim4_ARRAY-%A_JOB-%J_SIMNO-%4a.log",
    batchid = i,
    nsims = 1,
    nsteps = 3120,
    add_arrivals = 1.285,
    simdir = "~/scratch/sim4"
  )
}
