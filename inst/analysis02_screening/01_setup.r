################################################################################
## SETUP ##
################################################################################

library(pacman)
p_load(EpiModelHIV, rlecuyer)

selected_paramsets <-
  readRDS(here::here("inst", "cal", "main_analysis_inputs.rds"))

an02_path <- here::here("inst", "analysis02_screening")


################################################################################
## SAMPLE PARAMETER SETS ##
################################################################################

nsim_spec <- 1000

set.seed(5768933)
paramset_sampled <- sample(
  seq_len(length(selected_paramsets)),
  nsim_spec,
  replace = TRUE
)

saveRDS(
  paramset_sampled,
  here::here(an02_path, "paramset_sampled.rds")
)


################################################################################
## SETUP RANDOM NUMBER STREAMS ##
################################################################################

streams <- paste0("str", sprintf("%04d", 1:nsim_spec))
.lec.CreateStream(streams)

saveRDS(
  .lec.Random.seed.table,
  here::here(an02_path, "seeds.rds")
)


################################################################################
## MAKE BATCH SCRIPT ##
################################################################################

# This function writes a batch script to submit a job array.
make_batch_script <- function(jobname, walltime, partition, mem,
                              ncores, array, log_fullpath, batchid,
                              nsims, nsteps, add_arrivals, sti_screen_type,
                              kiss_exposure, simdir) {

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
      ",STI_SCREEN_TYPE=", sti_screen_type,
      ",STI_SCREEN_KISS_EXPOSURE=", kiss_exposure,
      ",SIMDIR=", simdir
    ),
    paste(  sb, "--mail-type=ALL"   ),
    paste(  sb, "--mail-user=jrgant@brown.edu"),
    "module load R/4.0.3",
    "cd ~/data/jgantenb/xgcmsm/",
    "Rscript ./inst/analysis02_screening/02.00_main.r --vanilla",
    sep = "\n"
  )

  sdir <- an02_path
  writeLines(
    text = specs,
    con = file.path(sdir, paste0(batchid, ".sh"))
  )

}

## NOTE: kiss_exposure is only looked for by the model in CDC scenarios.
##       Set to FALSE by default.
screen_specs <- list(
  sti_base = data.table(scenario = "base",        kiss_exposure = "FALSE"),
  sti_symp = data.table(scenario = "symptomatic", kiss_exposure = "FALSE"),
  sti_cdc1 = data.table(scenario = "cdc",         kiss_exposure = "FALSE"),
  sti_cdc2 = data.table(scenario = "cdc",         kiss_exposure = "TRUE"),
  sti_univ = data.table(scenario = "universal",   kiss_exposure = "FALSE")
)

lapply(names(screen_specs), function(.x) {
  index <- which(names(screen_specs) == .x)
  slug <- toupper(.x)
  tmp <- screen_specs[[.x]]
  make_batch_script(
    jobname = paste0("ScreenEpi-", slug),
    walltime = "3:00:00",
    partition = "batch",
    mem = "3GB",
    ncores = 1,
    array = paste0("1-", nsim_spec),
    log_fullpath = paste0("ScreenEpi-", slug, "_ARRAY-%A_JOB-%J_SIMNO-%4a.log"),
    batchid = paste0("02.", sprintf("%02d", index), "_", .x),
    nsims = 1,
    nsteps = 52 * 65,
    add_arrivals = 1.285,
    sti_screen_type = tmp[, scenario],
    kiss_exposure = tmp[, kiss_exposure],
    simdir = paste0("~/scratch/ScreenEpi-", slug)
  )
})
