################################################################################
## SETUP ##
################################################################################

pacman::p_load(data.table, rlecuyer)

sim3_dir         <- here::here("inst", "cal")
analysis_dir     <- here::here("inst", "analysis02_gc_epi")
lims             <- readRDS(here::here(sim3_dir, "sim3_sel_lhs_limits.rds"))
anlys_param_sets <- readRDS(here::here(sim3_dir, "analysis_param_sets.rds"))

nsims_spec <- readRDS(
  here::here("inst", "analysis01_epi", "nsims_per_scenario.rds")
)

set.seed(19670723)


################################################################################
## SENSITIVITY ANALYSIS PARAMETER RANGES ##
################################################################################

trans_probs <- c("U2RGC_PROB", "U2PGC_PROB", "R2UGC_PROB", "P2UGC_PROB")

gcprob_lims <-
  range(suppressWarnings(
    melt(lims[input %in% trans_probs, .(s3_ll, s3_ul)])[, value]
  ))

gcprob_lims_rnd <- round(gcprob_lims / 0.05) * 0.05
gcprob_vals <- seq(gcprob_lims_rnd[1], gcprob_lims_rnd[2], by = 0.1)


################################################################################
## SET UP KISSING SENSITIVITY ANALYSES ##
################################################################################

kiss_rate_main <- 7  # assume kissing 7 times/week in main partnerships
kiss_rate_casl <- 3  # assume kissing 3 times/week in casual partnerships
kiss_prob_oo   <- 1  # assume kissing always occurs in one-time partnerships

kiss_grid <- as.data.table(expand.grid(
  tprob = gcprob_vals,
  kiss_rate_main = kiss_rate_main,
  kiss_rate_casl = kiss_rate_casl,
  kiss_prob_oo   = kiss_prob_oo,
  kiss_flag      = 1
))

kiss_grid[, pretty_name := paste0("KISS", 1:.N)][]


################################################################################
## SET UP RIMMING SENSITIVITY ANALYSIS ##
################################################################################

## Background
## S. L. Reisner, M. J. Mimiaga, M. Skeer, and K. H. Mayer, “Beyond Anal Sex:
## Sexual Practices Associated with HIV Risk Reduction among Men Who Have Sex
## with Men in Boston, Massachusetts,” AIDS Patient Care and STDs, vol. 23,
## no. 7, Art. no. 7, 2009-07, doi: 10.1089/apc.2008.0249.
##
## NOTE
## Parameters here ARE NOT based on the reference above. Reisner et al.
## simply find that rimming is much less common than other sexual activities.
rim_rate_main <- 1 # assume level of rim rate likely above background
rim_rate_casl <- 1
rim_prob_oo   <- 0.01

rim_grid <- as.data.table(expand.grid(
  tprob = gcprob_vals,
  rim_rate_main = rim_rate_main,
  rim_rate_casl = rim_rate_casl,
  rim_prob_oo   = rim_prob_oo,
  rim_flag      = 1
))

rim_grid[, pretty_name := paste0("RIM", 1:.N)][]


################################################################################
## SET UP JOB SPECS ##
################################################################################

## The main analysis involves estimating epidemiologic parameters of gonorrhea:
## anatomic site specific measurements of prevalence and incidence, in addition
## to proportion of infections attributable to each transmission pathway.
## Kissing and rimming are turned off in the main analysis.

spec_table <- rbindlist(
  list(data.table(pretty_name = "MAIN"), kiss_grid, rim_grid),
  fill = TRUE
)

## Set all NAs to 0
cols <- names(spec_table)[-1] # omit pretty_name
spec_table[, (cols) := lapply(.SD, function(.x) fifelse(is.na(.x), 0, .x)),
             .SDcols = cols]

num_analyses <- spec_table[, .N]

## Parameter indexes to retrieve input parameter sets.
param_index <- lapply(
  seq_len(num_analyses),
  function(.x) {
    sample(
      seq_len(length(anlys_param_sets)),
      size = nsims_spec,
      replace = TRUE
    )
  }
)


## RNG streams
streams <- paste0("str", sprintf("%03d", seq_len(nsims_spec)))

## The process below warns that seeds are being overwritten during the loop.
## That should be fine, as we're simply creating a list to store RNG stream
## tables for use in analysis scripts, not calling on these seeds now.
seed_tables <- suppressWarnings(
  lapply(
    seq_len(num_analyses),
    function(.x) {
      .lec.CreateStream(streams)
      .lec.Random.seed.table
    }))

speclist <- lapply(
  setNames(seq_len(num_analyses), spec_table[, pretty_name]),
  function(.x) {
    list(
      specs = spec_table[.x, 2:length(names(spec_table))],
      param_index = param_index[[.x]],
      seed_table = seed_tables[[.x]]
    )
  })

saveRDS(speclist, here::here(analysis_dir, "analysis02_job_specs.rds"))


################################################################################
## MAKE BATCH SCRIPTS ##
################################################################################

# This function writes a batch script to submit a job array.
make_batch_script <- function(jobname, walltime, partition, mem,
                              ncores, array, log_fullpath, batchid,
                              nsims, nsteps, add_arrivals, simdir) {

  sb <- "#SBATCH"

  specs <- paste(
    "#!/bin/bash",
    paste(sb, "-J", jobname),
    paste0(sb, " --time=", walltime),
    paste(sb, "-p", partition),
    paste0(sb, " --mem=", mem),
    paste(sb, "-n", ncores),
    paste0(sb, " --array=", array),
    paste(sb, "-o", log_fullpath),
    paste0(
      sb,
      " --export=ALL,NSIMS=", nsims,
      ",NSTEPS=", nsteps,
      ",ARRIVE_RATE_ADD_PER20K=", add_arrivals,
      ",SIMDIR=", simdir,
      ",JOB_PRETTY_NAME=", jobname
    ),
    paste(sb, "--mail-type=ALL"),
    paste(sb, "--mail-user=jrgant@brown.edu"),
    "module load R/4.0.3",
    "cd ~/data/jgantenb/xgcmsm/",
    "Rscript ./inst/analysis02_gc_epi/02_episim.r --vanilla",
    sep = "\n"
  )

  sdir <- here::here("inst", "analysis02_gc_epi")
  writeLines(
    text = specs,
    con = file.path(sdir, paste0(jobname, "_episim", batchid, ".sh"))
  )
}

lapply(names(speclist), function(.x) {
  make_batch_script(
    jobname = .x,
    walltime = "3:00:00",
    partition = "batch",
    mem = "3GB",
    ncores = 1,
    array = paste0("1-", nsims_spec),
    log_fullpath = paste0(.x, "_Sim_ARRAY-%A_JOB-%J_SIMNO-%4a.log"),
    batchid = "",
    nsims = 1,
    nsteps = 52 * 65,
    add_arrivals = 1.285,
    simdir = sprintf("~/scratch/%s", tolower(.x))
  )
})
