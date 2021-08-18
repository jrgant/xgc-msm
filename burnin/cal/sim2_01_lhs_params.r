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

lhs_new_limits <- readRDS(here::here("inst", "cal", "sim1_sel_lhs_limits.rds"))


# CALIBRATED PRIOR RANGES ------------------------------------------------------

# convenience function to add a uniform prior to the prior list
pvec2 <- function(name) {
  tmp <- lhs_new_limits[input == name]
  c("unif", tmp[, s2_ll], tmp[, s2_ul], name)
}

priors <- list(
  # GC per-act transmission probability by pathways
  pvec2("U2RGC_PROB"), # u2rgc.prob
  pvec2("U2PGC_PROB"), # u2pgc.prob
  pvec2("R2UGC_PROB"), # r2ugc.prob
  pvec2("P2UGC_PROB"), # p2ugc.prob
  # Untreated infection durations
  pvec2("RECT_GC_DURAT_NOTX"), # rectal
  pvec2("URETH_GC_DURAT_NOTX"), # urethral
  pvec2("PHAR_GC_DURAT_NOTX"), # pharyngeal
  # Probability of still being GC-infected after x weeks of treatment
  pvec2("RECT_GC_RX_INFPR_WK1"), # treated rectal gc, 1st week still infected prop.
  pvec2("URETH_GC_RX_INFPR_WK1"), # treated urethral gc, 1st week still infected prop.
  pvec2("PHAR_GC_RX_INFPR_WK1"), # treated pharyngeal gc, 1st week still infected prop.
  # Symptom probability
  pvec2("RECT_GC_SYMPT_PROB"), # rectal
  pvec2("URETH_GC_SYMPT_PROB"), # urethral
  pvec2("PHAR_GC_SYMPT_PROB"), # pharyngeal
  # GC symptomatic, weekly probability of testing
  pvec2("STITEST_PROB_GC_SYMPT"),
  # Probability of testing at asymptomatic sites in clinic
  pvec2("RECT_ASYMP_STITEST_PROB"),  # rectal
  pvec2("URETH_ASYMP_STITEST_PROB"),  # urethral
  pvec2("PHAR_ASYMP_STITEST_PROB"),  # pharyngeal
  # Probability of testing at symptomatic sites in clinic
  pvec2("RECT_SYMP_STITEST_PROB"),  # rectal
  pvec2("URETH_SYMP_STITEST_PROB"),  # urethral
  pvec2("PHAR_SYMP_STITEST_PROB"),  # pharyngeal
  # HIV late-tester probability (priors from Singh 2017, MMWR)
  pvec2("HIV_LATE_TESTER_PROB_BLACK"),
  pvec2("HIV_LATE_TESTER_PROB_HISP"),
  pvec2("HIV_LATE_TESTER_PROB_OTHER"),
  pvec2("HIV_LATE_TESTER_PROB_WHITE"),
  # tx.init.prob
  pvec2("HIV_RX_INIT_PROB_BLACK"),
  pvec2("HIV_RX_INIT_PROB_HISP"),
  pvec2("HIV_RX_INIT_PROB_OTHER"),
  pvec2("HIV_RX_INIT_PROB_WHITE"),
  # tx.halt.part.prob (Singh 2017, MMWR)
  pvec2("HIV_RX_HALT_PROB_BLACK"),
  pvec2("HIV_RX_HALT_PROB_HISP"),
  pvec2("HIV_RX_HALT_PROB_OTHER"),
  pvec2("HIV_RX_HALT_PROB_WHITE"),
  # tx.reinit.part.prob
  pvec2("HIV_RX_REINIT_PROB_BLACK"),
  pvec2("HIV_RX_REINIT_PROB_HISP"),
  pvec2("HIV_RX_REINIT_PROB_OTHER"),
  pvec2("HIV_RX_REINIT_PROB_WHITE"),
  # sex act scalars
  pvec2("SCALAR_AI_ACT_RATE"), # ai.acts.scale
  pvec2("SCALAR_OI_ACT_RATE"), # oi.acts.scale
  # HIV transmission prob. scalars
  pvec2("SCALAR_HIV_TRANS_PROB_BLACK"), # black
  pvec2("SCALAR_HIV_TRANS_PROB_HISP"), # hispanic
  pvec2("SCALAR_HIV_TRANS_PROB_OTHER"), # other
  pvec2("SCALAR_HIV_TRANS_PROB_WHITE"), # white
  # condom efficacy
  pvec2("CONDOM_EFF_HIV"), # cond.eff (HIV), 1 - rel. risk, condom use vs. no
  pvec2("CONDOM_EFF_GC"),   # sti.cond.eff, relative risk, condom use vs. no
  # PrEP discontinuation rates
  pvec2("PREP_DISCONT_RATE_BLACK"),
  pvec2("PREP_DISCONT_RATE_HISP"),
  pvec2("PREP_DISCONT_RATE_OTHER"),
  pvec2("PREP_DISCONT_RATE_WHITE"),
  # probability that an agent is of the type to stop sexual activity
  # while symptomatic or on treatment
  pvec2("ACT_STOPPER_PROB"),
  pvec2("HIV_TRANS_RR_RGC"),
  pvec2("HIV_TRANS_RR_UGC")
)


# DRAW LATIN HYPERCUBE ---------------------------------------------------------

set.seed(1044)
lhs_unif <- randomLHS(1000, length(priors))

draw_param <- function(lhscol, lhsrow) {
  p_min <- as.numeric(priors[[lhscol]][2])
  p_max <- as.numeric(priors[[lhscol]][3])
  param <- priors[[lhscol]][4]
  draw <- (p_max - p_min) * lhsrow[lhscol] + p_min
  names(draw) <- param
  draw
}

lhs_real <-
  lapply(seq_len(nrow(lhs_unif)), function(.x) {
    lhsrow <- lhs_unif[.x,]
    sapply(
      seq_len(length(lhsrow)),
      function(.y) draw_param(.y, lhsrow = lhsrow)
    )
  })

# check that all priors are included in all parameter sets
priornames_ok <- sapply(priors, function(.x) {
  currname <- .x[4]
  sapply(lhs_real, function(.y) {
    currname %in% names(.y)
  })
}) %>% all

# check that all sampled values across all parameter sets
# fall within specified prior bounds
priorvals_ok <- sapply(priors, function(.x) {
  currname <- .x[4]
  p_min <- as.numeric(.x[2])
  p_max <- as.numeric(.x[3])
  sapply(lhs_real, function(.y) {
    (p_min <= .y[currname]) & (p_max >= .y[currname])
  })
}) %>% all

caldir <- here::here("burnin", "cal")

if (priornames_ok & priorvals_ok) {
  if (!file.exists(file.path(caldir, "sim2"))) {
    dir.create(file.path(caldir, "sim2"))
  }
  saveRDS(lhs_real, file = file.path(caldir, "sim2", "lhs_sim2.rds"))
} else {
  stop("Check prior specification.")
}


# CREATE RNG STREAMS FOR SIMULATIONS -------------------------------------------

## use the rlecuyer package to create independent random numbers streams
## for each simulation.

streams <- paste0("str", sprintf("%04d", seq_along(lhs_real)))
.lec.CreateStream(streams)

saveRDS(.lec.Random.seed.table, file.path(caldir, "sim2", "seeds_sim2.rds"))


# MAKE BATCH SCRIPTS -----------------------------------------------------------

## # make batch scripts to submit arrays
## # SLURM limit is 1200 job requests at a time
## nbatches <- 5
## start_index <- seq(1, length(lhs_real), length(lhs_real) / nbatches)

## arrays <- sapply(
##   start_index,
##   function(.x) paste0(.x, "-", .x + (length(lhs_real) / nbatches) - 1)
## )

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
    "Rscript ./burnin/cal/sim2_02_lhs.r --vanilla",
    sep = "\n"
  )

  sdir <- here::here("burnin", "cal", "sim2")
  writeLines(
    text = specs,
    con = file.path(sdir, paste0("sim2_batch", batchid, ".sh"))
  )

}

make_batch_script(
  jobname = "Sim2-LHS-XGC",
  walltime = "3:00:00",
  partition = "batch",
  mem = "3GB",
  ncores = 1,
  array = "1-1000",
  log_fullpath = "LHS-Sim2_ARRAY-%A_JOB-%J_SIMNO-%4a.log",
  batchid = 1,
  nsims = 1,
  nsteps = 3120,
  add_arrivals = 1.285,
  simdir = "~/scratch/sim2"
)
