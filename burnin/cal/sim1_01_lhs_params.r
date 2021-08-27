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


# INITIAL PRIOR RANGES ---------------------------------------------------------

# convenience function to add a uniform prior to the prior list
pvec <- function(ll, ul, name) c("unif", ll, ul, name)

priors <- list(
  # GC per-act transmission probability by pathways
  pvec(0, 0.35, "U2RGC_PROB"), # u2rgc.prob
  pvec(0, 0.45, "U2PGC_PROB"), # u2pgc.prob
  pvec(0, 0.55, "R2UGC_PROB"), # r2ugc.prob
  pvec(0, 0.25, "P2UGC_PROB"), # p2ugc.prob
  pvec(0, 0.45, "R2PGC_PROB"), # r2pgc.prob
  pvec(0, 0.20, "P2RGC_PROB"), # p2rgc.prob
  pvec(0, 0.10, "P2PGC_PROB"), # p2pgc.prob
  # Untreated infection durations
  pvec(2, 22, "RECT_GC_DURAT_NOTX"), # rectal
  pvec(1, 24, "URETH_GC_DURAT_NOTX"), # urethral
  pvec(3, 20, "PHAR_GC_DURAT_NOTX"), # pharyngeal
  # Probability of still being GC-infected after x weeks of treatment
  pvec(0.01, 0.11, "RECT_GC_RX_INFPR_WK1"), # treated rectal gc, 1st wk
  pvec(0.01, 0.20, "URETH_GC_RX_INFPR_WK1"), # treated urethral gc, 1st wk
  pvec(0.06, 0.20, "PHAR_GC_RX_INFPR_WK1"), # treated pharyngeal gc, 1st wk
  # Symptom probability
  pvec(0.06, 0.46, "RECT_GC_SYMPT_PROB"), # rectal
  pvec(0.46, 0.99, "URETH_GC_SYMPT_PROB"), # urethral
  pvec(0.00, 0.15, "PHAR_GC_SYMPT_PROB"), # pharyngeal
  # GC symptomatic, weekly probability of testing
  pvec(0.5, 1, "STITEST_PROB_UGC_SYMPT"),
  pvec(0.01, 0.99, "STITEST_RGC_RR_SYMPT"),
  pvec(0.01, 0.99, "STITEST_PGC_RR_SYMPT"),
  # Probability of testing at asymptomatic sites in clinic
  pvec(0.400, 0.800, "RECT_ASYMP_STITEST_PROB"),  # rectal
  pvec(0.497, 0.657, "URETH_ASYMP_STITEST_PROB"),  # urethral
  pvec(0.522, 0.749, "PHAR_ASYMP_STITEST_PROB"),  # pharyngeal
  # Probability of testing at symptomatic sites in clinic
  pvec(0.5, 1, "RECT_SYMP_STITEST_PROB"),  # rectal
  pvec(0.5, 1, "URETH_SYMP_STITEST_PROB"),  # urethral
  pvec(0.5, 1, "PHAR_SYMP_STITEST_PROB"),  # pharyngeal
  # HIV late-tester probability (priors from Singh 2017, MMWR)
  pvec(0.160, 0.419, "HIV_LATE_TESTER_PROB_BLACK"),
  pvec(0.202, 0.391, "HIV_LATE_TESTER_PROB_HISP"),
  pvec(0.207, 0.369, "HIV_LATE_TESTER_PROB_OTHER"),
  pvec(0.222, 0.377, "HIV_LATE_TESTER_PROB_WHITE"),
  # tx.init.prob
  pvec(0.001, 0.99, "HIV_RX_INIT_PROB_BLACK"),
  pvec(0.001, 0.99, "HIV_RX_INIT_PROB_HISP"),
  pvec(0.001, 0.99, "HIV_RX_INIT_PROB_OTHER"),
  pvec(0.001, 0.99, "HIV_RX_INIT_PROB_WHITE"),
  # tx.halt.part.prob (Singh 2017, MMWR)
  pvec(1 - (1 - 0.464)^(1/52), 1 - (1 - 0.464)^(1/52), "HIV_RX_HALT_PROB_BLACK"),
  pvec(1 - (1 - 0.416)^(1/52), 1 - (1 - 0.416)^(1/52), "HIV_RX_HALT_PROB_HISP"),
  pvec(1 - (1 - 0.366)^(1/52), 1 - (1 - 0.366)^(1/52), "HIV_RX_HALT_PROB_OTHER"),
  pvec(1 - (1 - 0.406)^(1/52), 1 - (1 - 0.406)^(1/52), "HIV_RX_HALT_PROB_WHITE"),
  # tx.reinit.part.prob
  pvec(0.001, 0.99, "HIV_RX_REINIT_PROB_BLACK"),
  pvec(0.001, 0.99, "HIV_RX_REINIT_PROB_HISP"),
  pvec(0.001, 0.99, "HIV_RX_REINIT_PROB_OTHER"),
  pvec(0.001, 0.99, "HIV_RX_REINIT_PROB_WHITE"),
  # sex act scalars
  pvec(1, 1, "SCALAR_AI_ACT_RATE"), # ai.acts.scale
  pvec(1, 1, "SCALAR_OI_ACT_RATE"), # oi.acts.scale
  # kissing rate/prob priors
  pvec(0, 10,  "KISS_RATE_MAIN"),
  pvec(0, 10,  "KISS_RATE_CASUAL"),
  pvec(0.5, 1, "KISS_PROB_ONETIME"),
  # rimming rate/prob priors
  pvec(0, 10,  "RIM_RATE_MAIN"),
  pvec(0, 10,  "RIM_RATE_CASUAL"),
  pvec(0, 1,   "RIM_PROB_ONETIME"),
  # HIV transmission prob. scalars
  pvec(0.75, 3, "SCALAR_HIV_TRANS_PROB_BLACK"), # black
  pvec(0.75, 3, "SCALAR_HIV_TRANS_PROB_HISP"), # hispanic
  pvec(0.75, 3, "SCALAR_HIV_TRANS_PROB_OTHER"), # other
  pvec(1, 1, "SCALAR_HIV_TRANS_PROB_WHITE"), # white
  # condom efficacy
  pvec(0.6, 1, "CONDOM_EFF_HIV"), # cond.eff (HIV), 1 - rel. risk, condom vs. no
  pvec(0.7, 1, "CONDOM_EFF_GC"),   # sti.cond.eff, rel. risk, condom vs. no
  # PrEP discontinuation rates
  pvec(0, 0.05, "PREP_DISCONT_RATE_BLACK"),
  pvec(0, 0.05, "PREP_DISCONT_RATE_HISP"),
  pvec(0, 0.05, "PREP_DISCONT_RATE_OTHER"),
  pvec(0, 0.05, "PREP_DISCONT_RATE_WHITE"),
  pvec(0.2, 1, "ACT_STOPPER_PROB"),
  pvec(1.2, 6.4, "HIV_TRANS_RR_RGC"), # Vaughan et al.
  pvec(1, 3, "HIV_TRANS_RR_UGC") # Sam's STI incidence paper
)

## Write initial ranges to a data set.
priors_dt <- as.data.table(Reduce(rbind, priors))[, 2:4]
sapply(priors_dt, class)

lim_labs <- c("s1_ll", "s1_ul")
setnames(priors_dt, c("V2", "V3", "V4"), c(lim_labs, "input"))

priors_dt[, (lim_labs) := lapply(.SD, as.numeric), .SDcols = lim_labs][]

saveRDS(priors_dt, here::here("burnin", "cal", "sim1", "sim1_priors.rds"))


# DRAW LATIN HYPERCUBE ---------------------------------------------------------

set.seed(1971)
lhs_unif <- randomLHS(5000, length(priors))

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
  if (!file.exists(file.path(caldir, "sim1"))) {
    dir.create(file.path(caldir, "sim1"))
  }
  saveRDS(lhs_real, file = file.path(caldir, "sim1", "lhs_sim1.rds"))
} else {
  stop("Check prior specification.")
}


# CREATE RNG STREAMS FOR SIMULATIONS -------------------------------------------

## use the rlecuyer package to create independent random numbers streams
## for each simulation.

streams <- paste0("str", sprintf("%04d", seq_along(lhs_real)))
.lec.CreateStream(streams)

saveRDS(.lec.Random.seed.table, file.path(caldir, "sim1", "seeds_sim1.rds"))


# MAKE BATCH SCRIPTS -----------------------------------------------------------

# make batch scripts to submit arrays
# SLURM limit is 1200 job requests at a time
nbatches <- 5
start_index <- seq(1, length(lhs_real), length(lhs_real) / nbatches)

arrays <- sapply(
  start_index,
  function(.x) paste0(.x, "-", .x + (length(lhs_real) / nbatches) - 1)
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
    "Rscript ./burnin/cal/sim1_02_lhs.r --vanilla",
    sep = "\n"
  )

  sdir <- here::here("burnin", "cal", "sim1")
  writeLines(
    text = specs,
    con = file.path(sdir, paste0("sim1_batch", batchid, ".sh"))
  )

}

for (i in seq_len(length(arrays))) {
  make_batch_script(
    jobname = "Sim1-LHS-XGC",
    walltime = "3:00:00",
    partition = "batch",
    mem = "3GB",
    ncores = 1,
    array = arrays[i],
    log_fullpath = "LHS-Sim1_ARRAY-%A_JOB-%J_SIMNO-%4a.log",
    batchid = i,
    nsims = 1,
    nsteps = 3120,
    add_arrivals = 1.285,
    simdir = "~/scratch/sim1"
  )
}
