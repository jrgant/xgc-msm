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

pretty_name <- Sys.getenv("JOB_PRETTY_NAME")
slurm_array_task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

job_specs <- readRDS(
  here::here("inst", "analysis02_gc_epi", "analysis02_job_specs.rds")
)

current_job <- job_specs[[pretty_name]]


# EPIMODEL SIM -----------------------------------------------------------------

## select seed for the current job and read it into the current environment
.lec.Random.seed.table <- current_job$seed_table

.lec.CurrentStream(
  .lec.Random.seed.table$name[slurm_array_task_id]
)

## read in network stats
netstats    <- get_est("netstats")
est         <- get_est("netest")
epistats    <- get_est("epistats")

param_set_sel    <- readRDS(
  here::here("inst", "cal", "analysis_param_sets.rds")
)

# Set environment variables (parameter sets) based on list position
# corresponding to SLURM_ARRAY_TASK_ID.
do.call(
  Sys.setenv,
  as.list(param_set_sel[[current_job$param_index[slurm_array_task_id]]])
)

# Check that manually set environment variables (set in sbatch)
if (Sys.getenv("SIMDIR") == "") {
  outpath <- "~/scratch/sim2"
} else {
  outpath <- Sys.getenv("SIMDIR")
}

sprintf("Output will be saved in %s", outpath)

## This function passes the environment variables to EpiModel functions
## in numeric format.
fmt_getenv <- function(x) as.numeric(Sys.getenv(x))

param <- param_msm(
  # external objects
  netstats          = netstats,
  epistats          = epistats,
  # demography
  arrival.age       = 18,
  # tweak marginal mortality to account for AIDS and achieve an arrival rate
  # that keeps population size at N = 20,000 (in expectation)
  a.rate            = netstats$demog$mortrate.marginal + fmt_getenv("ARRIVE_RATE_ADD_PER20K") / 20000,
  u2rgc.tprob       = fmt_getenv("U2RGC_PROB"),
  u2pgc.tprob       = fmt_getenv("U2PGC_PROB"),
  r2ugc.tprob       = fmt_getenv("R2UGC_PROB"),
  p2ugc.tprob       = fmt_getenv("P2UGC_PROB"),
  ## NOTE: Following tprobs used only if the kissing/rimming flags are active
  ## in control_msm.
  r2pgc.tprob       = current_job$specs[, tprob],
  p2rgc.tprob       = current_job$specs[, tprob],
  p2pgc.tprob       = current_job$specs[, tprob],
  rgc.ntx.int       = fmt_getenv("RECT_GC_DURAT_NOTX"),
  ugc.ntx.int       = fmt_getenv("URETH_GC_DURAT_NOTX"),
  pgc.ntx.int       = fmt_getenv("PHAR_GC_DURAT_NOTX"),
  ## Treated GC resolution probs c(after 1 week, after 2 weeks, after 3 weeks)
  rgc.tx.recov.pr   = c(1 - fmt_getenv("RECT_GC_RX_INFPR_WK1"), 0.5, 1),
  ugc.tx.recov.pr   = c(1 - fmt_getenv("URETH_GC_RX_INFPR_WK1"), 0.5, 1),
  pgc.tx.recov.pr   = c(1 - fmt_getenv("PHAR_GC_RX_INFPR_WK1"), 0.5, 1),
  rgc.sympt.prob    = fmt_getenv("RECT_GC_SYMPT_PROB"),
  ugc.sympt.prob    = fmt_getenv("URETH_GC_SYMPT_PROB"),
  pgc.sympt.prob    = fmt_getenv("PHAR_GC_SYMPT_PROB"),
  # STI testing
  gc.sympt.seek.test.prob = fmt_getenv("STITEST_PROB_GC_SYMPT"),
  ## NOTE: Changed the treatment probs to be conditional on someone's seeking
  ##       STI treatment. Repeat 3 times, one for each anatomic site.
  ##       Tune asymptomatic treatment probability to achieve the overall
  ##       probability of receiving an STI test.
  gc.sympt.prob.test  = rep(1, 3),
  ## NOTE: Order of probs: rectal, urethral, pharyngeal
  gc.asympt.prob.test = c(
    fmt_getenv("RECT_ASYMP_STITEST_PROB"),
    fmt_getenv("URETH_ASYMP_STITEST_PROB"),
    fmt_getenv("PHAR_ASYMP_STITEST_PROB")
  ),
  # HIV testing
  # @ORIG, prop. of MSM testing only at late-stage (AIDS)
  hiv.test.late.prob  = c(
    fmt_getenv("HIV_LATE_TESTER_PROB_BLACK"),
    fmt_getenv("HIV_LATE_TESTER_PROB_HISP"),
    fmt_getenv("HIV_LATE_TESTER_PROB_OTHER"),
    fmt_getenv("HIV_LATE_TESTER_PROB_WHITE")
  ),
  # HIV treatment parameters
  tt.part.supp        = rep(0, 4), # Partial VLS class, post ART
  tt.full.supp        = rep(1, 4), # Full VLS class, post ART
  tt.dur.supp         = rep(0, 4), # Durable VLS post ART
  tx.init.prob        = c(
    fmt_getenv("HIV_RX_INIT_PROB_BLACK"),
    fmt_getenv("HIV_RX_INIT_PROB_HISP"),
    fmt_getenv("HIV_RX_INIT_PROB_OTHER"),
    fmt_getenv("HIV_RX_INIT_PROB_WHITE")
  ),
  tx.halt.part.prob   = c(
    fmt_getenv("HIV_RX_HALT_PROB_BLACK"),
    fmt_getenv("HIV_RX_HALT_PROB_HISP"),
    fmt_getenv("HIV_RX_HALT_PROB_OTHER"),
    fmt_getenv("HIV_RX_HALT_PROB_WHITE")
  ),
  tx.halt.full.rr     = rep(1, 4), # not used
  tx.halt.dur.rr      = rep(1, 4), # not used
  tx.reinit.part.prob = c(
    fmt_getenv("HIV_RX_REINIT_PROB_BLACK"),
    fmt_getenv("HIV_RX_REINIT_PROB_HISP"),
    fmt_getenv("HIV_RX_REINIT_PROB_OTHER"),
    fmt_getenv("HIV_RX_REINIT_PROB_WHITE")
  ), # @ORIG
  tx.reinit.full.rr   = rep(1, 4),
  tx.reinit.dur.rr    = rep(1, 4),
  # Scaling parameters
  ai.acts.scale.mc    = fmt_getenv("SCALAR_AI_ACT_RATE"),
  oi.acts.scale.mc    = fmt_getenv("SCALAR_OI_ACT_RATE"),
  kiss.rate.main      = current_job$specs[, kiss_rate_main],
  kiss.rate.casl      = current_job$specs[, kiss_rate_casl],
  kiss.prob.oo        = current_job$specs[, kiss_prob_oo],
  rim.rate.main       = current_job$specs[, rim_rate_main],
  rim.rate.casl       = current_job$specs[, rim_rate_casl],
  rim.prob.oo         = current_job$specs[, rim_prob_oo],
  trans.scale         = c(
    fmt_getenv("SCALAR_HIV_TRANS_PROB_BLACK"),
    fmt_getenv("SCALAR_HIV_TRANS_PROB_HISP"),
    fmt_getenv("SCALAR_HIV_TRANS_PROB_OTHER"),
    fmt_getenv("SCALAR_HIV_TRANS_PROB_WHITE")
  ),
  cdc.sti.int         = 12,
  cdc.sti.hr.int      = 6,
  cond.eff            = fmt_getenv("CONDOM_EFF_HIV"),
  ## NOTE: Change cond.fail params to 0 to turn off (reflect uncertainty
  ##       in cond. effect)
  cond.fail           = rep(0, 4),
  sti.cond.fail       = rep(0, 4),
  sti.cond.eff        = fmt_getenv("CONDOM_EFF_GC"),
  circ.prob           = netstats$inputs$circ.probs,
  # PrEP parameters
  prep.discont.rate   = c(
    fmt_getenv("PREP_DISCONT_RATE_BLACK"),
    fmt_getenv("PREP_DISCONT_RATE_HISP"),
    fmt_getenv("PREP_DISCONT_RATE_OTHER"),
    fmt_getenv("PREP_DISCONT_RATE_WHITE")
  )
)

init <- init_msm(
  prev.ugc = 0.01,
  prev.rgc = 0.01,
  prev.pgc = 0.01
)

## Spec summary
nsim_env <- fmt_getenv("NSIMS")
nsteps_env <- fmt_getenv("NSTEPS")
ncores_env <- fmt_getenv("SLURM_NPROCS")

print(paste("NSIMS:", nsim_env))
print(paste("NSTEPS:", nsteps_env))
print(paste("SLURM_NPROCS:", ncores_env))

control <- control_msm(
  # Computing options (set in batch script or using sbatch on SLURM)
  nsims   = nsim_env,
  nsteps  = nsteps_env,
  ncores  = ncores_env,
  # Epidemic simulation Modules
  initialize.FUN    = initialize_msm,
  aging.FUN         = aging_msm,
  departure.FUN     = departure_msm,
  arrival.FUN       = arrival_msm,
  hivtest.FUN       = hivtest_msm,
  hivtx.FUN         = hivtx_msm,
  hivprogress.FUN   = hivprogress_msm,
  hivvl.FUN         = hivvl_msm,
  resim_nets.FUN    = simnet_msm,
  acts.FUN          = acts_msm,
  condoms.FUN       = condoms_msm, # NOTE All act lists are finalized here
  position.FUN      = position_msm,
  prep.FUN          = prep_msm,
  hivtrans.FUN      = hivtrans_msm,
  stirecov.FUN      = stirecov_msm,
  stitx.FUN         = stitx_msm,
  stitrans.FUN      = stitrans_msm_rand,
  prev.FUN          = prevalence_msm,
  # Epidemic simulation options
  transRoute_Kissing      = current_job$specs[, as.logical(kiss_flag)],
  transRoute_Rimming      = current_job$specs[, as.logical(rim_flag)],
  gcUntreatedRecovDist    = "geom",
  stiScreeningProtocol    = "base",
  skip.check              = TRUE,
  cdcExposureSite_Kissing = FALSE, # FLAG: Kissing considered exposure for CDC protocol?
  tergmLite               = TRUE, # NOTE Must set to avoid error from saveout.net()
  debug_stitx             = FALSE,
  save.network		        = FALSE,
  save.nwstats            = FALSE,
  verbose                 = FALSE
)

sim <- netsim(est, param, init, control)

.lec.CurrentStreamEnd()

sim[["seed.table.state"]] <- .lec.GetState(
  .lec.Random.seed.table$name[slurm_array_task_id]
)


################################################################################
## WRITE OUTPUT ##
################################################################################

if (!dir.exists(outpath)) dir.create(outpath)

saveRDS(
  object = sim,
  file = file.path(
    outpath,
    paste0(
      sprintf("episim_%04d", slurm_array_task_id),
      "_", Sys.getenv("SLURM_JOB_ID"), ".rds")
  )
)
