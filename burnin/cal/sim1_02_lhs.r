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

slurm_array_task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# EPIMODEL SIM -----------------------------------------------------------------

## select seed for the current job and read it into the current environment
.lec.Random.seed.table <- readRDS(
  here::here("burnin", "cal", "sim1", "seeds_sim1.rds")
)

.lec.CurrentStream(
  .lec.Random.seed.table$name[slurm_array_task_id]
)

## read in network stats
netstats    <- get_est("netstats")
est         <- get_est("netest")
epistats    <- get_est("epistats")

lhs_real    <- readRDS(
  here::here("burnin", "cal", "sim1", "lhs_sim1.rds")
)

# Set environment variables (parameter sets) based on list position
# corresponding to SLURM_ARRAY_TASK_ID.
do.call(
  Sys.setenv,
  as.list(lhs_real[[slurm_array_task_id]])
)

# Check that manually set environment variables (set in sbatch)
if (Sys.getenv("SIMDIR") == "") {
  outpath <- "~/scratch/sim1"
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
  u2rgc.tprob       = fmt_getenv("U2RGC_PROB"), # ureth-to-rect GC trans. prob.
  u2pgc.tprob       = fmt_getenv("U2PGC_PROB"), # ureth-to-phar GC trans. prob.
  r2ugc.tprob       = fmt_getenv("R2UGC_PROB"), # rect-to-ureth GC trans. prob.
  p2ugc.tprob       = fmt_getenv("P2UGC_PROB"), # phar-to-ureth GC trans. prob.
  ## NOTE: Following tprobs used only if the kissing/rimming flags are active
  ## in control_msm.
  r2pgc.tprob       = fmt_getenv("R2PGC_PROB"), # rect-to-phar GC trans. prob.
  p2rgc.tprob       = fmt_getenv("P2RGC_PROB"), # phar-to-rect GC trans. prob.
  p2pgc.tprob       = fmt_getenv("P2PGC_PROB"), # kissing GC trans. prob.
  ## GC infection duration
  rgc.ntx.int       = fmt_getenv("RECT_GC_DURAT_NOTX"),
  ugc.ntx.int       = fmt_getenv("URETH_GC_DURAT_NOTX"),
  pgc.ntx.int       = fmt_getenv("PHAR_GC_DURAT_NOTX"),
  ## Treated GC resolution probs c(after 1 week, after 2 weeks, after 3 weeks)
  rgc.tx.recov.pr   = c(1 - fmt_getenv("RECT_GC_RX_INFPR_WK1"), 0.5, 1),
  ugc.tx.recov.pr   = c(1 - fmt_getenv("URETH_GC_RX_INFPR_WK1"), 0.5, 1),
  pgc.tx.recov.pr   = c(1 - fmt_getenv("PHAR_GC_RX_INFPR_WK1"), 0.5, 1),
  rgc.sympt.prob    = fmt_getenv("RECT_GC_SYMPT_PROB"), # rectal symptom prob
  ugc.sympt.prob    = fmt_getenv("URETH_GC_SYMPT_PROB"), # ureth symptom prob
  pgc.sympt.prob    = fmt_getenv("PHAR_GC_SYMPT_PROB"), # phar symptom prob
  # STI testing
  ugc.sympt.seek.test.prob = fmt_getenv("STITEST_PROB_UGC_SYMPT"),
  rgc.sympt.seek.test.rr = fmt_getenv("STITEST_RGC_RR_SYMPT"),
  pgc.sympt.seek.test.rr = fmt_getenv("STITEST_PGC_RR_SYMPT"),
  ## NOTE: Changed the treatment probs to be conditional on someone's seeking
  ##       STI treatment. Repeat 3 times, one for each anatomic site.
  ##       Tune asymptomatic treatment probability to achieve the overall
  ##       probability of receiving an STI test.
  ## NOTE: Order of probs: rectal, urethral, pharyngeal
  gc.asympt.prob.test = c(
    fmt_getenv("RECT_ASYMP_STITEST_PROB"),
    fmt_getenv("URETH_ASYMP_STITEST_PROB"),
    fmt_getenv("PHAR_ASYMP_STITEST_PROB")
  ),
  gc.sympt.prob.test = c(
    fmt_getenv("RECT_SYMP_STITEST_PROB"),
    fmt_getenv("URETH_SYMP_STITEST_PROB"),
    fmt_getenv("PHAR_SYMP_STITEST_PROB")
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
  tt.part.supp        = rep(0, 4), # Partial VLS class, turned off
  tt.full.supp        = rep(1, 4), # Full VLS class, post ART
  tt.dur.supp         = rep(0, 4), # Durable VLS post ART, turned off
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
  # kissing rate/prob
  kiss.rate.main      = fmt_getenv("KISS_RATE_MAIN"),
  kiss.rate.casl      = fmt_getenv("KISS_RATE_CASUAL"),
  kiss.prob.oo        = fmt_getenv("KISS_PROB_ONETIME"),
  # rimming rate/prob
  rim.rate.main       = fmt_getenv("RIM_RATE_MAIN"),
  rim.rate.casl       = fmt_getenv("RIM_RATE_CASUAL"),
  rim.prob.oo         = fmt_getenv("RIM_PROB_ONETIME"),
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
  ),
  # Sexual activity reduction due to GC symptoms/treatment
  act.stopper.prob = fmt_getenv("ACT_STOPPER_PROB"),
  # HIV risk ratio due to GC
  hiv.rgc.rr = fmt_getenv("HIV_TRANS_RR_RGC"),
  hiv.ugc.rr = fmt_getenv("HIV_TRANS_RR_UGC")
)

init <- init_msm(
  prev.ugc = 0.05,
  prev.rgc = 0.05,
  prev.pgc = 0.05
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
  transRoute_Kissing      = TRUE,  # FLAG: Toggle kissing transmission
  transRoute_Rimming      = TRUE,  # FLAG: Toggle rimming transmission
  gcUntreatedRecovDist    = "geom",
  stiScreeningProtocol    = "base",
  skip.check              = TRUE,
  cdcExposureSite_Kissing = FALSE, # FLAG: Kissing considered exposure for CDC protocol?
  tergmLite               = TRUE, # NOTE Must set to avoid error from saveout.net()
  debug_stitx             = FALSE,
  save.network		        = FALSE,
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
