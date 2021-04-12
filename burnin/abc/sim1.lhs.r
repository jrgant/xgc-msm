# SETUP ------------------------------------------------------------------------

pacman::p_load(
  xgcmsm,
  EasyABC,
  methods,
  EpiModelHIV,
  data.table,
  magrittr,
  rms,
  stringr,
  pscl,
  lhs
)


# EPIMODEL SIM -----------------------------------------------------------------

netstats    <- get_est("netstats")
est         <- get_est("netest")
epistats    <- get_est("epistats")

## Start model at 1% HIV prevalence (reduce necessary burn-in time)
## The init_status_msm() function uses the initialized diag.status
## attribute to initialize HIV status.
netstats$attr$diag.status <- rbinom(
  length(netstats$attr$diag.status), 1, 0.01
)

param <- param_msm(
  # external objects
  netstats          = netstats,
  epistats          = epistats,
  # demography
  arrival.age       = 18,
  a.rate            = netstats$demog$mortrate.marginal,
  u2rgc.tprob       = Sys.getenv("U2RGC_PROB"), # urethral-to-rectal transmission probability
  u2pgc.tprob       = Sys.getenv("U2PGC_PROB"), # urethral-to-pharyngeal transmission probability
  r2ugc.tprob       = Sys.getenv("R2UGC_PROB"), # rectal-to-urethral transmission probability
  p2ugc.tprob       = Sys.getenv("P2UGC_PROB"), # pharyngeal-to-urethral transmission probability
  ## NOTE: Following tprobs used only if the kissing/rimming flags are active
  ## in control_msm.
  r2pgc.tprob       = 0, # rectal-to-pharyngeal transmission probability
  p2rgc.tprob       = 0, # pharyngeal-to-rectal transmission probability
  p2pgc.tprob       = 0, # kissing transmission probability
  rgc.ntx.int       = Sys.getenv("RECT_GC_DURAT_NOTX"),
  ugc.ntx.int       = Sys.getenv("URETH_GC_DURAT_NOTX"),
  pgc.ntx.int       = Sys.getenv("PHAR_GC_DURAT_NOTX"),
  ## Treated GC resolution probs c(after 1 week, after 2 weeks, after 3 weeks)
  rgc.tx.recov.pr   = c(1 - Sys.getenv("RECT_GC_RX_INFPR_WK1"), 0.5, 1),
  ugc.tx.recov.pr   = c(1 - Sys.getenv("URETH_GC_RX_INFPR_WK1"), 0.5, 1),
  pgc.tx.recov.pr   = c(1 - Sys.getenv("PHAR_GC_RX_INFPR_WK1"), 0.5, 1),
  rgc.sympt.prob    = Sys.getenv("RECT_GC_SYMPT_PROB"), # rectal symptom probability
  ugc.sympt.prob    = Sys.getenv("URETH_GC_SYMPT_PROB"), # urethral symptom probability
  pgc.sympt.prob    = Sys.getenv("PHAR_GC_SYMPT_PROB"), # pharyngeal symptom probability
  # STI testing
  gc.sympt.seek.test.prob = Sys.getenv("STITEST_PROB_GC_SYMPT"),
  ## NOTE: Changed the treatment probs to be conditional on someone's seeking
  ##       STI treatment. Repeat 3 times, one for each anatomic site.
  ##       Tune asymptomatic treatment probability to achieve the overall
  ##       probability of receiving an STI test.
  gc.sympt.prob.test  = rep(1, 3),
  ## NOTE: Order of probs: rectal, urethral, pharyngeal
  gc.asympt.prob.test = c(
    Sys.getenv("RECT_ASYMP_STITEST_PROB"),
    Sys.getenv("URETH_ASYMP_STITEST_PROB"),
    Sys.getenv("PHAR_ASYMP_STITEST_PROB")
  ),
  # HIV testing
  # @ORIG, prop. of MSM testing only at late-stage (AIDS)
  hiv.test.late.prob  = c(
    Sys.getenv("HIV_LATE_TESTER_PROB_BLACK"),
    Sys.getenv("HIV_LATE_TESTER_PROB_HISP"),
    Sys.getenv("HIV_LATE_TESTER_PROB_OTHER"),
    Sys.getenv("HIV_LATE_TESTER_PROB_WHITE")
  ),
  # HIV treatment parameters
  tt.part.supp        = rep(0, 4), # Partial VLS class, post ART
  tt.full.supp        = rep(1, 4), # Full VLS class, post ART
  tt.dur.supp         = rep(0, 4), # Durable VLS post ART
  tx.init.prob        = c(
    Sys.getenv("HIV_RX_INIT_PROB_BLACK"),
    Sys.getenv("HIV_RX_INIT_PROB_HISP"),
    Sys.getenv("HIV_RX_INIT_PROB_OTHER"),
    Sys.getenv("HIV_RX_INIT_PROB_WHITE")
  ),
  tx.halt.part.prob   = c(
    Sys.getenv("HIV_RX_HALT_PROB_BLACK"),
    Sys.getenv("HIV_RX_HALT_PROB_HISP"),
    Sys.getenv("HIV_RX_HALT_PROB_OTHER"),
    Sys.getenv("HIV_RX_HALT_PROB_WHITE")
  ),
  tx.halt.full.rr     = rep(1, 4), # not used
  tx.halt.dur.rr      = rep(1, 4), # not used
  tx.reinit.part.prob = c(
    Sys.getenv("HIV_RX_REINIT_PROB_BLACK"),
    Sys.getenv("HIV_RX_REINIT_PROB_HISP"),
    Sys.getenv("HIV_RX_REINIT_PROB_OTHER"),
    Sys.getenv("HIV_RX_REINIT_PROB_WHITE")
  ), # @ORIG
  tx.reinit.full.rr   = rep(1, 4),
  tx.reinit.dur.rr    = rep(1, 4),
  # Scaling parameters
  ai.acts.scale.mc    = Sys.getenv("SCALAR_AI_ACT_RATE"),
  oi.acts.scale.mc    = Sys.getenv("SCALAR_OI_ACT_RATE"),
  kiss.rate.main      = 0,
  kiss.rate.casl      = 0,
  kiss.prob.oo        = 0,
  rim.rate.main       = 0,
  rim.rate.casl       = 0,
  rim.prob.oo         = 0,
  trans.scale         = c(
    Sys.getenv("SCALAR_HIV_TRANS_PROB_BLACK"),
    Sys.getenv("SCALAR_HIV_TRANS_PROB_HISP"),
    Sys.getenv("SCALAR_HIV_TRANS_PROB_OTHER"),
    Sys.getenv("SCALAR_HIV_TRANS_PROB_WHITE")
  ),
  cdc.sti.int         = 12,
  cdc.sti.hr.int      = 6,
  cond.eff            = Sys.getenv("CONDOM_EFF_HIV"),
  ## NOTE: Change cond.fail params to 0 to turn off (reflect uncertainty
  ##       in cond. effect)
  cond.fail           = rep(0, 4),
  sti.cond.fail       = rep(0, 4),
  sti.cond.eff        = Sys.getenv("CONDOM_EFF_GC"),
  circ.prob           = netstats$inputs$circ.probs
)

init <- init_msm(
  prev.ugc = 0.01,
  prev.rgc = 0.01,
  prev.pgc = 0.01
)

control <- control_msm(
  # Computing options
  simno   = paste0("SIM1_LHS_", Sys.getenv("SIMNO")),
  nsteps  = as.numeric(Sys.getenv("NSTEPS")),
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
  transRoute_Kissing      = FALSE,  # FLAG: Toggle kissing transmission
  transRoute_Rimming      = FALSE,  # FLAG: Toggle rimming transmission
  gcUntreatedRecovDist    = "geom",
  stiScreeningProtocol    = "base",
  skip.check              = TRUE,
  cdcExposureSite_Kissing = FALSE, # FLAG: Kissing considered exposure for CDC protocol?
  tergmLite               = TRUE, # NOTE Must set to avoid error from saveout.net()
  debug_stitx             = FALSE
)

sim <- netsim(est, param, init, control)

