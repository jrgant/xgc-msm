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
  pscl
)


# CALIBRATION TARGETS ----------------------------------------------------------

targl <- readRDS(here::here("est", "caltargets.Rds"))
targl_sel <- targl[c(1:2,6:7)]

targets <- unlist(
  sapply(seq_along(targl_sel), function(x) unlist(targl_sel[x]))
)


# MAIN MODEL FUNCTION ----------------------------------------------------------

xgc <- function(x) {

  require(xgcmsm)

  # Source target stat helper functions.
  set.seed(x[1])

  # Load libraries in function so they're seen on all parallel cores.
  require(EpiModelHIV)
  require(data.table)
  require(magrittr)
  require(rms)
  require(stringr)
  require(pscl)

  tail_length <- as.numeric(Sys.getenv("TAIL_LENGTH"))
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
    u2rgc.tprob       = x[2], # urethral-to-rectal transmission probability
    u2pgc.tprob       = x[3], # urethral-to-pharyngeal transmission probability
    r2ugc.tprob       = x[4], # rectal-to-urethral transmission probability
    p2ugc.tprob       = x[5], # pharyngeal-to-urethral transmission probability
    ## NOTE: Following tprobs used only if the kissing/rimming flags are active
    ## in control_msm.
    r2pgc.tprob       = 0, # rectal-to-pharyngeal transmission probability
    p2rgc.tprob       = 0, # pharyngeal-to-rectal transmission probability
    p2pgc.tprob       = 0, # kissing transmission probability
    rgc.ntx.int       = x[6],
    ugc.ntx.int       = x[7],
    pgc.ntx.int       = x[8],
    ## Treated GC resolution probs c(after 1 week, after 2 weeks, after 3 weeks)
    rgc.tx.recov.pr   = c(1 - x[9], 0.5, 1),
    ugc.tx.recov.pr   = c(1 - x[10], 0.5, 1),
    pgc.tx.recov.pr   = c(1 - x[11], 0.5, 1),
    rgc.sympt.prob    = x[12], # rectal symptom probability
    ugc.sympt.prob    = x[13], # urethral symptom probability
    pgc.sympt.prob    = x[14], # pharyngeal symptom probability
    # STI testing
    gc.sympt.seek.test.prob = x[15],
    ## NOTE: Changed the treatment probs to be conditional on someone's seeking
    ##       STI treatment. Repeat 3 times, one for each anatomic site.
    ##       Tune asymptomatic treatment probability to achieve the overall
    ##       probability of receiving an STI test.
    gc.sympt.prob.test  = rep(1, 3),
    ## NOTE: Order of probs: rectal, urethral, pharyngeal
    gc.asympt.prob.test = c(x[16], x[17], x[18]),
    # HIV testing
    # @ORIG, prop. of MSM testing only at late-stage (AIDS)
    hiv.test.late.prob  = c(x[19], x[20], x[21], x[22]),
    # HIV treatment parameters
    tt.part.supp        = rep(0, 4), # Partial VLS class, post ART
    tt.full.supp        = rep(1, 4), # Full VLS class, post ART
    tt.dur.supp         = rep(0, 4), # Durable VLS post ART
    tx.init.prob        = c(x[23], x[24], x[25], x[26]),
    tx.halt.part.prob   = c(x[27], x[28], x[29], x[30]),
    tx.halt.full.rr     = rep(1, 4),
    tx.halt.dur.rr      = rep(1, 4),
    tx.reinit.part.prob = c(x[31], x[32], x[33], x[34]), # @ORIG
    tx.reinit.full.rr   = rep(1, 4),
    tx.reinit.dur.rr    = rep(1, 4),
    # Scaling parameters
    ai.acts.scale.mc    = x[35],
    oi.acts.scale.mc    = x[36],
    kiss.rate.main      = 0,
    kiss.rate.casl      = 0,
    kiss.prob.oo        = 0,
    rim.rate.main       = 0,
    rim.rate.casl       = 0,
    rim.prob.oo         = 0,
    trans.scale         = c(x[37], x[38], x[39], x[40]), # ORIGPARAM
    cdc.sti.int         = 12,
    cdc.sti.hr.int      = 6,
    cond.eff            = x[41],
    ## NOTE: Change cond.fail params to 0 to turn off (reflect uncertainty
    ##       in cond. effect)
    cond.fail           = rep(0, 4),
    sti.cond.fail       = rep(0, 4),
    sti.cond.eff        = x[42],
    circ.prob           = netstats$inputs$circ.probs
  )

  init <- init_msm(
    prev.ugc = 0.01,
    prev.rgc = 0.01,
    prev.pgc = 0.01
  )

  control <- control_msm(
    # Computing options
    simno   = 1,
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
    transRoute_Kissing    = FALSE,  # FLAG: Toggle kissing transmission
    transRoute_Rimming    = FALSE,  # FLAG: Toggle rimming transmission
    gcUntreatedRecovDist  = "geom",
    stiScreeningProtocol  = "base",
    skip.check = TRUE,
    cdcExposureSite_Kissing = FALSE, # FLAG: Kissing exposure for CDC protocol?
    tergmLite     = TRUE, # NOTE Must set to avoid error from saveout.net()
    debug_stitx   = FALSE
  )

  sim <- netsim(est, param, init, control)

  ## Calculate target stats.
  targ.hiv.prev <-
    mean(
      tail(unlist(sim$epi$i.num), n = tail_length) /
      tail(unlist(sim$epi$num), n = tail_length)
    )

  targ.hiv.incid <-
    mean(
      tail(unlist(sim$epi$incid), n = tail_length) /
      tail(unlist(sim$epi$s.num), n = tail_length)
    ) * 52 * 100000

  ## Return vector of target stats.
  targs.out <- c(
    hiv.incid.100k = targ.hiv.incid,
    hiv.prev = targ.hiv.prev,
    ## target_prob_hivdx(sim, tail_length, "race")
    ## target_prob_vls(sim, tail_length, "race")
    ## target_prob_hivdx(sim, tail_length, "age"),
    ## target_prob_vls(sim, tail_length, "age"),
    ## target_prob_prep_byrace(sim, tailn = tail_length),
    target_prop_anatsites_tested(sim, tailn = tail_length),
    target_prob_gcpos_tested_anatsites(sim, tailn = tail_length)
  )

  unname(targs.out)
}


# ABC PRIORS AND TARGET STATS --------------------------------------------------

# convenience function to add a uniform prior to the prior list
pvec <- function(ll, ul) c("unif", ll, ul)

# NOTE: When use_seed = TRUE abc in_test(), priors list starts with x[2].

# NOTE: This file contains a simplified set of prior inputs, where viral
#       suppression trajectory classes (tt.part, tt.full, tt.dur) are turned off
#       To achieved this, the following steps apply:
#         1) All agents are classified as tt.part
#         2) Treatment initiation, halting, and reinitiation probabilities are
#            fed to the corresponding parameters for the partial suppression
#            class, which would normally represent the baseline probability
#            in a model with trajectory classes.
#         3) Relative rates for the full and durable suppression classes are
#            set to 1, so that those assigned to tt.full (i.e., all MSM)
#            inherit the .part probabilities in step 2.
priors <- list(
  # GC per-act transmission probability by pathways
  pvec(0, 0.35), # u2rgc.prob
  pvec(0, 0.45), # u2pgc.prob
  pvec(0, 0.55), # r2ugc.prob
  pvec(0, 0.25), # p2ugc.prob
  # Untreated infection durations
  pvec(2, 22), # rectal
  pvec(1, 24), # urethral
  pvec(3, 20), # pharyngeal
  # Probability of still being GC-infected after x weeks of treatment
  pvec(0.01, 0.11), # treated rectal gc, 1st week still infected prop.
  pvec(0.01, 0.20), # treated urethral gc, 1st week still infected prop.
  pvec(0.06, 0.20), # treated pharyngeal gc, 1st week still infected prop.
  # Symptom probability
  pvec(0.06, 0.46), # rectal
  pvec(0.46, 0.99), # urethral
  pvec(0.00, 0.15), # pharyngeal
  # GC symptomatic, weekly probability of testing
  pvec(0.5, 1),
  # Probability of testing at asymptomatic sites in clinic
  pvec(0.400, 0.800),  # rectal
  pvec(0.497, 0.657),  # urethral
  pvec(0.522, 0.749),  # pharyngeal
  # HIV late-tester probability (priors from Singh 2017, MMWR)
  pvec(0.160, 0.419),
  pvec(0.202, 0.391),
  pvec(0.207, 0.369),
  pvec(0.222, 0.377),
  # tx.init.prob
  pvec(1 / 62, 1 / 27),
  pvec(1 / 62, 1 / 27),
  pvec(1 / 62, 1 / 27),
  pvec(1 / 62, 1 / 27),
  # tx.halt.part.prob (Singh 2017, MMWR)
  pvec(1 - (1 - 0.464)^(1/52), 1 - (1 - 0.464)^(1/52)),
  pvec(1 - (1 - 0.416)^(1/52), 1 - (1 - 0.416)^(1/52)),
  pvec(1 - (1 - 0.366)^(1/52), 1 - (1 - 0.366)^(1/52)),
  pvec(1 - (1 - 0.406)^(1/52), 1 - (1 - 0.406)^(1/52)),
  # tx.reinit.part.prob -- TODO: These are way different than Sam's CombPrev. Check plausibility.
  pvec(0.0001, 0.20),
  pvec(0.0001, 0.20),
  pvec(0.0001, 0.20),
  pvec(0.0001, 0.20),
  # ai.acts.scale
  pvec(1, 1),
  # oi.acts.scale
  pvec(1, 1),
  # HIV transmission prob. scalar
  pvec(0.01, 10), # black
  pvec(0.01, 10), # hispanic
  pvec(0.01, 10), # other
  pvec(0.01, 10), # white
  # cond.eff (HIV), 1 - relative risk, condom use vs. no
  pvec(0.6, 1),
  # sti.cond.eff, relative risk, condom use vs. no
  pvec(0.7, 1)
)


# RUN ABC ----------------------------------------------------------------------

abc <- ABC_sequential(
  method                = "Lenormand",
  model                 = xgc,
  prior                 = priors,
  nb_simul              = as.numeric(Sys.getenv("SIMS_BELOW_TOL")),
  n_cluster             = as.numeric(Sys.getenv("SLURM_NPROCS")),
  summary_stat_target   = targets,
  use_seed              = TRUE,
  inside_prior          = TRUE,
  progress_bar          = FALSE,
  verbose               = TRUE
)


# WRITE ABC RESULTS ------------------------------------------------------------
saveRDS(
  abc,
  paste0(
    "abc_posterior_",
    format(Sys.time(), "%Y-%m-%d_%H:%M"),
    "_JOBID-", Sys.getenv("SLURM_JOBID"),
    "_SIMS-", Sys.getenv("SIMS_BELOW_TOL"),
    "_TAIL-", Sys.getenv("TAIL_LENGTH"),
    "_STEPS-", Sys.getenv("NSTEPS"),
    ".Rds"
  )
)
