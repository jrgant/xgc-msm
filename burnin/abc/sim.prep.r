library(pacman)

p_load(
  methods,
  EasyABC,
  EpiModelHIV,
  data.table,
  magrittr,
  rms,
  stringr,
  pscl
)

# Source target stat helper functions.
source(here::here("R", "target-stats.R"))

# Main Model Function ----------------------------------------------------------

xgc <- function(x) {

  set.seed(x[1])

  require(EpiModelHIV)
  require(data.table)
  require(magrittr)
  require(rms)
  require(stringr)
  require(pscl)

  tail_length <- 12
  est_path  <- here::here("est")
  netstats  <- readRDS(file.path(est_path, "netstats.Rds"))
  est       <- readRDS(file.path(est_path, "netest.Rds"))
  epistats  <- readRDS(file.path(est_path, "epistats.Rds"))

  param <- param_msm(
    # external objects
    netstats          = netstats,
    epistats          = epistats,
    # demography
    arrival.age       = 18,
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
    ugc.ntx.int       = 2,
    pgc.ntx.int       = x[7],
    ## Treated GC resolution probs c(after 1 week, after 2 weeks, after 3 weeks)
    rgc.tx.recov.pr   = c(1 - x[8], x[9], 1),
    ugc.tx.recov.pr   = c(1 - x[10], x[11], 1),
    pgc.tx.recov.pr   = c(1 - x[12], x[13], 1),
    rgc.sympt.prob    = x[14], # rectal symptom probability
    ugc.sympt.prob    = x[15], # urethral symptom probability
    pgc.sympt.prob    = x[16], # pharyngeal symptom probability
    # STI testing
    gc.sympt.seek.test.scale  = x[17],
    ## NOTE: Changed the treatment probs to be conditional on someone's seeking
    ##       STI treatment. Repeat 3 times, one for each anatomic site.
    ##       Tune asymptomatic treatment probability to achieve the overall
    ##       probability of receiving an STI test.
    gc.sympt.prob.test  = rep(1, 3),
    ## NOTE: Order of probs: rectal, urethral, pharyngeal
    gc.asympt.prob.test = c(x[18], x[19], x[20]),
    # HIV testing
    # @ORIG, prop. of MSM testing only at late-stage (AIDS)
    hiv.test.late.prob  = rep(0.25, 4),
    # HIV treatment parameters
    tt.part.supp        = rep(0.2, 4), # ORIGPARAM, partial VLS post ART start
    tt.full.supp        = rep(0.4, 4), # ORIGPARAM, full VLS w/ post ART start
    tt.dur.supp         = rep(0.4, 4),  # ORIGPARAM, durable VLS post ART start
    tx.init.prob        = c(x[21], x[22], x[23], x[24]),
    tx.halt.part.prob   = c(x[25], x[26], x[27], x[28]),
    tx.halt.full.rr     = rep(0.9, 4), # ORIGPARAM
    tx.halt.dur.rr      = rep(0.5, 4),  # ORIGPARAM
    tx.reinit.part.prob = c(x[29], x[30], x[31], x[32]), # @ORIG
    tx.reinit.full.rr   = rep(1.0, 4), # ORIGPARAM
    tx.reinit.dur.rr    = rep(1.0, 4),  # ORIGPARAM
    # Scaling parameters
    ai.acts.scale       = 1,
    oi.acts.scale       = 1,
    kiss.rate.main      = 0,
    kiss.rate.casl      = 0,
    kiss.prob.oo        = 0,
    rim.rate.main       = 0,
    rim.rate.casl       = 0,
    rim.prob.oo         = 0,
    trans.scale         = rep(1.0, 4), # ORIGPARAM
    cdc.sti.int         = 12,
    cdc.sti.hr.int      = 6,
    cond.eff            = x[33],
    cond.fail           = rep(0, 4),
    # NOTE: Change to 0 to turn off (reflect uncertainty in cond. effect)
    sti.cond.fail       = rep(0, 4),
    sti.cond.eff        = x[34],
    circ.prob           = c(0.798, 0.435, 0.600, 0.927)
  )

  init <- init_msm(
    prev.ugc = 0.01,
    prev.rgc = 0.01,
    prev.pgc = 0.01
  )

  control <- control_msm(
    # Computing options
    simno   = 1,
    nsteps  = 20,
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
    verbose.FUN       = verbose.net,
    # Epidemic simulation options
    transRoute_Kissing    = FALSE,  # FLAG: Toggle kissing transmission
    transRoute_Rimming    = FALSE,  # FLAG: Toggle rimming transmission
    gcUntreatedRecovDist  = "geom",
    stiScreeningProtocol  = "base",
    cdcExposureSite_Kissing = FALSE, # FLAG: Kissing exposure for CDC protocol?
    tergmLite     = TRUE, # NOTE Must set to avoid error from saveout.net()
    debug_stitx   = FALSE
  )

  sim <<- netsim(est, param, init, control)

  ## Calculate target stats.
  targ.hiv.prev <- mean(tail(unlist(sim$epi$i.num)) / tail(unlist(sim$epi$num)))

  targ.hiv.incid <-
    mean(tail(unlist(sim$epi$incid)) / tail(unlist(sim$epi$num))) * 52 * 100000

  ## Return vector of target stats.
  targs.out <- c(
    targ.hiv.prev,
    targ.hiv.incid,
    target_prob_agedx_byrace(sim, tailn = tail_length),
    target_prob_hivdx_byraceage(sim, tailn = tail_length),
    target_prob_vls_byraceage(sim, tailn = tail_length),
    target_prop_anatsites_tested(sim, tailn = tail_length),
    target_prob_gcpos_tested_anatsites(sim, tailn = tail_length),
    target_prob_prep_byrace(sim, tailn = tail_length)
  )

  targs.out
}


# ABC Priors and Target Stats --------------------------------------------------

# NOTE: When use_seed = TRUE abc in_test(), priors list starts with x[2].
priors <- list(
  # GC per-act transmission probability by pathways
  c("unif", 0, 0.35), # u2rgc.prob
  c("unif", 0, 0.45), # u2pgc.prob
  c("unif", 0, 0.55), # r2ugc.prob
  c("unif", 0, 0.25), # p2ugc.prob
  # Untreated infection durations
  c("unif", 2, 22), # rectal
  c("unif", 3, 27), # pharyngeal
  # Symptom probability
  c("unif", 0.06, 0.46), # rectal
  c("unif", 0.46, 0.99), # urethral
  c("unif", 0, 0.15), # pharyngeal
  # GC symptomatic testing scalar
  c("unif", 1, 10),
  # Probability of testing at asymptomatic sites in clinic
  c("unif", 0.400, 0.800),  # rectal
  c("unif", 0.497, 0.657),  # urethral
  c("unif", 0.522, 0.749),  # pharyngeal
  c("unif", 0.01, 0.11), # treated rectal gc, 1st week still infected prop.
  c("unif", 0.05, 0.95), # treated rectal gc, 2nd week resolution prop.
  c("unif", 0.01, 0.20), # treated urethral gc, 1st week still infected prop.
  c("unif", 0.05, 0.95), # treated urethral gc, 2nd week resolution prop.
  c("unif", 0.06, 0.20), # treated pharyngeal gc, 1st week still infected prop.
  c("unif", 0.05, 0.95), # treated pharyngeal gc, 2nd week resolution prop.
  # tx.init.prob (from Sam's CombPrev paper)
  c("unif", 1 / 24, 1 / 1.1),
  c("unif", 1 / 24, 1 / 1.1),
  c("unif", 1 / 24, 1 / 1.1),
  c("unif", 1 / 24, 1 / 1.1),
  # tx.halt.part.prob (from Sam's CombPrev paper)
  c("unif", 1 / 36, 1 / 1.1),
  c("unif", 1 / 36, 1 / 1.1),
  c("unif", 1 / 36, 1 / 1.1),
  c("unif", 1 / 36, 1 / 1.1),
  # tx.reinit.part.prob (from Sam's CombPrev paper)
  c("unif", 1 / 36, 1 / 1.1),
  c("unif", 1 / 36, 1 / 1.1),
  c("unif", 1 / 36, 1 / 1.1),
  c("unif", 1 / 36, 1 / 1.1),
  # cond.eff (HIV)
  c("unif", 0.8, 1),
  # sti.cond.eff
  c("unif", 0.2, 0.8)
)

targets <- c(
  # HIV TARGETS
  0.124, 514, # HIV prevalence, HIV incidence
  ## Age distribution among HIV-diagnosed
  0.095, 0.340, # Age distro among Black HIV dx'd: ages 1, 2
  0.061, 0.183, # Age distro among Hispanic HIV dx'd: ages 1, 5
  0.047, # Age distro among Other HIV dx'd: ages 1
  0.021, 0.280, 0.425, # Age dist among White HIV dx'd: ages 1, 4, 5
  ## Age- and race-specific HIV diagnosed probabilities
  0.557, 0.724, 0.859, 0.937, 0.959, # Prob. HIV dx by age, among Black HIV+
  0.488, 0.651, 0.808, 0.909, 0.947, # Prob. HIV dx by age, among Hispanic HIV+
  0.577, 0.706, 0.835, 0.926, 0.962, # Prob. HIV dx by age, among Other HIV+
  0.589, 0.726, 0.846, 0.922, 0.952, # Prob. HIV dx by age, among White HIV+
  ## Viral suppression targets
  0.471, 0.549, # Prob. VLS for Black HIV dx: 18-34, 35-65
  0.588, 0.631, # Prob. VLS for White HIV dx: 18-44, 45-65
  0.603, 0.684, # Prob. VLS for Other HIV dx: 18-34, 35-65
  0.569, 0.667, # Prob. VLS for White HIV dx: 18-24, 25-65
  ## PrEP use among MSM with indications (Finlayson 2019, MMWR)
  0.262, 0.300, 0.398, 0.424, # Past 12-mo. PrEP use (Black, Hisp, Other, White)
  # GONORRHEA TARGETS
  ## targets refer to STI testing that's not part of routine testing as
  ## part of being on PrEP
  0.657, 0.657, 0.749, # Proportion of anat sites tested in clinic (R, U, P)
  0.148, 0.079, 0.129 # Proportion of tests positive in clinic (R, U, P)
)

abc <- ABC_sequential(
  method = "Lenormand",
  model = xgc,
  prior = priors,
  nb_simul = 500,
  n_cluster = 1,
  summary_stat_target = targets,
  use_seed = TRUE,
  progress_bar = TRUE
)

saveRDS(abc, "abc_posterior.Rds")
