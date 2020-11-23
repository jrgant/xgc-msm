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

# Main Model Function ----------------------------------------------------------

xgc <- function(x) {

  set.seed(x[1])

  require(EpiModelHIV)
  require(data.table)
  require(magrittr)
  require(rms)
  require(stringr)
  require(pscl)

  est_path  <- here::here("est")
  netstats  <- readRDS(file.path(est_path, "netstats.Rds"))
  est       <- readRDS(file.path(est_path, "netest.Rds"))
  epistats  <- readRDS(file.path(est_path, "epistats.Rds"))

  param <- param_msm(
    # external objects
    netstats          = netstats,
    epistats          = epistats,
    # demography
    a.rate            = 0.00052,
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
    rgc.sympt.prob    = x[8], # rectal symptom probability
    ugc.sympt.prob    = x[9], # urethral symptom probability
    pgc.sympt.prob    = x[10], # pharyngeal symptom probability
    # STI testing
    gc.sympt.seek.test.scale  = x[11],
    ## NOTE: Changed the treatment probs to be conditional on someone's seeking
    ##       STI treatment. Repeat 3 times, one for each anatomic site.
    ##       Tune asymptomatic treatment probability to achieve the overall
    ##       probability of receiving an STI test.
    gc.sympt.prob.test  = rep(1, 3),
    ## NOTE: Order of probs: rectal, urethral, pharyngeal
    gc.asympt.prob.test = c(x[12], x[13], x[14]),
    # HIV testing
    hiv.test.rate       = c(x[15], x[16], x[17], x[18]), # by race/eth
    # @ORIG, prop. of MSM testing only at late-stage (AIDS)
    hiv.test.late.prob  = rep(0.25, 4),
    # HIV treatment parameters
    tt.part.supp        = rep(0.2, 4), # ORIGPARAM, partial VLS post ART start
    tt.full.supp        = rep(0.4, 4), # ORIGPARAM, full VLS w/ post ART start
    tt.dur.supp         = rep(0.4, 4),  # ORIGPARAM, durable VLS post ART start
    tx.init.prob        = c(x[19], x[20], x[21], x[22]),
    tx.halt.part.prob   = c(x[23], x[24], x[25], x[26]),
    tx.halt.full.rr     = rep(0.9, 4), # ORIGPARAM
    tx.halt.dur.rr      = rep(0.5, 4),  # ORIGPARAM
    tx.reinit.part.prob = c(x[27], x[28], x[29], x[30]), # @ORIG
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
    cond.eff            = 0.95,
    cond.fail           = rep(0, 4),
    # NOTE: Change to 0 to turn off (reflect uncertainty in cond. effect)
    sti.cond.fail       = rep(0, 4),
    sti.cond.eff        = x[31],
    circ.prob           = c(0.798, 0.435, 0.600, 0.927)
  )

  init <- init_msm(
    prev.ugc = 0.0,
    prev.rgc = 0.0,
    prev.pgc = 0.0
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

  sim <- netsim(est, param, init, control)

  ## Calculate target stats.
  targ.hiv.prev <- mean(tail(unlist(sim$epi$i.num)) / tail(unlist(sim$epi$num)))

  targ.hiv.incid <-
    mean(tail(unlist(sim$epi$incid)) / tail(unlist(sim$epi$num))) * 52 * 100000

  i.num.dx.B <- tail(unlist(sim$epi$i.num.dx.B))

  targ.age1.among.B.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.B.age1)) / i.num.dx.B)

  targ.age2.among.B.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.B.age2)) / i.num.dx.B)

  targ.age3.among.B.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.B.age3)) / i.num.dx.B)

  targ.age4.among.B.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.B.age4)) / i.num.dx.B)

  targ.age5.among.B.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.B.age5)) / i.num.dx.B)


  i.num.dx.H <- tail(unlist(sim$epi$i.num.dx.H))

  targ.age1.among.H.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.H.age1)) / i.num.dx.H)

  targ.age2.among.H.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.H.age2)) / i.num.dx.H)

  targ.age3.among.H.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.H.age3)) / i.num.dx.H)

  targ.age4.among.H.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.H.age4)) / i.num.dx.H)

  targ.age5.among.H.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.H.age5)) / i.num.dx.H)


  i.num.dx.O <- tail(unlist(sim$epi$i.num.dx.O))

  targ.age1.among.O.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.O.age1)) / i.num.dx.O)

  targ.age2.among.O.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.O.age2)) / i.num.dx.O)

  targ.age3.among.O.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.O.age3)) / i.num.dx.O)

  targ.age4.among.O.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.O.age4)) / i.num.dx.O)

  targ.age5.among.O.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.O.age5)) / i.num.dx.O)


  i.num.dx.W <- tail(unlist(sim$epi$i.num.dx.W))

  targ.age1.among.W.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.W.age1)) / i.num.dx.W)

  targ.age2.among.W.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.W.age2)) / i.num.dx.W)

  targ.age3.among.W.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.W.age3)) / i.num.dx.W)

  targ.age4.among.W.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.W.age4)) / i.num.dx.W)

  targ.age5.among.W.dx <-
    mean(tail(unlist(sim$epi$i.num.dx.W.age5)) / i.num.dx.W)

  ## Return vector of target stats.
  targs.out <- c(
    targ.hiv.prev, targ.hiv.incid,
    targ.age1.among.B.dx,
    targ.age2.among.B.dx,
    targ.age3.among.B.dx,
    targ.age4.among.B.dx,
    targ.age5.among.B.dx,
    targ.age1.among.H.dx,
    targ.age2.among.H.dx,
    targ.age3.among.H.dx,
    targ.age4.among.H.dx,
    targ.age5.among.H.dx,
    targ.age1.among.O.dx,
    targ.age2.among.O.dx,
    targ.age3.among.O.dx,
    targ.age4.among.O.dx,
    targ.age5.among.O.dx,
    targ.age1.among.W.dx,
    targ.age2.among.W.dx,
    targ.age3.among.W.dx,
    targ.age4.among.W.dx,
    targ.age5.among.W.dx
  )

  targs.out

}


# ABC Priors and Target Stats --------------------------------------------------

# NOTE: When use_seed = TRUE abc in_test(), priors list starts with x[2].
priors <- list(
  # GC per-act transmission probability by pathways
  c("unif", 0.35, 0.35), # u2rgc.prob
  c("unif", 0.45, 0.45), # u2pgc.prob
  c("unif", 0.55, 0.55), # r2ugc.prob
  c("unif", 0.25, 0.25), # p2ugc.prob
  # Untreated infection durations
  c("unif", 2, 2), # rectal
  c("unif", 3, 3), # pharyngeal
  # Symptom probability
  c("unif", 0.06, 0.46), # rectal
  c("unif", 0.46, 0.99), # urethral
  c("unif", 0, 0.15), # pharyngeal
  # GC symptomatic testing scalar
  c("unif", 3, 3),
  # Probability of testing at asymptomatic sites in clinic
  c("unif", 0.85, 0.85),  # rectal
  c("unif", 0.85, 0.85),  # urethral
  c("unif", 0.85, 0.85),  # pharyngeal
  # HIV test rate (from Sam's CombPrev paper)
  c("unif", 0.001, 0.02),
  c("unif", 0.001, 0.02),
  c("unif", 0.001, 0.02),
  c("unif", 0.001, 0.02),
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
  # sti.cond.eff
  c("unif", 0.8, 0.8)
)

targets <- c(
  0.124, 514, # HIV prevalence, HIV incidence
  0.095, 0.340, 0.193, 0.191, 0.182, # Age distro among Black HIV-diagnosed
  0.061, 0.265, 0.244, 0.248, 0.183, # Age distro among Hispanic HIV-diagnosed
  0.047, 0.254, 0.220, 0.252, 0.226, # Age distro among Other HIV-diagnosed
  0.021, 0.125, 0.149, 0.280, 0.425, # Age distro among White HIV-diagnosed
  0.188, 0.212, 0.309, 0.271 # Past 12-mo. PrEP use (Black, Hisp, Other, White)
)

abc_test <- ABC_sequential(
  method = "Lenormand",
  model = xgc,
  prior = priors,
  nb_simul = 100,
  n_cluster = 6,
  summary_stat_target = targets,
  use_seed = TRUE,
  progress_bar = TRUE
)

saveRDS(abc_test, "C:/Users/jason/Downloads/abc_test.Rds")
