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


# INITIAL PRIOR RANGES ---------------------------------------------------------

# convenience function to add a uniform prior to the prior list
pvec <- function(ll, ul, name) c("unif", ll, ul, name)

priors <- list(
  # GC per-act transmission probability by pathways
  pvec(0, 0.35, "U2RGC_PROB"), # u2rgc.prob
  pvec(0, 0.45, "U2PGC_PROB"), # u2pgc.prob
  pvec(0, 0.55, "R2UGC_PROB"), # r2ugc.prob
  pvec(0, 0.25, "P2UGC_PROB"), # p2ugc.prob
  # Untreated infection durations
  pvec(2, 22, "RECT_GC_DURAT_NOTX"), # rectal
  pvec(1, 24, "URETH_GC_DURAT_NOTX"), # urethral
  pvec(3, 20, "PHAR_GC_DURAT_NOTX"), # pharyngeal
  # Probability of still being GC-infected after x weeks of treatment
  pvec(0.01, 0.11, "RECT_GC_RX_INFPR_WK1"), # treated rectal gc, 1st week still infected prop.
  pvec(0.01, 0.20, "URETH_GC_RX_INFPR_WK1"), # treated urethral gc, 1st week still infected prop.
  pvec(0.06, 0.20, "PHAR_GC_RX_INFPR_WK1"), # treated pharyngeal gc, 1st week still infected prop.
  # Symptom probability
  pvec(0.06, 0.46, "RECT_GC_SYMPT_PROB"), # rectal
  pvec(0.46, 0.99, "URETH_GC_SYMPT_PROB"), # urethral
  pvec(0.00, 0.15, "PHAR_GC_SYMPT_PROB"), # pharyngeal
  # GC symptomatic, weekly probability of testing
  pvec(0.5, 1, "STITEST_PROB_GC_SYMPT"),
  # Probability of testing at asymptomatic sites in clinic
  pvec(0.400, 0.800, "RECT_ASYMP_STITEST_PROB"),  # rectal
  pvec(0.497, 0.657, "URETH_ASYMP_STITEST_PROB"),  # urethral
  pvec(0.522, 0.749, "PHAR_ASYMP_STITEST_PROB"),  # pharyngeal
  # HIV late-tester probability (priors from Singh 2017, MMWR)
  pvec(0.160, 0.419, "HIV_LATE_TESTER_PROB_BLACK"),
  pvec(0.202, 0.391, "HIV_LATE_TESTER_PROB_HISP"),
  pvec(0.207, 0.369, "HIV_LATE_TESTER_PROB_OTHER"),
  pvec(0.222, 0.377, "HIV_LATE_TESTER_PROB_WHITE"),
  # tx.init.prob
  pvec(1 / 62, 1 / 27, "HIV_RX_INIT_PROB_BLACK"),
  pvec(1 / 62, 1 / 27, "HIV_RX_INIT_PROB_HISP"),
  pvec(1 / 62, 1 / 27, "HIV_RX_INIT_PROB_OTHER"),
  pvec(1 / 62, 1 / 27, "HIV_RX_INIT_PROB_WHITE"),
  # tx.halt.part.prob (Singh 2017, MMWR)
  pvec(1 - (1 - 0.464)^(1/52), 1 - (1 - 0.464)^(1/52), "HIV_RX_HALT_PROB_BLACK"),
  pvec(1 - (1 - 0.416)^(1/52), 1 - (1 - 0.416)^(1/52), "HIV_RX_HALT_PROB_HISP"),
  pvec(1 - (1 - 0.366)^(1/52), 1 - (1 - 0.366)^(1/52), "HIV_RX_HALT_PROB_OTHER"),
  pvec(1 - (1 - 0.406)^(1/52), 1 - (1 - 0.406)^(1/52), "HIV_RX_HALT_PROB_WHITE"),
  # tx.reinit.part.prob
  pvec(0.01, 0.20, "HIV_RX_REINIT_PROB_BLACK"),
  pvec(0.01, 0.20, "HIV_RX_REINIT_PROB_HISP"),
  pvec(0.01, 0.20, "HIV_RX_REINIT_PROB_OTHER"),
  pvec(0.01, 0.20, "HIV_RX_REINIT_PROB_WHITE"),
  # sex act scalars
  pvec(1, 1, "SCALAR_AI_ACT_RATE"), # ai.acts.scale
  pvec(1, 1, "SCALAR_OI_ACT_RATE"), # oi.acts.scale
  # HIV transmission prob. scalars
  pvec(0.01, 5, "SCALAR_HIV_TRANS_PROB_BLACK"), # black
  pvec(0.01, 5, "SCALAR_HIV_TRANS_PROB_HISP"), # hispanic
  pvec(0.01, 5, "SCALAR_HIV_TRANS_PROB_OTHER"), # other
  pvec(0.01, 5, "SCALAR_HIV_TRANS_PROB_WHITE"), # white
  # condom efficacy
  pvec(0.6, 1, "CONDOM_EFF_HIV"), # cond.eff (HIV), 1 - relative risk, condom use vs. no
  pvec(0.7, 1, "CONDOM_EFF_GC")   # sti.cond.eff, relative risk, condom use vs. no
)


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

abcdir <- here::here("burnin", "abc")

if (priornames_ok & priorvals_ok) {
  if (!file.exists(file.path(abcdir, "sim1"))) {
    dir.create(file.path(abcdir, "sim1"))
  }
  saveRDS(lhs_real, file = file.path(abcdir, "sim1", "lhs_sim1.rds"))
} else {
  stop("Check prior specification.")
}