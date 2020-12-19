raceage_grid <- function() {
  dt <- as.data.table(expand.grid(race = c("B", "H", "O", "W"), age.grp = 1:5))
  setkeyv(dt, c("race", "age.grp"))
  dt
}

# Probability of being diagnosed with HIV among those with HIV, by race and age.
target_prob_hivdx_byraceage <- function(dat, tailn) {
  df <- raceage_grid()

  sapply(seq_len(nrow(df)), function(x) {

    num <- tail(
      unlist(dat$epi[paste0("i.num.dx.", df$race[x], ".age", df$age.grp[x])]),
      n = tailn
    )

    den <- tail(
      unlist(dat$epi[paste0("i.num.", df$race[x], ".age", df$age.grp[x])]),
      n = tailn
    )

    out <- mean(num / den)
    out
  })
}

# Calculate age distribution among those diagnosed with HIV, by race/ethnicity.
target_prob_agedx_byrace <- function(dat, tailn) {

  dl <- list(
    list("B", 1, 2),
    list("H", 1, 5),
    list("O", 1),
    list("W", 1, 4, 5)
    )

  targs <- sapply(dl, function(x) {
    rslug <- x[[1]]
    n_ages <- length(x) - 1

    sapply(seq_len(n_ages), function(z) {

      num <- tail(
        unlist(dat$epi[paste0("i.num.dx.", rslug, ".age", x[[1 + z]])]),
        n = tailn
      )

      den <- tail(
        unlist(dat$epi[paste0("i.num.dx.", rslug)]),
        n = tailn
      )

      mean(num / den)
    })
  })

  out <- unlist(targs)
  out
}

# Calculate probability of viral suppression among those with diagnosed HIV,
# by race/ethnicity and age.
target_prob_vls_byraceage <- function(dat, tailn) {

  vsupp.targs <- names(dat$epi)[grepl("cct", names(dat$epi))]

  sapply(vsupp.targs, function(x) {
    vls.prob <- tail(unlist(dat$epi[x]), n = tailn)
    mean(vls.prob)
  })
}

# Calculate probability of PrEP use among indicated agents, by race/ethnicity.
target_prob_prep_byrace <- function(dat, tailn) {
  r <- unique(raceage_grid()$race)

  sapply(seq_along(r), function(x) {
    prep.cov <-
      tail(unlist(dat$epi[paste0("prepCurr.", r[x])]), n = tailn) /
      tail(unlist(dat$epi[paste0("prepElig.", r[x])]), n = tailn)

    ifelse(is.nan(prep.cov), 0, mean(prep.cov))
  })
}

# Return proportion of anatomic sites tested among MSM seeking testing at
# the clinic.
target_prop_anatsites_tested <- function(dat, tailn) {
  vars <- names(dat$epi)[grepl("^prop", names(dat$epi))]
  sapply(vars, function(x) {
    mean(tail(unlist(dat$epi[x]), n = tailn))
  })
}

# Return probability of being GC+ at tested anatomic sites.
target_prob_gcpos_tested_anatsites <- function(dat, tailn) {
  vars <- names(dat$epi)[grepl("^prob", names(dat$epi))]
  sapply(vars, function(x) {
    mean(tail(unlist(dat$epi[x]), n = tailn))
  })
}

# Helper function to use as input into the xgc() function
# defined in burnin/abc/sim.prep.r. Intended for quick tests.
test_xgc <- function(seed = 123, priorlist = priors) {
  pl <- sapply(priorlist, function(x) {
    l <- unlist(x)
    runif(1, as.numeric(l[2]), as.numeric(l[3]))
  })
  c(seed, pl)
}

# Helper function to set temporary environment variable.
setenv <- function(nsteps = 200,
                   tail_length = 52,
                   sims_below_tol = 50,
                   slurm_nprocs = Sys.getenv("SLURM_NPROCS")) {
  Sys.setenv(
    "NSTEPS"          = nsteps,
    "TAIL_LENGTH"     = tail_length,
    "SIMS_BELOW_TOL"  = sims_below_tol,
    "SLURM_NPROCS"    = slurm_nprocs
  )
}
