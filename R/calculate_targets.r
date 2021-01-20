#' Calculating Target Statistics
#'
#' These functions extract and calculate simulated outcomes
#' based on the calibration target statistics specified.
#'
#' @param dat An object of class netsim.
#' @param tailn The number of weeks (always at the tail) to extract and average over.

#' @describeIn calculate_targets Probability of being diagnosed with HIV among those with HIV, by race and age.
#' @export
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

    out1 <- mean(num / den)
    out <- ifelse(is.na(out1), 0, out1)
    names(out) <- paste0("prob.hivdx.", df$race[x], ".age", df$age.grp[x])
    out
  })
}

#' @describeIn calculate_targets Calculate age distribution among
#' those diagnosed with HIV, by race/ethnicity.
#' @export
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

      t <- mean(num / den)
      names(t) <- paste0("prob.agedx.hiv.", rslug, ".age", x[[1 + z]])
      t
    })
  })

  out <- unlist(targs)
  out
}

#' @describeIn calculate_targets Calculate probability of viral suppression among those with diagnosed HIV, by race/ethnicity and age.
#' @export
target_prob_vls_byraceage <- function(dat, tailn) {
  vsupp.targs <- names(dat$epi)[grepl("cct", names(dat$epi))]
  sapply(vsupp.targs, function(x) {
    vls.prob <- tail(unlist(dat$epi[x]), n = tailn)
    mean(vls.prob)
  })
}

#' @describeIn calculate_targets Calculate probability of PrEP use among indicated agents, by race/ethnicity.
#' @export
target_prob_prep_byrace <- function(dat, tailn) {
  r <- unique(raceage_grid()$race)

  sapply(seq_along(r), function(x) {
    prep.cov <-
      mean(
        tail(unlist(dat$epi[paste0("prepCurr.", r[x])]), n = tailn) /
        tail(unlist(dat$epi[paste0("prepElig.", r[x])]), n = tailn)
      )

    out <- ifelse(is.nan(prep.cov), 0, prep.cov)
    names(out) <- paste0("prep.cov.", r[x])
    out
  })
}

#' @describeIn calculate_targets Return proportion of anatomic sites tested among MSM seeking testing at the clinic.
#' @export
target_prop_anatsites_tested <- function(dat, tailn) {
  vars <- names(dat$epi)[grepl("^prop", names(dat$epi))]
  sapply(vars, function(x) {
    mean(tail(unlist(dat$epi[x]), n = tailn))
  })
}

#' @describeIn calculate_targets Return probability of being GC+ at tested anatomic sites.
#' @export
target_prob_gcpos_tested_anatsites <- function(dat, tailn) {
  vars <- names(dat$epi)[grepl("^prob", names(dat$epi))]
  sapply(vars, function(x) {
    mean(tail(unlist(dat$epi[x]), n = tailn))
  })
}
