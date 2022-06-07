# This script calculates outcome means over the final 5 years of simulation
# and saves the results to a dataset. The dataset is used in subsequent
# analysis scripts and is much smaller than the longitudinal simulation
# output which can't be loaded for all analyses at the same time (on my
# machine).

################################################################################
## SETUP ##
################################################################################

library(pacman)

p_load(
  data.table,
  stringr,
  foreach,
  doParallel
)

## path to current directory
simpath <- here::here("inst/analysis02_screening/")

## make sure length of main and sensitivity simulations are the same
timewarn <- function(maxstep_sens, maxstep_main = max(smain65$at)) {
  ## maxstep_X INT  last time step in a simulation file
  if (maxstep_sens != maxstep_main) {
    warning("Main and sensitivity timesteps don't match!")
  } else if (maxstep_sens == maxstep_main) {
    message("Timesteps match. Do a barrel roll!")
  }
}

## labels
rlabs <- c("B", "H", "O", "W")
agelabs <- paste0("age", 1:5)
anatlabs <- c("rect", "ureth", "phar")
gclabsl <- c("rgc", "ugc", "pgc")
gclabsh <- c("rGC", "uGC", "pGC")


################################################################################
## HELPERS FOR FORMATTING IMPORTED DATA ##
################################################################################

## Extract column names from a data table given a regex pattern
getnames <- function(pattern, dt) {
  names(dt)[names(dt) %like% pattern]
}

## Calculate incidence per 100 person-years at each timestep.
## Denominator is all men in a given strata.
calc_gc_ir100 <- function(incid, num) incid / num * 5200


################################################################################
## IMPORT AND CREATE 5-YEAR SUMMARIES BY SIMULATION ID ##
################################################################################

sims <- list.files(simpath, pattern = "^epi_02.*rds")

## incomplete epi files (in progress) will return filesizes of 43
keepsims <- which(file.size(file.path(simpath, sims)) != 43)

ncores <- detectCores() - 2 # avoid using all cores on a personal machine
registerDoParallel(ncores)

## loop
epil <- foreach(i = seq_along(sims[keepsims])) %dopar% {

  slug <- stringr::str_extract(sims[i], "(?<=epi_).*(?=\\.rds)")
  fullpath <- file.path(simpath, sims[i])
  dt <- readRDS(fullpath)

  ## calculate proportion of infections due to each transpath
  gcpath_vars <- names(dt)[names(dt) %like% "(r|u|p)2"]
  gcpath_props <- paste0("prop.", gcpath_vars)

  dt[, (gcpath_props) := lapply(.SD, function(.x) .x / incid.gc),
    .SDcols = gcpath_vars
    ]

  ## loop through these columns again and replace NaN with 0
  dt[, (gcpath_props) := lapply(.SD, function(.x) replace(.x, is.na(.x), 0)),
     .SDcols = gcpath_props]

  ## calculate race/ethnicity-specific and age-specific gonorrhea incidence
  dt[, ":=" (
    ir100.pop.B.rgc = calc_gc_ir100(incid.B.rgc, num.B),
    ir100.pop.H.rgc = calc_gc_ir100(incid.H.rgc, num.H),
    ir100.pop.O.rgc = calc_gc_ir100(incid.O.rgc, num.O),
    ir100.pop.W.rgc = calc_gc_ir100(incid.W.rgc, num.W),
    ir100.pop.B.ugc = calc_gc_ir100(incid.B.ugc, num.B),
    ir100.pop.H.ugc = calc_gc_ir100(incid.H.ugc, num.H),
    ir100.pop.O.ugc = calc_gc_ir100(incid.O.ugc, num.O),
    ir100.pop.W.ugc = calc_gc_ir100(incid.W.ugc, num.W),
    ir100.pop.B.pgc = calc_gc_ir100(incid.B.pgc, num.B),
    ir100.pop.H.pgc = calc_gc_ir100(incid.H.pgc, num.H),
    ir100.pop.O.pgc = calc_gc_ir100(incid.O.pgc, num.O),
    ir100.pop.W.pgc = calc_gc_ir100(incid.W.pgc, num.W),
    ir100.pop.age1.rgc = calc_gc_ir100(incid.age.1.rgc, num.age.1),
    ir100.pop.age2.rgc = calc_gc_ir100(incid.age.2.rgc, num.age.2),
    ir100.pop.age3.rgc = calc_gc_ir100(incid.age.3.rgc, num.age.3),
    ir100.pop.age4.rgc = calc_gc_ir100(incid.age.4.rgc, num.age.4),
    ir100.pop.age5.rgc = calc_gc_ir100(incid.age.5.rgc, num.age.5),
    ir100.pop.age1.ugc = calc_gc_ir100(incid.age.1.ugc, num.age.1),
    ir100.pop.age2.ugc = calc_gc_ir100(incid.age.2.ugc, num.age.2),
    ir100.pop.age3.ugc = calc_gc_ir100(incid.age.3.ugc, num.age.3),
    ir100.pop.age4.ugc = calc_gc_ir100(incid.age.4.ugc, num.age.4),
    ir100.pop.age5.ugc = calc_gc_ir100(incid.age.5.ugc, num.age.5),
    ir100.pop.age1.pgc = calc_gc_ir100(incid.age.1.pgc, num.age.1),
    ir100.pop.age2.pgc = calc_gc_ir100(incid.age.2.pgc, num.age.2),
    ir100.pop.age3.pgc = calc_gc_ir100(incid.age.3.pgc, num.age.3),
    ir100.pop.age4.pgc = calc_gc_ir100(incid.age.4.pgc, num.age.4),
    ir100.pop.age5.pgc = calc_gc_ir100(incid.age.5.pgc, num.age.5)
  )]

  ## calculate infection totals
  dt[, ":=" (
    sum_incid.gc = sum(incid.gc),
    sum_incid.rgc = sum(incid.rgc),
    sum_incid.ugc = sum(incid.ugc),
    sum_incid.pgc = sum(incid.pgc)
  ), by = simid]

  ## calculate means over last 5 years
  sumcols <- names(dt[, -c("simid", "at")])
  dts <- dt[, lapply(.SD, mean), .SDcols = sumcols, by = simid]
  dts[, analysis := slug]

  ## extract some metadata about the original data sets
  meta <- dt[, .(min_at = min(at), max_at = max(at)), by = simid]

  ## pull in the metadata for  each simid
  out <- merge(dts, meta, by = "simid")
  out
}


################################################################################
## CHECK BOUND DATA ##
################################################################################

epi <- rbindlist(epil, fill = T)
epi[, .N, analysis]
epi[, .N, .(analysis, min_at, max_at)]


################################################################################
## WRITE DATA ##
################################################################################

saveRDS(epi, file.path(simpath, "epi_5yr_BOUND.rds"))
