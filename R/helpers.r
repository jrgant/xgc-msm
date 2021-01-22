#' Helper functions
#'
#' These functions provide easy access to demographic distributions
#' in the target population as well as helper functions to generate lookup
#' tables for different combinations of parameter levels.
#'
#' @import data.table

#' @describeIn helpers Produce a lookup table of all possible combinations of
#' race/ethnicity and age group.
#' @export
raceage_grid <- function() {
  dt <- as.data.table(expand.grid(race = c("B", "H", "O", "W"), age.grp = 1:5))
  setkeyv(dt, c("race", "age.grp"))
  dt
}


#' @describeIn helpers Helper function to use as input into the xgc() function
# defined in burnin/abc/sim.prep.r. Intended for quick tests.
#' @export
test_xgc <- function(seed = 123, priorlist = priors) {
  pl <- sapply(priorlist, function(x) {
    l <- unlist(x)
    runif(1, as.numeric(l[2]), as.numeric(l[3]))
  })
  c(seed, pl)
}


#' @describeIn helpers Helper function to set temporary environment variable.
#' @export
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


#' @describeIn helpers Retrieve epistats object.
#' @export
#' @importFrom here here
get_est <- function(string = c("epistats", "netest", "netstats")) {
  readRDS(here::here("est", paste0(string, ".Rds")))
}
