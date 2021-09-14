################################################################################
## SETUP ##
################################################################################

library(data.table)

rates <- c(5, 10)
probs <- c(0.5, 1)


################################################################################
## KISSING AND RIMMING SENSITIVITIES ##
################################################################################

## k, kiss
## r, rim
## m, main
## c, casual
## o, one-time

kissX_rim0 <- expand.grid(
  km = rates, kc = rates, ko = probs,
  rm = 0, rc = 0, ro = 0
)

rimX_kiss0 <- expand.grid(
  km = 0, kc = 0, ko = 0,
  rm = rates, rc = rates, ro = probs
)

kisslabs <- c(
  paste("KISS_RATE", c("MAIN", "CASUAL"), "ALTPARAM", sep = "_"),
  "KISS_PROB_ONETIME_ALTPARAM"
)

rimlabs <- c(
  paste("RIM_RATE", c("MAIN", "CASUAL"), "ALTPARAM", sep = "_"),
  "RIM_PROB_ONETIME_ALTPARAM"
)

specit <- function(index, data) {
  vec <- data[index, ]
  capnames <- toupper(names(vec))
  jobnames <- paste(capnames, vec, sep = "_", collapse = "_")
  fullnames <- paste0("SENS_", jobnames)
  out <- list(jname = fullnames, params = vec)
  out
}

kissx_specs <- lapply(seq_len(nrow(kissX_rim0)), specit, data = kissX_rim0)
rimx_specs <- lapply(seq_len(nrow(rimX_kiss0)), specit, data = rimX_kiss0)


################################################################################
## ACT STOPPER PROB SENSITIVITY ##
################################################################################

actstop <- c(0.25, 0.5, 0.75)
actstop_names <- paste0("SENS_ActStopper_", format(actstop, digits = 2))

actstop_specs <- lapply(1:3, function(.x) {
  list(
    jname = actstop_names[.x],
    params = actstop[.x]
  )
})


################################################################################
## MAKE BATCHES ##
################################################################################

makescript <- function(x, type = c("kissrim", "actstop")) {
  if (type == "kissrim") {
    pline <- paste(c(kisslabs, rimlabs), x$params, sep = "=", collapse = ",")
  } else if (type == "actstop") {
    pline <- paste("ACT_STOPPER_PROB_ALTPARAM", x$params, sep = "=")
  }

  paste0(
    "sbatch -J ", x$jname, " --export=ALL,SIMDIR=~/scratch/", x$jname,
    paste0(",NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,", pline),
    " 02_main.sh"
  )
}

jobs <- c(
  sapply(kissx_specs, makescript, type = "kissrim"),
  sapply(rimx_specs, makescript, type = "kissrim")
)

for (i in seq_len(length(jobs))) {
  writeLines(
    jobs[[i]],
    here::here(
      "inst", "analysis01_epi",
      paste0(
        "07.", sprintf("%02d", i), "_",
        stringr::str_extract(jobs[[i]], "(?<=J ).*(?= \\-)"), ".sh"
      )
    )
  )
}

actindex <- c((length(jobs) + 1) : (length(jobs) + 3))
ajobs <- sapply(actstop_specs, makescript, type = "actstop")

for (i in seq_len(length(ajobs))) {
  writeLines(
    ajobs[[i]],
    here::here(
      "inst", "analysis01_epi",
      paste0(
        "07.", sprintf("%02d", actindex[i]), "_",
        stringr::str_extract(ajobs[[i]], "(?<=J ).*(?= \\-)"), ".sh"
      )
    )
  )
}
