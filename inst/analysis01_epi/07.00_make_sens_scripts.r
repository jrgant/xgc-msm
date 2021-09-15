################################################################################
## SETUP ##
################################################################################

library(data.table)
library(stringr)

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

specit <- function(index, data, ids) {
  vec <- data[index, ]
  capnames <- toupper(names(vec))
  jobnames <- paste(capnames, vec, sep = "_", collapse = "_")
  fullnames <- paste0("SENS_07.", sprintf("%02d", ids[index]), "_", jobnames)
  out <- list(jname = fullnames, params = vec, id = ids[index])
  out
}

kissx_specs <- lapply(
  seq_len(nrow(kissX_rim0)),
  specit,
  data = kissX_rim0,
  ids = seq_len(nrow(kissX_rim0))
)

rimx_specs <- lapply(
  seq_len(nrow(rimX_kiss0)),
  specit,
  data = rimX_kiss0,
  ids = nrow(kissX_rim0) + seq_len(nrow(rimX_kiss0))
)


################################################################################
## ACT STOPPER PROB SENSITIVITY ##
################################################################################

actstop <- c(0.25, 0.5, 0.75)

actstop_ids <-
  sum(c(nrow(kissX_rim0), nrow(rimX_kiss0))) + seq_len(length(actstop))

actstop_names <- paste0(
  "SENS_07.", sprintf("%02d", actstop_ids),
  "_ActStopper_", format(actstop, digits = 2)
)

actstop_specs <- lapply(1:3, function(.x) {
  list(
    jname = actstop_names[.x],
    params = actstop[.x],
    ids = actstop_ids[.x]
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
    "sbatch -J ", x$jname,
    " -o ", str_extract(x$jname, "SENS_07\\.[0-9]{2}"),
    "_ARRAY-%A_JOB-%J_SIMNO-%4a.log",
    " --export=ALL,SIMDIR=~/scratch/", x$jname,
    paste0(",NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,", pline),
    " 02_main.sh"
  )
}

jobs <- c(
  sapply(kissx_specs,   makescript,  type = "kissrim"),
  sapply(rimx_specs,    makescript,  type = "kissrim"),
  sapply(actstop_specs, makescript,  type = "actstop")
)

for (i in seq_len(length(jobs))) {
  writeLines(
    jobs[[i]],
    here::here(
      "inst", "analysis01_epi",
      paste0(
        str_extract(jobs[[i]], "(?<=J SENS_)07\\.[0-9]{2}(?=_)"), "_SENS_",
        str_extract(jobs[[i]], "(?<=J SENS_07\\.[0-9]{2}_).*(?= \\-o)"),
        ".sh"
      )
    )
  )
}
