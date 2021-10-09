################################################################################
## SETUP ##
################################################################################

library(data.table)
library(stringr)
library(lhs)

selparams <- readRDS(here::here("inst/cal/main_analysis_inputs.rds"))

pard <- rbindlist(
  lapply(selparams, as.data.table, keep.rownames = T),
  idcol = "pid"
)

setnames(pard, c("V1", "V2"), c("input", "value"))

pardlims <- pard[, .(min = min(value), max = max(value)), by = input]


################################################################################
## PREP SENSITIVITIES ##
################################################################################

prep_scalars <- seq(1.5, 3, 0.5)

prep_ids <- seq_along(prep_scalars)

prep_names <- paste0(
  "SENS_03.", sprintf("%02d", prep_ids),
  "_PrEP_", format(prep_scalars, digits = 2)
)

prep_specs <- lapply(prep_ids, function(.x) {
  list(
    jname = prep_names[.x],
    params = prep_scalars[.x],
    ids = prep_ids[.x]
  )
})


################################################################################
## SYMPTOM PROBABILITY SENSITIVITIES ##
################################################################################

gcsympt_ins <- pardlims[input %like% "GC_SYMPT_PROB"]

set.seed(7348938)

lhs_sample <- as.data.table(
  geneticLHS(n = 5, k = 3, pop = 10000, gen = 1000, pMut = 0.1)
)[, id := 1:.N]

names(lhs_sample)[1:3] <- gcsympt_ins$input

lhs_sample

lhsl <- melt(
  lhs_sample,
  id.vars = "id",
  measure.vars = names(lhs_sample)[-4],
  variable.name = "input",
  value.name = "lhsdraw"
)

lhsmatch <- merge(gcsympt_ins, lhsl, by = "input", all.x = T)
lhsmatch[, value := qunif(lhsdraw, min, max)]

library(plotly)

plot_ly(
  lhs_sample,
  x = ~RECT_GC_SYMPT_PROB,
  y = ~URETH_GC_SYMPT_PROB,
  z = ~PHAR_GC_SYMPT_PROB,
  type = "scatter3d",
  mode = "markers+text"
)

sympt_alt <- dcast(lhsmatch, id ~ input, value.var = "value")

sympt_labs <- paste0(names(sympt_alt)[-1], "_ALTPARAM")

## NOTE: kiss_exposure is only looked for by the model in CDC scenarios.
##       Set to FALSE by default.
screen_specs <- list(
  sti_base = data.table(scenario = "base",        kiss_exposure = "FALSE"),
  sti_symp = data.table(scenario = "symptomatic", kiss_exposure = "FALSE"),
  sti_cdc1 = data.table(scenario = "cdc",         kiss_exposure = "FALSE"),
  sti_cdc2 = data.table(scenario = "cdc",         kiss_exposure = "TRUE"),
  sti_univ = data.table(scenario = "universal",   kiss_exposure = "FALSE")
)

specit <- function(index, data, startnum = length(prep_ids)) {

  vec <- unlist(data[index, -c("id")])
  pnames <- str_extract(names(vec), "[A-Z]+(?=_)")

  out <- list(
    jname = paste0(
      "SENS_03.", sprintf("%02d", startnum + index),
      "_GCSYMPT_", paste(pnames, round(vec, 3), sep = "_", collapse = "_")
    ),
    params = vec,
    id = startnum + index
  )
  out
}

sympt_specs <- lapply(
  seq_len(nrow(sympt_alt)),
  specit,
  data = sympt_alt
)


################################################################################
## MAKE BATCHES ##
################################################################################

makescript <- function(x, type = c("prep", "sympt"), stiscreen) {
  if (type == "prep") {
    pline <- paste("PREP_SCALE_ALTPARAM", x$params, sep = "=")
  } else if (type == "sympt") {
    pline <- paste(sympt_labs, x$params, sep = "=", collapse = ",")
  }

  currscreen <- screen_specs[[stiscreen]]

  bfile <- paste0(names(match.arg(
    stiscreen,
    c("02.01" = "sti_base", "02.02" = "sti_symp", "02.03" = "sti_cdc1",
      "02.04" = "sti_cdc2.sh", "02.05" = "sti_univ")
  )), "_", stiscreen, ".sh")

  paste0(
    "sbatch -J ", paste0(x$jname, "_SCREENTYPE_", toupper(stiscreen)),
    " -o ", str_extract(x$jname, "SENS_03\\.[0-9]{2}"),
    "_ARRAY-%A_JOB-%J_SIMNO-%4a.log",
    " -t 5:00:00 ",
    " --export=ALL,SIMDIR=~/scratch/",
    paste0(x$jname, "_SCREEN_", toupper(stiscreen)),
    ",STI_SCREEN_TYPE=", currscreen$scenario,
    ",STI_SCREEN_KISS_EXPOSURE=", currscreen$kiss_exposure,
    paste0(",NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,", pline),
    " ", bfile
  )
}

jobs <- c(
  sapply(sympt_specs, makescript, type = "sympt", stiscreen = "sti_base"),
  sapply(sympt_specs, makescript, type = "sympt", stiscreen = "sti_symp"),
  sapply(sympt_specs, makescript, type = "sympt", stiscreen = "sti_cdc1"),
  sapply(sympt_specs, makescript, type = "sympt", stiscreen = "sti_cdc2"),
  sapply(sympt_specs, makescript, type = "sympt", stiscreen = "sti_univ"),
  sapply(prep_specs, makescript, type = "prep", stiscreen = "sti_base"),
  sapply(prep_specs, makescript, type = "prep", stiscreen = "sti_symp"),
  sapply(prep_specs, makescript, type = "prep", stiscreen = "sti_cdc1"),
  sapply(prep_specs, makescript, type = "prep", stiscreen = "sti_cdc2"),
  sapply(prep_specs, makescript, type = "prep", stiscreen = "sti_univ")
)


for (i in seq_len(length(jobs))) {
  writeLines(
    jobs[[i]],
    here::here(
      "inst", "analysis02_screening",
      paste0(
        str_extract(jobs[[i]], "(?<=J SENS_)03\\.[0-9]{2}(?=_)"), "_SENS_",
        str_extract(jobs[[i]], "(?<=J SENS_03\\.[0-9]{2}_).*(?= \\-o)"),
        ".sh"
      )
    )
  )
}


################################################################################
## SAVE SPECS ##
################################################################################

quickdt <- function(data) rbindlist(lapply(data, as.data.table))

speclist <- list(
  prep = quickdt(prep_specs),
  sympt = rbindlist(
    lapply(
      sympt_specs, function(.x) {
        dcast(
          as.data.table(.x[["params"]], keep.rownames = T),
          . ~ V1,
          value.var = "V2"
        )
      }
    )
  )[, -1]
)

saveRDS(
  speclist,
  here::here("inst/analysis02_screening", "specs_SENS03.rds")
)
