################################################################################
## SETUP ##
################################################################################

library(pacman)
p_load(EpiModelHIV, rlecuyer, data.table, stringr, lhs)

selected_paramsets <-
  readRDS(here::here("inst", "cal", "main_analysis_inputs.rds"))

an02_path <- here::here("inst", "analysis02_screening")


################################################################################
## SAMPLE PARAMETER SETS ##
################################################################################

nsim_spec <- 1000

set.seed(5768933)
paramset_sampled <- sample(
  seq_len(length(selected_paramsets)),
  nsim_spec,
  replace = TRUE
)

saveRDS(
  paramset_sampled,
  here::here(an02_path, "paramset_sampled.rds")
)


################################################################################
## SETUP RANDOM NUMBER STREAMS ##
################################################################################

streams <- paste0("str", sprintf("%04d", 1:nsim_spec))

## streams for burnin runs

.lec.CreateStream(streams)

saveRDS(
  .lec.Random.seed.table,
  here::here(an02_path, "seeds_burnin.rds")
)

.lec.exit()

## streams for scenario runs

.lec.CreateStream(streams)

saveRDS(
  .lec.Random.seed.table,
  here::here(an02_path, "seeds_scenario.rds")
)


################################################################################
## MAKE BATCH SCRIPT ##
################################################################################

# This function writes a batch script to submit a job array.
make_batch_script <- function(jobname, walltime, partition, mem,
                              ncores, array, log_fullpath, batchid,
                              nsims, nsteps, add_arrivals, runtype,
                              sti_screen_type, kiss_exposure, starttime,
                              isburnin, simdir) {

  sb <- "#SBATCH"

  specs <- paste(
    "#!/bin/bash",
    paste(sb, "-J", jobname),
    paste0(sb, " --time=", walltime),
    paste(sb, "-p", partition),
    paste0(sb, " --mem=", mem),
    paste(sb, "-n", ncores),
    paste0(sb, " --array=", array),
    paste(sb, "-o", log_fullpath),
    paste0(
      sb,
      " --export=ALL,NSIMS=", nsims,
      ",NSTEPS=", nsteps,
      ",ARRIVE_RATE_ADD_PER20K=", add_arrivals,
      ",EPI_RUN_TYPE=", runtype,
      ",STI_SCREEN_TYPE=", sti_screen_type,
      ",STI_SCREEN_KISS_EXPOSURE=", kiss_exposure,
      ",STARTTIME=", starttime,
      ",SIMDIR=", simdir
    ),
    paste(  sb, "--mail-type=ALL"   ),
    paste(  sb, "--mail-user=jrgant@brown.edu"),
    "module load R/4.0.3",
    "cd ~/data/jgantenb/xgcmsm/",
    "Rscript ./inst/analysis02_screening/02.00_epi.r --vanilla",
    sep = "\n"
  )

  sdir <- an02_path
  writeLines(
    text = specs,
    con = file.path(sdir, paste0(batchid, ".sh"))
  )

}


################################################################################
## MAIN ANALYSIS ##
################################################################################

## NOTE: kiss_exposure is only looked for by the model in CDC scenarios.
##       Set to FALSE by default.
screen_specs_main <- list(
  sti_burnin = data.table(
    runtype       = "burnin",
    scenario      = "base",
    kiss_exposure = "FALSE",
    starttime     = 1,
    num_steps     = 52 * 60
  ),
  sti_base = data.table(
    runtype       = "scenario",
    scenario      = "base",
    kiss_exposure = "FALSE",
    starttime     = 3121,
    num_steps     = 52 * 5
  ),
  sti_symp = data.table(
    runtype       = "scenario",
    scenario      = "symptomatic",
    kiss_exposure = "FALSE",
    starttime     = 3121,
    num_steps     = 52 * 5
  ),
  sti_cdc1 = data.table(
    runtype       = "scenario",
    scenario      = "cdc",
    kiss_exposure = "FALSE",
    starttime     = 3121,
    num_steps     = 52 * 5
  ),
  sti_cdc2 = data.table(
    runtype       = "scenario",
    scenario      = "cdc",
    kiss_exposure = "TRUE",
    starttime     = 3121,
    num_steps     = 52 * 5
  ),
  sti_univ = data.table(
    runtype       = "scenario",
    scenario      = "universal",
    kiss_exposure = "FALSE",
    starttime     = 3121,
    num_steps     = 52 * 5
  )
)

## Write scripts
lapply(names(screen_specs_main), function(.x) {
  index <- which(names(screen_specs_main) == .x)
  slug <- toupper(.x)
  tmp <- screen_specs_main[[.x]]
  make_batch_script(
    jobname = paste0("ScreenEpi-", slug),
    walltime = "3:00:00",
    partition = "batch",
    mem = "3GB",
    ncores = 1,
    array = paste0("1-", nsim_spec),
    log_fullpath = paste0("ScreenEpi-", slug, "_ARRAY-%A_JOB-%J_SIMNO-%4a.log"),
    batchid = paste0("02.", "01", "_", .x),
    nsims = 1,
    nsteps = tmp[, num_steps],
    add_arrivals = 1.285,
    runtype = tmp[, runtype],
    sti_screen_type = tmp[, scenario],
    kiss_exposure = tmp[, kiss_exposure],
    starttime = tmp[, starttime],
    simdir = paste0("~/scratch/ScreenEpi-", slug)
  )
})


################################################################################
## SENSITIVITY ANALYSES ##
################################################################################

pard <- rbindlist(
  lapply(selected_paramsets, as.data.table, keep.rownames = T),
  idcol = "pid"
)

setnames(pard, c("V1", "V2"), c("input", "value"))

pardlims <- pard[, .(min = min(value), max = max(value)), by = input]


################################################################################
## PREP SENSITIVITIES ##
################################################################################

prep_scalars <- seq(1.5, 3, 0.5)

prep_ids <- seq_along(prep_scalars) + 1

prep_names <- paste0(
  "SENS_02.", sprintf("%02d", prep_ids),
  "_PrEP_", format(prep_scalars, digits = 2)
)

prep_specs <- lapply(seq_along(prep_ids), function(.x) {
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

sympt_alt <- dcast(lhsmatch, id ~ input, value.var = "value")
sympt_labs <- paste0(names(sympt_alt)[-1], "_ALTPARAM")

## Specify sensitivity analysis specs.
specit <- function(index, data,
                   startnum = length(prep_ids) + 1) {

  vec <- unlist(data[index, -c("id")])
  pnames <- str_extract(names(vec), "[A-Z]+(?=_)")

  out <- list(
    jname = paste0(
      "SENS_02.", sprintf("%02d", startnum + index),
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

make_sens_script <- function(x, type = c("prep", "sympt"), stiscreen) {

  if (type == "prep") {
    pline <- paste("PREP_SCALE_ALTPARAM", x$params, sep = "=")
  } else if (type == "sympt") {
    pline <- paste(sympt_labs, x$params, sep = "=", collapse = ",")
  }

  currscreen <- screen_specs_main[[stiscreen]]
  bfile <- paste0("02.01_", stiscreen, ".sh")

  paste0(
    "sbatch -J ", paste0(x$jname, "_SCREENTYPE_", toupper(stiscreen)),
    " -o ", str_extract(x$jname, "SENS_02\\.[0-9]{2}"),
    "_ARRAY-%A_JOB-%J_SIMNO-%4a.log",
    " -t 5:00:00 ",
    " --export=ALL,SIMDIR=~/scratch/",
    paste0(x$jname, "_SCREEN_", toupper(stiscreen)),
    ",EPI_RUN_TYPE=", currscreen$runtype,
    ",STI_SCREEN_TYPE=", currscreen$scenario,
    ",STI_SCREEN_KISS_EXPOSURE=", currscreen$kiss_exposure,
    paste0(",NSIMS=1,NSTEPS=3380,ARRIVE_RATE_ADD_PER20K=1.285,", pline),
    " ", bfile
  )
}

jobs <- c(
  sapply(
    sympt_specs, make_sens_script, type = "sympt", stiscreen = "sti_burnin"
  ),
  sapply(sympt_specs, make_sens_script, type = "sympt", stiscreen = "sti_base"),
  sapply(sympt_specs, make_sens_script, type = "sympt", stiscreen = "sti_symp"),
  sapply(sympt_specs, make_sens_script, type = "sympt", stiscreen = "sti_cdc1"),
  sapply(sympt_specs, make_sens_script, type = "sympt", stiscreen = "sti_cdc2"),
  sapply(sympt_specs, make_sens_script, type = "sympt", stiscreen = "sti_univ"),
  sapply(prep_specs, make_sens_script, type = "prep", stiscreen = "sti_burnin"),
  sapply(prep_specs, make_sens_script, type = "prep", stiscreen = "sti_base"),
  sapply(prep_specs, make_sens_script, type = "prep", stiscreen = "sti_symp"),
  sapply(prep_specs, make_sens_script, type = "prep", stiscreen = "sti_cdc1"),
  sapply(prep_specs, make_sens_script, type = "prep", stiscreen = "sti_cdc2"),
  sapply(prep_specs, make_sens_script, type = "prep", stiscreen = "sti_univ")
)

for (i in seq_len(length(jobs))) {
  writeLines(
    jobs[[i]],
    here::here(
      "inst", "analysis02_screening",
      paste0(
        str_extract(jobs[[i]], "(?<=J SENS_)02\\.[0-9]{2}(?=_)"), "_SENS_",
        str_extract(jobs[[i]], "(?<=J SENS_02\\.[0-9]{2}_).*(?= \\-o)"),
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
  here::here("inst/analysis02_screening", "specs_SENS.rds")
)
