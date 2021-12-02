################################################################################
## SETUP ##
################################################################################

library(pacman)

p_load(
  data.table,
  EpiModel,
  ggplot2,
  ggthemes,
  foreach,
  doParallel,
  stringr
)


################################################################################
## FILES ##
################################################################################


simdir <- Sys.getenv("SIMDIR")
picksim <- Sys.getenv("PICKSIM")
first_timestep <- Sys.getenv("FIRST_TIME")

if (simdir == "" | picksim == "" | first_timestep == "") {
  stop("Must set the environment variables SIMDIR and PICKSIM in sbatch.")
} else if (
         str_extract(simdir, "sim[0-9]{1}") != str_extract(picksim, "sim[0-9]{1}")
       ) {
  stop(
    "SIMDIR and PICKSIM indicate different simulation batches."
  )
}

## Print some info to the batch logs
cat("PICKSIM=", picksim, "\n\n")
cat("Retrieving simulations from", simdir, "\n\n")

destdir <- here::here("burnin", "cal", picksim)
cat("Writing cleaned and combined simulations to", destdir, "\n\n")

## List and process files
cs <- list.files(
  simdir,
  pattern = "rds",
  full.names = TRUE
)

ncores <- detectCores()
registerDoParallel(ncores)

epi <- foreach(i = seq_along(cs), .combine = rbind) %dopar% {
  file_fullpath <- cs[i]
  cd <- readRDS(file_fullpath)

  # extract epi
  cde <- as.data.table(cd$epi)
  cde[, simid :=
          stringr::str_extract(file_fullpath, "(?<=sim_)[0-9]{4}")]

  cde[, at := 1:.N]
  cde[at >= as.numeric(first_timestep)] # keep last x years of 60-year burnin
}

## subset variables to reduce filesize
racelabs <- c(".B", ".H", ".O", ".W")
agelabs <- paste0(".age", 1:5)
agedots <- paste0(".age.", 1:5)

std_vec <- c("", racelabs, agelabs)
dot_vec <- c("", racelabs, agedots)

keepvars <- c(
  "simid", "at",
  paste0("num", dot_vec),
  paste0("i.num", std_vec),
  paste0("i.prev", dot_vec),
  paste0("i.prev.dx", dot_vec),
  paste0("incid", std_vec),
  paste0("cc.vsupp", std_vec),
  paste0("prepElig", std_vec[1:5]),
  paste0("prepCurr", std_vec[1:5]),
  paste0("i.num.", c("gc", "rgc", "ugc", "pgc")),
  "ir100",
  paste0("prop.", c("rect", "ureth", "phar"), ".tested"),
  paste0("prob.", c("rGC", "uGC", "pGC"), ".tested")
)


saveRDS(
  epi[, ..keepvars],
  here::here(destdir, paste0(picksim, "_epi_", as.numeric(Sys.time()), ".rds"))
)
