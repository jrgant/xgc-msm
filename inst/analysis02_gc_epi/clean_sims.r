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
  doParallel
)


################################################################################
## FILES ##
################################################################################


simdir <- Sys.getenv("SIMDIR")

if (Sys.getenv("SIMDIR") == "" | Sys.getenv("EPI_DEST") == "") {
  stop("Must set the environment variableS SIMDIR and EPI_DEST in sbatch.")
}

print(paste("Retrieving simulations from", simdir))

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
  cde[at > 3120] # keep last 5 years: 60-year burnin + 5-year analytic run
}

saveRDS(
  epi,
  here::here("inst", "analysis02_gc_epi", paste0(Sys.getenv("EPI_DEST"), ".rds"))
)
