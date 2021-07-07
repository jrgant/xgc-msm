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

if (Sys.getenv("SIMDIR") == "") {
  stop("Must set the environment variable SIMDIR in sbatch.")
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
  cde[at >= 3016] # keep last 2 years of 60-year burnin
}

if (Sys.getenv("EPI_DEST") != "") {
  saveRDS(
    epi,
    here::here("burnin", "cal", "sim2", paste0(Sys.getenv("EPI_DEST"), ".rds"))
  )
} else {
  saveRDS(epi, here::here("burnin", "cal", "sim2", "sim2_epi.rds"))
}
