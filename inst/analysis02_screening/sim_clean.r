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
  cd <- readRDS(cs[i])

  # extract epi
  cde <- as.data.table(cd$epi)
  cde[, simid := i]
  cde[, at := 1:.N]

  # choose columns (reduce filesize)
  selcols <- c(
    "simid", "at",
    names(cde)[grepl("gc|GC|prep|^num", names(cde))]
  )

  # return
  cde[at > 3120, ..selcols] # keep last 5 years after a 60-year burnin
}

an02_path <- here::here("inst", "analysis02_screening")

if (Sys.getenv("EPI_DEST") != "") {
  saveRDS(
    epi,
    here::here(an02_path, paste0(Sys.getenv("EPI_DEST"), ".rds"))
  )
} else {
  saveRDS(epi, here::here(an02_path, "main_epi.rds"))
}
