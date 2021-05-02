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

cs <- list.files(
  "~/scratch/sim1",
  pattern = "rds",
  full.names = TRUE
)

ncores <- detectCores()
registerDoParallel(ncores)

epi <- foreach(i = seq_along(cs), .combine = rbind) %dopar% {
  cd <- readRDS(cs[i])
  cde <- as.data.table(cd$epi)

  cde[, simid := i]
  cde[, at := 1:.N]

  cde[at >= 3016] # keep last 2 years of 60-year burnin
}

saveRDS(epi, here::here("burnin", "abc", "sim1", "sim1_epi.rds"))
