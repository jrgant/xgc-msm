library(pacman)
p_load(
  EpiModelHIV,
  EpiModel,
  xgcmsm,
  data.table,
  ggplot2,
  ggthemes,
  parallel,
  tictoc
)

nsl <- list.files(
  here::here("burnin", "abc", "sim1"),
  pattern = "rds",
  full.names = TRUE
)

ncores <- detectCores() - 2

fd <- as.data.table(list(file = nsl, id = seq_len(length(nsl))))[1:400]
split_int <- ceiling(nrow(fd) / ncores)

fd[, batch := rep(1:ncores, each = split_int, length.out = nrow(fd))]

batches <- split(fd, by = "batch")

tic()
cl <- makeCluster(ncores)

sims <- parLapply(cl = cl, batches, function(.x) {
  require(data.table)
  f <- .x$file
  fl <- lapply(f, function(.y) {
    s <- readRDS(.y)
    d <- as.data.table(s$epi)
    d[, at := 1:.N]
    d[, sim := stringr::str_extract(.y, "(?<=_)[0-9]+(?=)")]
    return(d)
  })
  return(fl)
})

stopCluster(cl)
toc()

dl <- rbindlist(lapply(sims, rbindlist))

fwrite(dl, here::here("burnin", "abc", "sim1", "sim1_epi.csv"))


