library(data.table)
library(stringr)

batch <- list.files(
  here::here("inst/analysis02_screening/"),
  pattern = "^02\\.[0-9]{2}.*(SENS|base|cdc1|cdc2|symp|univ)"
)

bpats <- str_extract_all(
  batch,
  "(^02\\.[0-9]{2})|(?<!BU)([A-Za-z1-2]{4}(?=\\.sh))"
)

bpats2 <- lapply(
  bpats[which(sapply(bpats, length) == 2)],
  function(.x) {
    paste0(.x[1], ".*", toupper(.x[2]), "$")
  }
)

scratchdirs <- lapply(bpats2, function(.x) {
  list.files("~/scratch/", .x)
})

sbatchlines <- sapply(
  scratchdirs[which(sapply(scratchdirs, length) > 0)],
  function(.x) {
    if (.x %like% "SENS") {
      paste0(
        "sbatch -n 64 --export=ALL,SIMDIR=~/scratch/", .x,
        ",EPI_DEST=epi_",
        str_extract(.x, "02\\.[0-1][0-9]"), "_sens_",
        tolower(str_extract(.x, "[A-Z1-2]{4}$")),
        " sim_clean.sh"
      )
    } else {
      paste0(
        "sbatch -n 64 --export=ALL,SIMDIR=~/scratch/", .x,
        ",EPI_DEST=epi_",
        str_extract(.x, "02\\.[0-1][0-9]"), "_main_",
        tolower(str_extract(.x, "[A-Z1-2]{4}$")),
        " sim_clean.sh"
      )
    }
})

writeLines(
  sbatchlines,
  file.path(
    here::here("inst/analysis02_screening/"),
    "03.01_submit_cleaning_scripts.sh"
  )
)
