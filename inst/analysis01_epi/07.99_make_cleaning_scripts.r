library(data.table)
library(stringr)

batch07 <- list.files("inst/analysis01_epi/", pattern = "07\\.(0[1-9]|1[0-9])")
batch07

scratchdirs <- paste0("~/scratch/SENS_", str_remove_all(batch07, "(SENS_)|(\\.sh$)"))

sbatchlines <- sapply(scratchdirs, function(.x) {
  paste0(
    "sbatch -n 64 --export=ALL,SIMDIR=", .x,
    ",EPI_DEST=epi_SENS_", str_extract(.x, "07\\.[0-1][0-9]"),
    " sim_clean.sh"
  )
})


writeLines(
  sbatchlines,
  file.path(here::here("inst/analysis01_epi/"), "07.99_submit_cleaning_scripts.sh")
)
