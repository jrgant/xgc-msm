library(futile.logger)
library(tryCatchLog)

options(keep.source = TRUE) # source code file name and line number tracking
options("tryCatchLog.write.error.dump.file" = TRUE) # dump for post-mortem analysis
options("tryCatchLog.write.error.folder" = here::here("burnin", "abc"))

flog.appender(
  appender.file(
    here::here(
      "burnin", "abc",
      paste0("abc_debug_gcsympt.", Sys.getenv("SLURM_JOBID"), ".log")
    )
  )
)  # to log into a file instead of console

flog.threshold(WARN)    # TRACE, DEBUG, INFO, WARN, ERROR, FATAL

tryCatchLog(source(here::here("burnin", "abc", Sys.getenv("SIMFILE"))))
