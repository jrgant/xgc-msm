################################################################################
## SETUP ##
################################################################################

picksim <- "sim4"

# NOTE The file dated 2021-7-27 contains entire burn-in period.
#      Needed to investigate extinction of GC in selected runs and lack of
#      equilibrium in HIV dynamics. We will narrow parameter sets based on
#      overall HIV prevalence in this script, and then take another
#      random Latin hypercube sample for sim3. Based on looking at 2021-07-27
#      it seems like we should get parameter sets that result in HIV
#      equilibrium and non-extinction of GC.
# Sys.setenv("BATCH_DATE" = "2021-07-26")

ic_dir <- here::here("inst", "cal")
source(here::here(ic_dir, "sim0.0_setup.r"))
ls()

sim_epis[sim_epis %like% batch_date]
pretty_batch_date


################################################################################
## VISUALIZE DISTRIBUTIONS OF TARGETS ACROSS ALL SIMS ##
################################################################################
source(here::here(ic_dir, "sim0.1_fmt_targetdat.r"))
ls()

out_vs_targ


################################################################################
## SEE SELECTED INPUTS ##
################################################################################

inputs <- readRDS(file.path(ic_dir, "main_analysis_inputs.rds"))

params <- rbindlist(
  lapply(inputs, as.data.table, keep.rownames = TRUE),
  idcol = "pset"
)[!V1 %like% "RX_HALT|OI_ACT|AI_ACT|TRANS_PROB_WHITE"]

setnames(params, c("V1", "V2"), c("param", "value"))

params %>%
  ggplot(aes(x = value)) +
  geom_histogram(color = "white") +
  geom_density() +
  facet_wrap(~ param, scales = "free")


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write summarized simulations to file.
saveRDS(
  out_vs_targ,
  here::here("inst", "cal", "sim4_calsum.rds")
)
