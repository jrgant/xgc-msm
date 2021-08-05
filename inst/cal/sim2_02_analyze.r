################################################################################
## SETUP ##
################################################################################

picksim <- "sim2"

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


################################################################################
## VISUALIZE DISTRIBUTIONS OF TARGETS ACROSS ALL SIMS ##
################################################################################
source(here::here(ic_dir, "sim0.1_fmt_targetdat.r"))
ls()


################################################################################
## SELECT SIMULATIONS ##
################################################################################
simid_sel_hivpr <-
  out_vs_targ[variable == "i.prev" & output_within_5pts == 1, simid]

## NOTE
## quicktargets() is a helper function and particular about its inputs,
## so we make some lists to get several targets one plot
set_jitter <- list(w = 0.1, h = 0)

cs_vsupp_reth <- lapply(1:4, function(x) simid_sel_hivpr)
names(cs_vsupp_reth) <- rlabs

quicktarget("i.prev", list(hiv = simid_sel_hivpr), set_jitter)
quicktarget("ir100.pop", list(hiv = simid_sel_hivpr), set_jitter)
quicktarget("cc.vsupp.%s", cs_vsupp_reth, set_jitter)


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim2.rds"))

hivpr_sel_inputs <- pull_params(list(simid_sel_hivpr))[, -c("selection_group")]

all(
  all(simid_sel_hivpr %in% hivpr_sel_inputs[, simid]),
  all(hivpr_sel_inputs[, simid] %in% simid_sel_hivpr)
)


## Since only a handful of simids were selected, get the min/max for each
## input parameter to use as the new sampling ranges in sim3.
narrowed_hiv_priors <- hivpr_sel_inputs[, .(
  q25 = quantile(value, 0.25),
  q75 = quantile(value, 0.75)
), input][input %like% "HIV"]

sim2_priors <- readRDS(here::here(ic_dir, "sim1_sel_lhs_limits.rds"))

new_priors <- merge(
  sim2_priors,
  narrowed_hiv_priors,
  by = "input",
  all.x = TRUE
)

new_priors <- new_priors[
  !is.na(q25) & !is.na(q75),
  ":=" (s2_ll = q25, s2_ul = q75)][, -c("q25", "q75")]

setnames(new_priors, c("s2_ll", "s2_ul"), c("s3_ll", "s3_ul"))


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(simid_sel_hivpr, here::here("inst", "cal", "sim2_simid_sel.rds"))

## Write pre-selection boxplots to files
lapply(
  ls(pattern = "^ps_"),
  function(x) {
    psave(
      f = sprintf("preselect_%s", stringr::str_extract(x, "(?<=ps_).*")),
      p = get(x)
    )
    invisible()
  }
)

## Write limits of input parameter values across selected simulations.
saveRDS(new_priors, here::here("inst", "cal", "sim2_sel_lhs_limits.rds"))
