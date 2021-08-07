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

simid_sel_vls <-
  lapply(
    setNames(rslugs, rslugs),
    function(.x) {
      out_vs_targ[
        variable == sprintf("cc.vsupp.%s", .x) & output_within_5pts == 1, simid]
    }
  )

simid_sel_prep <-
  lapply(
    setNames(rslugs, rslugs),
    function(.x) {
      out_vs_targ[
        variable == sprintf("prepCov.%s", .x) & output_within_5pts == 1, simid]
    }
  )

simid_sel_gcpos <-
  lapply(
    setNames(gcpos_slugs, gcpos_slugs),
    function(.x) {
      out_vs_targ[variable == sprintf("prob.%s.tested", .x)
                  ][output > 0 & output_within_5pts == 1, simid]
    }
  )
out_vs_targ[variable == "i.prev" & output_within_5pts == 1]
simid_sel_gcpos

## NOTE
## quicktargets() is a helper function and particular about its inputs,
## so we make some lists to get several targets one plot
set_jitter <- list(w = 0.1, h = 0)

cs_vsupp_reth <- lapply(1:4, function(x) simid_sel_gcpos)
names(cs_vsupp_reth) <- rslugs

quicktarget("cc.vsupp.%s", simid_sel_vls, set_jitter)

quicktarget("cc.vsupp.%s", simid_sel_vls, list(w = 0, h = 0))

quicktarget(
  "prob.%s.tested",
  simid_sel_gcpos,
  set_jitter
)


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim2.rds"))

vls_sel_inputs <- pull_params(simid_sel_vls)

prep_sel_inputs <- pull_params(simid_sel_prep)[input %like% "PREP"]

gcpos_sel_inputs <-
  pull_params(
    simid_sel_gcpos
  )[input %like% "DURAT_NOTX|RX_INFPR|SYMPT_PROB|ASYMP"]

gcpos_sel_oiscale_inputs <-
  pull_params(list(pGC = simid_sel_gcpos$pGC))[input %like% "OI_ACT"]

narrowed_vls_priors <- input_quantiles(vls_sel_inputs, "RX_|TESTER")

narrowed_prep_inputs <- input_quantiles(prep_sel_inputs, "DISCONT")

narrowed_gc_priors <- gcpos_sel_inputs[, .(
  q25 = quantile(value, 0.25),
  q75 = quantile(value, 0.75)
  ), .(selection_group, input)
  ][substring(input, 1, 1) == toupper(substring(selection_group, 1, 1))
  ][, -c("selection_group")]

sim2_priors <- readRDS(here::here(ic_dir, "sim1_sel_lhs_limits.rds"))

new_priors <- merge(
  sim2_priors,
  rbind(
    narrowed_prep_inputs,
    narrowed_vls_priors,
    narrowed_gc_priors
  ),
  by = "input",
  all.x = TRUE
)

new_priors[!is.na(q25), all(q25 >= s2_ll)]
new_priors[!is.na(q75), all(q75 <= s2_ul)]

new_priors <- new_priors[
  !is.na(q25) & !is.na(q75),
  ":=" (s2_ll = q25, s2_ul = q75)][, -c("q25", "q75")]

setnames(new_priors, c("s2_ll", "s2_ul"), c("s3_ll", "s3_ul"))

## widen HIV transmission prob. scalar ranges for sim3
new_priors[input %like% "SCALAR_HIV", ":=" (s3_ll = 0.5, s3_ul = 3)][]

## reset HIV condom efficacy prior
new_priors[input %like% "CONDOM_EFF_HIV", ":=" (s3_ll = 0.6, s3_ul = 1)][]


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(
  list(simid_sel_prep, simid_sel_vls, simid_sel_gcpos),
  here::here("inst", "cal", "sim2_simid_sel.rds")
)

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
