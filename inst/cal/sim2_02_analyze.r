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
pretty_batch_date

################################################################################
## VISUALIZE DISTRIBUTIONS OF TARGETS ACROSS ALL SIMS ##
################################################################################
source(here::here(ic_dir, "sim0.1_fmt_targetdat.r"))
ls()

ps_ct_hivprev
ps_ct_hivinc
ps_ct_vsupp_race
ps_ct_vsupp_age
ps_ct_prep
ps_ct_gctested
ps_ct_gcpos

################################################################################
## SELECT SIMULATIONS ##
################################################################################

simid_sel_hivpr <-
  out_vs_targ[variable == "i.prev" & output_within_5pts == 1, sort(simid)]

simid_sel_hivinc <- out_vs_targ[
  variable == "ir100.pop" & between(output, 0.4, 0.6), sort(simid)]

simid_sel_hivprinc <- intersect(simid_sel_hivpr, simid_sel_hivinc)
length(simid_sel_hivprinc)

simid_sel_dx_race_rel <- lapply(
  setNames(rslugs[1:3], rslugs[1:3]),
  function(.x) {
    out_vs_targ[
      variable == sprintf("i.prev.dx.inf.%s.rel.ref.W", .x) & output < 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_dx_race_rel, length)

simid_sel_dx_race_rel_intersect <- Reduce(intersect, simid_sel_dx_race_rel)
length(simid_sel_dx_race_rel_intersect)

ageslugs <- paste0("age", 1:5)
simid_sel_dx_age_rel <- lapply(
  setNames(age.rel.dx, age.rel.dx),
  function(.x) {
    out_vs_targ[
      variable == .x & output < 1,
      sort(simid)
    ]
  }
)

simid_sel_dx_age_rel_intersect <- Reduce(intersect, simid_sel_dx_age_rel)
length(simid_sel_dx_age_rel_intersect)

simid_sel_dp_intersect <- Reduce(
  intersect,
  list(
    simid_sel_hivprinc,
    simid_sel_dx_race_rel_intersect,
    simid_sel_dx_age_rel_intersect
  )
)

length(simid_sel_dp_intersect)


out_vs_targ[simid %in% simid_sel_dp_intersect] %>%
  ggplot(aes(x = variable, y = output)) +
  geom_point(position = position_jitter(width = 0.1)) +
  geom_point(
    aes(y = target_val, color = "target"),
    shape = 21,
    size = 7,
    stroke = 2,
    fill = "white"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_line(color = "gray80")
  )


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim2.rds"))
sim2_priors <- readRDS(here::here("inst/cal", "sim1_sel_lhs_limits.rds"))

sel_inputs <- pull_params(simid_sel_dp_intersect)[, -c("selection_group")]

new_priors <- sel_inputs[, .(
  ll = min(value),
  ul = max(value)
), by = input]

setnames(new_priors, c("ll", "ul"), c("s3_ll", "s3_ul"))


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(
  simid_sel_dp_intersect,
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
