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
  variable == "ir100.pop" & between(output, 0.444, 0.584), sort(simid)]

simid_sel_hivprinc <- intersect(simid_sel_hivpr, simid_sel_hivinc)
length(simid_sel_hivprinc)

simid_sel_vls_race <- lapply(
  setNames(rslugs, rslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("cc.vsupp.%s", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_vls_race, length)
simid_sel_vls_race_intersect <- Reduce(intersect, simid_sel_vls_race)

simid_sel_dx_race <- lapply(
  setNames(rslugs, rslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("i.prev.dx.inf.%s", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_dx_race, length)
simid_sel_dx_race_intersect <- Reduce(intersect, simid_sel_dx_race)

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

simid_sel_vls_age <- lapply(
  setNames(ageslugs, ageslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("cc.vsupp.%s", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_vls_age, length)
simid_sel_vls_age_intersect <- Reduce(intersect, simid_sel_vls_age)

simid_sel_dx_age <- lapply(
  setNames(ageslugs, ageslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("i.prev.dx.inf.%s", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

simid_sel_dx_age_intersect <- Reduce(intersect, simid_sel_dx_age)
length(simid_sel_dx_age_intersect)

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

simid_sel_gctest <- lapply(
  setNames(anatslugs, anatslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("prop.%s.tested", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_gctest, length)
simid_sel_gctest_intersect <- Reduce(intersect, simid_sel_gctest)

simid_sel_gcpos <- lapply(
  setNames(gcpos_slugs, gcpos_slugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("prob.%s.tested", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_gcpos, length)

simid_sel_gcpos_intersect <- Reduce(intersect, simid_sel_gcpos)
length(simid_sel_gcpos_intersect)

simid_sel_dp_intersect <- Reduce(
  intersect,
  list(
    simid_sel_hivprinc,
    simid_sel_vls_race_intersect,
    simid_sel_dx_age_intersect,
    simid_sel_dx_age_rel_intersect
  )
)

length(simid_sel_dp_intersect)

plot_targets <- function(simids) {
  out_vs_targ[simid %in% simids] %>%
    ggplot(aes(x = variable, y = output)) +
    geom_point(position = position_jitter(width = 0.1)) +
    geom_point(
      aes(y = target_val, color = "target"),
      shape = 21,
      size = 4,
      stroke = 2,
      fill = NA
    ) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major.x = element_line(color = "gray80")
    )
}

plot_targets(simid_sel_dp_intersect)


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim2.rds"))

lhsdt <- rbindlist(
  lapply(lhs_groups, as.data.table, keep.rownames = TRUE),
  idcol = "simid"
)

lhsdt[, simid := sprintf("%04d", simid)][]
setnames(lhsdt, c("V1", "V2"), c("input", "value"))

lhsdt[simid %in% simid_sel_dp_intersect] %>%
  ggplot(aes(x = input, y = value, group = simid)) +
  geom_line() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

sim2_priors <- readRDS(here::here("inst/cal", "sim1_sel_lhs_limits.rds"))

sel_hiv_inputs <-
  pull_params(simid_sel_dp_intersect)[, -c("selection_group")]

sel_hiv_inputs[!(input %like% "RX_HALT")] %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~ input, scales = "free_x")

sel_hiv_inputs_w <- dcast(
  sel_hiv_inputs[!(input %like% "RX_HALT")],
  simid ~ input,
  value.var = "value"
)

ggcorrplot(cor(sel_hiv_inputs_w[, -c(1:2)]), type = "upper")

narrowed_priors <- sel_hiv_inputs[, .(
  ll = min(value),
  ul = max(value)
), input]
##[!(input %like% "_GC|GC_|ASYMP|SYMP|OI_ACT")]

## plot_targets(intersect(simid_sel_gcpos$rGC, simid_sel_hivpr))
## plot_targets(intersect(simid_sel_gcpos$uGC, simid_sel_hivprinc))
## plot_targets(intersect(simid_sel_gcpos$pGC, simid_sel_hivprinc))

## lhsdt[simid %in% simid_sel_gcpos$pGC & !(input %like% "RX_HALT")] %>%
##   ggplot(aes(x = value)) +
##   geom_histogram() +
##   facet_wrap(~ input, scale = "free_x")

new_priors <- merge(
  sim2_priors,
  narrowed_priors,
  by = "input",
  all.x = TRUE
)

new_priors[
  !is.na(ll) & !is.na(ul),
  ":=" (
    s2_ll = ll,
    s2_ul = ul
  )][]

new_priors[, ":=" (ll = NULL, ul = NULL)][]

setnames(new_priors, c("s2_ll", "s2_ul"), c("s3_ll", "s3_ul"))


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
