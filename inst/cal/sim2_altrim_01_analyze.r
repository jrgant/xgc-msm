################################################################################
## SETUP ##
################################################################################

picksim <- "sim2_altrim"

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

# NOTE: These are the only simids used to select runs in this round.
simid_sel_gcpos_intersect <- Reduce(intersect, simid_sel_gcpos)
length(simid_sel_gcpos_intersect)

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

plot_targets(simid_sel_gcpos_intersect)


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim2_altrim.rds"))

lhsdt <- rbindlist(
  lapply(lhs_groups, as.data.table, keep.rownames = TRUE),
  idcol = "simid"
)

lhsdt[, simid := sprintf("%04d", simid)][]
setnames(lhsdt, c("V1", "V2"), c("input", "value"))

lhsdt[simid %in% simid_sel_gcpos_intersect] %>%
  ggplot(aes(x = input, y = value, group = simid)) +
  geom_line() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

sim2_priors <- readRDS(here::here("inst/cal", "sim1_altrim_sel_lhs_limits.rds"))

## selected params, narrowed by VLS targets
vls_sel_inputs <- pull_params(
  simid_sel_vls_race
)[input %like% "RX_INIT|RX_REINIT"]

vls_sel_inputs[,
  input_group := stringr::str_extract(input, "RX_[A-Z]+")][]

vls_sel_inputs_w <- dcast(
  vls_sel_inputs,
  selection_group + simid ~ input,
  value.var = "value"
)[simid %in% simid_sel_hivprinc][]

cbind(
  vls_sel_inputs_w[, .(N_input_pull = .N), selection_group],
  N_simid_list = sapply(simid_sel_vls_race, length)
)

## selected params, narrowed by GC targets
gcpos_sel_inputs <- pull_params(
  simid_sel_gcpos_intersect
)[, -c("selection_group")][
  !(input %like% "SCALAR_AI|SCALAR_OI|HIV_RX_|TRANS_PROB_WHITE")]

gcpos_sel_inputs %>%
  ggplot(aes(x = value)) +
  geom_histogram() +
  geom_density() +
  facet_wrap(~ input, scales = "free")

gcpos_sel_inputs_w <- dcast(
  gcpos_sel_inputs,
  simid ~ input,
  value.var = "value"
)

ggcorrplot(cor(gcpos_sel_inputs_w[, -c(1:2)]), type = "upper")


################################################################################
## NARROW PRIORS FOR NEXT ROUND ##
################################################################################

## Prior narrowing based on VLS targets
narrowed_priors_vls <- input_quantiles(vls_sel_inputs, "RX_")

## Prior narrowing based on GC targets
narrowed_priors_gc <- gcpos_sel_inputs[, .(
 ll = quantile(value, 0.25),
 ul = quantile(value, 0.75)
), input]
##[!(input %like% "_GC|GC_|ASYMP|SYMP|OI_ACT")]

narrowed_mid80_gc <- gcpos_sel_inputs[, .(
  ll = quantile(value, 0.1),
  ul = quantile(value, 0.9)
), input]

## Compare narrowing thresholds based on GC targets
lc <- merge(
  merge(
    narrowed_priors_gc,
    narrowed_mid80_gc,
    by = "input",
    suffixes = c("q", "m")
  ),
  sim2_priors,
  by = "input"
)

lcm <- melt(lc, id.vars = "input")

tlabs <- c("Sim2 Priors", "Mid80", "IQR")
lcm[, type := fcase(
        variable %like% "q", tlabs[3],
        variable %like% "m", tlabs[2],
        variable %like% "s2", tlabs[1]
      )][]

lcm[, type := factor(type, levels = tlabs)]
lcm[, limit := fifelse(variable %like% "ll", "lower", "upper")][]

lcm[!(input %like% "HIV_RX_|TRANS_PROB_WHITE|_AI|_OI")] %>%
  ggplot(aes(x = type, y = value)) +
  geom_line(aes(group = limit)) +
  geom_point() +
  facet_wrap(~input, scales = "free_y") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

narrowed_priors <- rbind(
  narrowed_mid80_gc,
  narrowed_priors_vls,
  use.names = FALSE
)

new_priors_merge <- merge(
  sim2_priors,
  narrowed_priors,
  by = "input",
  all.x = TRUE
)

npm <- melt(
  new_priors_merge,
  measure.vars = c("s2_ll", "s2_ul", "ll", "ul")
  )[, ":=" (
      round = fifelse(variable %like% "s2", "sim2", "sim3"),
      limit = fifelse(variable %like% "ll", "lower", "upper")
    )][]

npm[!(input %like% "SCALAR_AI|SCALAR_OI|RX_HALT|TRANS_PROB_WHITE")] %>%
  ggplot(aes(x = round, y = value)) +
  geom_line(aes(group = limit)) +
  geom_point(shape = 21, fill = "white", size = 2) +
  facet_wrap(~ input, scales = "free_y", ncol = 4)

## finalize new priors
new_priors <- copy(new_priors_merge)

new_priors[
  !is.na(ll) & !is.na(ul),
  ":=" (
    s2_ll = ll,
    s2_ul = ul
  )][]

head(new_priors)

new_priors[, ":=" (ll = NULL, ul = NULL)][]

setnames(new_priors, c("s2_ll", "s2_ul"), c("s3_ll", "s3_ul"))
new_priors


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(
  list(
    simid_sel_gcpos_intersect = simid_sel_gcpos_intersect,
    simid_sel_vls_race = simid_sel_vls_race
  ),
  here::here("inst", "cal", "sim2_altrim_simid_sel.rds")
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
saveRDS(new_priors, here::here("inst", "cal", "sim2_altrim_sel_lhs_limits.rds"))
