################################################################################
## SETUP ##
################################################################################

picksim <- "sim3"

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
ps_ct_gctested + theme(axis.text.x = element_text(angle = 90))
ps_ct_gcpos + theme(axis.text.x = element_text(angle = 90))


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
length(simid_sel_vls_age_intersect)

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
length(simid_sel_gctest_intersect)

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


plot_targets(simid_sel_gcpos_intersect)
plot_targets(intersect(simid_sel_gcpos_intersect, simid_sel_hivprinc))
plot_targets(intersect(simid_sel_gcpos_intersect, simid_sel_vls_race_intersect))

varpats <- paste0(c(
  "i.prev$",
  "ir100.pop",
  "cc.vsupp.[A-Z]{1}$",
  "dx.inf.[A-Z]{1}$",
  "i.prev.dx.inf.age[1-5]{1}$",
  "prob.[A-Za-z]{3}.tested"
), collapse = "|")

mase <- out_vs_targ[
  variable %like% varpats,
  .(mase = sum(abs_std_diff)),
  simid
]

mase %>%
  ggplot(aes(x = mase)) +
  geom_histogram(color = "white", fill = "gray80")

mase_fitsqt <- lapply(
  c(0.05, 0.01, 0.0025),
  function(.x) {
    mase[mase <= quantile(mase, .x), simid]
})

best_fits10 <- mase[order(mase), simid][1:10]

plot_targets(mase_fitsqt[[1]]) # smallest 5% of mase
plot_targets(mase_fitsqt[[2]]) # smallest 1% of mase
plot_targets(mase_fitsqt[[3]]) # smallest 0.25% of mase
plot_targets(best_fits10) # smallest 10 mases


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim3.rds"))

lhsdt <- rbindlist(
  lapply(lhs_groups, as.data.table, keep.rownames = TRUE),
  idcol = "simid"
)

lhsdt[, simid := sprintf("%04d", simid)][]
setnames(lhsdt, c("V1", "V2"), c("input", "value"))

input_dist <- function(simids) {
  lhsdt[simid %in% simids] %>%
    ggplot(aes(x = value)) +
    geom_histogram(fill = "gray70", color = "white") +
    geom_density() +
    facet_wrap(~input, scales = "free")
}

input_dist(mase_fitsqt[[1]])
input_dist(mase_fitsqt[[2]])
input_dist(mase_fitsqt[[3]])

main_analysis_inputs <- lapply(
  setNames(best_fits10, paste0("s3_", best_fits10)),
  function(.x) {
    lhs_groups[[as.numeric(.x)]]
  })


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(
  main_analysis_inputs,
  here::here("inst", "cal", "main_analysis_inputs.rds")
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
