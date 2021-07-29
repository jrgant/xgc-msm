################################################################################
## DEPENDENCY: VISUALIZE DISTRIBUTIONS OF TARGETS ACROSS ALL SIMS ##
################################################################################

## Prior to any selection of simulations (except for those chosen based on
## overall population size), see whether there is overlap with calibration
## targets.

epi5k <- melt(
  epi_mn_selnum,
  id.vars = "simid",
  value.name = "mn_lastyr"
)

precal_plot_title <- "Pre-calibration distribution"

## ... HIV prevalence
ps_ct_hivprev <- pbox(
  "i.prev",
  targets[target == "ct_hiv_prev", value],
  targets[target == "ct_hiv_prev", ll95],
  targets[target == "ct_hiv_prev", ul95],
  data = epi5k,
  plot_title = precal_plot_title
)

## ... HIV incidence
ps_ct_hivinc <- pbox(
  "ir100.pop",
  targets[target == "ct_hiv_incid_per100pop", value],
  targets[target == "ct_hiv_incid_per100pop", ll95],
  targets[target == "ct_hiv_incid_per100pop", ul95],
  data = epi5k,
  plot_title = precal_plot_title
)

## ... Viral suppression, by age group
ps_ct_vsupp_age <- pbox(
  paste0("cc.vsupp.age", 1:5),
  targets[target == "ct_vls_pr_byage", value],
  targets[target == "ct_vls_pr_byage", ll95],
  targets[target == "ct_vls_pr_byage", ul95],
  data = epi5k,
  plot_title = precal_plot_title
)

## ... Viral suppression, by race/ethnicity
ps_ct_vsupp_race <- pbox(
  paste0("cc.vsupp.", c("B", "H", "O", "W")),
  targets[target == "ct_vls_pr_byrace", value],
  targets[target == "ct_vls_pr_byrace", ll95],
  targets[target == "ct_vls_pr_byrace", ul95],
  data = epi5k,
  plot_title = precal_plot_title
)

## ... PrEP coverage
ps_ct_prep <- pbox(
  race.prep.cov[-1],
  targets[target == "ct_prep", value],
  targets[target == "ct_prep", ll95],
  targets[target == "ct_prep", ul95],
  data = epi5k,
  plot_title = precal_plot_title
)

## ... Proportion anatomic sites tested for GC at STI clinic
ps_ct_gctested <- pbox(
  sprintf("prop.%s.tested", c("rect", "ureth", "phar")),
  targets[target == "ct_prop_anatsite_tested", value],
  targets[target == "ct_prop_anatsite_tested", ll95],
  targets[target == "ct_prop_anatsite_tested", ul95],
  data = epi5k,
  plot_title = precal_plot_title
)

## ... GC test positivity among tested
ps_ct_gcpos <- pbox(
  sprintf("prob.%s.tested", c("rGC", "uGC", "pGC")),
  targets[target == "ct_prop_anatsite_pos", value],
  targets[target == "ct_prop_anatsite_pos", ll95],
  targets[target == "ct_prop_anatsite_pos", ul95],
  data = epi5k,
  plot_title = precal_plot_title
)


################################################################################
## ABSOLUTE DIFFERENCES BETWEEN MODEL OUTPUT AND TARGETS ##
################################################################################

t_hiv_inc <- get_absdiff("ct_hiv_incid_per100pop", "ir100.pop")
t_hiv_prv <- get_absdiff("ct_hiv_prev", "i.prev")

# HIV diagnosis prevalence among infected, by race/ethnicity
t_hivdx_race <- get_absdiffv(
  "ct_hivdx_pr_byrace",
  rlabs,
  "i.prev.dx.inf.%s"
)

# ratio of HIV diagnosis probability among infected, by race/ethnicity
# (ref. White)
t_hivdx_race_rel <- get_absdiffv(
  "ct_hivdx_byrace_relative",
  rlabs[1:3],
  "i.prev.dx.inf.%s.rel.ref.W"
)
# HIV diagnosis prevalence among infected, by age
t_hivdx_age <- get_absdiffv(
  "ct_hivdx_pr_byage",
  paste0("age", 1:5),
  "i.prev.dx.inf.%s"
)

# ratio of HIV diagnosis probability among infected, by race/ethnicity
# (ref. White)
t_hivdx_age_rel <- get_absdiffv(
  "ct_hivdx_byage_relative",
  paste0("age", 1:4),
  "i.prev.dx.inf.%s.rel.ref.age5"
)

# viral suppression among HIV-diagnosed, by race/ethnicity
t_vls_race <- get_absdiffv(
  "ct_vls_pr_byrace",
  rlabs,
  "cc.vsupp.%s"
)

# ratio of HIV viral suppression probability, by race/ethnicity
t_vls_race_rel <- get_absdiffv(
  "ct_vls_byrace_relative",
  rlabs[1:3],
  "cc.vsupp.%s.rel.ref.W"
)

# viral suppression among HIV-diagnosed, by age group
t_vls_age <- get_absdiffv(
  "ct_vls_pr_byage",
  paste0("age", 1:5),
  "cc.vsupp.%s"
)

# ratio of HIV viral suppression probability, by age group
# (ref. age group 5)
t_vls_age_rel <- get_absdiffv(
  "ct_vls_byage_relative",
  paste0("age", 1:4),
  "cc.vsupp.%s.rel.ref.age5"
)

# PrEP coverage among indicated, by race/ethnicity
t_prepcov <- get_absdiffv("ct_prep", rlabs, "prepCov.%s")

# gonorrhea testing in clinic, by anatomic site
t_gc_anat_test <- get_absdiffv(
  "ct_prop_anatsite_tested",
  c("rect", "ureth", "phar"),
  "prop.%s.tested"
)

# gonorrhea positivity in clinic, among tested anatomic sites
t_gc_pos <- get_absdiffv(
  "ct_prop_anatsite_pos",
  c("rGC", "uGC", "pGC"),
  "prob.%s.tested"
)

## Concatenate all t_ objects into a single list for input into
## rbindlist. The for loop exists to unnest nested lists.
tobjs <- ls(pattern = "^t_")

tsum <- data.table(
  obj = tobjs,
  class = sapply(tobjs, function(.x) class(get(.x))[1]),
  length = sapply(tobjs, function(.x) {
    if (class(get(.x))[1] == "data.table") {
      1
    } else {
      length(get(.x))
    }
  })
)

tmp_t <- list()
slot <- 1

for (i in seq_len(tsum[, .N])) {
  tmp <- get(tobjs[i])
  if (class(tmp)[1] == "data.table") {
    vn <- names(tmp)[2]
    names(tmp)[2] <- "output"
    tmp_t[[slot]] <- tmp
    names(tmp_t)[slot] <- vn
    slot <- slot + 1
  } else if (is.list(tmp)) {
    for (j in seq_len(length(tmp))) {
      vn <- names(tmp[[j]])[2]
      names(tmp[[j]])[2] <- "output"
      tmp_t[[slot]] <- tmp[[j]]
      names(tmp_t)[slot] <- vn
      slot <- slot + 1
    }
  } else {
    stop("Burn it all down!")
  }
}

length(tmp_t)

out_vs_targ_pre <- rbindlist(tmp_t, idcol = "variable")

out_iqr <-
  out_vs_targ_pre[,
    .(
      out_q25 = quantile(output, 0.25),
      out_q75 = quantile(output, 0.75)
    ),
    by = variable
  ]

out_vs_targ <- merge(out_vs_targ_pre, out_iqr, by = "variable")
out_vs_targ[, target_within_iqr :=
  as.numeric(between(target_val, out_q25, out_q75))][]

iqr_capture <- out_vs_targ[, .N, .(variable, target_within_iqr)]

out_vs_targ[, .(
  N = .N,
  min = min(pct_diff),
  q25 = quantile(pct_diff, 0.25),
  median = median(pct_diff),
  mean = mean(pct_diff),
  q75 = quantile(pct_diff, 0.75),
  max = max(pct_diff)
), by = c("variable", "subgroup")]

