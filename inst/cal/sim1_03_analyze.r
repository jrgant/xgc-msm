###############################################################################
## SETUP ##
################################################################################

picksim <- "sim1"
source(here::here("inst", "cal", "sim0.0_setup.r"))
ls()
sim_epis[sim_epis %like% batch_date]


################################################################################
## VISUALIZE DISTRIBUTIONS OF TARGETS ACROSS ALL SIMS ##
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

out_vs_targ %>%
  ggplot(aes(x = diff)) +
  geom_histogram(color = "white") +
  geom_vline(aes(xintercept = 0), color = "red") +
  facet_wrap(~variable, scales = "free_x") +
  theme_base()

targsum_limits <- out_vs_targ[, .(
  N_within_10pct = sum(output_within_10pct),
  N_within_5pts = sum(output_within_5pts),
  N_within_10pts = sum(output_within_10pts),
  N_within_cl = sum(output_in_targcl)
),
by = .(variable, subgroup)
]

targsum_limits

sad <- function(pattern, showpairs = FALSE) {
  tmp <- out_vs_targ[variable %like% pattern]

  cat(
    "Variables selected:",
    paste(tmp[, unique(variable)], collapse = ", "),
    "\n"
  )

  ## wide version of output data
  wide <- dcast(
    tmp[, .(simid, variable, output)],
    simid ~ variable,
    value.var = "output"
  )

  ## stratum-specific standardized differences (target vs. output)
  strata_diffs <- dcast(
    tmp[, .(simid, variable, abs_std_diff)],
    simid ~ variable,
    value.var = "abs_std_diff"
  )

  ## summaries of absolute standardized differences across all selected strata
  ## (target vs. output)
  sum_diffs <- tmp[strata, .(
    sum_abs_stdiff = sum(abs_std_diff),
    range_abs_stdiff = max(abs_std_diff) - min(abs_std_diff)
  ), by = simid][order(sum_abs_stdiff)]

  if (showpairs) {
    varn <- length(wide)
    pairs(as.data.frame(wide)[2:varn], upper.panel = NULL)
  }

  out <- list(sum_diffs = sum_diffs, strata_diffs = strata_diffs, wide = wide)
}

simid_sel_vls <- lapply(
  setNames(rlabs, rlabs),
  function(.x) {
    out_vs_targ[variable == paste0("cc.vsupp.", .x) & output_within_5pts, simid]
  }
)

simid_sel_prep <- lapply(
  setNames(rlabs, rlabs),
  function(.x) {
    out_vs_targ[variable == paste0("prepCov.", .x) & output_within_5pts, simid]
  }
)

episel <- epi_noNA[simid %in% unlist(list(simid_sel_vls, simid_sel_prep))]

out_vs_targ[simid %in% unlist(simid_sel_vls), .(
  min = min(output),
  max = max(output),
  target = mean(target_val)
), variable]

quicktarget("cc.vsupp.%s", simid_sel_vls)
quicktarget("prepCov.%s", simid_sel_prep)


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim1.rds"))

vls_sel_inputs <- pull_params(simid_sel_vls)[input %like% "RX_INIT|RX_REINIT"]
vls_sel_inputs[, input_group := stringr::str_extract(input, "RX_[A-Z]+")][]

vls_sel_inputs_w <- dcast(
  vls_sel_inputs,
  selection_group + simid ~ input,
  value.var = "value"
)

prep_sel_inputs <- pull_params(simid_sel_prep)[input %like% "PREP"]

prep_sel_inputs_w <- dcast(
  prep_sel_inputs,
  selection_group + simid ~ input,
  value.var = "value"
)

# Function to check to make sure we've captured all selected simids.
check_input_simids <- function(x, simid_list, simid_inputs) {
  data.table(
    captured = sort(simid_inputs[selection_group == x, unique(simid)]),
    sel = sort(simid_list[[x]])
  )[, .N, .(captured, sel)][, all(N) == 1]
}

# check selected VLS inputeter sets
sapply(
  rlabs,
  check_input_simids,
  simid_list = simid_sel_vls,
  simid_inputs = vls_sel_inputs
)

# check selected PrEP inputeter sets
sapply(
  rlabs,
  check_input_simids,
  simid_list = simid_sel_prep,
  simid_inputs = prep_sel_inputs
)

## View distributions of RX_INIT and RX_REINIT for each selection group,
## meaning simid sets chosen for matching race/ethnicity-specific viral
## suppression. We're interested in seeing the how the corresponding
## input parameter distribution (e.g., HIV_RX_INIT_PROB_BLACK) within
## the B selection group. Here, the REINIT parameters are much stronger
## drivers of VLS.
ggplot(vls_sel_inputs, aes(input, value)) +
  geom_boxplot() +
  facet_grid(selection_group ~ input_group, scale = "free_x") +
  theme_base(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90))

corplots_vls <- lapply(rlabs, function(.x) {
  tmp <- vls_sel_inputs_w[selection_group == .x, 3:ncol(vls_sel_inputs_w)]
  ggcorrplot(
    cor(tmp, method = "spearman"),
    type = "lower",
    lab = TRUE,
    title = sprintf("Selection Group: %s", .x),
    p.mat = cor_pmat(tmp, method = "spearman"),
    sig.level = 0.05
  )
})

gridExtra::grid.arrange(
  corplots_vls[[1]],
  corplots_vls[[2]],
  corplots_vls[[3]],
  corplots_vls[[4]]
)

corplots_prep <- lapply(rlabs, function(.x) {
  tmp <- prep_sel_inputs_w[selection_group == .x, 3:ncol(prep_sel_inputs_w)]
  ggcorrplot(
    cor(tmp, method = "spearman"),
    type = "lower",
    lab = TRUE,
    title = sprintf("Selection Group: %s", .x),
    p.mat = cor_pmat(tmp, method = "spearman"),
    sig.level = 0.05
  )
})

gridExtra::grid.arrange(
  corplots_prep[[1]],
  corplots_prep[[2]],
  corplots_prep[[3]],
  corplots_prep[[4]]
)

## Function that takes the long version of inputeter subsets from
## simid selections and extracts the 25% and 75% quantiles for use as prior
## limits in the next round of calibration.
input_quantiles <- function(data, inputstring, ql = 0.25, qu = 0.75) {
  inputsub <- data[
    input %like% inputstring,
    .(
      q25 = quantile(value, ql),
      q75 = quantile(value, qu)
    ),
    .(selection_group, input)
    ## this step matches the race/eth selection group to the corresponding
    ## input inputeter
  ][selection_group == substring(str_extract(input, "(?<=_)[A-Z]+$"), 1, 1)]

  inputsub[, -c("selection_group")]
}

narrow_rx_init <-   input_quantiles(vls_sel_inputs, "RX_INIT")
narrow_rx_reinit <- input_quantiles(vls_sel_inputs, "RX_REINIT")
narrow_prep_discont <- input_quantiles(prep_sel_inputs, "DISCONT")

narrowed <- rbindlist(
  list(narrow_rx_init, narrow_rx_reinit, narrow_prep_discont)
)

## Create a new set of priors for sim2
sim1_priors <- readRDS(here::here("burnin", "cal", "sim1", "sim1_priors.rds"))

new_priors <- merge(sim1_priors, narrowed, by = "input", all.x = TRUE)
new_priors <- new_priors[
  !is.na(q25) & !is.na(q75), ":="(s1_ll = q25, s1_ul = q75)][, -c("q25", "q75")]

setnames(new_priors, c("s1_ll", "s1_ul"), c("s2_ll", "s2_ul"))

## Add a prior for the act.stopper.prob parameter, which governs the proportion
## of GC-symptomatic or -treated patients cease sexual activity while these
## conditions hold. The limits are somewhat arbitrary, but the empirical
## estimates from Beck/Kramer of 0.8 are from 1980 and applied only to GC
## symptoms and reduction of sex.
new_priors <- rbindlist(
  list(
    new_priors,
    data.table(s2_ll = 0.2, s2_ul = 1, input = "ACT_STOPPER_PROB")
  ),
  use.name = TRUE
)


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(
  list(
    vls = simid_sel_vls,
    prep = simid_sel_prep
  ),
  here::here("inst", "cal", "sim1_simid_sel.rds")
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
saveRDS(new_priors, here::here("inst", "cal", "sim1_sel_lhs_limits.rds"))
