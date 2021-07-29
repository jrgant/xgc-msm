################################################################################
## SETUP ##
################################################################################

picksim <- "sim3"
source(here::here("inst", "cal", "sim0.0_setup.r"))
ls()


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

## ... HIV prevalence
ps_ct_hivprev <- pbox(
  "i.prev",
  targets[target == "ct_hiv_prev", value],
  targets[target == "ct_hiv_prev", ll95],
  targets[target == "ct_hiv_prev", ul95],
  data = epi5k
)

## ... HIV incidence
ps_ct_hivinc <- pbox(
  "ir100.pop",
  targets[target == "ct_hiv_incid_per100pop", value],
  targets[target == "ct_hiv_incid_per100pop",  ll95],
  targets[target == "ct_hiv_incid_per100pop",  ul95],
  data = epi5k
)

## ... Viral suppression, by age group
ps_ct_vsupp_age <- pbox(
  paste0("cc.vsupp.age", 1:5),
  targets[target == "ct_vls_pr_byage", value],
  targets[target == "ct_vls_pr_byage", ll95],
  targets[target == "ct_vls_pr_byage", ul95],
  data = epi5k
)

## ... Viral suppression, by race/ethnicity
ps_ct_vsupp_race <- pbox(
  paste0("cc.vsupp.", c("B", "H", "O", "W")),
  targets[target == "ct_vls_pr_byrace", value],
  targets[target == "ct_vls_pr_byrace", ll95],
  targets[target == "ct_vls_pr_byrace", ul95],
  data = epi5k
)

## ... PrEP coverage
ps_ct_prep <- pbox(
  race.prep.cov[-1],
  targets[target == "ct_prep", value],
  targets[target == "ct_prep", ll95],
  targets[target == "ct_prep", ul95],
  data = epi5k
)

## ... Proportion anatomic sites tested for GC at STI clinic
ps_ct_gctested <- pbox(
  sprintf("prop.%s.tested", c("rect", "ureth", "phar")),
  targets[target == "ct_prop_anatsite_tested", value],
  targets[target == "ct_prop_anatsite_tested", ll95],
  targets[target == "ct_prop_anatsite_tested", ul95],
  epi5k
)

## ... GC test positivity among tested
ps_ct_gcpos <- pbox(
  sprintf("prob.%s.tested", c("rGC", "uGC", "pGC")),
  targets[target == "ct_prop_anatsite_pos", value],
  targets[target == "ct_prop_anatsite_pos", ll95],
  targets[target == "ct_prop_anatsite_pos", ul95],
  epi5k
)


################################################################################
## ABSOLUTE DIFFERENCES BETWEEN MODEL OUTPUT AND TARGETS ##
################################################################################

get_absdiff <- function(targname, output,
                        target_dt = targets, data = epi_mn_selnum) {

  cols <- c("simid", output)
  d <- data[, ..cols]

  d[, diff := get(output) - target_dt[target == targname, value]]
  d[, abs_diff := abs(diff)]

  d[, pct_diff :=
        round(abs_diff / target_dt[target == targname, value] * 100, 2)]

  targsub <- target_dt[target == targname]

  d[, output_in_targcl := as.numeric(
        between(
          get(output),
          targsub[target == targname, ll95],
          targsub[target == targname, ul95],
          NAbounds = NA
        )
      )]

  d[, output_within_10pct := as.numeric(pct_diff <= 10)]

  # NOTE
  ## For HIV incidence per 100 pop/yr, we want the threshold to correspond to
  ## the modeled scale. All other measures are expressed as probabilities and
  ## can use the 5/10 percentage point threshold.
  if (targname == "ct_hiv_incid_per100pop") {
    d[, output_within_5pts  := as.numeric(abs_diff <= 0.5)]
    d[, output_within_10pts := as.numeric(abs_diff <= 1)]
  } else {
    d[, output_within_5pts  := as.numeric(abs_diff <= 0.05)]
    d[, output_within_10pts  := as.numeric(abs_diff <= 0.10)]
  }

  d[, ":=" (
    target_val = targsub[, value],
    ll95 = targsub[, ll95],
    ul95 = targsub[, ul95],
    subgroup = targsub[, subgroups]
  )]
  dout <- d[order(abs_diff)]
  return(dout)
}

t_hiv_inc <- get_absdiff("ct_hiv_incid_per100pop", "ir100.pop")
t_hiv_prv <- get_absdiff("ct_hiv_prev", "i.prev")

# diagnosis prevalence among infected, by race
t_hiv_dx_race <- lapply(
  targets[target == "ct_hivdx_pr_byrace", subgroups],
  function(.x) {
    tmp <- targets[target == "ct_hivdx_pr_byrace" & subgroups == .x]
    get_absdiff(
      "ct_hivdx_pr_byrace",
      target_dt = tmp,
      paste0("i.prev.dx.inf.", .x)
    )
  }
)

# diagnosis prevalence among infected, by age
t_hiv_dx_age <- lapply(
  targets[target == "ct_hivdx_pr_byage", subgroups],
  function(.x) {
    tmp <- targets[target == "ct_hivdx_pr_byage" & subgroups == .x]
    get_absdiff(
      "ct_hivdx_pr_byage",
      target_dt = tmp,
      paste0("i.prev.dx.inf.", .x)
    )
  }
)

# viral suppression among HIV-diagnosed, by age group
t_vls_age <- lapply(
  targets[target == "ct_vls_pr_byage", subgroups],
  function(.x) {
    tmp <- targets[target == "ct_vls_pr_byage" & subgroups == .x]
    get_absdiff(
      "ct_vls_pr_byage",
      target_dt = tmp,
      paste0("cc.vsupp.", .x)
    )
  }
)

# viral suppression among HIV-diagnosed, by race
t_vls_race <- lapply(
  targets[target == "ct_vls_pr_byrace", subgroups],
  function(.x) {
    tmp <- targets[target == "ct_vls_pr_byrace" & subgroups == .x]
    get_absdiff(
      "ct_vls_pr_byrace",
      target_dt = tmp,
      paste0("cc.vsupp.", .x)
    )
  }
)

# PrEP coverage among indicated, by race/ethnicity
t_prepcov <- lapply(
  c("B", "H", "O", "W"),
  function(.x) {
    tmp <- targets[target == "ct_prep" & subgroups == .x]
    get_absdiff("ct_prep", target_dt = tmp, paste0("prepCov.", .x))
  }
)

# gonorrhea testing in clinic, by anatomic site
t_gc_anat_test <- lapply(
  targets[target == "ct_prop_anatsite_tested", subgroups],
  function(.x) {
    var <- sprintf("prop.%s.tested", .x)
    tmp <- targets[target == "ct_prop_anatsite_tested" & subgroups == .x]
    get_absdiff(
      "ct_prop_anatsite_tested",
      target_dt = tmp,
      var
    )
   }
)

# gonorrhea positivity in clinic, among tested anatomic sites
t_gc_pos <- lapply(
  targets[target == "ct_prop_anatsite_pos", subgroups],
  function(.x) {
    tmp  <- targets[target == "ct_prop_anatsite_pos" & subgroups == .x]
    lugc <- c("uGC" = "ureth", "rGC" = "rect", "pGC" = "phar")
    slug <- names(match.arg(.x, lugc))
    var  <- sprintf("prob.%s.tested", slug)
    get_absdiff(
      "ct_prop_anatsite_pos",
      target_dt = tmp,
      var
    )
  }
)

t_hiv_inc
t_hiv_prv
t_hiv_dx_race
t_hiv_dx_age
t_vls_race
t_vls_age
t_prepcov
t_gc_anat_test
t_gc_pos

tlist <- list(
  t_hiv_inc,
  t_hiv_prv,
  t_hiv_dx_race[[1]],
  t_hiv_dx_race[[2]],
  t_hiv_dx_race[[3]],
  t_hiv_dx_race[[4]],
  t_hiv_dx_age[[1]],
  t_hiv_dx_age[[2]],
  t_hiv_dx_age[[3]],
  t_hiv_dx_age[[4]],
  t_hiv_dx_age[[5]],
  t_vls_race[[1]],
  t_vls_race[[2]],
  t_vls_race[[3]],
  t_vls_race[[4]],
  t_vls_age[[1]],
  t_vls_age[[2]],
  t_vls_age[[3]],
  t_vls_age[[4]],
  t_vls_age[[5]],
  t_prepcov[[1]],
  t_prepcov[[2]],
  t_prepcov[[3]],
  t_prepcov[[4]],
  t_gc_anat_test[[1]],
  t_gc_anat_test[[2]],
  t_gc_anat_test[[3]],
  t_gc_pos[[1]],
  t_gc_pos[[2]],
  t_gc_pos[[3]]
)

names(tlist) <- sapply(tlist, function(.x) names(.x)[2])
for (i in seq_along(tlist)) {
  names(tlist[[i]])[2] <- "output"
}

tlist

## Compare simulation outputs to calibration targets (summary).
out_vs_targ <- rbindlist(tlist, idcol = "variable")

out_vs_targ[, .(
  N       = .N,
  min     = min(pct_diff),
  q25     = quantile(pct_diff, 0.25),
  median  = median(pct_diff),
  mean    = mean(pct_diff),
  q75     = quantile(pct_diff, 0.75),
  max     = max(pct_diff)
), by = c("variable", "subgroup")]

out_vs_targ %>%
  ggplot(aes(x = diff)) +
  geom_histogram(color = "white") +
  geom_vline(aes(xintercept = 0), color = "red") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_base()

targsum_limits <- out_vs_targ[, .(
  N_within_10pct  = sum(output_within_10pct),
  N_within_5pts   = sum(output_within_5pts),
  N_within_10pts  = sum(output_within_10pts),
  N_within_cl     = sum(output_in_targcl)
  ),
  by = .(variable, subgroup)]

targsum_limits


## Specify explicit targets and select simulation IDs
prop_tested_labs <- sprintf("prop.%s.tested", c("rect", "ureth", "phar"))
rlabs <- c("B", "H", "O", "W")

caltargs <- c(
  paste0("i.prev.dx.inf.age", 1:2),
  paste0("i.prev.dx.inf.", rlabs[c(1, 4)])
)

match_simids <- lapply(
  setNames(caltargs, caltargs),
  function(.x) {
    tmp <- out_vs_targ[variable == .x]
    if (.x == "i.prev") {
      tmp[output_within_5pts == 1, sort(simid)]
    } else {
      tmp[output_within_10pts == 1, sort(simid)]
    }
  })


names(match_simids)
match_simids

simid_sel <- Reduce(intersect, match_simids)
simid_sel

episel <- epi_noNA[simid %in% simid_sel]
episel

out_vs_targ[simid %in% simid_sel, .(
  min = min(output),
  max = max(output),
  target = mean(target_val)
), variable]

out_vs_targ[simid %in% simid_sel] %>%
  ggplot(aes(x = 0, y = output)) +
  geom_point(aes(color = "Sim"), alpha = 0.3, size = 5) +
  geom_pointrange(
    aes(x = 0, y = target_val,
        ymin = ll95, ymax = ul95,
        fill = "Target"),
    size  = 0.75,
    shape = 21,
    color = "red"
  ) +
  scale_color_scico_d() +
  facet_wrap(variable ~ subgroup) +
  theme_tufte(base_size = 20)


################################################################################
## Time Series Plots ##
################################################################################

plotepi <- function(var, target = NULL, type = "line",
                    data = epi, varname = NULL, line_alpha = 1) {

  keepv <- c("simid", "at", var)
  data <- data[, ..keepv]

  if (length(var) > 1) {
    data <- melt(
      data,
      id.vars = c("simid", "at"),
      measure.vars = var,
      variable.name = varname
    )

    p <- ggplot(
      data,
      aes(x = at, y = value, color = get(varname))
    ) +
      facet_wrap(~ get(varname), nrow = 1)
  } else {
    p <- ggplot(data, aes(x = at, y = get(var))) + ylab(var)
  }

  if (type == "line") p <- p + geom_line(aes(group = simid), alpha = line_alpha)
  if (type == "boxplot") p <- p + geom_boxplot()

  if (!is.null(target)) {
    p <- p +
      geom_hline(
        aes(yintercept = target),
        linetype = "dashed",
        color = "salmon"
      )
  }

  return(p)
}


spec_demog_num_insel <- plotepi("num", 20000, data = episel)

spec_demog_reth_insel <- plotepi(
  c(paste0("prop.", rlabs)),
  varname = "prop",
  line_alpha = 0.4,
  data = episel
  ) +
  geom_hline(
    data =
      data.table(
        prop = paste0("prop.", rlabs),
        targ = with(netstats$demog, c(num.B, num.H, num.O, num.W) / num)
      ),
    aes(yintercept = targ),
    color = "black"
  )

targ_hivprev_insel <- plotepi(
  c("i.prev", paste0("i.prev.", rlabs)),
  varname = "HIV prevalence",
  line_alpha = 0.4,
  data = episel
)

targ_hivdx_byreth_insel <- plotepi(
  c(paste0("i.prev.dx.inf.", rlabs)),
  varname = "hivdx_race",
  line_alpha = 0.4,
  data = episel
)

targ_hivdx_byage_insel <- plotepi(
  paste0("i.prev.dx.inf.age", 1:5),
  varname = "hivdx_age",
  line_alpha = 0.4,
  data = episel
  ) +
  geom_hline(
    data = data.table(
      hivdx_age = c(paste0("i.prev.dx.inf.age", 1:5)),
      value = targets[target == "ct_hivdx_pr_byage", value]
    ), aes(yintercept = value)
  )

targ_ir100pop_insel <- plotepi(
  "ir100.pop",
  targets[target == "ct_hiv_incid_per100pop", value],
  data = episel
  ) +
  geom_hline(
    aes(yintercept = targets[target == "ct_hiv_incid_per100pop", value]),
    color = "red"
  )


################################################################################
## AVERAGES OVER LAST YEAR ##
################################################################################

episel_avg <- melt(
  episel,
  id.vars = "simid"
)[, .(mn_lastyr = mean(value)), .(simid, variable)][]

episel_avg[, length(unique(simid))]

pbox(
  paste0("i.prev.dx.inf.age", 1:5),
  targets[target %like% "hivdx.*age", value],
  data = episel_avg
)

pbox(
  paste0("i.prev"),
  targets[target %like% "prev", value],
  targets[target %like% "prev", ll95],
  targets[target %like% "prev", ul95],
  data = episel_avg
)

pbox(
  paste0("ir100.pop"),
  targets[target %like% "inc", value],
  targets[target %like% "inc", ll95],
  targets[target %like% "inc", ul95],
  data = episel_avg
)

pbox(
  paste0("prepCov", rdxlabs[-1]),
  targets[target %like% "prep", value],
  targets[target %like% "prep", ll95],
  targets[target %like% "prep", ul95],
  data = episel_avg
)

pbox(
  paste0("cc.vsupp", rdxlabs[-1]),
  targets[target %like% "vls" & subgroups %in% c("B", "H", "O", "W"), value],
  targets[target %like% "vls" & subgroups %in% c("B", "H", "O", "W"), ll95],
  targets[target %like% "vls" & subgroups %in% c("B", "H", "O", "W"), ul95],
  data = episel_avg
)

pbox(
  paste0("cc.vsupp.age", 1:5),
  targets[target %like% "vls" & subgroups %like% "age", value],
  targets[target %like% "vls" & subgroups %like% "age", ll95],
  targets[target %like% "vls" & subgroups %like% "age", ul95],
  data = episel_avg
)



################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

sel_simid <- episel_avg[, unique(simid)]
lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim3.rds"))

sel_lhs <- lapply(
  setNames(as.numeric(sel_simid), paste0("simid_", sel_simid)),
  function(.x) as.data.table(lhs_groups[[.x]], keep.rownames = TRUE)
)

length(sel_lhs) == length(sel_simid)

sel_lhs_d <- rbindlist(sel_lhs, idcol = "simid")
sel_lhs_d[, simid := stringr::str_remove(simid, "simid_")]
setnames(sel_lhs_d, c("V1", "V2"), c("input", "value"))

sel_lhs_lims <-
  sel_lhs_d[, .(s3_ll = min(value), s3_ul = max(value), n = .N), by = input]

sim1_priors <- readRDS(here::here("burnin", "cal", "sim1", "sim1_priors.rds"))

## subsetting step here drops scalars or other inputs that were kept static
s13_lims <- merge(sel_lhs_lims, sim1_priors, by = "input")[s1_ul - s1_ll > 0]
s13_lims[input %like% "CONDOM", group := "CONDOM_EFF"][]
s13_lims[input %like% "LATE", group := "HIV_LATE_TESTER"][]
s13_lims[input %like% "HIV_RX_INIT", group := "HIV_RX_INIT"][]
s13_lims[input %like% "HIV_RX_REINIT", group := "HIV_RX_REINIT"][]
s13_lims[input %like% "2", group := "GC_TRANS_PROB"][]
s13_lims[input %like% "DURAT_NOTX", group := "GC_DURAT_NOTX"][]
s13_lims[input %like% "ASYMP_STITEST", group := "ASYMP_STITEST_PROB"][]
s13_lims[input %like% "PREP", group := "PREP_DISCONT"][]
s13_lims[input %like% "SYMPT_PROB", group := "GC_SYMPT_PROB"][]
s13_lims[input %like% "INFPR", group := "GC_TX_CURE_1WK_PROB"][]
s13_lims[input %like% "PROB_GC_SYMPT", group := "GC_SYMPT_STITEST_PROB"][]

s13_lims %>%
  ggplot(aes(x = 0, xend = 0)) +
  geom_segment(
    aes(y = s1_ll, yend = s1_ul, size = "sim1"),
    width = 0.1
  ) +
  geom_segment(
    aes(y = s3_ll, yend = s3_ul, size = picksim),
    width = 0.1
  ) +
  facet_wrap(~ group + input, scales  = "free_y", as.table = TRUE) +
  scale_size_manual(values = c(1, 5)) +
#  scale_color_scico_d(direction = -1, end = 0.5, palette = "bilbao") +
  xlim(c(-1, 1))


## Visualize correlations between input parameters among the selected model
## outputs.
sel_lhs_wide <- dcast(sel_lhs_d, simid ~ input)

### drop static parameters (no prior ranges submitted)
keeplhs_vars <-
  names(sel_lhs_wide)[!(names(sel_lhs_wide) %like% "SCALAR|RX_HALT|simid")]

lhs_cor <- cor(sel_lhs_wide[, ..keeplhs_vars], method = "spearman")

incorr <- ggcorrplot(
  lhs_cor,
  type = "upper",
  outline.color = "black"
  ) +
  scale_fill_scico(
    name = "Spearman\ncorrelation\n",
    limits = c(-1, 1),
    palette = "vik"
  ) +
  labs(
    caption = sprintf(
      "\nCorrelations between input values among %s parameter sets selected.",
      sel_lhs_d[, format(length(unique(simid)), big.mark = ",")]
    )
  ) +
  theme(
    text                = element_text(size = 18, family = "sans"),
    axis.text.x         = element_text(size = 12),
    axis.text.y         = element_text(size = 12),
    legend.key.height   = unit(1, "in"),
    legend.title        = element_text(face = "bold"),
    plot.title          = element_text(face = "bold"),
    plot.caption        = element_text(face = "italic"),
    panel.grid.major.y  = element_blank(),
    panel.grid.major.x  = element_line("gray90")
)

incorr


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(simid_sel, here::here("inst", "cal", paste0(picksim, "_simid_sel.rds")))

saveRDS(
  lhs_groups[as.numeric(simid_sel)],
  here::here("inst", "cal", "analysis_param_sets.rds")
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

## Write model output from selected simulations to files
lapply(
  ls(pattern = "insel$"),
  function(x) psave(x, get(x))
)

## Write plot showing correlations between inputs in selected simulations.
psave("incorr", incorr)

## Write limits of input parameter values across selected simualtions.
saveRDS(sel_lhs_lims, here::here("inst", "cal", "sim3_sel_lhs_limits.rds"))
