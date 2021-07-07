################################################################################
## SETUP ##
################################################################################

picksim <- "sim2"
source(here::here("inst", "cal", "sim0_setup.r"))
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


################################################################################
## ABSOLUTE DIFFERENCES BETWEEN MODEL OUTPUT AND TARGETS ##
################################################################################

get_absdiff <- function(targname, output,
                        target_dt = targets, data = epi_mn_selnum) {

  cols <- c("simid", output)
  d <- data[, ..cols]

  d[, abs_diff := abs(get(output) - target_dt[target == targname, value])]

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
  min     = min(pct_diff),
  q25     = quantile(pct_diff, 0.25),
  median  = median(pct_diff),
  mean    = mean(pct_diff),
  q75     = quantile(pct_diff, 0.75),
  max     = max(pct_diff)
), by = c("variable", "subgroup")]

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

caltargs <- c(
  "i.prev",
  paste0("prepCov.", c("B", "W")),
  paste0("i.prev.dx.inf.age", 1:2)
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
  geom_point(aes(color = "Sim"), alpha = 0.3) +
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

## Write selected simulation ids to file.
saveRDS(simid_sel, here::here("inst", "cal", paste0(picksim, "_simid_sel.rds")))


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
      facet_wrap(~ get(varname))
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

rlabs <- c("B", "H", "O", "W")

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

targ_hivdx_btyreth_insel <- plotepi(
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

lapply(
  ls(pattern = "insel$"),
  function(x) psave(x, get(x))
)


################################################################################
## AVERAGES OVER LAST YEAR ##
################################################################################

episel_avg <- melt(
  episel,
  id.vars = "simid"
)[, .(mn_lastyr = mean(value)), .(simid, variable)][]

episel_avg[, length(unique(simid))]


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

sel_simid <- episel_avg[, unique(simid)]
lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim1.rds"))

sel_lhs <- lapply(
  setNames(sel_simid, paste0("simid_", sel_simid)),
  function(.x) as.data.table(lhs_groups[[.x]], keep.rownames = TRUE)
)

length(sel_lhs) == length(sel_simid)

sel_lhs_d <- rbindlist(sel_lhs, idcol = "simid")
sel_lhs_d[, simid := stringr::str_remove(simid, "simid_")]
setnames(sel_lhs_d, c("V1", "V2"), c("param", "value"))

sel_lhs_lims <-
  sel_lhs_d[, .(min = min(value), max = max(value), n = .N), by = param][]

saveRDS(sel_lhs_lims, here::here("inst", "cal", "sim1_sel_lhs_limits.rds"))

## Visualize correlations between input parameters among the selected model
## outputs.
sel_lhs_wide <- dcast(sel_lhs_d, simid ~ param)

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
    panel.grid.major.x  = element_line("gray60")
)

incorr

psave("incorr", incorr)
