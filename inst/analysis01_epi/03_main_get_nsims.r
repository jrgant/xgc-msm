################################################################################
## SETUP ##
################################################################################

library(pacman)

p_load(
  EpiModelHIV,
  xgcmsm,
  ggplot2,
  ggthemes,
  scico
)

netstats <- get_est("netstats")
epistats <- get_est("epistats")
targets  <- readRDS(here::here("est", "caltargets.Rds"))

sim <- readRDS(
  here::here("inst", "analysis01_epi", "main_epi_2021-07-08.rds")
)

theme_set(theme_base(base_size = 22))


################################################################################
## CALCULATE DERIVED MEASURES ##
################################################################################

# HIV diagnosis prevalence among the infected.
# Need this to compare to calibration target.

# ... by race/ethnicity
rdxlabs <- c("", ".B", ".H", ".O", ".W")
race.dx <- paste0("i.prev.dx.inf", rdxlabs)

sim[, (race.dx) := lapply(1:5, function(x) {
  get(paste0("i.prev.dx", rdxlabs[x])) / get(paste0("i.prev", rdxlabs[x]))
})][]

# ... by age group
age.dx <- paste0("i.prev.dx.inf.age", 1:5)

sim[, (age.dx) := lapply(1:5, function(x) {
  get(paste0("i.prev.dx.age.", x)) / get(paste0("i.prev.age.", x))
})][]

# PrEP coverage
race.prep.cov <- paste0("prepCov", rdxlabs)

sim[, (race.prep.cov) := lapply(1:5, function(x) {
  get(paste0("prepCurr", rdxlabs[x])) / get(paste0("prepElig", rdxlabs[x]))
})][]

# HIV incidence per 100,000 MSM/year, per population (not at-risk).
# See NOTE in inst/cal/calval_targets.r about aligning this output with the
# the target source (Singh et al., 2018, Annals of Internal Medicine).
#
# Here, we need to calculate that value within every week, hence the
# multiplier 5200.

# ... overall
sim[, ir100.pop := incid / num * 5200][]

# ... by race/ethnicity
race.ir100.pop <- paste0("ir100.pop", rdxlabs)[-1]

sim[, (race.ir100.pop) := lapply(1:4, function(x) {
  get(paste0("incid", rdxlabs[x])) / get(paste0("num", rdxlabs[x]))
})][]

# ... demographics
sim[, paste0("prop", rdxlabs[-1]) := lapply(1:4, function(x) {
  get(paste0("num", rdxlabs[x + 1])) / num
})][]

sim[, paste0("prop.age.", 1:5) := lapply(1:5, function(x) {
  get(paste0("num.age.", x)) / num
})][]


################################################################################
## GET VARIANCES at VARIOUS NSIMS ##
################################################################################

set.seed(2001)
v1000_samp <- seq_len(sim[, max(simid)])
v50_samp    <- sample(v1000_samp,  50)
v100_samp   <- sample(v1000_samp, 100)
v250_samp   <- sample(v1000_samp, 250)
v500_samp   <- sample(v1000_samp, 500)
v750_samp   <- sample(v1000_samp, 750)

vlist <- list(
  v50   = v50_samp,
  v100  = v100_samp,
  v250  = v250_samp,
  v500  = v500_samp,
  v750  = v750_samp,
  v1000 = v1000_samp
)

outputcols <- names(sim)[!names(sim) %in% c("at", "simid")]

sigma_dt <- rbindlist(
  lapply(setNames(vlist, names(vlist)), function(.x) {
    melt(
      sim[simid %in% .x,
          lapply(.SD, sd),
          .SDcols = outputcols,
          by = at],
      id.vars = "at"
    )
  }), idcol = "nsims"
)

sigma_dt[, nsims := factor(nsims, levels = names(vlist))]

## Plot variance in model outputs by time (at) and number of sims used to
## calculate variance
plot_sigma <- function(varsel, data = sigma_dt) {

  p <- ggplot(data[variable %in% varsel], aes(x = at, y = value)) +
    geom_line(aes(color = nsims)) +
    scale_color_scico_d(palette = "roma", end = 0.75) +
    ylab("mean time-specific sigma") +
    labs(caption = "v1000, across all 1,000 test simulations.\nv500 and v750 are single realizations of a sample of indicated size (w/out replacement) from 1,000 test simulations.\nTherefore, this plot does not display sampling variability.")
  theme_base(base_size = 24)

  if (length(varsel) > 1) {
    p <- p +
      facet_wrap(~variable)
  }


  return(p + theme(plot.caption = element_text(size = 10, face = "italic")))
}

sigma_dt_sub <- sigma_dt[nsims %in% c("v250", "v500", "v750", "v1000")]

plot_sigma("num", data = sigma_dt_sub)
plot_sigma(paste0("ir100", rdxlabs), sigma_dt_sub)
plot_sigma(paste0("i.prev", rdxlabs), sigma_dt_sub)
plot_sigma(paste0("ir100.gc", rdxlabs), sigma_dt_sub)
plot_sigma(paste0("ir100.rgc", rdxlabs), sigma_dt_sub)
plot_sigma(paste0("ir100.ugc", rdxlabs), sigma_dt_sub)
plot_sigma(paste0("ir100.pgc", rdxlabs), sigma_dt_sub)
plot_sigma(paste0("i.prev.dx.inf.age", 1:5), sigma_dt_sub)
plot_sigma(paste0("incid.age.", 1:5, ".rgc"), sigma_dt_sub)
plot_sigma(paste0("incid.age.", 1:5, ".ugc"), sigma_dt_sub)
plot_sigma(paste0("incid.age.", 1:5, ".pgc"), sigma_dt_sub)

sigma_dt

# NOTE
# Everything over v500 seems to provide similar standard deviation estimates.
# n = 500, good tradeoff between computation time and output variance.

nsim_per_scenario <- 500

saveRDS(
  nsim_per_scenario,
  here::here("inst", "analysis01_epi", "nsims_per_scenario.rds")
)
