################################################################################
## SETUP ##
################################################################################

## Setup script used in subsequent sim analyses:
##   sim1_03_analyze.r

library(pacman)

p_load(
  data.table,
  EpiModel,
  ggplot2,
  ggthemes,
  ggcorrplot,
  scico,
  xgcmsm,
  magrittr,
  lubridate,
  stringr
)

extrafont::loadfonts()
theme_set(theme_tufte(base_size = 22))

## Get date from system environment if set. Otherwise, find date of most recent
## sim1_epi run.
sys_batch_date <- Sys.getenv("BATCH_DATE")

if (!exists("picksim")) {
  stop("Must specify picksim before sourcing sim1_02_setup.r")
}

sim_epis <- list.files(
  here::here("burnin", "cal", picksim),
  pattern = "sim1_epi",
  full.names = TRUE
)

if (sys_batch_date != "") {
  batch_date <- Sys.getenv("BATCH_DATE")
} else {
  batch_date <- max(ymd(str_extract(sim_epis, "[0-9]{4}-[0-9]{2}-[0-9]{2}")))
}

## Load calibration targets, netstats, and formatted simulated output.
epi      <- readRDS(sim_epis[sim_epis %like% batch_date])
netstats <- get_est("netstats")
targets  <- readRDS(here::here("est", "caltargets.Rds"))

## Create directory for plot files (if necessary)
simpath <- here::here("inst", "cal", "calfigs", picksim)
if (!dir.exists(simpath)) dir.create(simpath)

figpath <- here::here(simpath, batch_date)
if (!dir.exists(figpath)) dir.create(figpath)

cat("Figures are being saved to:\n", figpath, "\n")

## Helper function to save plots.
## Just a thin wrapper around ggsave() with some dimension presets
psave <- function(f, p, w = 10, h = 10, caldir = figpath) {
  # f = filename (no extension)
  # p = plot object
  # w = width
  # h = height
  # caldir = directory in which to save calibration figures
  ggsave(
    here::here(caldir, paste0(f, ".pdf")),
    p, width = w, height = h, units = "in"
  )
}


## Plot label sets.
racelabs <- c("", ", Black", ", Hispanic", ", Other", ", White")
anatlabs <- c("", ", Rectal", ", Urethral", ", Pharyngeal")


## This function produces boxplots of selected variables.
## varnames: varnames output by the episims
## targets: optional set of calibration targets
## data: longform dataset of simulated outputs; the value variable must be named
##       mn_lastyr
pbox <- function(varnames, targets = NULL,
                 targets_ll = NULL, targets_ul = NULL, data = epi_hivgc) {

  cdat <- data[variable %in% varnames]
  mdat <- cdat[, .(mean = mean(mn_lastyr)), by = variable]

  basep <-
    ggplot(data = cdat, aes(variable, mn_lastyr)) +
    geom_boxplot(color = "gray", outlier.size = 0.5) +
    geom_point(
      data = mdat, aes(variable, mean, fill = "Simulated mean"),
      shape = 21, color = "black", size = 6
    )

  if (!is.null(targets)) {
    tdat <- data.table(
      target = mdat[, variable],
      targval = targets,
      targnames = names(targets)
    )

    if (!is.null(targets_ll) & !is.null(targets_ul)) {
      tdat[, ":=" (
        targval_ll = targets_ll,
        targval_ul = targets_ul
      )]
    }

    basep <- basep +
      geom_pointrange(
        data = tdat,
        aes(
          x = target, y = targval,
          ymin = targval_ll, ymax = targval_ul,
          fill = "Target"
        ),
        shape = 21, color = "black", size = 1
      )

  }

  basep +
    scale_fill_scico_d(name = "Value", begin = 0.7) +
    theme_tufte(base_family = "sans", base_size = 22)
}


################################################################################
## CALCULATE DERIVED MEASURES ##
################################################################################

# HIV diagnosis prevalence among the infected.
# Need this to compare to calibration target.

# ... by race/ethnicity
rdxlabs <- c("", ".B", ".H", ".O", ".W")
race.dx <- paste0("i.prev.dx.inf", rdxlabs)

epi[, (race.dx) := lapply(1:5, function(x) {
  get(paste0("i.prev.dx", rdxlabs[x])) / get(paste0("i.prev", rdxlabs[x]))
})][]

# ... by age group
age.dx <- paste0("i.prev.dx.inf.age", 1:5)

epi[, (age.dx) := lapply(1:5, function(x) {
  get(paste0("i.prev.dx.age.", x)) / get(paste0("i.prev.age.", x))
})][]

# PrEP coverage
race.prep.cov <- paste0("prepCov", rdxlabs)

epi[, (race.prep.cov) := lapply(1:5, function(x) {
  get(paste0("prepCurr", rdxlabs[x])) / get(paste0("prepElig", rdxlabs[x]))
})][]

# HIV incidence per 100,000 MSM/year, per population (not at-risk).
# See NOTE in inst/cal/calval_targets.r about aligning this output with the
# the target source (Singh et al., 2018, Annals of Internal Medicine).
#
# Here, we need to calculate that value within every week, hence the
# multiplier 5200.

# ... overall
epi[, ir100.pop := incid / num * 5200][]

# ... by race/ethnicity
race.ir100.pop <- paste0("ir100.pop", rdxlabs)[-1]

epi[, (race.ir100.pop) := lapply(1:4, function(x) {
  get(paste0("incid", rdxlabs[x])) / get(paste0("num", rdxlabs[x]))
})][]


# ... demographics
epi[, paste0("prop", rdxlabs[-1]) := lapply(1:4, function(x) {
  get(paste0("num", rdxlabs[x + 1])) / num
})][]

epi[, paste0("prop.age.", 1:5) := lapply(1:5, function(x) {
  get(paste0("num.age.", x)) / num
})][]


################################################################################
## CLEAN UP OUTPUT AND CALCULATE LAST-YEAR AVERAGES ##
################################################################################

sumcols <- names(epi)[!names(epi) %like% "at|simid"]

## Mark NAs as 0
epi_noNA <- epi[, lapply(.SD, function(.x) fifelse(is.na(.x), 0, .x)),
                  .SDcols = c(sumcols, "at", "simid")]

## Get output averages over the last 52 weeks
epi_mn <- epi_noNA[
  at >= 3068, lapply(.SD, mean),
  by = simid, .SDcols = sumcols]


################################################################################
## SELECT INITIAL SIMS BASED ON POP. SIZE ##
################################################################################

# Pick parameter sets that keep us within 5% of population size N = 20,000
# (average over final burnin-in year)
pop_N <- 20000
epi_mn_selnum <- epi_mn[abs(pop_N - num) / pop_N <= 0.05]
epi_mn_selnum[, .(N = .N, P = .N / 5000)]
