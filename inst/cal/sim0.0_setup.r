################################################################################
## SETUP ##
################################################################################

## Setup script used in subsequent sim analyses:
##   sim1_02_analyze.r
##   sim2_02_analyze.r
##   sim3_01_analyze.r
##   sim4_01_analyze.r

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

suppressMessages(extrafont::loadfonts())
theme_set(theme_tufte(base_size = 22))

## Get date from system environment if set. Otherwise, find date of most recent
## simXX_epi run.
sys_batch_date <- Sys.getenv("BATCH_DATE")

if (!exists("picksim")) {
  stop("Must specify picksim before sourcing sim0.0_setup.r")
}

sim_epis <- list.files(
  here::here("burnin", "cal", picksim),
  pattern = paste0(picksim, "_epi"),
  full.names = TRUE
)

if (sys_batch_date != "") {
  batch_date <- sys_batch_date
} else {

  # Standard dates will show as missing after switch to Unix date-time
  # labeling of files. To retrieve an older run labeled with standard
  # YYYY-MM-DD, set BATCH_DATE environment variable or sys_batch_date.
  batch_date <- max(
    as.numeric(str_extract(sim_epis, "(?<=_)[0-9]+(?=\\.)")),
    na.rm = TRUE
  )

  # save a human readable version of the batch date-time
  # January 1, 1970 is Unix time origin
  pretty_batch_date <- as.POSIXct(batch_date, origin = "1970-01-01")
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
psave <- function(f, p, w = 14, h = 14, ext = ".png", caldir = figpath) {
  # f = filename (no extension)
  # p = plot object
  # w = width
  # h = height
  # caldir = directory in which to save calibration figures
  ggsave(
    here::here(caldir, paste0(f, ext)),
    p, width = w, height = h, units = "in"
  )
}


## Label sets.
racelabs <- c("", ", Black", ", Hispanic", ", Other", ", White")
anatlabs <- c("", ", Rectal", ", Urethral", ", Pharyngeal")
rslugs <- c("B", "H", "O", "W")
ageslugs <- paste0("age", 1:5)
anatslugs <- c("rect", "ureth", "phar")
gcpos_slugs <- c("rGC", "uGC", "pGC")
gcprev_slugs <- gcpos_slugs[gcpos_slugs %in% c("rGC", "uGC")]

## This function produces boxplots of selected variables.
## varnames: varnames output by the episims
## targets: optional set of calibration targets
## data: longform dataset of simulated outputs; the value variable must be named
##       mn_lastyr
pbox <- function(varnames, targets = NULL, targets_ll = NULL, targets_ul = NULL,
                 data = episel, plot_title = "") {

  cdat <- data[variable %in% varnames]
  mdat <- cdat[, .(mean = mean(mn_lastyr)), by = variable]

  basep <-
    ggplot(data = cdat, aes(variable, mn_lastyr)) +
    geom_boxplot(outlier.size = 0, outlier.stroke = 0) +
    geom_line(aes(group = simid), alpha = 0.1) +
    geom_point(
      alpha = 0.3,
      size = 3,
      position = position_jitter(width = 0.2, height = 0)
    ) +
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
    } else {
      tdat[, ":=" (
        targval_ll = NA_real_,
        targval_ul = NA_real_
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
    ggtitle(plot_title) +
    scale_fill_scico_d(name = "Value", begin = 0.7) +
    theme_tufte(base_family = "sans", base_size = 22)
}

################################################################################
## CALCULATE RELATIVE TARGETS ##
################################################################################

calc_relative_target <- function(targ_pattern, targ_vec,
                                 ref_group, newtarg_slug = "") {
  tcalc <-
    targets[target %like% targ_pattern,
            value[subgroups %in% targ_vec] / value[subgroups == ref_group]]

  data.table(
    target = paste0(targ_pattern, "_", newtarg_slug, "_relative"),
    value = tcalc,
    subgroups = targ_vec
  )
}

trel_vls_byrace <- calc_relative_target("ct_vls", rslugs[1:3], "W", "byrace")
trel_vls_byage <-
  calc_relative_target("ct_vls", paste0("age", 1:4), "age5", "byage")

trel_hivdx_byrace <- calc_relative_target("ct_hivdx", rslugs[1:3], "W", "byrace")
trel_hivdx_byage <-
  calc_relative_target("ct_hivdx", paste0("age", 1:4), "age5", "byage")

trel_prep_byrace <- calc_relative_target("ct_prep", rslugs[1:3], "W", "byrace")

targets <- rbindlist(
  list(
    targets,
    trel_vls_byrace,
    trel_vls_byage,
    trel_hivdx_byrace,
    trel_hivdx_byage
  ),
  fill = TRUE
)

################################################################################
## CALCULATE DERIVED MEASURES in SIMULATION OUTPUT ##
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

# ... PrEP coverage
race.prep.cov <- paste0("prepCov", rdxlabs)

epi[, (race.prep.cov) := lapply(1:5, function(x) {
  get(paste0("prepCurr", rdxlabs[x])) / get(paste0("prepElig", rdxlabs[x]))
})][]

# ... relative HIV diagnosis, by race/ethicity
race.rel.dx <- paste0(race.dx[2:4], ".rel.ref.W")

epi[, (race.rel.dx) := lapply(race.dx[2:4], function(x) {
  get(paste0(x)) / i.prev.dx.inf.W
})][]

# ... relative HIV diagnosis, by age group
age.rel.dx <- paste0(age.dx[1:4], ".rel.ref.age5")

epi[, (age.rel.dx) := lapply(age.dx[1:4], function(x) {
  get(paste(x)) / i.prev.dx.inf.age5
})][]

# ... relative viral load suppression, by race/ethnicity
race.vls <- paste0("cc.vsupp", rdxlabs[2:4])
race.rel.vls <- paste0(race.vls, ".rel.ref.W")

epi[, (race.rel.vls) := lapply(race.vls, function(x) {
  get(x) / cc.vsupp.W
})][]

# ... relative viral load suppression, by age group
age.vls <- paste0("cc.vsupp.age", 1:4)
age.rel.vls <- paste0(age.vls[1:4], ".rel.ref.age5")

epi[, (age.rel.vls) := lapply(age.vls, function(x) {
  get(x) / cc.vsupp.age5
})][]

# ... relative PrEP coverage, by race/ethnicity
race.rel.prep <- paste0(race.prep.cov[2:4], ".rel.ref.W")

epi[, (race.rel.prep) := lapply(race.prep.cov[2:4], function(x) {
  get(x) / prepCov.W
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
  get(paste0("incid", rdxlabs[x])) / get(paste0("num", rdxlabs[x])) * 5200
})][]

# ... demographics
epi[, paste0("prop", rdxlabs[-1]) := lapply(1:4, function(x) {
  get(paste0("num", rdxlabs[x + 1])) / num
})][]

epi[, paste0("prop.age.", 1:5) := lapply(1:5, function(x) {
  get(paste0("num.age.", x)) / num
})][]

## ... population-level rectal and urethral GC
epi[, paste0("prev.", c("rGC", "uGC")) := lapply(c("rgc", "ugc"), function(x) {
          get(paste0("i.num.", x)) / num
        })][]

################################################################################
## CLEAN UP OUTPUT AND CALCULATE LAST-YEAR AVERAGES ##
################################################################################

sumcols <- names(epi)[!names(epi) %like% "at|simid"]

## Mark NAs as 0
epi_noNA <-
  epi[, (sumcols) := lapply(.SD, function(.x) fifelse(is.na(.x), 0, .x)),
      .SDcols = sumcols]

names(epi_noNA)

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
epi_mn_selnum[, .(N = .N, P = .N / length(unique(epi$simid)))]


################################################################################
## FUNCTION TO CALCULATE DEVIATION FROM TARGETS ##
################################################################################

## This function calculates the discrepancy between simulated output and
## a target statistic of choice. The following metrics are calculated
##   - Difference between target and simulated output
##   - Absolute difference between target and simulated output
##   - Absolute standardized difference between target and simulated output
##   - Percent difference between target and simulated output
##   - Output within following thresholds (5%, and 5/10 percentage points)
get_absdiff <- function(targname, output,
                        target_dt = targets, data = epi_mn_selnum) {

  cols <- c("simid", output)
  d <- data[, ..cols]

  d[, diff := get(output) - target_dt[target == targname, value]]
  d[, abs_diff := abs(diff)]
  d[, abs_std_diff := abs_diff / sd(get(output))]

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

## Run get_absdiff() on multiple targets.
##   target_str = value of target variable in targets object
##   subgroups  = slugs for race/ethnicity, age group, etc.
##   output_pattern = pattern to retrieve variable from simulation output,
##                    passed to sprintf() internally
get_absdiffv <- function(target_str, subgroups, output_pattern) {
  lapply(subgroups, function(.x) {
    get_absdiff(
      target_str,
      sprintf(output_pattern, .x),
      target_dt = targets[target == target_str & subgroups == .x])
  })
}


################################################################################
## FUNCTION TO PLOT TIME SERIES ##
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


################################################################################
## FUNCTION TO QUICKLY VISUALIZE TARGET DISTRIBUTIONS AMONG SIMID SETS ##
################################################################################

quicktarget <- function(outcome_pattern, simid_list,
                        custom_jitter = list(w = 0.1, h = 0)) {
  tmp <- rbindlist(
    lapply(
      names(simid_list),
      function(.x) {
        out_vs_targ[
          simid %in% simid_list[[.x]] &
            variable == sprintf(outcome_pattern, .x)
        ][, selection_group := .x][]
      }
    )
  )

  tmp %>%
    ggplot(aes(x = variable, y = output)) +
    geom_line(aes(group = simid), color = "gray70") +
    geom_point(
      aes(fill = selection_group),
      color = "white", shape = 21, size = 3,
      position = position_jitter(
        width = custom_jitter$w,
        height = custom_jitter$h
      )
    ) +
    geom_boxplot(alpha = 0.5, outlier.size = 0) +
    geom_point(
      aes(y = target_val),
      size = 7, shape = 21, color = "red"
    ) +
    scale_fill_scico_d() +
    theme_minimal(base_size = 25) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0))
}


################################################################################
## FUNCTION TO RETRIEVE THE INPUT PARAMETERS SETS ASSOCIATED ##
## WITH PARTICULAR LISTS OF SIMIDS                           ##
################################################################################

pull_params <- function(simid_list, pattern) {
  out <- rbindlist(
    lapply(
      simid_list,
      function(.x) {
        rbindlist(
          lapply(
            setNames(as.numeric(.x), .x),
            function(.y) as.data.table(lhs_groups[[.y]], keep.rownames = TRUE)
          ),
          idcol = "simid"
        )
      }
    ),
    idcol = "selection_group"
  )
  setnames(out, c("V1", "V2"), c("input", "value"))
  out
}


################################################################################
## GET QUANTILES FOR PRIORS ##
################################################################################

## Function that takes the long version of input parameter subsets from
## simid selections and extracts the 25% and 75% quantiles for use as prior
## limits in the next round of calibration.
## Default pattern extracts the race/ethnicity slug
input_quantiles <- function(data, inputstring,
                            extract_pat = "(?<=_)[A-Z]+$",
                            ql = 0.25, qu = 0.75) {
  inputsub <- data[
    input %like% inputstring,
    .(
      q25 = quantile(value, ql),
      q75 = quantile(value, qu)
    ),
    .(selection_group, input)
    ## this step matches the selection group to the corresponding
    ## input input parameter
  ][selection_group == substring(str_extract(input, extract_pat), 1, 1)]

  inputsub[, -c("selection_group")]
}


################################################################################
## PLOT TARGETS AMONG SELECTION OF SIMIDS ##
################################################################################
plot_targets <- function(simids, filter_out = " ") {
  out_vs_targ[simid %in% simids][!(variable %like% filter_out)] %>%
    ggplot(aes(x = variable, y = output)) +
    geom_point(position = position_jitter(width = 0.1)) +
    geom_point(
      aes(y = target_val, color = "target"),
      shape = 21,
      size = 4,
      stroke = 2,
      fill = NA
    ) +
    ggtitle(paste("Number of simids selected:", length(unique(simids)))) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major.x = element_line(color = "gray80"),
      plot.title = element_text(face = "bold")
    )
}
