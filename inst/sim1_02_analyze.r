################################################################################
## SETUP ##
################################################################################

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
theme_set(theme_tufte(base_size = 20))

## Get date from system environment if set. Otherwise, find date of most recent
## sim1_epi run.
sys_batch_date <- Sys.getenv("BATCH_DATE")

sim_epis <- list.files(
    here::here("burnin", "abc", "sim1"),
    pattern = "sim1_epi",
    full.names = TRUE
)

if (sys_batch_date != "") {
  batch_date <- Sys.getenv("BATCH_DATE")
} else {
  batch_date <- max(ymd(str_extract(sim_epis, "[0-9]{4}-[0-9]{2}-[0-9]{2}")))
}

## Load calibration targets, netstats, and formatted simulated output.
epi <- readRDS(sim_epis[sim_epis %like% batch_date])
targets <- readRDS(here::here("est", "caltargets.Rds"))
netstats <- get_est("netstats")

## Helper function to save plots.
psave <- function(f, p, w = 22, h = 22) {
  caldir <- here::here("inst", "calfigs", batch_date)
  if (!dir.exists(caldir)) dir.create(caldir)
  ggsave(
    here::here(caldir, paste0(f, ".pdf")),
    p, width = w, height = h, units = "in"
  )
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


################################################################################
## SELECT INITIAL SIMS BASED ON POP. SIZE ##
################################################################################

pop_N <- 20000
t_eval <- 3068

sumcols <- names(epi)[!names(epi) %like% "at|simid"]

## Mark NAs as 0
epi_noNA <- epi[, lapply(.SD, function(.x) fifelse(is.na(.x), 0, .x)),
                  .SDcols = c(sumcols, "at", "simid")]

## Get output averages over the last 52 weeks
epi_mn <- epi_noNA[
  at >= t_eval, lapply(.SD, mean),
  by = simid, .SDcols = sumcols]

# Pick parameter sets that keep us within 5% of population size N = 20,000
# (average over final burnin-in year)
epi_mn_selnum <- epi_mn[abs(pop_N - num) / pop_N <= 0.05]
epi_mn_selnum[, .N]

get_absdiff <- function(target, output, data = epi_mn_selnum) {
  cols <- c("simid", output)
  d <- data[, ..cols]
  d[, absdiff := abs(get(output) - target)]
  dout <- d[order(absdiff)]
  dout[, rank := 1:.N]
  return(dout)
}

t_hiv_inc <- get_absdiff(targets$ct_hiv_incid_100k / 1000, "ir100")
t_hiv_prv <- get_absdiff(targets$ct_hiv_prev, "i.prev")

# diagnosis prevalence among infected
t_hiv_dx_race <- lapply(
  names(targets$ct_hivdx_pr_byrace),
  function(.x) {
    get_absdiff(targets$ct_hivdx_pr_byrace[[.x]], paste0("i.prev.dx.inf.", .x))
  }
)

# viral suppression among HIV-diagnosed, by age group
t_vls_age <- lapply(
  names(targets$ct_vls_pr_byage),
  function(.x) {
    get_absdiff(targets$ct_vls_pr_byage[[.x]], paste0("cc.vsupp.", .x))
  }
)

# gonorrhea testing in clinic, by anatomic site
t_gc_anat_test <- lapply(
  names(targets$ct_prop_anatsite_tested),
  function(.x) {
    var <- sprintf("prop.%s.tested", .x)
    pass1 <- get_absdiff(
      targets$ct_prop_anatsite_tested[[.x]],
      var
    )
    # choose only runs where men received tests at the anat site
    pass2 <- pass1[get(var) > 0][, rank := 1:.N][]
    pass2
  }
)

# gonorrhea positivity in clinic, among tested anatomic sites
t_gc_pos <- lapply(
  names(targets$ct_prop_anatsite_pos),
  function(.x) {
    lugc <- c("uGC" = "ureth", "rGC" = "rect", "pGC" = "phar")
    slug <- names(match.arg(.x, lugc))
    var <- sprintf("prob.%s.tested", slug)
    get_absdiff(
      targets$ct_prop_anatsite_pos[[.x]],
      var
    )[get(var) > 0][, rank := 1:.N]
  }
)

t_hiv_inc
t_hiv_prv
t_hiv_dx_race
t_vls_age
t_gc_anat_test
t_gc_pos

tlist <- list(
  t_hiv_inc,
  t_hiv_prv,
  t_hiv_dx_race[[1]],
  t_hiv_dx_race[[2]],
  t_hiv_dx_race[[3]],
  t_hiv_dx_race[[4]],
  t_vls_age[[1]],
  t_vls_age[[2]],
  t_vls_age[[3]],
  t_vls_age[[4]],
  t_vls_age[[5]],
  t_gc_anat_test[[1]],
  t_gc_anat_test[[2]],
  t_gc_anat_test[[3]],
  t_gc_pos[[1]],
  t_gc_pos[[2]],
  t_gc_pos[[3]]
)

names(tlist) <- sapply(tlist, function(.x) names(.x)[2])

simid_sel <- lapply(tlist, function(.x) .x[rank <= .N / 2, simid])

intsel <- Reduce(intersect, simid_sel)

lapply(tlist, function(.x) .x[simid %in% intsel])

episel <- epi[simid %in% intsel]
episel


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

    p <- ggplot(data, aes(x = at, y = value, color = get(varname))) +
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

plotepi("num", 20000, data = episel)

plotepi(
  c(paste0("num.", rlabs)),
  varname = "N",
  line_alpha = 0.4,
  data = episel
)

plotepi(
  "i.prev",
  targets$ct_hiv_prev,
  data = episel
)

plotepi(
  c("i.prev", paste0("i.prev.", rlabs)),
  varname = "HIV prevalence",
  line_alpha = 0.4,
  data = episel
)

plotepi("ir100", targets$ct_hiv_incid_100k / 1000, data = episel)


melt(
  epi_mn_selnum[, lapply(.SD, mean), .SDcols = names(tlist)]
)

melt(
  epi_mn_selnum[simid %in% intsel, lapply(.SD, mean), .SDcols = names(tlist)]
)

tnames <- names(tlist)
epi_mn_intsel <- epi_mn_selnum[simid %in% intsel, ..tnames]

ggcorrplot::ggcorrplot(cor(epi_mn_intsel))
