################################################################################
## SETUP ##
################################################################################

library(pacman)
p_load(
  data.table,
  ggplot2,
  ggthemes,
  magrittr,
  scico,
  stringr
)

## path to current directory
an_path <- here::here("inst/analysis01_epi/")

## make sure length of main and sensitivity simulations are the same
timewarn <- function(maxstep_sens, maxstep_main = max(smain65$at)) {
  ## maxstep_X INT  last time step in a simulation file
  if (maxstep_sens != maxstep_main) {
    warning("Main and sensitivity timesteps don't match!")
  } else if (maxstep_sens == maxstep_main) {
    message("Timesteps match. Do a barrel roll!")
  }
}

## labels
rlabs    <- c("B", "H", "O", "W")
agelabs  <- paste0("age", 1:5)
anatlabs <- c("rect", "ureth", "phar")
gclabsl  <- c("rgc", "ugc", "pgc")
gclabsh  <- c("rGC", "uGC", "pGC")


################################################################################
## DATA ##
################################################################################

sims <- list.files(an_path, pattern = "(MAIN|SENS).*rds", full.names = TRUE)

## All 65 years
smain65 <- readRDS(sims[sims %like% "MAIN"])
spois65 <- readRDS(sims[sims %like% "poisdurat"])

names(smain65) <- str_replace(names(smain65), "age\\.(?=[1-5])", "age")
names(spois65) <- str_replace(names(spois65), "age\\.(?=[1-5])", "age")

## Last 5 years
l5start <- max(smain65$at) - (52 * 5)
smain <- smain65[at >= l5start]

timewarn(max(spois65$at))
spois <- spois65[at >= l5start]

## Create some derivative variables
gcpath_vars <- names(smain)[names(smain) %like% "(r|u|p)2"]
gcpath_props <- paste0("prop.", gcpath_vars)

smain[, (gcpath_props) := lapply(.SD, function(.x) .x / incid.gc),
      .SDcols = gcpath_vars]

spois[, (gcpath_props) := lapply(.SD, function(.x) .x / incid.gc),
      .SDcols = gcpath_vars]

## Mean over last 5 years
sumcols <- names(smain[, -c("simid", "at")])

smains <- smain[, lapply(.SD, mean), .SDcols = sumcols, by = simid]
spoiss <- spois[, lapply(.SD, mean), .SDcols = sumcols, by = simid]


################################################################################
## VISUALIZATION FUNCTIONS ##
################################################################################

head(smain)
head(spois)

## Function: Plot time series of simulations.
##  - data must include the 'at' variable (i.e., time step)
plotts <- function(var, data) {

  ## encase var in a quosure
  var <- enquo(var)

  ## plot
  data %>%
    ggplot(aes(x = at, y = !!var, group = simid)) +
    geom_line(color = "gray80") +
    labs(x = var) +
    theme_few()
}

## Function: Plot distribution of last-5-year averages
##  - data must be summarized by simid (simulation ID)
plotavg <- function(var, data) {

  if (length(var) == 1) {
    ## encase var in quosure
    var <- enquo(var)

    ## jitter seed
    set.seed(31415)

    ## plot
    data %>%
      ggplot(aes(x = factor(1), y = !!var)) +
      geom_point(alpha = 0.3, position = position_jitter(width = 0.03)) +
      geom_boxplot(outlier.size = 0.5, width = 0.1, alpha = 0.5) +
      labs(x = var) +
      theme_few()
  }

  if (length(var) > 1) {
    tmp <- melt(data[, ..var], measure.vars = var)
    tmp %>%
      ggplot(aes(x = variable, y = value)) +
      geom_point(alpha = 0.1, position = position_jitter(width = 0.03)) +
      geom_boxplot(
        outlier.size = 0,
        outlier.alpha = 0,
        width = 0.1,
        alpha = 0.5
      ) +
      theme_few()
  }
}



################################################################################
## PLOTS ##
################################################################################

## variable names
## tags:
##  i = incidence
##  p = prevalence
##  v = variable
##  r = race/ethnicity
##  a = age group
##  s = anatomic site
##  g = anatomic site-specific gonorrhea
##  p = pathway
hivp_vr <- paste0("i.prev.", rlabs)
hivi_vr <- paste0("ir100.", rlabs)
hivp_ar <- paste0("i.prev.", agelabs)
hivi_ar <- paste0("ir100.", agelabs)

gcp_vg <- paste0("prev.", c("gc", gclabsl))
gci_vg <- paste0("ir100.", c("gc", gclabsl))
gci_vp <- names(smain)[names(smain) %like% "(r|u|p)2"]

plotts(i.prev, smain)
plotts(num.B, smain)
plotts(incid.rgc, smain)
plotts(incid.r2pgc, smain)

plotavg(hivp_vr, smains)
plotavg(hivi_vr, smains)

plotavg(hivp_ar, smains)

plotavg(gcp_vg, smains)
plotavg(gci_vg, smains)
plotavg(gci_vp, smains)
plotavg(gcpath_props, smains)


################################################################################
## TABLES ##
################################################################################

simb <- rbindlist(
  list(
    Main = smains,
    Poisson = spoiss
  ),
  idcol = "analysis"
)

simbm <- melt(simb, id.vars = c("simid", "analysis"))

droppat <-
  paste("cc\\.linked1m", "cc\\.dx\\.delay", "mean\\.tx\\.", sep = "|")

simbm <- simbm[!(variable %like% droppat)]

sresult <- simbm[, .(
  median = median(value),
  mean = mean(value),
  sd = sd(value),
  si025 = quantile(value, 0.025),
  si975 = quantile(value, 0.975),
  q25 = quantile(value, 0.25),
  q75 = quantile(value, 0.75)
), keyby = .(analysis, variable)]

sresult_gcpop <-
  sresult[variable %in% c(gcp_vg, gci_vg, gci_vp, "ir100.anatsite.yrs.gc")]

sumquants <- c("median", "mean", "sd", "si025", "si975", "q25", "q75")

sresult_gcpop[,
  (sumquants) := lapply(.SD, round, digits = 2),
  .SDcols = sumquants
  ]

## No by processing is actually done here. The 'by' statement simply preserves
## the analysis and variable columns.
sresult_gcpop[, .(
  median = median,
  iqr = paste0(" [", q25, ", ", q75, "]"),
  si95 = paste0("(", si025, ", ", si975, ")")
), by = .(analysis, variable)]


getnames <- function(pattern, dt) {
  names(dt)[names(dt) %like% pattern]
}


################################################################################
## GC INCIDENCE BY RACE/ETHNICITY ##
################################################################################

## Bind the weekly datasets so we can calculate incidence rates by
## rac/ethnicity and age group.
simb_at <- rbindlist(
  list(
    Main = smain,
    Poisson = spois
  ),
  idcol = "analysis"
)

## Variable lists
idvars <- c("analysis", "simid", "at")

rethvars <- c(
  idvars,
  getnames("(incid.[B,H,O,W].[r,u,p]gc)|(^num.[B,H,O,W])$", smain)
)

simb_at_reth <- simb_at[, ..rethvars]

reth_incid_names <- getnames("incid.[B,H,O,W]", simb_at_reth)
reth_incnames <- str_replace(reth_incid_names, "incid", "ir100.pop")

calc_gc_ir100 <- function(incid, num) incid / num * 5200

simb_at_reth[, ":=" (
  ir100.pop.B.rgc = calc_gc_ir100(incid.B.rgc, num.B),
  ir100.pop.H.rgc = calc_gc_ir100(incid.H.rgc, num.H),
  ir100.pop.O.rgc = calc_gc_ir100(incid.O.rgc, num.O),
  ir100.pop.W.rgc = calc_gc_ir100(incid.W.rgc, num.W),
  ir100.pop.B.ugc = calc_gc_ir100(incid.B.ugc, num.B),
  ir100.pop.H.ugc = calc_gc_ir100(incid.H.ugc, num.H),
  ir100.pop.O.ugc = calc_gc_ir100(incid.O.ugc, num.O),
  ir100.pop.W.ugc = calc_gc_ir100(incid.W.ugc, num.W),
  ir100.pop.B.pgc = calc_gc_ir100(incid.B.pgc, num.B),
  ir100.pop.H.pgc = calc_gc_ir100(incid.H.pgc, num.H),
  ir100.pop.O.pgc = calc_gc_ir100(incid.O.pgc, num.O),
  ir100.pop.W.pgc = calc_gc_ir100(incid.W.pgc, num.W)
)][]

reths <- simb_at_reth[,
  lapply(.SD, mean),
  .SDcols = reth_incnames,
  by = .(analysis, simid)
]

rethsm <- melt(reths, id.vars = c("analysis", "simid"))

rethsum <- rethsm[, .(
    median = median(value),
    mean = mean(value),
    sd = sd(value),
    si025 = quantile(value, 0.025),
    si975 = quantile(value, 0.975),
    q25 = quantile(value, 0.25),
    q75 = quantile(value, 0.75)
  ), keyby = .(analysis, variable)]

rethsm[, ":=" (
  anatsite = str_extract(variable, "[r,u,p]gc$"),
  race = str_extract(variable, "[B,H,O,W]")
)][]

rethsm %>%
  ggplot(aes(x = race, y = value)) +
  geom_point(
    alpha = 0.1,
    position = position_jitter(width = 0.2)
  ) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.6,
    outlier.alpha = 0
  ) +
  facet_grid(analysis ~ anatsite) +
  theme_few()


################################################################################
## GC INCIDENCE BY AGE GROUP ##
################################################################################

simb_at_age <- simb_at[, ..agevars]

age_incid_names <- getnames("incid.age", simb_at_age)
age_incnames <- str_replace(age_incid_names, "incid", "ir100.pop")

simb_at_age[, ":=" (
  ir100.pop.age1.rgc = calc_gc_ir100(incid.age1.rgc, num.age1),
  ir100.pop.age2.rgc = calc_gc_ir100(incid.age2.rgc, num.age2),
  ir100.pop.age3.rgc = calc_gc_ir100(incid.age3.rgc, num.age3),
  ir100.pop.age4.rgc = calc_gc_ir100(incid.age4.rgc, num.age4),
  ir100.pop.age5.rgc = calc_gc_ir100(incid.age5.rgc, num.age5),
  ir100.pop.age1.ugc = calc_gc_ir100(incid.age1.ugc, num.age1),
  ir100.pop.age2.ugc = calc_gc_ir100(incid.age2.ugc, num.age2),
  ir100.pop.age3.ugc = calc_gc_ir100(incid.age3.ugc, num.age3),
  ir100.pop.age4.ugc = calc_gc_ir100(incid.age4.ugc, num.age4),
  ir100.pop.age5.ugc = calc_gc_ir100(incid.age5.ugc, num.age5),
  ir100.pop.age1.pgc = calc_gc_ir100(incid.age1.pgc, num.age1),
  ir100.pop.age2.pgc = calc_gc_ir100(incid.age2.pgc, num.age2),
  ir100.pop.age3.pgc = calc_gc_ir100(incid.age3.pgc, num.age3),
  ir100.pop.age4.pgc = calc_gc_ir100(incid.age4.pgc, num.age4),
  ir100.pop.age5.pgc = calc_gc_ir100(incid.age5.pgc, num.age5)
)][]

ages <- simb_at_age[,
  lapply(.SD, mean),
  .SDcols = age_incnames,
  by = .(analysis, simid)
]

agesm <- melt(ages, id.vars = c("analysis", "simid"))

agesum <- agesm[, .(
    median = median(value),
    mean = mean(value),
    sd = sd(value),
    si025 = quantile(value, 0.025),
    si975 = quantile(value, 0.975),
    q25 = quantile(value, 0.25),
    q75 = quantile(value, 0.75)
  ), keyby = .(analysis, variable)]

agesm[, ":=" (
  anatsite = str_extract(variable, "[r,u,p]gc$"),
  agegrp = str_extract(variable, "age[1-5]")
)][]

agesm %>%
  ggplot(aes(x = agegrp, y = value)) +
  geom_point(
    alpha = 0.1,
    position = position_jitter(width = 0.2)
  ) +
  geom_boxplot(
    width = 0.6,
    alpha = 0.6,
    outlier.alpha = 0
  ) +
  facet_grid(analysis ~ anatsite) +
  theme_few()
