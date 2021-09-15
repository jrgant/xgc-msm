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
