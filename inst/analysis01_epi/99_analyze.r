################################################################################
## SETUP ##
################################################################################

library(pacman)

p_load(
  data.table,
  ggplot2,
  ggthemes,
  ggridges,
  magrittr,
  scico,
  stringr
)

simpath <- here::here("inst", "analysis01_epi")
epi <- readRDS(file.path(simpath, "epi_5yr_BOUND.rds"))

## clean up age column names
epnm <- names(epi)
names(epi) <- str_replace(epnm, "age\\.", "age")

## labels
rlabs <- c("B", "H", "O", "W")
rlabslong <- c("Black", "Hispanic", "Other", "White")
agelabs <- paste0("age", 1:5)
agelabs_pretty <- c("18-24", "25-34", "35-44", "45-54", "55+")
anatlabs <- c("rect", "ureth", "phar")
gclabsl <- c("rgc", "ugc", "pgc")
gclabsh <- c("rGC", "uGC", "pGC")

## sensitivity analysis parameter lookup
sens07_specs <- readRDS(file.path(simpath, "specs_SENS07.rds"))


################################################################################
## VISUALIZATION FUNCTIONS ##
################################################################################

## Function: Plot distribution of last-5-year averages
##  - data must be summarized by simid (simulation ID)
plotepi <- function(var, data) {

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
gci_vp <- names(epi)[names(epi) %like% "(r|u|p)2"]

## HIV
plotepi(hivp_vr, epi[analysis == "MAIN"])
plotepi(hivi_vr, epi[analysis == "MAIN"])
plotepi(hivp_ar, epi[analysis == "MAIN"])

## Gonorrhea
plotepi(gcp_vg, epi[analysis == "MAIN"])
plotepi(gci_vg, epi[analysis == "MAIN"])
plotepi(gci_vp, epi[analysis == "MAIN"])


################################################################################
## TABLES ##
################################################################################

droppat <- paste(
  "cc\\.linked1m", "cc\\.dx\\.delay", "mean\\.tx\\.",
  "min_at", "max_at",
  sep = "|"
)

## NOTE Suppress warnings about coercing integers to doubles.
epim <- suppressWarnings(
  melt(
    epi,
    id.vars = c("simid", "analysis"),
    variable.factor = FALSE
  )[!(variable %like% droppat)]
)

## create grouping variables
epim[, ":=" (
  anatsite = str_extract(variable, "[r,u,p]gc$"),
  race = str_extract(variable, "[B,H,O,W]"),
  agegrp = str_extract(variable, "age(\\.[1-5]|[1-5])")
)]

epim[, agegrp := str_replace(agegrp, "\\.", "")]
epim[, .N, agegrp]

## replace NA/NaN with 0
epim[is.na(value), value := 0]

sresult <- epim[, .(
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
  (sumquants) := lapply(.SD, round, digits = 3),
  .SDcols = sumquants
  ]


################################################################################
## GC INCIDENCE by RACE/ETHNICITY ##
################################################################################

rethsm <- epim[!is.na(race) & !is.na(anatsite) &
               variable %like% "ir100" & analysis %like% "MAIN|pois"]

rethsm[, race := fcase(race == "B", "Black",
                       race == "H", "Hispanic",
                       race == "O", "Other",
                       race == "W", "White")]

rethsm[, race := factor(race, levels = rlabslong)]
rethsm[, race := factor(race, levels = rev(levels(factor(race))))]
str(rethsm)

## plot labelers
## NOTE analysis labeler is reused for age plots below
anat_lblr <- c(rgc = "Rectal", ugc = "Urethral", pgc = "Pharyngeal")
analysis_lblr <- c(
  MAIN = "Main",
  MAIN_altrim = "Alternative rimming",
  SENS_poisdurat = "Poisson"
)

rethsm %>%
  ggplot(aes(x = value, y = race, fill = stat(x))) +
  geom_density_ridges_gradient(
    color = "white",
    scale = 0.7, alpha = 1,
    jittered_points = TRUE,
    point_size = 0.5, point_alpha = 1, point_color = "#F9E480",
    vline_size = 0.5, vline_color = "white",
    quantile_lines = TRUE,
    position = position_raincloud(height = 0.1)
  ) +
  geom_boxplot(
    aes(x = value, y = race, fill = NULL),
    alpha = 0,
    width = 0.2,
    outlier.size = 0,
    position = position_nudge(y = -0.1)
  ) +
  scale_fill_scico("Incidence rate", palette = "lajolla") +
  facet_grid(
    analysis ~ anatsite,
    labeller = labeller(anatsite = anat_lblr, analysis = analysis_lblr)
  ) +
  labs(x = "Incidence rate", y = "Race/ethnicity") +
#  coord_cartesian(clip = "off") +
  theme_few(base_size = 18) +
  theme(
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "top",
    legend.title = element_text(vjust = 1),
    legend.key.width = unit(2, "cm")
  )


################################################################################
## GC INCIDENCE by AGE GROUP ##
################################################################################

agesm <- epim[!is.na(agegrp) & !is.na(anatsite) &
              variable %like% "ir100" & analysis %like% "MAIN|pois"]

agesm[, agegrp := agelabs_pretty[match(agegrp, agelabs)]]
agesm[, agegrp := factor(agegrp, levels = agelabs_pretty)]
agesm[, agegrp := factor(agegrp, levels = rev(levels(factor(agegrp))))]
str(agesm)

agesm %>%
  ggplot(aes(x = value, y = agegrp, fill = stat(x))) +
  geom_density_ridges_gradient(
    color = "white",
    scale = 0.7, alpha = 1,
    jittered_points = TRUE,
    point_size = 0.5, point_alpha = 1, point_color = "#F9E480",
    vline_size = 0.5, vline_color = "white",
    quantile_lines = TRUE,
    position = position_raincloud(height = 0.1)
  ) +
  geom_boxplot(
    aes(x = value, y = agegrp, fill = NULL),
    alpha = 0,
    width = 0.1,
    outlier.size = 0,
    position = position_nudge(y = -0.1)
  ) +
  scale_fill_scico("Incidence rate", palette = "lajolla") +
  facet_grid(
    analysis ~ anatsite,
    labeller = labeller(anatsite = anat_lblr, analysis = analysis_lblr)
  ) +
  labs(x = "Incidence rate", y = "Age group") +
  theme_few(base_size = 18) +
  theme(
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "top",
    legend.title = element_text(vjust = 1),
    legend.key.width = unit(2, "cm")
  )


################################################################################
## WRITE FILES ##
################################################################################

an01_dir <- here::here("inst/analysis01_epi")

saveRDS(
  sresult_gcpop,
  file.path(an01_dir, "result_gc_summary.rds")
)

saveRDS(
  rethsm,
  file.path(an01_dir, "result_gcinc_byreth_bypois.rds")
)

saveRDS(
  agesm,
  file.path(an01_dir, "result_gcinc_byage_bypois.rds")
)
