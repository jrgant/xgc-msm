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

simpath <- here::here("inst", "analysis02_screening")
epi <- readRDS(file.path(simpath, "epi_5yr_BOUND.rds"))

# NOTE kissing exposure tracking not implemented
epi <- epi[!(analysis %like% "CDC2")]

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
sens03_specs <- readRDS(file.path(simpath, "specs_SENS03.rds"))


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
##  i, incidence
##  p, prevalence
##  v, variable
##  r, race/ethnicity
##  a, age group
##  s, anatomic site
##  g, anatomic site-specific gonorrhea
##  p, pathway
hivp_vr <- paste0("i.prev.", rlabs)
hivi_vr <- paste0("ir100.", rlabs)
hivp_ar <- paste0("i.prev.", agelabs)
hivi_ar <- paste0("ir100.", agelabs)

gcp_vg <- paste0("prev.", c("gc", gclabsl))
gci_vg <- paste0("ir100.", c("gc", gclabsl))
gci_vp <- names(epi)[names(epi) %like% "(r|u|p)2"]

## HIV
plotepi(hivp_vr, epi[analysis == "MAIN_STI_BASE"])
plotepi(hivi_vr, epi[analysis == "MAIN_STI_BASE"])
plotepi(hivp_ar, epi[analysis == "MAIN_STI_BASE"])

## Gonorrhea
plotepi(gcp_vg, epi[analysis == "MAIN_STI_BASE"])
plotepi(gci_vg, epi[analysis == "MAIN_STI_BASE"])
plotepi(gci_vp, epi[analysis == "MAIN_STI_BASE"])

plotepi(gcp_vg, epi[analysis == "MAIN_STI_CDC1"])
plotepi(gci_vg, epi[analysis == "MAIN_STI_CDC1"])
plotepi(gci_vp, epi[analysis == "MAIN_STI_CDC1"])

plotepi(gcp_vg, epi[analysis == "MAIN_STI_CDC2"])
plotepi(gci_vg, epi[analysis == "MAIN_STI_CDC2"])
plotepi(gci_vp, epi[analysis == "MAIN_STI_CDC2"])

plotepi(gcp_vg, epi[analysis == "MAIN_STI_SYMP"])
plotepi(gci_vg, epi[analysis == "MAIN_STI_SYMP"])
plotepi(gci_vp, epi[analysis == "MAIN_STI_SYMP"])

plotepi(gcp_vg, epi[analysis == "MAIN_STI_UNIV"])
plotepi(gci_vg, epi[analysis == "MAIN_STI_UNIV"])
plotepi(gci_vp, epi[analysis == "MAIN_STI_UNIV"])


################################################################################
## MATCH SIMULATION RUNS ##
################################################################################

## TODO After thinking about the current approach, I don't think burning all
##      scenarios in makes a ton of sense. All the models reach similar
##      steady states, except for CDC1. Calculating NIA, PIA, etc. based on
##      models that were all burned in is not comparing what I originally
##      intended. Rather, for each set of parameters, 1,000 BASE simulations
##      should be burned in + 5 years. Then, for each alternative screening
##      protocol, the model should be rewound so that each screening strategy
##      is simulated as a switch to that strategy for 5 years.

## make row IDs
epi[, rowid := 1:.N]

## make analysis abbreviations
epi[analysis %like% "MAIN",
  aslug := str_flatten(
    str_extract_all(
      analysis,
      "^[A-Z]{1}(?=.*_.*_)|[A-Z0-9]{4}$",
      simplify = T
    )
  ),
  by = rowid
]

epi[analysis %like% "SENS",
  aslug := paste0(
    "S",
    str_remove_all(substring(analysis, 6, 10), "0|\\."),
    str_extract(analysis, "[A-Z0-9]{4}$")
  ),
  by = rowid
]

aslugtab <- epi[, .N, aslug]
aslugtab[, all(N == 1000)] # check analysis slugs

## calculate total number of tests administered during each simulation run
epi[, num.gc.test := num.rgc.test + num.ugc.test + num.pgc.test]

epi[, ":="(
  sum_num.gc.test = num.gc.test * length(min_at:max_at),
  sum_num.rgc.test = num.rgc.test * length(min_at:max_at),
  sum_num.ugc.test = num.ugc.test * length(min_at:max_at),
  sum_num.pgc.test = num.pgc.test * length(min_at:max_at)
), by = rowid]

## long-to-wide
relcols <- c(
  "simid", "aslug",
  "sum_incid.gc", "sum_incid.rgc", "sum_incid.ugc", "sum_incid.pgc"
)

getnames <- function(pattern, data) {
  names(data)[names(data) %like% pattern]
}

calc_relmeas <- function(sens_slug, alabs = c("gc", gclabsl), data = epi) {
  tmp_base <- data[aslug == "MBASE"]
  tmp_comp <- data[aslug == sens_slug]
  tmp_merge <- merge(tmp_base, tmp_comp, by = "simid")

  ## NIA, number of infections averted
  tmp_merge[, ":="(
    nia_gc = sum_incid.gc.x - sum_incid.gc.y,
    nia_rgc = sum_incid.rgc.x - sum_incid.rgc.y,
    nia_ugc = sum_incid.ugc.x - sum_incid.ugc.y,
    nia_pgc = sum_incid.pgc.x - sum_incid.pgc.y
  )]

  ## PIA, proportion of infections averted
  tmp_merge[, ":="(
    pia_gc = nia_gc / sum_incid.gc.x,
    pia_rgc = nia_rgc / sum_incid.rgc.x,
    pia_ugc = nia_ugc / sum_incid.ugc.x,
    pia_pgc = nia_ugc / sum_incid.pgc.x
  )]

  ## NTIA, number tested per infection averted
  tmp_merge[, ":="(
    ntia_gc = sum_num.gc.test.y / nia_gc,
    ntia_rgc = sum_num.rgc.test.y / nia_rgc,
    ntia_ugc = sum_num.ugc.test.y / nia_ugc,
    ntia_pgc = sum_num.pgc.test.y / nia_pgc
  )]

  tmp_merge[, analysis := sens_slug]

  savecols <- c("simid", "analysis", getnames("nia|pia|ntia", tmp_merge))
  tmp_merge[, ..savecols]
}

nobase_aslugs <- aslugtab[aslug %like% "^S" | !(aslug %like% "^M.*BASE"), aslug]

epil <- rbindlist(lapply(nobase_aslugs, calc_relmeas))


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

sresult_gcpop


################################################################################
## WRITE FILES ##
################################################################################

saveRDS(
  sresult_gcpop,
  file.path(simpath, "result_gc_summary.rds")
)

saveRDS(
  rethsm,
  file.path(simpath, "result_gcinc_byreth_bypois.rds")
)

saveRDS(
  agesm,
  file.path(simpath, "result_gcinc_byage_bypois.rds")
)
