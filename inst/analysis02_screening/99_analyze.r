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

## 5-year averages
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

## time series
epitime <- rbindlist(
  lapply(
    file.path(simpath, list.files(simpath, pattern = "^epi_02.*rds")),
    function(.x) {
      sname <- str_extract(.x, "(?<=epi_).*(?=\\.rds)")
      tmp <- readRDS(.x)
      tmp[, .(
        simid, at, aid = sname,
        incid.gc, incid.rgc, incid.ugc, incid.pgc,
        prev.gc, prev.rgc, prev.ugc, prev.pgc,
        ir100.gc, ir100.rgc, ir100.ugc, ir100.pgc
      )]
    }
  )
)

format(object.size(epitime), units = "GB")

## sensitivity analysis parameter lookup
sens_specs <- readRDS(file.path(simpath, "specs_SENS.rds"))


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
## MATCH SIMULATION RUNS ##
################################################################################

## make row IDs
epi[, rowid := 1:.N]

## make analysis abbreviations
## aid = analysis ID
## groupid = analysis group ID (main vs. sensitivity batches)
## scenid = scenario ID
epi[,
  aid := str_flatten(
    str_extract_all(
      analysis,
      "([0-9]{2}\\.[0-9]{2})|((?<=[0-9]{1}_)[a-z]{1})|([a-z0-9]{4}$)",
      simplify = T
    )
  ),
  by = rowid
]

epi[, groupid := str_extract(aid, "[0-9]{2}\\.[0-9]{2}")]
epi[, scenid := str_extract(aid, "[a-z0-9]{4}$")]

aidtab <- epi[, .N, aid]
aidtab[, all(N == 1000)] # check analysis slugs

## calculate total number of tests administered during each simulation run
epi[, num.gc.test := num.rgc.test + num.ugc.test + num.pgc.test]

epi[, ":=" (
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

calc_relmeas <- function(group, alabs = c("gc", gclabsl), data = epi) {

  tmp_base <- data[groupid == group & scenid == "base"]
  tmp_comp <- data[groupid == group & scenid != "base"]
  tmp_merge <- merge(tmp_comp, tmp_base, by = "simid", all.x = T)

  ## NIA, number of infections averted
  tmp_merge[, ":=" (
    nia_gc = sum_incid.gc.y - sum_incid.gc.x,
    nia_rgc = sum_incid.rgc.y - sum_incid.rgc.x,
    nia_ugc = sum_incid.ugc.y - sum_incid.ugc.x,
    nia_pgc = sum_incid.pgc.y - sum_incid.pgc.x
  )]

  ## PIA, proportion of infections averted
  tmp_merge[, ":=" (
    pia_gc = nia_gc / sum_incid.gc.y,
    pia_rgc = nia_rgc / sum_incid.rgc.y,
    pia_ugc = nia_ugc / sum_incid.ugc.y,
    pia_pgc = nia_pgc / sum_incid.pgc.y
  )]

  ## NTIA, number tested per infection averted
  tmp_merge[, ":=" (
    ntia_gc = sum_num.gc.test.x / nia_gc,
    ntia_rgc = sum_num.rgc.test.x / nia_rgc,
    ntia_ugc = sum_num.ugc.test.x / nia_ugc,
    ntia_pgc = sum_num.pgc.test.x / nia_pgc
  )]

  tmp_merge[, ":=" (
    groupid = groupid.x,
    scenid = scenid.x,
    aid = aid.x
  )]

  #tmp_merge[, analysis := sens_slug]
  savecols <- c(
    "simid", "groupid", "aid", "scenid",
    getnames("nia|pia|ntia|.x|.y", tmp_merge)
  )

  tmp_merge[, ..savecols]
}

epirel <- rbindlist(lapply(unique(epi[, groupid]), calc_relmeas))


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
    epi[, -c("aid", "groupid", "scenid")],
    id.vars = c("simid", "analysis"),
    variable.factor = FALSE
  )[!(variable %like% droppat)]
)

nrow(epim)

## create grouping variables
epim[, ":=" (
  anatsite = str_extract(variable, "[r,u,p]gc$"),
  race = str_extract(variable, "[B,H,O,W]"),
  agegrp = str_extract(variable, "age(\\.[1-5]|[1-5])")
)]

epim[, agegrp := str_replace(agegrp, "\\.", "")]
epim[, .N, agegrp]

## replace epi measures marked NA/NaN with 0
epim[is.na(value), value := 0]

sresult <- epim[, .(
  median  = median(value),
  mean    = mean(value),
  sd      = sd(value),
  si025   = quantile(value, 0.025),
  si975   = quantile(value, 0.975),
  q25     = quantile(value, 0.25),
  q75     = quantile(value, 0.75)
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
## GC EPI DISTRIBUTIONS ##
################################################################################

epim[variable == "prev.gc"] %>%
  ggplot(aes(x = analysis, y = value)) +
  geom_point() +
  geom_boxplot() +
  theme_base(base_size = 40) +
  theme(axis.text.x = element_text(angle = 90))

epim[variable == "incid.gc",
     .(mean = mean(value),
       median = median(value),
       q25 = quantile(value, 0.25),
       q75 = quantile(value, 0.75)
       ),
     by = analysis]


################################################################################
## EPI TIME SERIES ##
################################################################################

timesum <- epitime[aid %like% "main", .(mn = mean(prev.gc),
                                        md = median(prev.gc)), .(aid, at)]

convert_cols <- names(epitime)[!(names(epitime) %in% c("aid", "simid", "at"))]
epitime[, (convert_cols) := lapply(.SD, as.numeric), .SDcols = convert_cols]

epitime_mainl <- melt(
  epitime[aid %like% "main"],
  id.vars = c("aid", "simid", "at")
)

epitime_mainl[, ":=" (
  measure = str_extract(variable, ".*(?=\\.)"),
  anatsite = str_extract(variable, "(?<=\\.).*")
)][]

epitime_mainl[measure == "prev" & aid %like% "base"] %>%
  ggplot(aes(x = at, y = value, color = anatsite)) +
  geom_line(alpha = 0.05) +
  facet_wrap(~ aid, nrow = 1) +
  theme_base()


epitime[aid %like% "main"] %>%
  ggplot(aes(x = at, y = prev.gc)) +
  geom_line(aes(group = simid), alpha = 0.05) +
  geom_line(
    data = timesum,
    aes(y = md),
    color = "cyan",
    size = 1
  ) +
  facet_wrap(~ aid, nrow = 1) +
  theme_base(base_size = 30) +
  theme(axis.text.x = element_text(angle = 90))

## epitime %>%
##   ggplot(aes(x = at, y = ir100.gc)) +
##   geom_line(aes(group = simid), alpha = 0.3) +
##   facet_wrap(~ aid, ncol = 5) +
##   theme_base(base_size = 30) +
##   theme(axis.text.x = element_text(angle = 90))


################################################################################
## GC SCREENING EFFECT DISTRIBUTIONS ##
################################################################################

epirel[groupid == "02.01"] %>%
  ggplot(aes(x = aid, y = pia_gc)) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_point(alpha = 0.1, position = position_jitter(width = 0.05)) +
  geom_boxplot(width = 0.2, fill = NA, outlier.alpha = 0) +
  theme_tufte(base_size = 35) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))


################################################################################
## PLOTS ##
################################################################################

## HIV
plotepi(hivp_vr, epi[analysis %like% "main_base"])
plotepi(hivi_vr, epi[analysis %like% "main_base"])
plotepi(hivp_ar, epi[analysis %like% "main_base"])

## Gonorrhea
plotepi(gcp_vg, epi[analysis %like% "main_base"])
plotepi(gci_vg, epi[analysis %like% "main_base"])
plotepi(gci_vp, epi[analysis %like% "main_base"])

plotepi(gcp_vg, epi[analysis %like% "main_cdc1"])
plotepi(gci_vg, epi[analysis %like% "main_cdc1"])
plotepi(gci_vp, epi[analysis %like% "main_cdc1"])

plotepi(gcp_vg, epi[analysis %like% "main_cdc2"])
plotepi(gci_vg, epi[analysis %like% "main_cdc2"])
plotepi(gci_vp, epi[analysis %like% "main_cdc2"])

plotepi(gcp_vg, epi[analysis %like% "main_symp"])
plotepi(gci_vg, epi[analysis %like% "main_symp"])
plotepi(gci_vp, epi[analysis %like% "main_symp"])

plotepi(gcp_vg, epi[analysis %like% "main_univ"])
plotepi(gci_vg, epi[analysis %like% "main_univ"])
plotepi(gci_vp, epi[analysis %like% "main_univ"])


################################################################################
## WRITE FILES ##
################################################################################

saveRDS(
  sresult_gcpop,
  file.path(simpath, "result_gc_summary.rds")
)

saveRDS(
  epirel,
  file.path(simpath, "result_screening_fx.rds")
)

saveRDS(
  epitime[aid %like% "main"],
  file.path(simpath, "result_epitime_main.rds")
)
