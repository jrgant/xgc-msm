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


################################################################################
## CORRELATIONS BETWEEN MODEL OUTPUTS ##
################################################################################

sel_regex <- paste(
  c("i.prev$", "i.prev.dx.inf.*[A-Z]{1}$", "cc.vsupp$",
    "cc.vsupp.[A-Z]{1}$", "prob.*GC", "prop.*tested", "ir100$"),
  collapse = "|"
)

selepi <- names(epi_mn)[names(epi_mn) %like% sel_regex]

outcorr <-
  ggcorrplot(
    cor(epi_mn[, ..selepi]),
    outline.color = "black",
    type = "upper"
  ) +
  scale_fill_scico(
    name = "",
    limits = c(-1, 1),
    palette = "vik"
  ) +
  ggtitle("Correlations between Selected Model Outputs") +
  theme(
    text = element_text(size = 18, family = "sans"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.key.height = unit(1, "in"),
    plot.title = element_text(face = "bold")
  )

outcorr

ggsave(
  here::here("inst", "calfigs", batch_date, "outcorr.pdf"),
  outcorr,
  width = 10, height = 10, units = "in"
)


################################################################################
## CORRELATIONS BETWEEN INPUTS AND OUTPUTS ##
################################################################################

sel_popn <- epi_mn[simid %in% epi_mn_selnum[, simid], simid]

lhs_groups <- readRDS(here::here("burnin/abc/sim1", "lhs_sim1.rds"))

lhs_inputs <- rbindlist(
  lapply(
    lhs_groups,
    function(.x) {
      as.data.table(.x, keep.rownames = TRUE)
    }),
  idcol = "simid"
)[simid %in% sel_popn]

setnames(lhs_inputs, c("rn", ".x"), c("param", "inval"))

io_join <- epi_mn[lhs_inputs, on = "simid"]

drop_static <- c(
  paste0("HIV_RX_HALT_PROB_", c("BLACK", "HISP", "OTHER", "WHITE")),
  "SCALAR_AI_ACT_RATE",
  "SCALAR_OI_ACT_RATE",
  paste0("SCALAR_HIV_TRANS_PROB_", c("BLACK", "HISP", "OTHER", "WHITE"))
)

io_joinfil <- io_join[!param %in% drop_static]

io_scatter <- function(outcome, ylab = NULL, data = io_joinfil) {

  vars <- data[, unique(param)]
  cors <- lapply(
    setNames(vars, vars),
    function(.x) data[param == .x, .(corr = cor(inval, get(outcome)))]
  )

  yrange <- data[, get(outcome)]
  label_locs <- data[, .(xl = min(inval), xu = max(inval)), by = param]
  label_locs[, xloc := xl + ((xu - xl) / 2)]
  label_locs[, yloc := min(yrange) + ((max(yrange) - min(yrange)) / 1.11)]

  corl <- rbindlist(cors, idcol = "param")

  # print(label_locs)
  # print(corl)
  ifelse(!is.null(ylab), yl <- ylab, yl <- outcome)

  plotout <-
    ggplot(
      corl[data, on = "param"], # join correlation data
      aes(inval, get(outcome), color = corr)
    ) +
      geom_point(size = 0.5) +
      geom_smooth(se = FALSE, color = "black") +
      geom_text(
        data = corl[label_locs, on = "param"],
        aes(xloc, yloc, label = format(round(corr, 3), digit = 3)),
        color = "black",
        size = 6
      ) +
      facet_wrap(~param, scales = "free_x") +
      scale_color_scico(
        name = "Pearson correlation",
        limits = c(-1, 1),
        palette = "vik"
      ) +
      labs(
        y = yl, x = "Input parameter value",
        caption = sprintf(
          paste(
            "\n Among %s simulations that produced overall population sizes",
            "within 5 percent of the target (N = 20,000)"
          ),
          format(data[, length(unique(simid))], big.mark = ",")
        )
      ) +
    theme_base(base_size = 18) +
    theme(
      axis.title = element_text(size = 28),
      plot.caption = element_text(size = 18, hjust = 0, face = "italic"),
      legend.position = "top",
      legend.spacing.x = unit(1, "cm"),
      legend.key.width = unit(1, "in"),
      legend.title = element_text(face = "bold", vjust = 1)
    )

  return(plotout)
}

## Plot label sets.
racelabs <- c("", ", Black", ", Hispanic", ", Other", ", White")
anatlabs <- c("", ", Rectal", ", Urethral", ", Pharyngeal")

hivp_labs <- paste0("Simulated HIV prevalence", racelabs)
hivd_labs <- paste0(
  "Simulated HIV diagnosis prevalence among infected",
  racelabs
)
hivi_labs <- paste0("Simulated HIV incidence (per 100 person-years)", racelabs)
hivv_labs <- paste0("Simulated HIV viral suppression among diagnosed", racelabs)
gcpr_labs <- paste0("Simulated GC prevalence", anatlabs)
gcin_labs <- paste0("Simulated GC incidence", anatlabs)
gcpt_labs <- paste0("Simulated proportion tested for GC in clinic", anatlabs)
gcpp_labs <- paste0(
  "Simulated proportion positive for GC in clinic among tested",
  anatlabs
)


## HIV prevalence
psave("i.prev", io_scatter("i.prev",   hivp_labs[1]))
psave("i.prev.B", io_scatter("i.prev.B", hivp_labs[2]))
psave("i.prev.H", io_scatter("i.prev.H", hivp_labs[3]))
psave("i.prev.O", io_scatter("i.prev.O", hivp_labs[4]))
psave("i.prev.W", io_scatter("i.prev.W", hivp_labs[5]))

## HIV diagnosis prevelance among infected
psave("i.prev.dx.inf", io_scatter("i.prev.dx.inf",   hivd_labs[1]))
psave("i.prev.dx.inf.B", io_scatter("i.prev.dx.inf.B", hivd_labs[2]))
psave("i.prev.dx.inf.H", io_scatter("i.prev.dx.inf.H", hivd_labs[3]))
psave("i.prev.dx.inf.O", io_scatter("i.prev.dx.inf.O", hivd_labs[4]))
psave("i.prev.dx.inf.W", io_scatter("i.prev.dx.inf.W", hivd_labs[5]))

## HIV incidence rate
psave("ir100", io_scatter("ir100",   hivi_labs[1]))
psave("ir100.B", io_scatter("ir100.B", hivi_labs[2]))
psave("ir100.H", io_scatter("ir100.H", hivi_labs[3]))
psave("ir100.O", io_scatter("ir100.O", hivi_labs[4]))
psave("ir100.W", io_scatter("ir100.W", hivi_labs[5]))

## HIV viral suppression among diagnosed
psave("cc.vsupp", io_scatter("cc.vsupp",   hivv_labs[1]))
psave("cc.vsupp.B", io_scatter("cc.vsupp.B", hivv_labs[2]))
psave("cc.vsupp.H", io_scatter("cc.vsupp.H", hivv_labs[3]))
psave("cc.vsupp.O", io_scatter("cc.vsupp.O", hivv_labs[4]))
psave("cc.vsupp.W", io_scatter("cc.vsupp.W", hivv_labs[5]))

## GC prevalence
psave("prev.gc", io_scatter("prev.gc",  gcpr_labs[1]))
psave("prev.rgc", io_scatter("prev.rgc", gcpr_labs[2]))
psave("prev.ugc", io_scatter("prev.ugc", gcpr_labs[3]))
psave("prev.pgc", io_scatter("prev.pgc", gcpr_labs[4]))

## GC incidence rate, by anatomic site
psave("ir100.gc", io_scatter("ir100.gc",  gcin_labs[1]))
psave("ir100.rgc", io_scatter("ir100.rgc", gcin_labs[2]))
psave("ir100.ugc", io_scatter("ir100.ugc", gcin_labs[3]))
psave("ir100.pgc", io_scatter("ir100.pgc", gcin_labs[4]))

## Proportion tested for GC in clinic, by anatomic site
psave("prop.rect.tested", io_scatter("prop.rect.tested",  gcpt_labs[2]))
psave("prop.ureth.tested", io_scatter("prop.ureth.tested", gcpt_labs[3]))
psave("prop.phar.tested", io_scatter("prop.phar.tested",  gcpt_labs[4]))

## Proportion GC-positive among tested in clinic, by anatomic site
psave("prob.rGC.tested", io_scatter("prob.rGC.tested", gcpp_labs[2]))
psave("prob.uGC.tested", io_scatter("prob.uGC.tested", gcpp_labs[3]))
psave("prob.pGC.tested", io_scatter("prob.pGC.tested", gcpp_labs[4]))


################################################################################
## ABSOLUTE DIFFERENCES BETWEEN MODEL OUTPUT AND TARGETS ##
################################################################################

get_abs_stdiff <- function(target, output, data = epi_mn_selnum) {
  cols <- c("simid", output)
  d <- data[, ..cols]
  d[, abs_std_diff := abs((get(output) - target) / sd(get(output)))]
  dout <- d[order(abs_std_diff)]
  return(dout)
}

t_hiv_inc <- get_abs_stdiff(targets$ct_hiv_incid_100k / 1000, "ir100")
t_hiv_prv <- get_abs_stdiff(targets$ct_hiv_prev, "i.prev")

# diagnosis prevalence among infected, by race
t_hiv_dx_race <- lapply(
  names(targets$ct_hivdx_pr_byrace),
  function(.x) {
    get_abs_stdiff(
      targets$ct_hivdx_pr_byrace[[.x]],
      paste0("i.prev.dx.inf.", .x)
    )
  }
)

# diagnosis prevalence among infected, by age
t_hiv_dx_age <- lapply(
  names(targets$ct_hivdx_pr_byage),
  function(.x) {
    get_abs_stdiff(targets$ct_hivdx_pr_byage[[.x]], paste0("i.prev.dx.inf.", .x))
  }
)

# viral suppression among HIV-diagnosed, by age group
## t_vls_age <- lapply(
##   names(targets$ct_vls_pr_byage),
##   function(.x) {
##     get_abs_stdiff(targets$ct_vls_pr_byage[[.x]], paste0("cc.vsupp.", .x))
##   }
## )

# viral suppression among HIV-diagnosed, by race
t_vls_race <- lapply(
  names(targets$ct_vls_pr_byrace),
  function(.x) {
    get_abs_stdiff(targets$ct_vls_pr_byrace[[.x]], paste0("cc.vsupp.", .x))
  }
)


# gonorrhea testing in clinic, by anatomic site
t_gc_anat_test <- lapply(
  names(targets$ct_prop_anatsite_tested),
  function(.x) {
    var <- sprintf("prop.%s.tested", .x)
    pass1 <- get_abs_stdiff(
      targets$ct_prop_anatsite_tested[[.x]],
      var
    )
    pass1
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
