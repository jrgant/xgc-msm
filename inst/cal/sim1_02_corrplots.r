################################################################################
## SETUP ##
################################################################################

picksim <- "sim1"
source(here::here("inst", "cal", "sim0.0_setup.r"))
ls()


################################################################################
## CORRELATIONS BETWEEN MODEL OUTPUTS ##
################################################################################

 sel_regex <- paste(
  c("i.prev$", "i.prev.dx.inf\\.[A-Z]{1}$", "cc.vsupp$", "prepCov\\.[A-Z]{1}$",
    "cc.vsupp.[A-Z]{1}$", "prob.*GC", "prop.*tested", "ir100.pop"),
  collapse = "|"
)

selepi <- names(epi_mn)[names(epi_mn) %like% sel_regex]

outcorr <-
  ggcorrplot(
    cor(epi_mn[, ..selepi], method = "spearman"),
    outline.color = "black",
    type = "upper"
  ) +
  scale_fill_scico(
    name = "Spearman\ncorrelation\n",
    limits = c(-1, 1),
    palette = "vik"
  ) +
#  labs(caption = "\nSpearman correlations between selected model outputs") +
  theme(
    text              = element_text(size = 20, family = "sans"),
    axis.text.x       = element_text(size = 20),
    axis.text.y       = element_text(size = 20),
    legend.key.height = unit(1, "cm"),
    legend.title      = element_text(face = "bold"),
    plot.title        = element_text(face = "bold"),
    plot.caption      = element_text(face = "italic")
  )

outcorr

## Write plot to file.
psave("outcorr", outcorr)


################################################################################
## CORRELATIONS BETWEEN INPUTS AND OUTPUTS (PRE-SELECTION) ##
################################################################################

sel_popn <- epi_mn[simid %in% epi_mn_selnum[, simid], simid]
lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim1.rds"))

lhs_inputs <- rbindlist(
  lapply(
    lhs_groups,
    function(.x) {
      as.data.table(.x, keep.rownames = TRUE)
    }),
  idcol = "simid"
)[, simid := sprintf("%04d", simid)][simid %in% sel_popn]

setnames(lhs_inputs, c("rn", ".x"), c("param", "inval"))

## check to make sure all sel_popn simids included properly
all(lhs_inputs[, simid] %in% sel_popn & sel_popn %in% lhs_inputs[, simid])

io_join <- epi_mn[lhs_inputs, on = "simid"]

drop_static <- c(
  paste0("HIV_RX_HALT_PROB_", c("BLACK", "HISP", "OTHER", "WHITE"))
)

io_joinfil <- io_join[!param %in% drop_static]

## Function to create scatterplots of inputs vs. outputs and correlations.
io_scatter <- function(outcome, ylab = NULL, data = io_joinfil,
                       base_font_size = 7, axis_title_font_size = 14,
                       caption_font_size = 12, numcols = 9) {

  vars <- data[, unique(param)]
  cors <- lapply(
    setNames(vars, vars),
    function(.x) {
      data[param == .x, .(corr = cor(inval, get(outcome), method = "spearman"))]
    }
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
      geom_smooth(se = FALSE, color = "black", size = 0.5) +
      geom_text(
        data = corl[label_locs, on = "param"],
        aes(xloc, yloc, label = format(round(corr, 3), digit = 3)),
        color = "black",
        size = 2
      ) +
      facet_wrap(~param, scales = "free_x", ncol = numcols) +
      scale_color_scico(
        name = "Spearman correlation",
        limits = c(-1, 1),
        palette = "vik"
      ) +
      labs(
        y = paste0(yl, "\n"), x = "\nInput parameter value",
        caption = sprintf(
          paste(
            "\nAmong %s simulations that produced overall population sizes",
            "within 5 percent of the target (N = 20,000)"
          ),
          format(data[, length(unique(simid))], big.mark = ",")
        )
      ) +
    theme_base(base_size = base_font_size) +
    theme(
      axis.title = element_text(size = axis_title_font_size, face = "bold"),
      plot.caption = element_text(
        size = caption_font_size,
        hjust = 0,
        face = "italic"
      ),
      plot.margin = margin(1, 1, 1, 1, "cm"),
      legend.position = "top",
      legend.spacing.x = unit(1, "cm"),
      legend.key.width = unit(1, "cm"),
      legend.title = element_text(
        face = "bold", vjust = 1, size = axis_title_font_size
      )
    )

  return(plotout)
}

hivp_labs <- paste0("Simulated HIV prevalence", racelabs)

hivd_labs <- paste0(
  "Simulated HIV diagnosis prevalence among infected",
  racelabs
)

hivi_labs <- paste0("Simulated HIV incidence (per 100 person-years)", racelabs)
hivv_labs <- paste0("Simulated HIV viral suppression among diagnosed", racelabs)
prep_labs <- paste0("PrEP coverage among eligible MSM", racelabs)
gcpr_labs <- paste0("Simulated GC prevalence", anatlabs)
gcin_labs <- paste0("Simulated GC incidence", anatlabs)
gcpt_labs <- paste0("Simulated proportion tested for GC in clinic", anatlabs)

gcpp_labs <- paste0(
  "Simulated proportion positive for GC in clinic among tested",
  anatlabs
)

## HIV prevalence
psave("output_i.prev",   io_scatter("i.prev", hivp_labs[1]))
psave("output_i.prev.B", io_scatter("i.prev.B", hivp_labs[2]))
psave("output_i.prev.H", io_scatter("i.prev.H", hivp_labs[3]))
psave("output_i.prev.O", io_scatter("i.prev.O", hivp_labs[4]))
psave("output_i.prev.W", io_scatter("i.prev.W", hivp_labs[5]))

## HIV diagnosis prevalence among infected
psave("output_i.prev.dx.inf",   io_scatter("i.prev.dx.inf",   hivd_labs[1]))
psave("output_i.prev.dx.inf.B", io_scatter("i.prev.dx.inf.B", hivd_labs[2]))
psave("output_i.prev.dx.inf.H", io_scatter("i.prev.dx.inf.H", hivd_labs[3]))
psave("output_i.prev.dx.inf.O", io_scatter("i.prev.dx.inf.O", hivd_labs[4]))
psave("output_i.prev.dx.inf.W", io_scatter("i.prev.dx.inf.W", hivd_labs[5]))

## HIV incidence rate
psave("output_ir100.pop",   io_scatter("ir100.pop",   hivi_labs[1]))
psave("output_ir100.pop.B", io_scatter("ir100.pop.B", hivi_labs[2]))
psave("output_ir100.pop.H", io_scatter("ir100.pop.H", hivi_labs[3]))
psave("output_ir100.pop.O", io_scatter("ir100.pop.O", hivi_labs[4]))
psave("output_ir100.pop.W", io_scatter("ir100.pop.W", hivi_labs[5]))

## HIV viral suppression among diagnosed (by race/ethnicity)
psave("output_cc.vsupp",   io_scatter("cc.vsupp",   hivv_labs[1]))
psave("output_cc.vsupp.B", io_scatter("cc.vsupp.B", hivv_labs[2]))
psave("output_cc.vsupp.H", io_scatter("cc.vsupp.H", hivv_labs[3]))
psave("output_cc.vsupp.O", io_scatter("cc.vsupp.O", hivv_labs[4]))
psave("output_cc.vsupp.W", io_scatter("cc.vsupp.W", hivv_labs[5]))

## HIV viral suppression among diagnosed (by age group)
psave("output_cc.vsupp.age1", io_scatter("cc.vsupp", "VLS, Age 1"))
psave("output_cc.vsupp.age2", io_scatter("cc.vsupp", "VLS, Age 2"))
psave("output_cc.vsupp.age3", io_scatter("cc.vsupp", "VLS, Age 3"))
psave("output_cc.vsupp.age4", io_scatter("cc.vsupp", "VLS, Age 4"))
psave("output_cc.vsupp.age5", io_scatter("cc.vsupp", "VLS, Age 5"))

## Proportion tested for GC in clinic, by anatomic site
psave("output_prop.rect.tested",  io_scatter("prop.rect.tested",  gcpt_labs[2]))
psave("output_prop.ureth.tested", io_scatter("prop.ureth.tested", gcpt_labs[3]))
psave("output_prop.phar.tested",  io_scatter("prop.phar.tested",  gcpt_labs[4]))

## Proportion GC-positive among tested in clinic, by anatomic site
psave("output_prob.rGC.tested", io_scatter("prob.rGC.tested", gcpp_labs[2]))
psave("output_prob.uGC.tested", io_scatter("prob.uGC.tested", gcpp_labs[3]))
psave("output_prob.pGC.tested", io_scatter("prob.pGC.tested", gcpp_labs[4]))

## Proportion eligible for PrEP on PrEP
psave("output_prepCov",   io_scatter("prepCov", prep_labs[1]))
psave("output_prepCov.B", io_scatter("prepCov.B", prep_labs[2]))
psave("output_prepCov.H", io_scatter("prepCov.H", prep_labs[3]))
psave("output_prepCov.O", io_scatter("prepCov.O", prep_labs[4]))
psave("output_prepCov.W", io_scatter("prepCov.W", prep_labs[5]))
