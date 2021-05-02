################################################################################
## SETUP ##
################################################################################

library(pacman)

p_load(
  data.table,
  EpiModel,
  ggplot2,
  ggthemes,
  xgcmsm
)

theme_set(theme_tufte(base_size = 20))

targets <- readRDS(here::here("est", "caltargets.Rds"))

epi <- readRDS(here::here("burnin", "abc", "sim1", "sim1_epi.rds"))
epi

plotepi <- function(var, target = NULL, type = "line", data = epi) {
  p <- ggplot(data, aes(x = at, y = get(var))) +
    ylab(var)

  if (type == "line") p <- p + geom_line(aes(group = simid), alpha = 0.01)
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

plotepi("num", 20000)

plotepi("i.prev", targets$ct_hiv_prev)

plotepi("ir100", targets$ct_hiv_incid_100k / 1000)
