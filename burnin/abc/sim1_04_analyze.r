# Plot -------------------------------------------------------------------------
pacman::p_load(data.table, magrittr, ggplot2, ggthemes)

dl <- fread(here::here("burnin", "abc", "sim1", "sim1_epi.csv"))

theme_set(theme_tufte(base_size = 20))

plotg <- function(outcomes, data = dl) {
  melt(data,
     id.vars = c("sim", "at"),
     measure.vars = outcomes
     ) %>%
  ggplot(aes(x = at, y = value)) +
  geom_line(aes(group = sim), alpha = 0.1) +
  facet_wrap(~ variable)
}

racelabs <- c("B", "H", "O", "W")
anatlabs <- c("rgc", "ugc", "pgc")

plotg(c("num", paste0("num.", racelabs)))

plotg(paste0("ir100.", racelabs))
plotg(paste0("incid.", racelabs))
plotg(c("i.prev", paste0("i.prev.", racelabs)))
plotg(paste0("i.prev.dx.", racelabs))
plotg(paste0("prepElig.", racelabs))
plotg(paste0("tot.tests.", racelabs))
plotg(paste0("hstage.", c("acute", "chronic", "aids")))
plotg(paste0("cc.dx.", racelabs))

plotg(paste0("incid.", anatlabs))
plotg(paste0("prop.", c("rect", "ureth", "phar"), ".tested"))
plotg(paste0("prob.", c("rGC", "uGC", "pGC"), ".tested"))
