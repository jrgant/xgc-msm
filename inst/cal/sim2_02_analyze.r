################################################################################
## SETUP ##
################################################################################

picksim <- "sim2"

# NOTE The file dated 2021-7-27 contains entire burn-in period.
#      Needed to investigate extinction of GC in selected runs and lack of
#      equilibrium in HIV dynamics. We will narrow parameter sets based on
#      overall HIV prevalence in this script, and then take another
#      random Latin hypercube sample for sim3. Based on looking at 2021-07-27
#      it seems like we should get parameter sets that result in HIV
#      equilibrium and non-extinction of GC.
# Sys.setenv("BATCH_DATE" = "2021-07-26")

ic_dir <- here::here("inst", "cal")
source(here::here(ic_dir, "sim0.0_setup.r"))
ls()
sim_epis[sim_epis %like% batch_date]


################################################################################
## VISUALIZE DISTRIBUTIONS OF TARGETS ACROSS ALL SIMS ##
################################################################################
source(here::here(ic_dir, "sim0.1_fmt_targetdat.r"))
ls()


################################################################################
## SELECT SIMULATIONS ##
################################################################################
simid_sel_hivpr <-
  out_vs_targ[variable == "i.prev" & output_within_5pts == 1, simid]

## NOTE
## quicktargets() is a helper function and particular about its inputs,
## so we make some lists to get several targets one plot
set_jitter <- list(w = 0.1, h = 0)
coerced_simid_gcpos <- list(
  "rGC" = simid_sel_hivpr,
  "uGC" = simid_sel_hivpr,
  "pGC" = simid_sel_hivpr
)

cs_gctest <- coerced_simid_gcpos
names(cs_gctest) <- c("rect", "ureth", "phar")

cs_vsupp_reth <- lapply(1:4, function(x) simid_sel_hivpr)
names(cs_vsupp_reth) <- rlabs

quicktarget("i.prev", list(hiv = simid_sel_hivpr), set_jitter)
quicktarget("ir100.pop", list(hiv = simid_sel_hivpr), set_jitter)
quicktarget("cc.vsupp.%s", cs_vsupp_reth, set_jitter)
quicktarget("prop.%s.tested", cs_gctest, set_jitter)
quicktarget("prob.%s.tested", coerced_simid_gcpos, set_jitter)


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim2.rds"))

hivpr_sel_inputs <- pull_params(list(simid_sel_hivpr))[, -c("selection_group")]

all(
  all(simid_sel_hivpr %in% hivpr_sel_inputs[, simid]),
  all(hivpr_sel_inputs[, simid] %in% simid_sel_hivpr)
) == TRUE

## Since only a handful of simids were selected, get the min/max for each
## input parameter to use as the new sampling ranges in sim3.
new_priors <- hivpr_sel_inputs[, .(
  s3_ll = min(value),
  s3_ul = max(value)
), input]

sim2_priors <- readRDS(here::here(ic_dir, "sim1_sel_lhs_limits.rds"))

prior_bind <- rbindlist(
  list(new_priors[, batch := "s3"][], sim2_priors[, batch := "s2"][]),
  use.names = FALSE
)

setnames(prior_bind, c("s3_ll", "s3_ul"), c("ll", "ul"))

ulcoord <- dcast(prior_bind, input ~ batch, value.var = c("batch", "ul"))
llcoord <- dcast(prior_bind, input ~ batch, value.var = c("batch", "ll"))

pcoord <- rbindlist(
  list(
    ulcoord[, .(input, x = batch.1_s2, y = ul_s2)],
    ulcoord[, .(input, x = batch.1_s3, y = ul_s3)],
    llcoord[, .(input, x = batch.1_s2, y = ll_s2)],
    llcoord[, .(input, x = batch.1_s3, y = ll_s3)]
  )
)[order(input)]

pcoord[, x := fcase(x == "s2", 1, x == "s3", 2)]
pcoord[, rowid := rowid(input)][]

sort_ref <- pcoord[, -c("x", "y", "rowid")]
sort_ref[, rowid := c(3, 4, 2, 1), by = input] # put rowids in the order I want
merge(sort_ref, pcoord, by = c("input", "rowid"))

pcoord_fin <- pcoord[sort_ref, on = c("input", "rowid")][!(input %like% "HALT")]
pcoord_fin[, batch := fifelse(x == 1, "s2", "s3")]

geom_polyside <- function(xvals, lt, ...) {
  geom_line(
    data = pcoord_fin[rowid %in% xvals],
    aes(x = x, y = y),
    linetype = lt,
    ...
  )
}

# TODO Consider a figure that shows the whole history of the calibration
#      Would probably want to move into its own file
ggplot(pcoord_fin, aes(x = x, y = y, fill = batch, color = batch)) +
  geom_rect(
    data = dcast(
      pcoord_fin[x == 1], input ~ rowid, value.var = "y"
    )[, batch := "s2"][],
    aes(
      x = NULL, y = NULL, color = NULL, fill = NULL,
      xmin = 1, xmax = 2, ymin = `3`, ymax = `1`
    ),
    fill = "gray90",
    color = "white"
  ) +
  geom_polygon(fill = "white", color = "white") +
  geom_segment(
    data = pcoord_fin[x == 2][, batch := "s3"][],
    aes(x = 1, xend = 2, y = y, yend = y),
    linetype = "dashed", color = "black"
  ) +
  xlim(c(1, 2)) +
  ## geom_polyside(3:4, "dotted", color = "black") +
  ## geom_polyside(1:2, "dotted", color = "black") +
  geom_polyside(c(2, 4), "solid", color = "black") +
  geom_polyside(c(1, 3), "solid", color = "black") +
  geom_point(aes(fill = batch), size = 5, shape = 21, color = "black") +
  facet_wrap(~ input, scales = "free_y") +
  scale_fill_scico_d(
    name = "Prior set", palette = "berlin", begin = 0.1, end = 0.2
  ) +
  scale_color_scico_d(
    name = "Prior set", palette = "berlin", begin = 0.75, end = 0.75
  ) +
  scale_x_continuous(breaks = c(1, 2), labels = c("s2", "s3")) +
  guides(color = FALSE) +
  theme_tufte(base_size = 8)


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(simid_sel_hivpr, here::here("inst", "cal", "sim2_simid_sel.rds"))

## Write pre-selection boxplots to files
lapply(
  ls(pattern = "^ps_"),
  function(x) {
    psave(
      f = sprintf("preselect_%s", stringr::str_extract(x, "(?<=ps_).*")),
      p = get(x)
    )
    invisible()
  }
)

## Write limits of input parameter values across selected simulations.
saveRDS(new_priors, here::here("inst", "cal", "sim2_sel_lhs_limits.rds"))
