###############################################################################
## SETUP ##
################################################################################

picksim <- "sim1"
ic_dir <- here::here("inst", "cal")
source(here::here(ic_dir, "sim0.0_setup.r"))
ls()

sim_epis[sim_epis %like% batch_date]
pretty_batch_date


################################################################################
## VISUALIZE DISTRIBUTIONS OF TARGETS ACROSS ALL SIMS ##
################################################################################
source(here::here(ic_dir, "sim0.1_fmt_targetdat.r"))
ls()


################################################################################
## SELECT SIMULATIONS ##
################################################################################

simid_sel_vls <- lapply(
  setNames(rlabs, rlabs),
  function(.x) {
    out_vs_targ[variable == paste0("cc.vsupp.", .x) & output_within_5pts, simid]
  }
)

simid_sel_prep <- lapply(
  setNames(rlabs, rlabs),
  function(.x) {
    out_vs_targ[variable == paste0("prepCov.", .x) & output_within_5pts, simid]
  }
)

episel <- epi_noNA[simid %in% unlist(list(simid_sel_vls, simid_sel_prep))]

out_vs_targ[simid %in% unlist(simid_sel_vls), .(
  min = min(output),
  max = max(output),
  target = mean(target_val)
), variable]

quicktarget("cc.vsupp.%s", simid_sel_vls)
quicktarget("prepCov.%s", simid_sel_prep)


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim1.rds"))

vls_sel_inputs <- pull_params(simid_sel_vls)[input %like% "RX_INIT|RX_REINIT"]
vls_sel_inputs[, input_group := stringr::str_extract(input, "RX_[A-Z]+")][]

vls_sel_inputs_w <- dcast(
  vls_sel_inputs,
  selection_group + simid ~ input,
  value.var = "value"
)

prep_sel_inputs <- pull_params(simid_sel_prep)[input %like% "PREP"]

prep_sel_inputs_w <- dcast(
  prep_sel_inputs,
  selection_group + simid ~ input,
  value.var = "value"
)

# Function to check to make sure we've captured all selected simids.
check_input_simids <- function(x, simid_list, simid_inputs) {
  data.table(
    captured = sort(simid_inputs[selection_group == x, unique(simid)]),
    sel = sort(simid_list[[x]])
  )[, .N, .(captured, sel)][, all(N) == 1]
}

# check selected VLS inputeter sets
sapply(
  rlabs,
  check_input_simids,
  simid_list = simid_sel_vls,
  simid_inputs = vls_sel_inputs
)

# check selected PrEP inputeter sets
sapply(
  rlabs,
  check_input_simids,
  simid_list = simid_sel_prep,
  simid_inputs = prep_sel_inputs
)

## View distributions of RX_INIT and RX_REINIT for each selection group,
## meaning simid sets chosen for matching race/ethnicity-specific viral
## suppression. We're interested in seeing the how the corresponding
## input parameter distribution (e.g., HIV_RX_INIT_PROB_BLACK) within
## the B selection group. Here, the REINIT parameters are much stronger
## drivers of VLS.
ggplot(vls_sel_inputs, aes(input, value)) +
  geom_boxplot() +
  facet_grid(selection_group ~ input_group, scale = "free_x") +
  theme_base(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90))

corplots_vls <- lapply(rlabs, function(.x) {
  tmp <- vls_sel_inputs_w[selection_group == .x, 3:ncol(vls_sel_inputs_w)]
  ggcorrplot(
    cor(tmp, method = "spearman"),
    type = "lower",
    lab = TRUE,
    title = sprintf("Selection Group: %s", .x),
    p.mat = cor_pmat(tmp, method = "spearman"),
    sig.level = 0.05
  )
})

gridExtra::grid.arrange(
  corplots_vls[[1]],
  corplots_vls[[2]],
  corplots_vls[[3]],
  corplots_vls[[4]]
)

corplots_prep <- lapply(rlabs, function(.x) {
  tmp <- prep_sel_inputs_w[selection_group == .x, 3:ncol(prep_sel_inputs_w)]
  ggcorrplot(
    cor(tmp, method = "spearman"),
    type = "lower",
    lab = TRUE,
    title = sprintf("Selection Group: %s", .x),
    p.mat = cor_pmat(tmp, method = "spearman"),
    sig.level = 0.05
  )
})

gridExtra::grid.arrange(
  corplots_prep[[1]],
  corplots_prep[[2]],
  corplots_prep[[3]],
  corplots_prep[[4]]
)

## Function that takes the long version of inputeter subsets from
## simid selections and extracts the 25% and 75% quantiles for use as prior
## limits in the next round of calibration.
input_quantiles <- function(data, inputstring, ql = 0.25, qu = 0.75) {
  inputsub <- data[
    input %like% inputstring,
    .(
      q25 = quantile(value, ql),
      q75 = quantile(value, qu)
    ),
    .(selection_group, input)
    ## this step matches the race/eth selection group to the corresponding
    ## input inputeter
  ][selection_group == substring(str_extract(input, "(?<=_)[A-Z]+$"), 1, 1)]

  inputsub[, -c("selection_group")]
}

narrow_rx_init <-   input_quantiles(vls_sel_inputs, "RX_INIT")
narrow_rx_reinit <- input_quantiles(vls_sel_inputs, "RX_REINIT")
narrow_prep_discont <- input_quantiles(prep_sel_inputs, "DISCONT")

narrowed <- rbindlist(
  list(narrow_rx_init, narrow_rx_reinit, narrow_prep_discont)
)

## Create a new set of priors for sim2
sim1_priors <- readRDS(here::here("burnin", "cal", "sim1", "sim1_priors.rds"))

new_priors <- merge(sim1_priors, narrowed, by = "input", all.x = TRUE)
new_priors <- new_priors[
  !is.na(q25) & !is.na(q75), ":="(s1_ll = q25, s1_ul = q75)][, -c("q25", "q75")]

setnames(new_priors, c("s1_ll", "s1_ul"), c("s2_ll", "s2_ul"))

## Add a prior for the act.stopper.prob parameter, which governs the proportion
## of GC-symptomatic or -treated patients cease sexual activity while these
## conditions hold. The limits are somewhat arbitrary, but the empirical
## estimates from Beck/Kramer of 0.8 are from 1980 and applied only to GC
## symptoms and reduction of sex.
new_priors <- rbindlist(
  list(
    new_priors,
    data.table(s2_ll = 0.2, s2_ul = 1, input = "ACT_STOPPER_PROB")
  ),
  use.name = TRUE
)


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(
  list(
    vls = simid_sel_vls,
    prep = simid_sel_prep
  ),
  here::here("inst", "cal", "sim1_simid_sel.rds")
)

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
saveRDS(new_priors, here::here("inst", "cal", "sim1_sel_lhs_limits.rds"))
