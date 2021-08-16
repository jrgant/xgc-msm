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
  setNames(rslugs, rslugs),
  function(.x) {
    out_vs_targ[
      variable == paste0("cc.vsupp.", .x) &
      output_within_5pts == 1, sort(simid)]
  }
)

simid_sel_prep <- lapply(
  setNames(rslugs, rslugs),
  function(.x) {
    out_vs_targ[
      variable == paste0("prepCov.", .x) &
      output_within_5pts == 1, sort(simid)]
  }
)

simid_sel_gcpos <- lapply(
  setNames(gcpos_slugs, gcpos_slugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("prob.%s.tested", .x) &
      output_within_5pts, sort(simid)]
  }
)

simid_sel_hivpr <-
  out_vs_targ[variable == "i.prev" & output_within_5pts == 1, sort(simid)]

simid_sel_hivinc <-
  out_vs_targ[variable == "ir100.pop" & between(output, 0.4, 0.6), sort(simid)]

simid_sel_hivprinc <- intersect(simid_sel_hivpr, simid_sel_hivinc)

fmt_list <- function(labvec, simids) {
  lapply(setNames(labvec, labvec), function(.x) simids)
}

quicktarget("i.prev",       list(hiv = simid_sel_hivpr))
quicktarget("ir100.pop",    list(hiv = simid_sel_hivpr))
quicktarget("cc.vsupp.%s",  fmt_list(rslugs, simid_sel_hivpr))
quicktarget("prepCov.%s",   fmt_list(rslugs, simid_sel_hivpr))
quicktarget("prob.%s.tested", fmt_list(gcpos_slugs, simid_sel_hivpr))

quicktarget("i.prev",       list(hiv = simid_sel_hivinc))
quicktarget("ir100.pop",    list(hiv = simid_sel_hivinc))
quicktarget("cc.vsupp.%s",  fmt_list(rslugs, simid_sel_hivinc))
quicktarget("prepCov.%s",   fmt_list(rslugs, simid_sel_hivinc))
quicktarget("prob.%s.tested", fmt_list(gcpos_slugs, simid_sel_hivinc))

quicktarget("i.prev",       list(hiv = simid_sel_hivprinc))
quicktarget("ir100.pop",    list(hiv = simid_sel_hivprinc))
quicktarget("cc.vsupp.%s",  fmt_list(rslugs, simid_sel_hivprinc))
quicktarget("prepCov.%s",   fmt_list(rslugs, simid_sel_hivprinc))
quicktarget("prob.%s.tested", fmt_list(gcpos_slugs, simid_sel_hivprinc))

quicktarget("i.prev",       simid_sel_vls) + facet_wrap(~selection_group)
quicktarget("ir100.pop",    simid_sel_vls) + facet_wrap(~selection_group)
quicktarget("cc.vsupp.%s",  simid_sel_vls)
quicktarget("prepCov.%s",   simid_sel_vls)
quicktarget("prepCov.%s",   simid_sel_vls)

quicktarget("i.prev",       simid_sel_prep) + facet_wrap(~selection_group)
quicktarget("ir100.pop",    simid_sel_prep) + facet_wrap(~selection_group)
quicktarget("cc.vsupp.%s",  simid_sel_prep)
quicktarget("prepCov.%s",   simid_sel_prep)

quicktarget("i.prev",       simid_sel_gcpos) + facet_wrap(~selection_group)
quicktarget("ir100.pop",    simid_sel_gcpos) + facet_wrap(~selection_group)
quicktarget("cc.vsupp.%s",  fmt_list(rslugs, simid_sel_gcpos[[1]]))
quicktarget("cc.vsupp.%s",  fmt_list(rslugs, simid_sel_gcpos[[2]]))
quicktarget("cc.vsupp.%s",  fmt_list(rslugs, simid_sel_gcpos[[3]]))
quicktarget("prepCov.%s",   fmt_list(rslugs, simid_sel_gcpos[[1]]))
quicktarget("prepCov.%s",   fmt_list(rslugs, simid_sel_gcpos[[2]]))
quicktarget("prepCov.%s",   fmt_list(rslugs, simid_sel_gcpos[[3]]))


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim1.rds"))

vls_sel_inputs <- pull_params(
  simid_sel_vls
)[input %like% "RX_INIT|RX_REINIT"]

vls_sel_inputs[,
  input_group := stringr::str_extract(input, "RX_[A-Z]+")][]

vls_sel_inputs_w <- dcast(
  vls_sel_inputs,
  selection_group + simid ~ input,
  value.var = "value"
)

#[simid %in% simid_sel_hivprinc][]

cbind(
  vls_sel_inputs_w[, .(N_input_pull = .N), selection_group],
  N_simid_list = sapply(simid_sel_vls, length)
)

prep_sel_inputs <- pull_params(simid_sel_prep)[input %like% "PREP"]

prep_sel_inputs_w <- dcast(
  prep_sel_inputs,
  selection_group + simid ~ input,
  value.var = "value"
)

#[simid %in% simid_sel_hivprinc]

cbind(
  prep_sel_inputs_w[, .(N_input_pull = .N), selection_group],
  N_simid_list = sapply(simid_sel_prep, length)
)

hivprinc_sel_inputs <-
  pull_params(
    simid_sel_hivprinc
  )[input %like% "EFF_HIV|AI_ACT"][, -c("selection_group")]

hivprinc_sel_inputs_w <- dcast(
  hivprinc_sel_inputs,
  simid ~ input,
  value.var = "value"
)

#[simid %in% simid_sel_hivprinc]

## gcpos_sel_inputs <-
##   pull_params(
##     simid_sel_gcpos
##   )[input %like% "DURAT_NOTX|RX_INFPR|SYMPT_PROB|ASYMP"]

## gcpos_sel_inputs[,
##   input_group := stringr::str_extract(input, "(?<=_)[A-Z0-9]+_[A-Z0-9]+$")]

## gcpos_sel_inputs_w <- dcast(
##   gcpos_sel_inputs,
##   selection_group + simid ~ input,
##   value.var = "value"
## )[simid %in% simid_sel_hivprinc]

## gcpos_sel_inputs_w[, .N, selection_group]

# Function to check to make sure we've captured all selected simids.
check_input_simids <- function(x, simid_list, simid_inputs) {
  data.table(
    captured = sort(simid_inputs[selection_group == x, unique(simid)]),
    sel = sort(simid_list[[x]])
  )[, .N, .(captured, sel)][, all(N) == 1]
}

# check selected VLS input parameter sets
sapply(
  rslugs,
  check_input_simids,
  simid_list = simid_sel_vls,
  simid_inputs = vls_sel_inputs
)

# check selected PrEP input parameter sets
sapply(
  rslugs,
  check_input_simids,
  simid_list = simid_sel_prep,
  simid_inputs = prep_sel_inputs
)

# check select HIV input parameter sets
data.table(
  capture = sort(hivprinc_sel_inputs[, unique(simid)]),
  sel = sort(simid_sel_hivprinc)
)[, all(capture %in% sel) & all(sel %in% capture)]

# check selected GC input parameter sets
## sapply(
##   gcpos_slugs,
##   check_input_simids,
##   simid_list = simid_sel_gcpos,
##   simid_inputs = gcpos_sel_inputs
## )

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

## corplots_vls <- lapply(rslugs, function(.x) {
##   tmp <- vls_sel_inputs_w[selection_group == .x, 3:ncol(vls_sel_inputs_w)]
##   ggcorrplot(
##     cor(tmp, method = "spearman"),
##     type = "lower",
##     lab = TRUE,
##     title = sprintf("Selection Group: %s", .x),
##     p.mat = cor_pmat(tmp, method = "spearman"),
##     sig.level = 0.05
##   )
## })

## gridExtra::grid.arrange(
##   corplots_vls[[1]],
##   corplots_vls[[2]],
##   corplots_vls[[3]],
##   corplots_vls[[4]]
## )

## corplots_prep <- lapply(rslugs, function(.x) {
##   tmp <- prep_sel_inputs_w[selection_group == .x, 3:ncol(prep_sel_inputs_w)]
##   ggcorrplot(
##     cor(tmp, method = "spearman"),
##     type = "lower",
##     lab = TRUE,
##     title = sprintf("Selection Group: %s", .x),
##     p.mat = cor_pmat(tmp, method = "spearman"),
##     sig.level = 0.05
##   )
## })

## gridExtra::grid.arrange(
##   corplots_prep[[1]],
##   corplots_prep[[2]],
##   corplots_prep[[3]],
##   corplots_prep[[4]]
## )

narrow_rx_init      <- input_quantiles(vls_sel_inputs, "RX_INIT")
narrow_rx_reinit    <- input_quantiles(vls_sel_inputs, "RX_REINIT")
narrow_late_tester  <- input_quantiles(vls_sel_inputs, "LATE_TESTER")
narrow_trans_prob   <- input_quantiles(vls_sel_inputs, "TRANS")

narrow_prep_discont <- input_quantiles(prep_sel_inputs, "DISCONT")


narrow_aiscale_hivcond <- hivprinc_sel_inputs[,
  .(q25 = quantile(value, 0.25), q75 = quantile(value, 0.75)),
  input
]

## gcpos_sel_inputs_labeled <- copy(gcpos_sel_inputs)

## gcpos_sel_inputs_labeled[,
##   selection_group := fcase(
##     selection_group == "rGC", "R",
##     selection_group == "uGC", "U",
##     selection_group == "pGC", "P"
##   )]

anatsite_pat <- "^[A-Z]+(?=_)"

## ## Narrow anatomic site-specific parameters
## narrow_gc_durat <- input_quantiles(
##   gcpos_sel_inputs_labeled,
##   "DURAT_NOTX",
##   anatsite_pat
## )

## narrow_gc_rxinf <- input_quantiles(
##   gcpos_sel_inputs_labeled,
##   "RX_INFPR",
##   anatsite_pat
## )

## narrow_gc_sympt <- input_quantiles(
##   gcpos_sel_inputs_labeled,
##   "SYMPT_PROB",
##   anatsite_pat
## )

## narrow_gc_asymp_test <- input_quantiles(
##   gcpos_sel_inputs_labeled,
##   "ASYMP",
##   anatsite_pat
## )

## ## Narrow transmission pathway, universal prob. of STI testing
## ## due to symptoms, and act stopper prob.
## gcpos_union <- gcpos_sel_inputs[, sort(unique(simid))]

## lhs_gcpos_union <- rbindlist(
##   lapply(
##     gcpos_union,
##     function(.x) {
##       as.data.table(
##         lhs_groups[[as.numeric(.x)]], keep.rownames = TRUE
##       )[V1 %like% "CONDOM_EFF_GC|GC_PROB|^STITEST|STOPPER"
##         ][, .(input = V1, value = V2)][, simid := .x][]
##     }
##   )
## )

## narrow_gctrans_common <-
##   lhs_gcpos_union[
##     simid %in% simid_sel_hivprinc][, .(
##                    q25 = quantile(value, 0.25),
##                    q75 = quantile(value, 0.75)
##                  ), by = input]

narrowed <- rbindlist(
  list(
    narrow_rx_init,
    narrow_rx_reinit,
    narrow_late_tester,
    narrow_trans_prob,
    narrow_prep_discont,
    narrow_aiscale_hivcond
    ## narrow_gc_durat,
    ## narrow_gc_rxinf,
    ## narrow_gc_sympt,
    ## narrow_gc_asymp_test,
    ## narrow_gctrans_common
  ),
  use.name = FALSE
)

## Create a new set of priors for sim2
sim1_priors <- readRDS(here::here("burnin", "cal", "sim1", "sim1_priors.rds"))
new_priors <- merge(sim1_priors, narrowed, by = "input", all.x = TRUE)

new_priors <- new_priors[
  !is.na(q25) & !is.na(q75),
  ":=" (s1_ll = q25, s1_ul = q75)][, -c("q25", "q75")]

setnames(new_priors, c("s1_ll", "s1_ul"), c("s2_ll", "s2_ul"))

new_priors

################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(
  list(
    vls = simid_sel_vls,
    prep = simid_sel_prep,
    hivp = simid_sel_hivpr,
    hivi = simid_sel_hivinc,
    hivpi = simid_sel_hivprinc
#    gcpos = simid_sel_gcpos
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
