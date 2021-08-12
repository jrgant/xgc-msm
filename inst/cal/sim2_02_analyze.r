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
pretty_batch_date

################################################################################
## VISUALIZE DISTRIBUTIONS OF TARGETS ACROSS ALL SIMS ##
################################################################################
source(here::here(ic_dir, "sim0.1_fmt_targetdat.r"))
ls()


################################################################################
## SELECT SIMULATIONS ##
################################################################################

simid_sel_hivpr <-
  out_vs_targ[variable == "i.prev" & output_within_5pts == 1, sort(simid)]

simid_sel_hivinc <- out_vs_targ[
  variable == "ir100.pop" & between(output, 0.4, 0.6), sort(simid)]

simid_sel_hivprinc <- intersect(simid_sel_hivpr, simid_sel_hivinc)

simid_sel_vls_reth <-
  lapply(
    setNames(rslugs, rslugs),
    function(.x) {
      out_vs_targ[
        variable == sprintf("cc.vsupp.%s", .x) &
          output_within_5pts == 1, sort(simid)
      ]
    }
  )

sapply(simid_sel_vls_reth, length)
simid_sel_vls_reth_reduce <- Reduce(intersect, simid_sel_vls_reth)

simid_sel_vls_reth_rel <- lapply(
  setNames(rslugs[-4], rslugs[-4]),
  function(.x) {
    out_vs_targ[
      variable == sprintf("cc.vsupp.%s.rel.ref.W", .x) &
        output_within_5pts == 1, sort(simid)
    ]
  }
)

simid_sel_vls_reth_rel_reduce <- Reduce(intersect, simid_sel_vls_reth_rel)

simid_sel_vls_reth_intersect <-
  intersect(simid_sel_vls_reth_reduce, simid_sel_vls_reth_rel_reduce)

simid_sel_dx_age <- lapply(
  setNames(age.rel.dx, age.rel.dx),
  function(.x) {
    out_vs_targ[variable == .x &
                output_within_5pts == 1, sort(simid)]
  }
)

simid_sel_dx_age_reduce <- Reduce(intersect, simid_sel_dx_age)

simid_sel_vdp_intersect <- Reduce(
  intersect,
  list(
    simid_sel_hivprinc,
    simid_sel_vls_reth_intersect,
    simid_sel_dx_age_reduce
  )
)

## ageslugs <- paste0("age", 1:5)

## simid_sel_vls_age <-
##   lapply(
##     setNames(ageslugs, ageslugs),
##     function(.x) {
##       out_vs_targ[
##         variable == sprintf("cc.vsupp.%s", .x) &
##           output_within_5pts == 1, sort(simid)
##       ]
##     }
##   )

## sapply(simid_sel_vls_age, head)
## simid_sel_vls_age345 <- Reduce(intersect, simid_sel_vls_age[1:2])

## simid_sel_vls_age_rel <- lapply(
##   setNames(age.rel.vls, age.rel.vls),
##   function(.x) {
##     out_vs_targ[variable == .x & output_within_5pts == 1, sort(simid)]
##   }
## )

## simid_sel_vls_age_intersect <- Reduce(intersect, simid_sel_vls_age_rel[2:4])

## intersect(simid_sel_vls_reth_intersect, simid_sel_vls_age_intersect)


## simid_sel_prep <-
##   lapply(
##     setNames(rslugs, rslugs),
##     function(.x) {
##       out_vs_targ[
##         variable == sprintf("prepCov.%s", .x) & output_within_5pts == 1, simid]
##     }
##   )

## sapply(simid_sel_prep, length)

## simid_sel_gcpos <-
##   lapply(
##     setNames(gcpos_slugs, gcpos_slugs),
##     function(.x) {
##       out_vs_targ[variable == sprintf("prob.%s.tested", .x)
##                   ][output > 0 & output_within_5pts == 1, simid]
##     }
##   )

## sapply(simid_sel_gcpos, length)

out_vs_targ[simid %in% simid_sel_vdp_intersect
            ][, rel := grepl("rel", variable)][] %>%
  ggplot(aes(x = variable, y = output)) +
  geom_point(position = position_jitter(width = 0.1)) +
  geom_point(
    aes(y = target_val, color = "target"),
    shape = 21,
    size = 7,
    stroke = 2,
    fill = "white"
  ) +
  facet_wrap(~ rel, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.major.x = element_line(color = "gray80"))


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim2.rds"))
sim2_priors <- readRDS(here::here("inst/cal", "sim1_sel_lhs_limits.rds"))

sel_inputs <- pull_params(simid_sel_vdp_intersect)[, -c("selection_group")]


narrowed_vrp_priors <- sel_inputs[, .(
  q25 = quantile(value, 0.25),
  q75 = quantile(value, 0.75)
), by = input][!(input %like% "U2|P2|R2|STOPPER")]

## narrowed_vls_priors <- input_quantiles(vls_sel_inputs, "RX_|TESTER")

## narrowed_prep_inputs <- input_quantiles(prep_sel_inputs, "DISCONT")

## narrowed_gc_priors <- gcpos_sel_inputs[, .(
##   q25 = quantile(value, 0.25),
##   q75 = quantile(value, 0.75)
##   ), .(selection_group, input)
##   ][substring(input, 1, 1) == toupper(substring(selection_group, 1, 1))
##   ][, -c("selection_group")]

## sim2_priors <- readRDS(here::here(ic_dir, "sim1_sel_lhs_limits.rds"))

## new_priors <- merge(
##   sim2_priors,
##   rbind(
##     narrowed_prep_inputs,
##     narrowed_vls_priors,
##     narrowed_gc_priors
##   ),
##   by = "input",
##   all.x = TRUE
## )

## new_priors[!is.na(q25), all(q25 >= s2_ll)]
## new_priors[!is.na(q75), all(q75 <= s2_ul)]

new_priors <- merge(
  sim2_priors,
  narrowed_vrp_priors,
  by = "input",
  all.x = TRUE
)

new_priors[!is.na(q25) & !is.na(q75), ":=" (s2_ll = q25, s2_ul = q75)]

new_priors <- new_priors[
  !is.na(q25) & !is.na(q75),
  ":=" (s2_ll = q25, s2_ul = q75)][, -c("q25", "q75")]

setnames(new_priors, c("s2_ll", "s2_ul"), c("s3_ll", "s3_ul"))

## ## widen HIV transmission prob. scalar ranges for sim3
## new_priors[input %like% "SCALAR_HIV", ":=" (s3_ll = 0.5, s3_ul = 3)][]

## ## reset HIV condom efficacy prior
## new_priors[input %like% "CONDOM_EFF_HIV", ":=" (s3_ll = 0.6, s3_ul = 1)][]


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(
  simid_sel_vdp_intersect,
  here::here("inst", "cal", "sim2_simid_sel.rds")
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
saveRDS(new_priors, here::here("inst", "cal", "sim2_sel_lhs_limits.rds"))
