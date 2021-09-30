################################################################################
## SETUP ##
################################################################################

picksim <- "sim3"

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

ps_ct_hivprev
ps_ct_hivinc
ps_ct_vsupp_race
ps_ct_vsupp_age
ps_ct_prep
ps_ct_gctested + theme(axis.text.x = element_text(angle = 90))
ps_ct_gcpos + theme(axis.text.x = element_text(angle = 90))


################################################################################
## SELECT SIMULATIONS ##
################################################################################

simid_sel_hivpr <-
  out_vs_targ[variable == "i.prev" & output_within_5pts == 1, sort(simid)]

simid_sel_hivinc <- out_vs_targ[
  variable == "ir100.pop" & between(output, 0.444, 0.584), sort(simid)]

simid_sel_hivprinc <- intersect(simid_sel_hivpr, simid_sel_hivinc)
length(simid_sel_hivprinc)

simid_sel_vls_race <- lapply(
  setNames(rslugs, rslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("cc.vsupp.%s", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_vls_race, length)
simid_sel_vls_race_intersect <- Reduce(intersect, simid_sel_vls_race)

simid_sel_dx_race <- lapply(
  setNames(rslugs, rslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("i.prev.dx.inf.%s", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_dx_race, length)
simid_sel_dx_race_intersect <- Reduce(intersect, simid_sel_dx_race)

simid_sel_dx_race_rel <- lapply(
  setNames(rslugs[1:3], rslugs[1:3]),
  function(.x) {
    out_vs_targ[
      variable == sprintf("i.prev.dx.inf.%s.rel.ref.W", .x) & output < 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_dx_race_rel, length)

simid_sel_dx_race_rel_intersect <- Reduce(intersect, simid_sel_dx_race_rel)
length(simid_sel_dx_race_rel_intersect)

ageslugs <- paste0("age", 1:5)

simid_sel_vls_age <- lapply(
  setNames(ageslugs, ageslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("cc.vsupp.%s", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_vls_age, length)
simid_sel_vls_age_intersect <- Reduce(intersect, simid_sel_vls_age)
length(simid_sel_vls_age_intersect)

simid_sel_dx_age <- lapply(
  setNames(ageslugs, ageslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("i.prev.dx.inf.%s", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

simid_sel_dx_age_intersect <- Reduce(intersect, simid_sel_dx_age)
length(simid_sel_dx_age_intersect)

simid_sel_dx_age_rel <- lapply(
  setNames(age.rel.dx, age.rel.dx),
  function(.x) {
    out_vs_targ[
      variable == .x & output < 1,
      sort(simid)
    ]
  }
)

simid_sel_dx_age_rel_intersect <- Reduce(intersect, simid_sel_dx_age_rel)
length(simid_sel_dx_age_rel_intersect)

simid_sel_gctest <- lapply(
  setNames(anatslugs, anatslugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("prop.%s.tested", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_gctest, length)
simid_sel_gctest_intersect <- Reduce(intersect, simid_sel_gctest)
length(simid_sel_gctest_intersect)

simid_sel_gcpos <- lapply(
  setNames(gcpos_slugs, gcpos_slugs),
  function(.x) {
    out_vs_targ[
      variable == sprintf("prob.%s.tested", .x) & output_within_5pts == 1,
      sort(simid)
    ]
  }
)

sapply(simid_sel_gcpos, length)

simid_sel_gcpos_intersect <- Reduce(intersect, simid_sel_gcpos)
length(simid_sel_gcpos_intersect)


plot_targets(simid_sel_gcpos_intersect)
plot_targets(intersect(simid_sel_gcpos_intersect, simid_sel_hivprinc))
plot_targets(intersect(simid_sel_gcpos_intersect, simid_sel_vls_race_intersect))

varpats <- paste0(c(
  "i.prev$",
  "ir100.pop",
  "cc.vsupp.[A-Z]{1}$",
  "dx.inf.[A-Z]{1}$",
  "i.prev.dx.inf.age[1-5]{1}$",
  "prob.[A-Za-z]{3}.tested"
), collapse = "|")

mase <- out_vs_targ[
  variable %like% varpats,
  .(mase = sum(abs_std_diff)),
  simid
]

mase %>%
  ggplot(aes(x = mase)) +
  geom_histogram(color = "white", fill = "gray80")

mase_fitsqt <- lapply(
  c(0.05, 0.01, 0.0025),
  function(.x) {
    mase[mase <= quantile(mase, .x), simid]
})

best_fits10 <- mase[order(mase), simid][1:10]

gridExtra::grid.arrange(
             plot_targets(mase_fitsqt[[1]]) +
             coord_cartesian(ylim = c(0, 1.25)), # smallest 5% of mase

             plot_targets(mase_fitsqt[[2]]) +
             coord_cartesian(ylim = c(0, 1.25)), # smallest 1% of mase

             plot_targets(mase_fitsqt[[3]]) +
             coord_cartesian(ylim = c(0, 1.25)), # smallest 0.25% of mase

             plot_targets(best_fits10) +
             coord_cartesian(ylim = c(0, 1.25)) # smallest 10 mases
           )

mase[simid %in% mase_fitsqt[[1]], summary(mase)]
mase[simid %in% mase_fitsqt[[2]], summary(mase)]
mase[simid %in% mase_fitsqt[[3]], summary(mase)]
mase[simid %in% best_fits10, summary(mase)]


################################################################################
## GET PARAMETER SETS FROM SELECTED SIMULATIONS ##
################################################################################

lhs_groups <- readRDS(here::here("burnin/cal/", picksim, "lhs_sim3.rds"))

lhsdt <- rbindlist(
  lapply(lhs_groups, as.data.table, keep.rownames = TRUE),
  idcol = "simid"
)

lhsdt[, simid := sprintf("%04d", simid)][]
setnames(lhsdt, c("V1", "V2"), c("input", "value"))

input_dist <- function(simids) {
  lhsdt[simid %in% simids] %>%
    ggplot(aes(x = value)) +
    geom_histogram(fill = "gray70", color = "white") +
    geom_density() +
    facet_wrap(~input, scales = "free")
}

input_dist(mase_fitsqt[[1]])
input_dist(mase_fitsqt[[2]])
input_dist(mase_fitsqt[[3]])

cor_sel <- cor(
  dcast(
    lhsdt[simid %in% mase_fitsqt[[2]] &
          !(input %like% "HIV_RX_HALT|SCALAR_(AI|OI)|HIV_TRANS_PROB_WHITE")],
    simid ~ input,
    value.var = "value"
  )[, -1],
  method = "spearman"
)

ggcorrplot(cor_sel, type = "upper")

cordat <- as.data.table(cor_sel, keep.rownames = TRUE)
names(cordat)

cordatl <- melt(cordat, id.vars = "rn", value.name = "spearman")

## create coarse variable groups (cvg) to split up plots for readability
## identify variables using regular expressions

cvg <- list(
  stibehave     = "ACT_STOPPER|RIM|KISS",
  clingc_phar   = "PHAR.*(DURAT|INFPR|GC_SYMPT_PROB)",
  clingc_rect   = "RECT.*(DURAT|INFPR|GC_SYMPT_PROB)",
  clingc_ureth  = "URETH.*(DURAT|INFPR|GC_SYMPT_PROB)",
  stitest_phar  = "^PHAR.*STITEST|PGC_RR",
  stitest_rect  = "^RECT.*STITEST|RGC_RR",
  stitest_ureth = "^URETH.*STITEST|UGC_SYMPT",
  transpath     = "((R|P|U)2)|CONDOM_EFF_GC",
  hivlate       = "HIV_LATE",
  hivart        = "HIV_RX",
  hivtrans      = "HIV_TRANS",
  prep          = "PREP"
)

cortiles <- function(subpattern, numrows = 1) {
  cordatl[rn != variable & rn %like% subpattern] %>%
    ggplot(aes(y = variable, x = factor(1))) +
    geom_tile(
      aes(fill = spearman),
      color = "white",
      size = 0,
      width = 0.75,
      height = 0.9
    ) +
    facet_wrap(~rn, nrow = numrows) +
    scale_fill_scico(name = "Spearman correlation", palette = "vikO") +
    labs(
      x = "",
      caption = "Correlations between input parameters among selected parameter sets."
    ) +
    theme(
      panel.spacing = unit(0, "in"),
      panel.grid.major.y = element_line(color = "black", size = 0.1),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "top",
      legend.key.width = unit(2, "cm"),
      legend.title = element_text(vjust = 1),
      plot.caption = element_text(face = "italic", size = 28)
    )
}


## Save main analysis inputs
main_analysis_inputs <- lapply(
  setNames(mase_fitsqt[[2]], paste0("s3_", mase_fitsqt[[2]])),
  function(.x) {
    lhs_groups[[as.numeric(.x)]]
  })


################################################################################
## WRITE OBJECTS TO FILEs ##
################################################################################

## Write selected simulation ids to file.
saveRDS(
  main_analysis_inputs,
  here::here("inst", "cal", "main_analysis_inputs.rds")
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

## Save correlation plots
plabs <- paste0("calibrated_corr_", names(cvg))

psave(plabs[1],   cortiles(cvg$stibehave),      w = 25, h = 20)
psave(plabs[2],   cortiles(cvg$clingc_phar),    w = 16, h = 20)
psave(plabs[3],   cortiles(cvg$clingc_rect),    w = 16, h = 20)
psave(plabs[4],   cortiles(cvg$clingc_ureth),   w = 16, h = 20)
psave(plabs[5],   cortiles(cvg$stitest_phar),   w = 20, h = 20)
psave(plabs[6],   cortiles(cvg$stitest_rect),   w = 20, h = 20)
psave(plabs[7],   cortiles(cvg$stitest_ureth),  w = 20, h = 20)
psave(plabs[8],   cortiles(cvg$transpath),      w = 30, h = 20)
psave(plabs[9],   cortiles(cvg$hivlate),        w = 25, h = 20)
psave(plabs[10],  cortiles(cvg$hivart),         w = 40, h = 20)
psave(plabs[11],  cortiles(cvg$hivtrans),       w = 35, h = 20)
psave(plabs[12],  cortiles(cvg$prep),           w = 25, h = 20)
