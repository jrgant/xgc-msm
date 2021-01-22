# Calculate the calibration targets against which to compare the simulated
# results.

library(xgcmsm)
library(data.table)

racelabs <- c("B", "H", "O", "W")
anatlabs <- c("rect", "ureth", "phar")


# CALIBRATION TARGETS ----------------------------------------------------------

## Age groups: 13-24, 35-34, 35-44, 45-54, 55+

# SOURCE:
# Centers for Disease Control and Prevention, “Estimated HIV Incidence and
# Prevalence in the United States, 2014–2018,” HIV Surveillance Report,
# vol. 25, no. 1, Art. no. 1, 2020-05, [Online]. Available:
# http://www.cdc.gov/hiv/library/reports/hiv-surveillance.html.
# Table 12. HIV Surveillance Report, Vol. 25, No. 1 -- 2018

base_prop_hivdx <- data.table(
  race = rep(c(racelabs[c(1:2, 4)], "Total"), each = 5),
  age.grp = rep(1:5, 4),
  hiv_prev_N = c(
    20700,  74300,  42200,  41700,  39700,
    10600,  45800,  42200,  42900,  31600,
    5000,   30200,  36100,  67800,  102700,
    38500,  162100, 130700, 164100, 184500
  ),
  hiv_dx_N = c(
    11522,  53772,  36247,  39073,  38071,
    5169,   29838,  34094,  38988,  29934,
    2886,   21307,  30142,  62802,  98819,
    20873,  113483, 109109, 151648, 176821
  )
)

base_prop_hivdx[, sum(hiv_dx_N) / sum(hiv_prev_N), .(race)]

hiv_prev_N.O <-
  base_prop_hivdx[race == "Total", sum(hiv_prev_N), age.grp][, V1] -
  base_prop_hivdx[race != "Total", sum(hiv_prev_N), age.grp][, V1]

hiv_dx_N.O <-
  base_prop_hivdx[race == "Total", sum(hiv_dx_N), age.grp][, V1] -
  base_prop_hivdx[race != "Total", sum(hiv_dx_N), age.grp][, V1]

bp_hivdx <- rbind(
  base_prop_hivdx,
  data.table(
    race        = rep("O", 5),
    age.grp     = 1:5,
    hiv_prev_N  = hiv_prev_N.O,
    hiv_dx_N    = hiv_dx_N.O
  )
)[race != "Total"]

setkeyv(bp_hivdx, c("race", "age.grp"))

ct_hivdx_pr_byrace <- bp_hivdx[, sum(hiv_dx_N) / sum(hiv_prev_N), race][, V1]
names(ct_hivdx_pr_byrace) <- racelabs

ct_hivdx_pr_byage  <- bp_hivdx[, sum(hiv_dx_N) / sum(hiv_prev_N), age.grp][, V1]
names(ct_hivdx_pr_byage) <- paste0("age", 1:5)

# SOURCE:
# S. Singh, A. Mitsch, and B. Wu, "HIV Care Outcomes among Men Who
# Have Sex with Men with Diagnosed HIV Infection - United States, 2015,"
# Morbidity and Mortality Weekly Report, vol. 66, no. 37, Art. no. 37,
# 2017-09, doi: 10.15585/mmwr.mm6637a2. Table 3.

base_prop_vls <- data.table(
  race = rep(racelabs, each = 5),
  age.grp = rep(1:5, 4),
  total_n = c(
    773 + 9381, 27792, 24205, 31274, 16437,
    316 + 3242, 16715, 22581, 24927, 11364,
    88 + 1046,  4264,  5832,  7134,  3896,
    115 + 2107, 13797, 26550, 58134, 46177
  ),
  vls_n = c(
    360 + 4246, 13379,  12712,  17378,  9315,
    191 + 1891, 9637,   13419,  15703,  7186,
    52 + 605,   2642,   3793,   4981,   2740,
    61 + 1279,  8850,   17404,  39492,  31731
  )
)

setkeyv(base_prop_vls, c("race", "age.grp"))

ct_vls_pr_byrace <- base_prop_vls[, sum(vls_n) / sum(total_n), race][, V1]
names(ct_vls_pr_byrace) <- racelabs

ct_vls_pr_byage  <- base_prop_vls[, sum(vls_n) / sum(total_n), age.grp][, V1]
names(ct_vls_pr_byage) <- paste0("age", 1:5)

# HIV: PrEP, HIV incidence, HIV prevalence
ct_prep <- c(0.262, 0.300, 0.398, 0.424)
ct_hiv_incid_100k <- 514 # per 100,000 MSM per year
ct_hiv_prev <- 0.124


# Gonorrhea targets (in STD clinic)

## Tests received by anatomic site
ct_prop_anatsite_tested <- c(0.657, 0826, 0.749)
names(ct_prop_anatsite_tested) <- anatlabs

## Proportion of tested anatomic sites positive
ct_prop_anatsite_pos <- c(0.18, 0.079, 0.129)
names(ct_prop_anatsite_pos) <- anatlabs


# VALIDATION TARGETS -----------------------------------------------------------

## Age distribution among HIV diagnosed
vt_age_among_hivdx <- raceage_grid()[, prop := c(
  0.095, 0.340, 0.193, 0.191, 0.182,
  0.061, 0.265, 0.244, 0.248, 0.183,
  0.047, 0.254, 0.220, 0.252, 0.226,
  0.021, 0.125, 0.149, 0.280, 0.425
)]

vt_hivdx_pr_byraceage <-
  base_prop_hivdx[, sum(hiv_dx_N) / sum(hiv_prev_N), .(race, age.grp)]

vt_vls_pr_byraceage <-
  base_prop_vls[, sum(vls_n) / sum(total_n), .(race, age.grp)]

# Donna Hubbard McCree et al. "Changes in Disparities in Estimated Hiv Incidence
# Rates among Black, Hispanic/Latino, and White Men Who Have Sex with Men (MSM)
# in the United States, 2010-2015". In: Journal of Acquired Immune Deﬁciency
# Syndromes 81.1 (May 2019), pp. 57–62. doi: 10.1097/qai.0000000000001977.

# Black, Hispanic, White
vt_hiv_incid_100k <- c(2459, 1140, 234)


# WRITE CALIBRATION TARGETS ----------------------------------------------------

caltargets <- ls(pattern = "ct")
caltargets_l <- lapply(setNames(caltargets, caltargets), function(x) get(x))

saveRDS(caltargets_l, here::here("est", "caltargets.Rds"))


# WRITE VALIDATION TARGETS -----------------------------------------------------

valtargets <- ls(pattern = "vt")
valtargets_l <- lapply(setNames(valtargets, valtargets), function(x) get(x))

saveRDS(valtargets_l, here::here("est", "valtargets.Rds"))
