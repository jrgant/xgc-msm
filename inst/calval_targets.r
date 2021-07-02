# Calculate the calibration targets against which to compare the simulated
# results.

pacman::p_load(xgcmsm, data.table)

racelabs <- c("B", "H", "O", "W")
anatlabs <- c("rect", "ureth", "phar")

################################################################################
## CALIBRATION TARGETS ##
################################################################################

# HIV DIAGNOSIS TARGETS --------------------------------------------------------

## Age groups: 13-24, 35-34, 35-44, 45-54, 55+

# SOURCE:
# Centers for Disease Control and Prevention, “Estimated HIV Incidence and
# Prevalence in the United States, 2014–2018,” HIV Surveillance Report,
# vol. 25, no. 1, Art. no. 1, 2020-05, [Online]. Available:
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

## ... by race/ethnicity
hivdx_byrace <- bp_hivdx[, .(
  n       = sum(hiv_prev_N),
  hivdx_n = sum(hiv_dx_N),
  prop    = sum(hiv_dx_N) / sum(hiv_prev_N)
), by = race]

ct_hivdx_pr_byrace_dt <- data.table(
  target    = "ct_hivdx_pr_byrace",
  value     = hivdx_byrace[, prop],
  subgroups = racelabs,
  ll95      = NA,
  ul95      = NA
)

## ... by age group
hivdx_byage <- bp_hivdx[, .(
  n       = sum(hiv_prev_N),
  hivdx_n = sum(hiv_dx_N),
  prop    = sum(hiv_dx_N) / sum(hiv_prev_N)
), by = age.grp]

ct_hivdx_pr_byage_dt <- data.table(
  target    = "ct_hivdx_pr_byage",
  value     = hivdx_byage[, prop],
  subgroups = paste0("age", hivdx_byage[, age.grp]),
  ll95      = NA,
  ul95      = NA
)


# VIRAL SUPPRESSION TARGETS ----------------------------------------------------

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

# ... by race/ethnicity
vls_byrace <- base_prop_vls[, .(
  n     = sum(total_n),
  vls_n = sum(vls_n),
  prop  = sum(vls_n) / sum(total_n)
), by = race
  ][, ":=" (
     ll95 = prop - qnorm(0.975) * sqrt(prop * (1 - prop) / n),
     ul95 = prop + qnorm(0.975) * sqrt(prop * (1 - prop) / n)
  )][]

ct_vls_pr_byrace_dt <- data.table(
  target    = "ct_vls_pr_byrace",
  value     = vls_byrace[, prop],
  subgroups = racelabs,
  ll95      = vls_byrace[, ll95],
  ul95      = vls_byrace[, ul95]
)

# ... by age group
vls_byage <- base_prop_vls[, .(
  n     = sum(total_n),
  vls_n = sum(vls_n),
  prop  = sum(vls_n) / sum(total_n)
), by = age.grp
  ][, ":=" (
     ll95 = prop - qnorm(0.975) * sqrt(prop * (1 - prop) / n),
     ul95 = prop + qnorm(0.975) * sqrt(prop * (1 - prop) / n)
  )][]

ct_vls_pr_byage_dt <- data.table(
  target    = "ct_vls_pr_byage",
  value     = vls_byage[, prop],
  subgroups = paste0("age", 1:5),
  ll95      = vls_byage[, ll95],
  ul95      = vls_byage[, ul95]
)

ct_vls_pr_byage_dt


# HIV PrEP COVERAGE TARGETS ----------------------------------------------------

## SOURCE:
## Finlayson et al. Changes in HIV Preexposure Prophylaxis Awareness and Use
## among Men Who Have Sex with Men - 20 Urban Areas, 2014 and 2017. 2019;07:
## 597-603. Morbidity and Mortality Weekly Report.
##
## Values from Table 2 (year 2017)
prep_n  <- c(222,  357,  137,  697)
tot_n   <- c(846, 1191,  344, 1645)
ct_prep <- prep_n / tot_n

cl <- sapply(seq_len(4), function(.x) {
  bt <- binom.test(prep_n[.x], tot_n[.x])
  ll <- bt$conf.int[[1]]
  ul <- bt$conf.int[[2]]
  c(ll, ul)
})

ct_prep_dt <- data.table(
  target    = "ct_prep",
  value     = ct_prep,
  subgroups = racelabs,
  ll95      = cl[1,],
  ul95      = cl[2,]
)


# HIV INCIDENCE & PREVALENCE TARGETS -------------------------------------------

# SOURCE:
# Singh et al. HIV Incidence, HIV Prevalence, and Undiagnosed HIV Infections in
# Men Who Have Sex with Men, United States 2018;03:1-10. Annals of Internal
# Medicine.

## Both parameters from Table 1 (year 2015)
ct_hiv_incid_100k_dt <- data.table(
  target    = "ct_hiv_incid_100k",
  value     = 514, # per 100,000 MSM per year in 2015 (took most recent year)
  subgroups = NA_character_,
  ll95      = 444,
  ul95      = 584
)

# HIV PREVALENCE
ct_hiv_prev_dt <- data.table(
  target    = "ct_hiv_prev",
  value     = 0.124,
  subgroups = NA_character_,
  ll95      = 0.109,
  ul95      = 0.138
)


# GONORRHEA TARGETS (IN STD CLINIC) --------------------------------------------

## SOURCE:
## Abara, Winston E. et al. Extragenital Gonorrhea and Chlamydia
## Positivity and the Potential for Missed Extragenital Gonorrhea with
## Concurrent Urethral Chlamydia among Men Who Have Sex with Men Attending
## Sexually Transmitted Disease Clinics-sexually Transmitted Disease
## Surveillance Network, 2015-2019, 2020-06. Sexually Transmitted Diseases.

## NOTE
## Tests received by anatomic site (see Table 1 in Abara et al.)
## Calculating weighted probability of receiving test given that GC&CT testing
## occured during visit. The analytic sample was limited to visits where
## BOTH GC and CT testing occurred, so this N, if used to calculate the
## probability of receiving a GC test given attending an STI visit,
## would underestimates GC test receipt probability. See NOTE later in this
## section regarding Pr(Receive GC test at anat site | STD visit).
ssun_N_gcanalytic_visits <- 139718

site <- c(
  "Baltimore", "LA County", "Miami", "Boston", "Multnomah", "Minneapolis",
  "New York City", "Philadelphia", "San Francisco", "Seattle"
)

site_n_gctest <- c(
  2214, 7083, 4065, 1066, 8623, 17714,
  46495, 9989, 28448, 14021
)

site_n_rect <- c(
  1315, 4863, 2844, 895, 6279, 12762,
  33842, 6511, 20872, 11183
)

site_n_uret <- c(
  2128, 6992, 3989, 500, 8391, 17361,
  44371, 9767, 27568, 5905
)

# NA = Miami didn't do pharyngeal testing
site_n_phar <- c(
  1718, 6363, NA, 1028, 7802, 16029,
  41576, 9434, 25904, 13472
)

# calculate proportion of visits with a test at anatomic site
site_p_rect <- site_n_rect / site_n_gctest
site_p_uret <- site_n_uret / site_n_gctest
site_p_phar <- site_n_phar / site_n_gctest

# Calculate overall weighted proportion (using inverse variance weighting,
# as in Abara et al.)
ssun_weighted_prop <- function(props, nums) {
  q <- 1 - props
  weighted.mean(props, 1 / (props * q / nums))
}

site_pweighted_rect <- ssun_weighted_prop(site_p_rect, site_n_gctest)
site_pweighted_uret <- ssun_weighted_prop(site_p_uret, site_n_gctest)
site_pweighted_phar <- ssun_weighted_prop(site_p_phar[-3], site_n_gctest[-3])

## NOTE
## Among all STD clinic visits, where
## Y = indicator of receiving GC test at given anatomic site
## X = receiving BOTH GC and CT tests at an STI clinic visit
## Pr(Y = y) = Pr(Y = y|X)Pr(X)
##
## This calculation will underestimate the probability of receiving a
## GC test given visiting the clinic, but it's as close as possible
## with the numbers presented in Abara et al. CIs are are big enough
## to allow some flexibility when choosing runs that fit the targets.
gctests_N_anat <- c(101466, 126972, 123326)
ct_prop_anatsite_tested <-
  gctests_N_anat / ssun_N_gcanalytic_visits * 0.732

ct_prop_anatsite_tested_ll95 <-
  gctests_N_anat / ssun_N_gcanalytic_visits * 0.648

ct_prop_anatsite_tested_ul95 <-
  gctests_N_anat / ssun_N_gcanalytic_visits * 0.817

ct_prop_anatsite_tested_dt <- data.table(
  target    = "ct_prop_anatsite_tested",
  value     = ct_prop_anatsite_tested,
  subgroups = anatlabs,
  ll95      = ct_prop_anatsite_tested_ll95,
  ul95      = ct_prop_anatsite_tested_ul95
)

## Proportion of tested anatomic sites positive
## 95% CLs: [10.4, 13.2], [5.7, 9.3], [7.9, 10.3]
## Weighted prevalence across SSuN sites (see Abara et al.)
ct_prop_anatsite_pos <- c(0.118, 0.075, 0.091)
names(ct_prop_anatsite_pos) <- anatlabs

ct_prop_anatsite_pos_dt <- data.table(
  target    = "ct_prop_anatsite_pos",
  value     = ct_prop_anatsite_pos,
  subgroups = anatlabs,
  ll95      = c(0.104, 0.057, 0.079),
  ul95      = c(0.132, 0.093, 0.103)
)


################################################################################
## VALIDATION CHECKS ##
################################################################################

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
names(vt_hiv_incid_100k) <- racelabs[c(1:2, 4)]


################################################################################
## WRITE FILES ##
################################################################################

# CALIBRATION TARGETS ----------------------------------------------------------

caltargets    <- ls(pattern = "_dt$")
caltargets_l  <- lapply(setNames(caltargets, caltargets), function(x) get(x))
caltargets_dt <- rbindlist(caltargets_l, fill = TRUE)

saveRDS(caltargets_dt, here::here("est", "caltargets.Rds"))


# VALIDATION TARGETS -----------------------------------------------------------

valtargets <- ls(pattern = "vt")
valtargets_l <- lapply(setNames(valtargets, valtargets), function(x) get(x))

saveRDS(valtargets_l, here::here("est", "valtargets.Rds"))
