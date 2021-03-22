# Calculate selected input parameters

pacman::p_load(data.table)


# HIV TESTING RATE -------------------------------------------------------------

# Source:
# Control, C. f. D., & Prevention, (). HIV Infection Risk, Prevention, and
# Testing Behaviors Among Men Who Have Sex with Men—National HIV Behavioral
# Surveillance, 23 U.S. Cities, 2017.

calc_hivtest_weekrate <- function(x) 1 - (1 - x)^(1 / 52)

# Age groups: 18-24, 25-29, 30-39, 40-49, >= 50
hivtest_past12mo_byage <- c(0.788, 0.823, 0.785, 0.726, 0.639)
hivtest_wk_byage <- calc_hivtest_weekrate(hivtest_past12mo_byage)
