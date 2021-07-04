library(data.table)
suppressMessages(library(EpiModelHIV))
library(ggplot2)
library(ggthemes)

arate_dir <- here::here("burnin", "cal", "sim1")
af <- list.files(arate_dir, "Arrive_", full.names = TRUE)
an <- list.files(arate_dir, "Arrive_")

al <- lapply(setNames(af, an), function(.x) {
  as.data.table(readRDS(.x))[at > 3120 - 52][, .(simid, at, num)]
})

ad <- rbindlist(al, idcol = "rate_adj")
ad

adm <- ad[, .(mn_num52 = mean(num)), by = .(rate_adj, simid)]
adm

ggplot(adm, aes(x = mn_num52)) +
  geom_histogram(color = "white") +
  geom_vline(aes(xintercept = 20000), color = "salmon") +
  facet_wrap(~ rate_adj) +
  theme_tufte(base_size = 20)

add <- dcast(simid ~ rate_adj, data = adm, value.var = "mn_num52")
add_sum <- sapply(add[, -c("simid")], summary)

add_sum

ar <- data.table(
  adj = c(1.1, 1.25, 1.5, 1, 2),
  mn_num = add_sum[rownames(add_sum) == "Mean"]
)

fit <- lm(mn_num ~ adj, data = ar)
fit

## solve for N = 20,000
## 20,000 = 17170 + 2203x
## 2203x = 20000 - 17170
## x = (20000 - 17170) / 2203

(20000 - 17170) / 2203

pd <- data.table(adj = seq(1.25, 1.5, 0.01))
pd[, pred := predict(fit, newdata = pd)]

pd
