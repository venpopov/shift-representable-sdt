library(here)
source(here("R/readexp.R"))
source(here("R/fitSR.R"))
source(here("R/standardize.R"))

set.seed(124)

y1 <- readexp(1, "data/dunn2025/")



N <- 1800


gen_sdt_conf_counts <- function(N, dprime, sd, crits) {
  ps <- rev(1 - pnorm(crits, dprime, sd))
  props <- c(ps[1], diff(ps), 1 - ps[length(ps)])
  choices <- sample(seq_along(props), size = N, prob = props, replace = T)
  choices <- factor(choices, levels = seq_along(props))
  table(choices)
}

ds <- rev((0:6) / 4)
dat <- t(sapply(ds, function(d) gen_sdt_conf_counts(N, dprime = d, sd = d / max(ds) / 1.5 + 1, crits = c((-1):3) - 0.25)))

out <- fitSR(dat, nstep = 100)
str(out, 1)

plot(out$z, out$data, pch = rep(c(1, 17), each = 35))


out0 <- fitSR(list(y1$total), nstep = 100)
str(out0, 1)

out2 <- fitSR(list(y1$total, dat), nstep = 400)
str(out2, 1)
plot(out2$z, out2$data, pch = rep(c(1, 17), each = 35))


par(mfrow = c(1, 2))
plot(out0$z, out2$z[1:length(out0$z)], xlab = "Exp 1 recovered linear component", ylab = "Exp 1 estimate from joint analysis")
plot(out$z, out2$z[(length(out$z) + 1):(2 * length(out$z))], xlab = "UVSD recovered linear component", ylab = "UVSD estimate from joint analysis")
