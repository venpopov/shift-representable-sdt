# Original script by John Dunn, retrieved from https://osf.io/n62y4/files/osfstorage
# Minor edits for my workflow by Ven Popov

library(here)
source(here("R/readexp.R"))
source(here("R/fitSR.R"))
source(here("R/standardize.R"))

# combine the three experiments

y1 <- readexp(1, "data/dunn2025/") # acquire response data from Experiment 1
y2 <- readexp(2, "data/dunn2025/") # acquire response data from Experiment 2
y3 <- readexp(3, "data/dunn2025/") # acquire response data from Experiment 3


out <- fitSR(list(y1$total, y2$total, y3$total), nstep = 60)
out <- fitSR(list(y1$total, y2$total, y3$total), nstep = 100, input = out)

# plot changes in fits and minfits
plot(1:out$nstep, out$fits,
  main = "Experiment 1 total", cex.main = 1.5, col = "red", type = "l",
  xlab = "Steps", ylab = "Fit", cex.lab = 1.5
)
points(1:out$nstep, out$minfits, col = "blue", type = "l", lty = 2)
legend("topright", legend = c("fits", "minfits"), col = c("red", "blue"), cex = 1.5, lty = c(1, 2))

out$py # show predicted counts (y1$total)
out$g2 # show G-squared value

cp <- t(matrix(t(1 - out$data), 5, 7))
cp # show data as matrix of cumulative proportions
# plot data (one minus cumulative proportions) against recovered linear component
plot(out$z, out$data,
  pch = 16, main = "Experiment 1 total", cex.main = 1.5,
  xlab = "Recovered Linear Component (z*)", ylab = "One minus Cumulative Proportions", cex.lab = 1.5
)

outm1 <- fitSR(list(y1$mean, y1$participants[[1]]), nstep = 100)
# joint analysis of response data for mean and participant 1

# plot changes in fits and minfits
plot(1:outm1$nstep, outm1$fits,
  main = "Joint analysis: Exp 1 mean + Participant 1", cex.main = 1.5,
  col = "red", type = "l", xlab = "Steps", ylab = "Fit", cex.lab = 1.5
)
points(1:outm1$nstep, outm1$minfits, col = "blue", type = "l", lty = 2)
legend("topright", legend = c("fits", "minfits"), col = c("red", "blue"), cex = 1.5, lty = c(1, 2))

# plot data (one minus cumulative proportions) against recovered linear component
plot(outm1$z[1:35], outm1$data[1:35],
  pch = 16, col = "blue", main = "Joint analysis: Exp 1 mean + Participant 1", cex.main = 1.5,
  xlab = "Joint Recovered Linear Component", ylab = "Data", cex.lab = 1.5
)
points(outm1$z[36:70], outm1$data[36:70], pch = 16, col = "red")
legend("topleft", legend = c("Experiment 1 mean", "Participant 1"), col = c("blue", "red"), pch = 16, cex = 1.5)

outm <- fitSR(y1$mean, nstep = 60)
out1 <- fitSR(y1$participants[[1]], nstep = 25) # separate analyses of mean and participant 1
c(outm$g2, out1$g2)
outm1$g2 # show separate and joint G-squared

c(cor(outm$z, outm1$z[1:35]), cor(out1$z, outm1$z[36:70])) # correlations

sj <- standardize(outm, out1, outm1) # standard joint analysis
# plot standardized recovered linear components against recovered linear components for standard and target
par(mfrow = c(2, 2))
plot(outm$z, outm1$z[1:35],
  pch = 16, col = "blue", main = "Experiment 1 mean", cex.main = 1.5,
  xlab = "Recovered Linear Component (z)", ylab = "Joint z", cex.lab = 1.25
)
plot(out1$z, outm1$z[36:70],
  pch = 16, col = "red", main = "Participant 1", cex.main = 1.5,
  xlab = "Recovered Linear Component (z)", ylab = "Joint z", cex.lab = 1.25
)
plot(outm$z, sj$zs,
  pch = 16, col = "blue", main = "Experiment 1 mean", cex.main = 1.5,
  xlab = "Recovered Linear Component (z)", ylab = "Standardized Joint z", cex.lab = 1.25
)
plot(out1$z, sj$zt,
  pch = 16, col = "red", main = "Participant 1", cex.main = 1.5,
  xlab = "Recovered Linear Component (z)", ylab = "Standardized Joint z", cex.lab = 1.25
)
par(mfrow = c(1, 1))

# plot parameter differences
dm <- sj$xs[2:7] - sj$xs[1]
d1 <- sj$xt[2:7] - sj$xt[1]
plot(c(1, 2, 4), dm[1:3], type = "b", pch = 16, col = "blue", xlab = "No. of presentations")
par(mfrow = c(1, 2))
