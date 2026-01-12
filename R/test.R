library(profvis)
library(here)
source(here("R/readexp.R"))
source(here("R/fitSR.R"))
source(here("R/standardize.R"))


y1 <- readexp(1, "data/dunn2025/")
y2 <- readexp(2, "data/dunn2025/")
y3 <- readexp(3, "data/dunn2025/")

# toy <- y1$total[c(3, 6), c(3, 5, 6)]
# out <- fitSR(list(toy))
# str(out,1)
#
# o1 <- fitSR(list(y1$total))
# o2 <- fitSR(list(y1$total, y2$total))
# o3 <- fitSR(list(y1$total, y2$total, y3$total))
#
#
# o1a <- fitSR(list(y1$total), lp_solver = "highs")
# o2a <- fitSR(list(y1$total, y2$total), lp_solver = "highs")
# o3a <- fitSR(list(y1$total, y2$total, y3$total), lp_solver = "highs")


o3a <- fitSR(list(y1$total, y2$total, y3$total), lp_solver = "highs", nstep = 100)
o3b <- fitSR(list(y1$total, y2$total, y3$total), lp_solver = "lpSolveAPI", nstep = 100)

str(o3a, 1)
str(o3b, 1)
