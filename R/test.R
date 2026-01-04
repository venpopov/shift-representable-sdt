library(profvis)
library(here)
source(here("R/readexp.R"))
source(here("R/fitSR.R"))
source(here("R/standardize.R"))


y1 <- readexp(1, "data/dunn2025/")
y2 <- readexp(2, "data/dunn2025/")
y3 <- readexp(3, "data/dunn2025/")

# times with one dataset
# 3.473 sec elapsed
# 2.529 sec elapsed
# 1.784 sec elapsed
# 0.846 sec elapsed
o1 = fitSR(list(y1$total), solveLP_fun = solveLP, simplematrix_fun = simplematrix)
o2 = fitSR(list(y1$total), solveLP_fun = solveLP, simplematrix_fun = simplematrix_fast)
o3 = fitSR(list(y1$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix)
o4 = fitSR(list(y1$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix_fast)

# times with 2 datasets joint analysis
# 55.779 sec elapsed
# 37.268 sec elapsed
# 34.022 sec elapsed
# 14.925 sec elapsed
o1_2 = fitSR(list(y1$total, y2$total), solveLP_fun = solveLP, simplematrix_fun = simplematrix)
o2_2 = fitSR(list(y1$total, y2$total), solveLP_fun = solveLP, simplematrix_fun = simplematrix_fast)
o3_2 = fitSR(list(y1$total, y2$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix)
o4_2 = fitSR(list(y1$total, y2$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix_fast)

# times with 3 datasets joint analysis
# 289.155 sec
# 178.219 sec
# 191.703 sec
# 80.092 sec
o1_3 = fitSR(list(y1$total, y2$total, y3$total), solveLP_fun = solveLP, simplematrix_fun = simplematrix)
o2_3 = fitSR(list(y1$total, y2$total, y3$total), solveLP_fun = solveLP, simplematrix_fun = simplematrix_fast)
o3_3 = fitSR(list(y1$total, y2$total, y3$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix)
o4_3 = fitSR(list(y1$total, y2$total, y3$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix_fast)

# times with 3 datasets joint analysis
o1_3a = fitSR(list(y1$total, y2$total, y3$total), solveLP_fun = solveLP, simplematrix_fun = simplematrix, nstep=100)
o2_3a = fitSR(list(y1$total, y2$total, y3$total), solveLP_fun = solveLP, simplematrix_fun = simplematrix_fast, nstep=100)
o3_3a = fitSR(list(y1$total, y2$total, y3$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix, nstep=100)
o4_3a = fitSR(list(y1$total, y2$total, y3$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix_fast, nstep=100)



profvis({
  out = fitSR(list(y1$total, y2$total, y3$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix_fast, nstep=1000)
})






o4 = fitSR(list(y1$total), solveLP_fun = solveLP_col, simplematrix_fun = simplematrix_fast, nstep=100)
o5 = fitSR(list(y1$total), solveLP_fun = solveLP_col_ind, simplematrix_fun = simplematrix_fast, nstep=100)

profvis({
  out = fitSR(list(y1$total), solveLP_fun = solveLP_col_ind, simplematrix_fun = simplematrix_fast, nstep=1000)
})