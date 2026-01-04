if (require("lpSolveAPI") == F) {
  install.packages("lpSolveAPI")
  library(lpSolveAPI)
}

if (require("bench") == F) {
  install.packages("bench")
}

library(profvis)

solveLP_row <- function(y = NULL, a = NULL) {
  # JD's original function
  # determines if sign(y) is a covector of matrix a
  # solves LP problem
  # if no solution then returns NULL

  n <- nrow(a)
  m <- ncol(a)
  lprec <- make.lp(n, m)
  set.objfn(lprec, rep(0, m))
  C <- a
  for (i in 1:nrow(C)) {
    vec = C[i, ]
    if (y[i] > 0) {
      add.constraint(lprec, vec, ">=", 1)
    } else if (y[i] < 0) {
      add.constraint(lprec, vec, "<=", -1)
    } else {
      add.constraint(lprec, vec, "=", 0)
    }
  }
  set.bounds(lprec, lower = rep(-Inf, m), upper = rep(Inf, m), columns = 1:m)
  solve(lprec)
  exitflag <- 1
  if (get.solutioncount(lprec) == 0) {
    exitflag <- 0
  }
  if (exitflag == 1) {
    x <- get.variables(lprec)
  } else {
    x <- NULL
  }
  return(x)
}


solveLP_col <- function(y = NULL, a = NULL) {
  # determines if sign(y) is a covector of matrix a
  # solves LP problem
  # if no solution then returns NULL

  n <- nrow(a)
  m <- ncol(a)
  sign_y <- sign(y)
  lprec <- make.lp(n, m)

  for (i in seq.int(m)) {
    set.column(lprec, i, a[,i])
  }

  set.objfn(lprec, rep(0, m))
  constraint_types <- c("<=", "=", ">=")[sign_y + 2L]
  set.constr.type(lprec, constraint_types)
  set.rhs(lprec, sign_y)


  set.bounds(lprec, lower = rep(-Inf, m), upper = rep(Inf, m), columns = 1:m)
  solve(lprec)
  exitflag <- 1
  if (get.solutioncount(lprec) == 0) {
    exitflag <- 0
  }
  if (exitflag == 1) {
    x <- get.variables(lprec)
  } else {
    x <- NULL
  }
  return(x)
}

lp_input_1s <- readRDS("data/lp_benchmark_input_1s.rds")
lp_input_3s <- readRDS("data/lp_benchmark_input_3s.rds")

bench::mark(
  solveLP_row(lp_input_1s$y, lp_input_1s$a),
  solveLP_col(lp_input_1s$y, lp_input_1s$a),
  min_iterations = 100
)

bench::mark(
  solveLP_row(lp_input_3s$y, lp_input_3s$a),
  solveLP_col(lp_input_3s$y, lp_input_3s$a),
  min_iterations = 100
)


profvis({
  for (k in seq.int(1000)) {
    y = lp_input_3s$y
    a = lp_input_3s$a
    n <- nrow(a)
    m <- ncol(a)
    sign_y <- sign(y)
    lprec <- make.lp(n, m)

    for (i in seq.int(m)) {
      x = a[,i]
      set.column(lprec, i, x)
    }

    set.objfn(lprec, rep(0, m))
    constraint_types <- c("<=", "=", ">=")[sign_y + 2L]
    set.constr.type(lprec, constraint_types)
    set.rhs(lprec, sign_y)


    set.bounds(lprec, lower = rep(-Inf, m), upper = rep(Inf, m), columns = 1:m)
    solve(lprec)
    exitflag <- 1
    if (get.solutioncount(lprec) == 0) {
      exitflag <- 0
    }
    if (exitflag == 1) {
      x <- get.variables(lprec)
    } else {
      x <- NULL
    }
  }
})


