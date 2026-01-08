# Original script by John Dunn, retrieved from https://osf.io/n62y4/files/osfstorage

fitSR <- function(
  data = NULL,
  design = NULL,
  nstep = 20,
  init = NULL,
  input = NULL,
  solveLP_fun = solveLP_col,
  simplematrix_fun = simplematrix_fast
) {
  # takes as input a nr x nc matrix of counts
  # columns correspond to rating categories with confidence that the item is old decreasing from left to right
  # rows correspond to conditions
  # although there are options to change them, usually the call will be:
  # out <- fitSR (y), where y is the data matrix or a list of matrices in which case a joint analysis is performed

  weights <- NULL
  y <- data
  if (!is.list(y)) {
    y <- list(y)
  }

  yy <- 1 - makedatavector(y)

  if (is.null(input)) {
    for (i in seq_along(y)) {
      nr <- nrow(y[[i]])
      nc <- ncol(y[[i]]) - 1
      a <- sdtdesign(c(nr, nc))
      if (is.null(design)) {
        design <- a
      } else {
        design <- magic::adiag(design, a)
      }
    }

    weights <- makeweights(y)
  }

  out <- mLR(
    data = yy,
    design = design,
    weights = weights,
    input = input,
    nstep = nstep,
    solveLP_fun = solveLP_fun,
    simplematrix_fun = simplematrix_fun
  )

  if (length(y) == 1) {
    y <- y[[1]]
  }
  out$y <- y
  out$py <- makedatamatrix(1 - out$predicted, out$y)
  out$g2 <- getgsq(out$py, out$y)
  out$g2[out$g2 < 0] <- 0

  out
}

makedatavector <- function(y) {
  # concatenates a list of data matrices to a single vector of cumulative proportions
  if (!is.list(y)) {
    y <- list(y)
  } # Ensure y is a list

  yy <- numeric(0) # Initialize an empty vector

  for (i in 1:length(y)) {
    cy <- t(apply(y[[i]], 1, cumsum)) # Cumulative sum over columns
    t <- matrix(rep(cy[, ncol(cy)], ncol(cy)), ncol = ncol(cy)) # Replicate the last column
    py <- cy / t # Normalize by dividing by the last column
    py <- py[, -ncol(py)] # Remove the last column

    # Reshape py and append it to yy
    yy <- c(yy, as.vector(t(py)))
  }

  return(yy)
}

makedatamatrix <- function(p, y) {
  # makedatamatrix (p, y)
  # p: numeric vector of predicted values (cumulative probabilities, last missing)
  # y: list (cell array equivalent) of matrices of observed counts

  # Ensure y is a list
  if (!is.list(y)) {
    y <- list(y)
  }

  py <- vector("list", length(y))

  for (i in 1:length(y)) {
    yi <- y[[i]]

    # number of parameters corresponding to y{i}
    n <- length(yi) - nrow(yi)

    # extract and remove those values from p
    z <- p[1:n]
    p <- p[-(1:n)]

    # reshape z to matrix (transpose for MATLAB's column-major behavior)
    z <- matrix(z, nrow = ncol(yi) - 1, ncol = nrow(yi))
    z <- t(z)

    # add a column of ones
    z <- cbind(z, 1)

    # ensure cumulative probabilities don't decrease
    for (j in seq_len(nrow(z))) {
      for (k in 2:ncol(z)) {
        if (z[j, k] < z[j, k - 1]) {
          z[j, k] <- z[j, k - 1]
        }
      }
    }

    # compute row totals and scale
    sy <- rowSums(yi)
    tmat <- matrix(rep(sy, ncol(yi)), nrow = length(sy))

    # differences of cumulative probabilities
    zz <- cbind(z[, 1], t(diff(t(z)))) * tmat

    py[[i]] <- zz
  }

  # simplify if only one matrix
  if (length(py) == 1) {
    py <- py[[1]]
  }

  return(py)
}

makeweights <- function(y) {
  # Ensure input is a list (cell array equivalent)
  if (!is.list(y)) {
    y <- list(y)
  }

  w <- NULL # Initialize empty weight matrix

  for (i in seq_along(y)) {
    cy <- t(apply(y[[i]], 1, cumsum)) # cumulative sum along rows

    for (j in seq_len(nrow(cy))) {
      # Each block is (size(cy,2)-1) x (size(cy,2)-1) identity * cy[j,end]
      block <- diag(ncol(cy) - 1) * cy[j, ncol(cy)]

      # Append block diagonally to existing w
      if (is.null(w)) {
        w <- block
      } else {
        w <- Matrix::bdiag(w, block)
      }
    }
  }

  # Convert sparse Matrix (from bdiag) to standard matrix
  w <- as.matrix(w)

  return(w)
}

getgsq <- function(py, y) {
  if (!is.list(y)) {
    y <- list(y)
  }
  if (!is.list(py)) {
    py <- list(py)
  }

  g2 <- numeric(length(y)) # Initialize the result vector

  for (i in seq_along(y)) {
    # Get the sum of the matrix y[[i]]
    n <- sum(y[[i]], na.rm = TRUE)
    p <- py[[i]]

    # Prevent zeroes in p
    p[p == 0] <- 1 / n
    py[[i]] <- p

    # Compute the z matrix (element-wise log and multiplication)
    z <- 2 * y[[i]] * log(y[[i]] / py[[i]])
    z[y[[i]] == 0] <- 0 # Set z to 0 where y[[i]] is 0

    # Sum up the values of z and assign to g2[i]
    g2[i] <- sum(z, na.rm = TRUE)
  }

  return(g2)
}

sdtdesign <- function(levels, zeroflag = 0) {
  # Returns design matrix for 2-dimensional conjoint SDT model

  a <- factdesign(levels)

  # Remove additive constant
  a <- a[, -1, drop = FALSE]

  # Remove final columns (from levels[1] to end)
  if (levels[1] <= ncol(a)) {
    a <- a[, 1:(levels[1] - 1), drop = FALSE]
  } else {
    a <- a[, 0, drop = FALSE]
  }

  A <- NULL
  for (i in 1:levels[1]) {
    i1 <- (i - 1) * levels[2] + 1
    i2 <- i1 + levels[2] - 1
    block <- cbind(-a[i1:i2, , drop = FALSE], diag(levels[2]))
    A <- rbind(A, block)
  }

  # Adjust for zeroflag
  if (zeroflag == 0) {
    A <- cbind(0, A)
    A[1:levels[2], 1] <- -1
  }

  return(A)
}

factdesign <- function(levels, nrep = 0) {
  # Returns design matrix for factorial combination defined by levels
  # levels: vector of factor levels
  # nrep: number of replications of each row
  # if nrep==0 (default) additive design matrix is returned

  nfact <- length(levels)
  c <- 1:prod(levels)

  # d(i) = product of levels(i+1:end)
  d <- rep(1, nfact)
  if (nfact > 1) {
    for (i in 1:(nfact - 1)) {
      d[i] <- prod(levels[(i + 1):nfact])
    }
  }

  # Compute base index matrix b (combinations of levels)
  b <- matrix(0, nrow = length(c), ncol = nfact)
  for (i in seq_along(c)) {
    j <- i
    for (k in 1:nfact) {
      z <- floor((j - 1) / d[k])
      b[i, k] <- z + 1
      j <- j - z * d[k]
    }
  }

  # Create main effect coding matrices (dummy-coded)
  m <- vector("list", nfact)
  for (k in 1:nfact) {
    z <- diag(levels[k])
    z <- z[, -1, drop = FALSE] # remove first column
    m[[k]] <- z
  }

  # Construct additive design
  A <- NULL
  for (i in 1:nrow(b)) {
    c_row <- NULL
    for (k in 1:ncol(b)) {
      c_row <- cbind(c_row, m[[k]][b[i, k], , drop = FALSE])
    }
    A <- rbind(A, c_row)
  }

  # Add interaction terms if nrep > 0
  if (nrep > 0) {
    for (ifact in 2:nfact) {
      kset <- combn(1:nfact, ifact, simplify = FALSE)
      for (k in kset) {
        b_list <- vector("list", length(k))
        for (i in seq_along(k)) {
          jj <- k[i]
          if (jj == 1) {
            i1 <- 1
          } else {
            i1 <- sum(levels[1:(jj - 1)] - 1) + 1
          }
          i2 <- i1 + levels[jj] - 2
          b_list[[i]] <- A[, i1:i2, drop = FALSE]
        }
        B <- b_list[[1]]
        if (length(b_list) > 1) {
          for (i in 2:length(b_list)) {
            bb <- NULL
            for (ib in 1:ncol(b_list[[i]])) {
              bb <- cbind(bb, B * b_list[[i]][, ib])
            }
            B <- bb
          }
        }
        A <- cbind(A, B)
      }
    }
  }

  # Add intercept
  A <- cbind(1, A)

  # Add replications
  if (nrep > 1) {
    A <- A[rep(1:nrow(A), each = nrep), , drop = FALSE]
  }

  return(A)
}

mLR <- function(
  data = NULL,
  weights = NULL,
  design = NULL,
  nstep = Inf,
  parallel = NULL,
  init = NULL,
  tol = 1e-5,
  input = NULL,
  solveLP_fun = solveLP_col,
  simplematrix_fun = simplematrix_fast
) {
  # solves the monotonic Linear Regression problem
  # data is a n-vector of observations
  # weights is an n x n matrix of weights (default=identity matrix)
  # design is a design matrix (e.g., for main effect, additive model, etc.)
  # (default=intercept model)
  # nstep is number of topes to search (default=exhaustive search)
  # parallel islist of parallelism classes of the difference matrix of design (calculated if NULL)
  # init is optional initial permutation to initiate search
  # tol is tolerance (absolute values less than tol set to zero)
  # input is optional output from a previous call to mLR (allows solution to proceed in stages)
  # returns: optimal fit, fitted values (predictions), and best-fitting permutation,
  # and vectors of fits and best fits for each step
  # checked against matlab code on 23 Oct 2025

  if (!is.null(input)) {
    y <- input$data
    w <- input$weights
    a <- input$design
    init <- input$permutation
    parallel <- input$parallel
    tol <- input$tol
  } else {
    y <- data
    if (is.data.frame(y)) {
      y <- as.vector(y)
    }
    n <- length(y)
    w <- weights %||% diag(rep(1, n))
    a <- design %||% matrix(rep(1, n), n, 1) # default intercept model

    if (is.data.frame(w)) {
      w <- data.matrix(w)
    }

    if (is.vector(w)) {
      w <- diag(w)
    }

    if (is.data.frame(a)) {
      a <- data.matrix(a)
    }
  }

  # initialize some matrices etc.
  itr <- 0 # step counter
  tictoc::tic() # start timer
  fits <- c()
  minfits <- c() # initialize arrays to contain current fit and best fit on each step
  if (is.vector(a)) a <- as.matrix(a)
  da <- mat2diff(a, 1) # difference matrix of a# re-order open by fit# re-order open by fit
  dy <- mat2diff(y) # difference vector of y
  d <- mat2diff(diag(nrow(a))) # difference matrix

  # deal with some special cases
  r <- qr(da)$rank
  if (r == 0) {
    type <- "DA rank 0"
    soltope <- rep(0, nrow(da)) # soltope is the null vector
    mr <- lsqisotonic1(y, w, soltope, d = d)
    pred <- mr$pred
    minfit <- mr$fit
  } else if (r == 1) {
    type <- "DA rank 1"
    k <- which.max(colSums(abs(da)))
    t <- da[, k]
    mr1 <- lsqisotonic1(y, w, t, d = d)
    mr2 <- lsqisotonic1(y, w, -t, d = d)
    if (mr1$fit < mr2$fit) {
      minfit <- mr1$fit
      soltope <- t
      pred <- mr1$pred
    } else {
      minfit <- mr2$fit
      soltope <- -t
      pred <- mr2$pred
    }
  } else if (!is.null(solveLP_fun(dy, da))) {
    type <- "Data conform to permitted permutation"
    soltope <- sign(dy)
    minfit <- 0
    pred <- y
  } else {
    # search topes for a solution
    type <- "Search"
    cond <- condense(a)
    dac <- mat2diff(cond$matrix, 1) # difference matrix of condensation of a
    dc <- mat2diff(diag(nrow(cond$matrix))) # difference matrix of order nrow(ac)

    parallel <- parallel %||% parallelclass(dac, simplematrix_fun)
    rankdc <- sapply(parallel, function(x) qr(dc[x, ])$rank) # store dc sub-matrix ranks
    u <- a %*% MASS::ginv(t(a) %*% w %*% a) %*% t(a) %*% w # calculate useful matrix (???)

    # calculate initial tope if required
    init <- init %||% (u %*% y) # weighted projection of y onto col(a)
    dyp <- mat2diff(init[cond$id])
    init <- as.vector(sign(dyp))

    # repair initial vector if not tope or not feasible
    if (any(init == 0) || is.null(solveLP_fun(init, dac))) {
      md <- mean(abs(dyp))
      my <- y + .1 * runif(length(y), -md, md) # add some random jitter to y
      yp <- u %*% my # weighted projection of my onto col(a)
      dyp <- mat2diff(yp[cond$id])
      init <- as.vector(sign(dyp))
    }

    # warn if initial tope is not in permitted set
    if (is.null(solveLP_fun(init, dac))) {
      warning("Initial tope not in permitted set")
    }

    # set up search
    tc <- init
    pred <- NULL
    t <- decondense(tc, cond)
    mr <- lsqisotonic1(y, w, t, d = d) # decondense initial tope and fit to data

    # precompute expensive empty diff matrix
    t2p_n <- (1 + sqrt(1 + 8 * length(tc))) / 2
    t2perm_d <- mat2diff(diag(t2p_n))
    ptc <- tope2perm(tc, d = t2perm_d) # convert to permutation

    open <- list(list("fit" = mr$fit, "pred" = mr$pred, "perm" = ptc)) # initialize open list
    closed <- new.env(hash = TRUE, parent = emptyenv(), size = 2000L) # store visited permutations to avoid checking

    # initialize optimal fit and solution tope
    if (is.null(input)) {
      minfit <- Inf
      soltope <- init
    } else {
      minfit <- input$fit
      soltope <- sign(mat2diff(input$permutation))
    }

    H <- dac %*% MASS::ginv(t(dac) %*% dac) %*% t(dac) # pre-calculate projection matrix on to col(dac)

    # start search
    while (length(open) > 0 && itr < nstep && minfit > 0) {
      itr <- itr + 1 # increment step count
      tc <- sign(mat2diff(open[[1]]$perm))
      f <- open[[1]]$fit
      yp <- open[[1]]$pred # retrieve values from open
      open <- open[-1] # remove first element from open

      # update if improved fit
      if (f < minfit) {
        minfit <- f
        soltope <- decondense(tc, cond)
        pred <- yp
      }
      fits <- c(fits, f)
      minfits <- c(minfits, minfit)

      # find adjacent topes
      ds <- abs(dc %*% t(dc) %*% tc) / 2 # indicator of adjacent topes in dc

      # precompute q values
      q_base <- as.vector(t(tc) %*% H)

      for (i in 1:length(parallel)) {
        s <- parallel[[i]]
        if (sum(ds[s] == 1) != rankdc[i]) {
          next
        }

        # parallelism class defines a face of the cone of tc in col(d)
        xc <- tc
        xc[s] <- -xc[s] # create candidate adjacent tope
        pxc <- tope2perm(xc, d = t2perm_d)

        # check if candidate already searched
        permutation_key <- paste(pxc, collapse = ",")
        if (exists(permutation_key, envir = closed, inherits = FALSE)) {
          next
        }
        closed[[permutation_key]] <- TRUE

        q_diff <- t(xc[s]) %*% H[s, ]
        q <- as.vector(q_base + 2 * q_diff)
        q[abs(q) < tol] <- 0 # projection test

        # solve the LP problem
        is_adjacent_tope <- all(xc == sign(q)) || !is.null(solveLP_fun(xc, dac))

        if (is_adjacent_tope) {
          x <- decondense(xc, cond)
          mr <- lsqisotonic1(y, w, x, d = d) # calculate fit
          open <- c(
            open,
            list(list(fit = mr$fit, pred = mr$pred, perm = pxc))
          )
        }
      }

      if (length(open) > 0) {
        # re-order open by fit
        f <- sapply(open, function(x) x[[1]])
        open <- open[order(f)]
      }
    }
  }

  # arrange output
  if (length(pred) == 0) {
    mr <- lsqisotonic1(y, w, soltope, d = d)
    pred <- mr$pred
    minfit <- mr$fit
  }

  # append to previous results if input specified
  if (!is.null(input)) {
    itr <- itr + input$nstep
    fits <- c(input$fits, fits)
    minfits <- c(input$minfits, minfits)
  }

  # calculate regression weights and calculate linear estimates
  x <- if (qr(da)$rank > 0) solveLP_fun(soltope, da) else rep(0, ncol(a))
  xp <- if (!is.null(x)) as.vector(a %*% x) else NULL

  list(
    "fit" = minfit,
    "data" = y,
    "predicted" = pred,
    "permutation" = tope2perm(as.vector(soltope)),
    "x" = x,
    "z" = xp,
    "nstep" = itr,
    "mstep" = which(minfits == minfit)[1],
    "fits" = fits,
    "minfits" = minfits,
    "design" = a,
    "weights" = w,
    "parallel" = parallel,
    "tol" = tol,
    "type" = type,
    "init" = init,
    "duration" = tictoc::toc()$callback_msg
  )
}

factdesign <- function(levels = NULL, addflag = 0) {
  # returns design matrix for factorial combination defined by levels,
  # a vector consisting of the number of levels of each factor
  # if addflag is 0 then it returns additive model
  # otherwise it returns all the interaction terms as well
  nfact <- length(levels)
  c <- 1:prod(levels)
  d <- c(rep(1, nfact))
  for (i in 1:nfact - 1) {
    d[i] <- prod(levels[(i + 1):length(levels)])
  }
  b <- matrix(0, length(c), nfact)
  for (i in 1:length(c)) {
    j <- i
    for (k in 1:nfact) {
      z <- floor((j - 1) / d[k])
      b[i, k] <- z + 1
      j <- j - z * d[k]
    }
  }
  m <- vector(mode = "list", length = nfact)
  for (k in 1:nfact) {
    z <- diag(levels[k])
    m[[k]] <- z[, -1]
    if (!is.matrix(m[[k]])) {
      m[[k]] <- matrix(m[[k]], length(m[[k]]), 1)
    }
  }
  A <- matrix()
  A <- A[-1]
  for (i in 1:nrow(b)) {
    aa <- c()
    for (k in 1:ncol(b)) {
      aa <- c(aa, m[[k]][b[i, k], ])
    }
    A <- rbind(A, aa)
  }
  rownames(A) <- NULL # remove annoying row names

  if (addflag != 0) {
    # add interaction terms if asked for
    for (ifact in 2:nfact) {
      k <- combn(nfact, ifact)
      k <- t(k)
      for (j in 1:nrow(k)) {
        b <- vector(mode = "list", length = ncol(k))
        for (i in 1:ncol(k)) {
          jj <- k[j, i]
          if (jj == 1) {
            i1 <- 1
          } else {
            i1 <- sum(levels[1:jj - 1] - 1) + 1
          }
          i2 <- i1 + levels[jj] - 2 # delimits main effect block
          aa <- A[, i1:i2]
          if (is.matrix(aa)) {
            b[[i]] <- aa
          } else {
            b[[i]] <- matrix(aa, length(aa), 1)
          }
        }
        B <- b[[1]]
        for (i in 2:length(b)) {
          bb <- matrix()
          bb <- bb[-1]
          for (ib in 1:ncol(b[[i]])) {
            bb <- cbind(bb, B * b[[i]][, ib])
          }
          B <- bb
        }
        A <- cbind(A, B)
      }
    }
  }
  cbind(c(rep(1, nrow(A))), A) # add intercept
}

mat2diff <- function(A = NULL, flag = 0) {
  # returns difference matrix DA = D*A
  # where A is a matrix - usually a design matrix
  # if flag == 0 then all null columns are removed from DA
  A <- as.matrix(A)
  n <- nrow(A)
  if (n > 1) {
    DA <- matrix(0, choose(n, 2), ncol(A))
    k <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        k <- k + 1
        DA[k, ] <- A[i, ] - A[j, ]
      }
    }
    # remove null columns
    if (flag == 0) {
      if (ncol(DA) > 1) {
        s <- apply(DA != 0, 2, sum)
        s <- which(s == 0)
        if (length(s) > 0) {
          DA <- DA[, -s]
          if (!is.matrix(DA)) {
            DA <- matrix(DA, length(DA), 1)
          }
        }
      }
    }
  } else {
    DA <- matrix(0, 1, length(A))
  }
  return(DA)
}

lsqisotonic1 <- function(y = NULL, w = NULL, t = NULL, d = NULL) {
  # calculates isotonic regression on vector y with weights w
  # t is a sign vector specifying the required order of fitted values
  # d is optional difference matrix
  # output is fitted values ($pred) and weighted least-squares fit ($fit)
  #
  # based on Matlab function lsqisotonic
  # Copyright 2003-2004 The MathWorks, Inc.
  # Revision: 1.1.6.3 Date: 2004/02/01 22:10:40

  n <- length(y)
  if (is.null(w)) {
    w <- diag(rep(1, n))
  }
  if (is.null(d)) {
    d <- mat2diff(diag(n))
  }
  if (!is.null(t)) {
    x <- as.vector(t) %*% d

    # Sort points ascending in x, break ties with y.
    # force secondary approach to ties
    iord <- 1:n
    xy <- cbind(t(x), -y, iord)
    xyord <- xy[order(xy[, 1], xy[, 2]), ]
    ord <- xyord[, 3]
    iord[ord] <- 1:n
    xyord[, 2] <- -xyord[, 2]

    # Initialize fitted values to the given values.
    yhat <- xyord[, 2]
    # convert vector of weights to matrix form or matrix of weights to vector
    # form
    if (is.null(w)) {
      w <- diag(rep(1, length(y)))
    } # default w = identity matrix
    if (is.vector(w)) {
      W <- diag(w)
    } else {
      W <- w
      w <- diag(W)
    }
    block <- 1:n
    w <- w[ord] # reorder w as a column

    # Merge zero-weight points with preceding pos-weighted point (or
    # with the following pos-weighted point if at start).
    posWgts <- (w > 0)
    if (any(!posWgts)) {
      idx <- cumsum(posWgts)
      idx[idx == 0] <- 1
      w <- w[posWgts]
      yhat <- yhat[posWgts]
      block <- idx[block]
    }

    diffs <- diff(yhat)
    while (any(diffs < 0)) {
      # If all blocks are monotonic, then we're done.
      # Otherwise, merge blocks of non-increasing fitted values, and set the
      # fitted value within each block equal to a constant, the weighted mean
      # of values in that block.
      idx <- cumsum(c(1, diffs > 0))
      sumyhat <- tapply(w * yhat, idx, sum)
      w <- tapply(w, idx, sum)
      yhat <- sumyhat / w
      block <- idx[block]
      diffs <- diff(yhat)
    }

    # Broadcast merged blocks out to original points, and put back in original order.
    yhat <- yhat[block]
    yhat <- yhat[iord]
    # calculate fit
    v <- as.vector(y - yhat)
    names(yhat) <- NULL
    fit <- v %*% W %*% v
    fit <- as.numeric(fit)
  } else {
    out <- NULL
  }

  list("pred" = yhat, "fit" = fit)
}

d2t <- function(x = NULL, n = NULL) {
  # Converts a decimal number x into a trinary array of (optional) length n
  # Based on the Matlab function d2b.m
  t <- 3
  if (is.null(n)) {
    n <- ceiling(log(max(x)) / log(t)) + 1
  }
  y <- matrix(0, length(x), n)
  for (j in 1:length(x)) {
    z <- x[j]
    if (z < 0) {
      y[j, n] <- -1
    } else if (z == 0) {
      y[j, n] <- 0
    } else if (z == 1) {
      y[j, n] <- 1
    } else {
      c <- ceiling(log(z) / log(t)) + 1 # Number of divisions necessary ( rounding up the log2(x) )
      yy <- rep(0, c) # Initialize output array
      for (i in 1:c) {
        r <- floor(z / t)
        yy[c + 1 - i] <- z - t * r
        z <- r
      }
      y[j, (n - c + 1):n] <- yy
    }
  }
  if (all(y[, 1] == 0)) {
    y <- y[, -1]
  } # delete leading zeros if not required
  return(y)
}

parallelclass <- function(a, simplematrix_fun = simplematrix_fast) {
  # returns list of parallel classes of matrix a
  out <- simplematrix_fun(a)
  p <- list()
  for (i in 1:length(out$id)) {
    p <- append(p, list(which(out$index == out$id[i])))
  }
  return(p)
}

simplematrix <- function(a) {
  # simplifies matrix a by removing null and parallel rows
  # s is the simple matrix
  # ia is a list of numbers of rows of a included in s
  # ic is a list of numbers of rows of s corresponding to each row of a
  # if any ic==0 then the corresponding row of a is null
  #
  ia <- matrix()
  ia <- ia[-1]
  ic <- rep(0, nrow(a))
  j <- apply(a != 0, 1, sum)
  ic[j == 0] <- -1
  # find parallel rows
  for (i in 1:nrow(a)) {
    if (ic[i] == 0) {
      ic[i] <- i
      ia <- c(ia, i)
      if (i < nrow(a)) {
        for (j in (i + 1):nrow(a)) {
          if (qr(a[c(i, j), ])$rank == 1 && ic[j] == 0) {
            ic[j] <- i
          }
        }
      }
    }
  }
  s <- a[ia, ]
  ic[ic < 0] <- 0
  list("matrix" = s, "id" = ia, "index" = ic)
}

simplematrix_fast <- function(a) {
  n <- nrow(a)
  if (is.null(n) || n == 0L) {
    return(list(
      matrix = a[FALSE, , drop = FALSE],
      id = integer(0),
      index = integer(0)
    ))
  }

  # Make it integer if it is numeric-but-whole (so we can be exact)
  if (!is.integer(a)) {
    if (all(is.finite(a)) && all(a == round(a))) {
      storage.mode(a) <- "integer"
    } else {
      stop(
        "simplematrix_fast: to guarantee identical output, a must be integer (or numeric with all whole numbers)."
      )
    }
  }

  # Zero-row detection (faster than apply)
  nnz <- rowSums(a != 0L)
  nonzero <- which(nnz != 0L)

  ic <- integer(n) # output index; 0 means null row (same as original final output)

  if (length(nonzero) == 0L) {
    return(list(matrix = a[FALSE, , drop = FALSE], id = integer(0), index = ic))
  }

  gcd2 <- function(x, y) {
    # x,y are nonnegative integers
    while (y != 0L) {
      tmp <- x %% y
      x <- y
      y <- tmp
    }
    x
  }

  row_gcd <- function(r) {
    z <- abs(r[r != 0L])
    g <- z[1L]
    if (length(z) > 1L) {
      for (k in 2L:length(z)) {
        g <- gcd2(g, z[k])
      }
    }
    g
  }

  keys <- character(length(nonzero))
  for (k in seq_along(nonzero)) {
    i <- nonzero[k]
    r <- a[i, ]

    g <- row_gcd(r)
    v <- r / g

    # fix sign: make first nonzero positive
    first <- which(v != 0L)[1L]
    if (v[first] < 0L) {
      v <- -v
    }

    keys[k] <- paste(v, collapse = ",")
  }

  # representatives are first occurrences (matching original behavior)
  first_of_key <- !duplicated(keys)
  ia <- nonzero[first_of_key]

  # map every nonzero row to its representative's ORIGINAL row index
  rep_for_key <- ia[match(keys, keys[first_of_key])]
  ic[nonzero] <- rep_for_key

  s <- a[ia, , drop = FALSE]
  list(matrix = s, id = ia, index = ic)
}

condense <- function(a) {
  # returns condensed version of a (redundant rows removed)
  # id is vector of retained rows of a
  # index is vector of rows of a indexed by id
  #
  a <- as.matrix(a)
  ic <- rep(0, nrow(a))
  id <- ic
  k <- 0
  # find identical rows
  for (i in 1:nrow(a)) {
    if (ic[i] == 0) {
      k <- k + 1
      ic[i] <- k
      id[k] <- i
      if (i < nrow(a)) {
        for (j in (i + 1):nrow(a)) {
          if (all(a[i, ] == a[j, ]) & ic[j] == 0) {
            ic[j] <- k
          }
        }
      }
    }
  }
  if (k < length(id)) {
    j <- (k + 1):length(id)
    id <- id[-j]
  }
  ac <- a[id, ]
  list("matrix" = ac, "id" = id, "index" = ic)
}

decondense <- function(x = NULL, cond = NULL) {
  # decondenses a condensed sign difference vector according to cond
  # cond is output from condense(a) where a is the design matrix
  if (is.null(x)) {
    stop("Error: Tope not specified.")
  } 
  ic <- cond$index
  u <- unique(ic)
  if (length(ic) > length(u)) {
    if (length(ic) == 1) {
      y <- rep(x, length(ic))
    } else {
      p <- tope2perm(x)
      p <- p[ic]
      y <- sign(mat2diff(p))
    }
  } else {
    y <- x
  }
  as.vector(y)
}

solveLP <- function(y = NULL, a = NULL) {
  # determines if sign(y) is a covector of matrix a
  # solves LP problem
  # if no solution then returns NULL
  n <- nrow(a)
  m <- ncol(a)
  lprec <- lpSolveAPI::make.lp(n, m)
  lpSolveAPI::set.objfn(lprec, rep(0, m))
  C <- a
  for (i in 1:nrow(C)) {
    vec = C[i, ]
    if (y[i] > 0) {
      lpSolveAPI::add.constraint(lprec, vec, ">=", 1)
    } else if (y[i] < 0) {
      lpSolveAPI::add.constraint(lprec, vec, "<=", -1)
    } else {
      lpSolveAPI::add.constraint(lprec, vec, "=", 0)
    }
  }
  lpSolveAPI::set.bounds(
    lprec,
    lower = rep(-Inf, m),
    upper = rep(Inf, m),
    columns = 1:m
  )
  solve(lprec)
  exitflag <- 1
  if (lpSolveAPI::get.solutioncount(lprec) == 0) {
    exitflag <- 0
  }
  if (exitflag == 1) {
    x <- lpSolveAPI::get.variables(lprec)
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
  lprec <- lpSolveAPI::make.lp(n, m)

  for (i in seq.int(m)) {
    lpSolveAPI::set.column(lprec, i, a[, i])
  }

  lpSolveAPI::set.objfn(lprec, rep(0, m))
  constraint_types <- c("<=", "=", ">=")[sign_y + 2L]
  lpSolveAPI::set.constr.type(lprec, constraint_types)
  lpSolveAPI::set.rhs(lprec, sign_y)

  lpSolveAPI::set.bounds(
    lprec,
    lower = rep(-Inf, m),
    upper = rep(Inf, m),
    columns = 1:m
  )
  solve(lprec) 
  if (lpSolveAPI::get.solutioncount(lprec) == 0) {
    return(NULL)
  }
  lpSolveAPI::get.variables(lprec)
}

tope2perm <- function(x = NULL, d = NULL) {
  # converts tope x to corresponding permutation
  n <- (1 + sqrt(1 + 8 * length(x))) / 2
  if (is.null(d)) {
    d <- mat2diff(diag(n))
  }
  p <- t(d) %*% x
  p[] <- (p + n + 1) / 2
  as.vector(p)
}