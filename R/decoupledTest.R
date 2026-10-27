# Decoupled invariance test used by PTC.TestInvariance.
#
# This is the part of seqICP::seqICP.s() (seqICP 1.1, by Niklas Pfister and
# Jonas Peters, GPL-3) that PTC uses: test = "decoupled", model = "iid",
# link = sum and a parametric bootstrap (permutation = FALSE). The computations
# and the random draws are the same as in seqICP.s. The bootstrap residuals are
# computed without the n x n matrix M, which is faster and gives the same p-values
# for the published analyses (checked with identical() on the TCGA and single-cell runs).
#
# Pfister, N., P. Buhlmann and J. Peters (2019). Invariant Causal Prediction for
# Sequential Data. Journal of the American Statistical Association.
# https://doi.org/10.1080/01621459.2018.1491403
#
# X: n x p matrix of candidate predictors, Y: n x 1 response,
# S: columns of X to test, grid: change points that define the environments,
# complements: compare each environment against its complement (TRUE) or
# against all non-overlapping environments (FALSE), B: bootstrap samples.
# Returns the Bonferroni-adjusted p-value of H0: "S is invariant".
decoupledTest <- function(X, Y, S, grid, complements, B = 100) {
  n <- nrow(X)
  Y <- as.matrix(Y)
  # iid model: intercept plus the predictors in S
  X <- cbind(matrix(1, n, 1), X[, S, drop = FALSE])
  d <- ncol(X)

  ###
  # Environments defined by the grid, and the matrices Z to fit each of them
  ###
  L <- length(grid)
  if (L <= 1) {
    stop("not enough grid points given, specify correct grid")
  }
  if (grid[1] == 1) {
    grid[1] <- 0
  }
  L.tot <- choose(L, 2)
  cuts <- combn(L, 2, simplify = TRUE)
  if (complements) {
    # each environment and its complement
    Z <- vector("list", 2 * L.tot)
    env <- vector("list", 2 * L.tot)
    fullenv <- L - 1
    for (i in 1:L.tot) {
      cp1 <- grid[cuts[1, i]]
      cp2 <- grid[cuts[2, i]]
      env[[i]] <- (cp1 + 1):cp2
      if (cp1 == 0 & cp2 != n) {
        env[[i + L.tot]] <- (cp2 + 1):n
      } else if (cp2 == n & cp1 != 0) {
        env[[i + L.tot]] <- 1:cp1
      } else if (cp1 == 0 & cp2 == n) {
        env[[i + L.tot]] <- 1:n
      } else {
        env[[i + L.tot]] <- c(1:cp1, (cp2 + 1):n)
      }
      X.comp <- X[env[[i + L.tot]], , drop = FALSE]
      X.env <- X[env[[i]], , drop = FALSE]
      Z[[i]] <- solve(t(X.env) %*% X.env, t(X.env))
      Z[[i + L.tot]] <- solve(t(X.comp) %*% X.comp, t(X.comp))
    }
    pairwise <- cbind((1:L.tot)[-fullenv], (1:L.tot)[-fullenv] + L.tot)
    L.tot <- 2 * L.tot
  } else {
    # every environment against every non-overlapping environment
    Z <- vector("list", L.tot)
    env <- vector("list", L.tot)
    for (i in 1:L.tot) {
      cp1 <- grid[cuts[1, i]]
      cp2 <- grid[cuts[2, i]]
      env[[i]] <- (cp1 + 1):cp2
      X.env <- X[env[[i]], , drop = FALSE]
      Z[[i]] <- solve(t(X.env) %*% X.env, t(X.env))
    }
    index <- cuts[1, ] >= cuts[2, 1]
    pairwise <- rbind(rep(1, sum(index)), (1:L.tot)[index])
    for (i in 2:L.tot) {
      index <- cuts[1, ] >= cuts[2, i]
      pairwise <- cbind(pairwise, rbind(rep(i, sum(index)), (1:L.tot)[index]))
    }
    pairwise <- t(cbind(pairwise, rbind(pairwise[2, ], pairwise[1, ])))
  }

  ###
  # Test statistics: differences in coefficients (row 1) and in
  # residual variances (row 2) between pairs of environments
  ###
  test.stat.fun <- function(R.scaled) {
    var.env <- vector("numeric", L.tot)
    gamma.env <- vector("list", L.tot)
    for (i in 1:L.tot) {
      R.env <- R.scaled[env[[i]], , drop = FALSE]
      X.env <- X[env[[i]], , drop = FALSE]
      gamma.env[[i]] <- Z[[i]] %*% R.env
      resid <- R.env - X.env %*% gamma.env[[i]]
      var.env[i] <- sum(resid^2) / (length(env[[i]]) - d)
    }
    statfun <- function(i, j) {
      T1 <- sum((gamma.env[[i]] - gamma.env[[j]])^2)
      T2 <- abs(var.env[i] / var.env[j] - 1)
      c(T1, T2)
    }
    T.pairs <- mapply(statfun, pairwise[, 1], pairwise[, 2])
    rbind(sum(T.pairs[1, ]), sum(T.pairs[2, ]))
  }

  # Scaled residuals of the pooled linear fit. M = I - X (X'X)^-1 X' is the
  # residual-maker matrix (n x n)
  XtX.Xt <- solve(t(X) %*% X, t(X))
  M <- diag(n) - X %*% XtX.Xt
  R <- M %*% Y
  R.scaled <- R / sqrt(sum(R^2))
  test.stat <- test.stat.fun(R.scaled)

  ###
  # Parametric bootstrap of the scaled residuals for the null distribution
  ###
  boot.fun <- function(i) {
    # Same as M %*% e, computed without the n x n matrix (much faster for large n)
    e <- matrix(rnorm(n), n, 1)
    tmp <- e - X %*% (XtX.Xt %*% e)
    R.scaled.boot <- tmp / sqrt(sum(tmp^2))
    test.stat.fun(R.scaled.boot)
  }
  test.stat.boot <- matrix(sapply(1:B, boot.fun), 2, B)
  p.value <- (rowSums(test.stat.boot >= matrix(rep(test.stat, B), 2, B)) + 1) / (B + 1)

  # Bonferroni adjustment for the two statistics
  min(2 * min(p.value), 1)
}
