#' Generate data for integrative analysis
#' @param seed random seed for data generation
#' @param beta0 the coefficients of the linear function for generating Y(0)
#' @param psi0 the coefficients of the linear function for generating treatment effect
#' @param n the size of finite population
#' @param mean.x the mean of X for the finite population
#' @param m the size of the real-world data (RWD)
#' @param tlocalpar the local bias parameter for local alternative
#' @param n.t the size of the randomized clinical trial (RCT)
#' @param equal.weight the indicator for whether to assign the treatment equally or not
#' @param overlap the indicator for assessing the violation of overlap condition
#' @export
GenData <- function(seed = NULL,
                    beta0 = c(0, 1, 1, 1), # for the mu0 function
                    psi0 = c(0, 1, 1), # for contrast function
                    n = 1e5, mean.x = 1, # the whole population
                    m = 2000, tlocalpar = 0, # size and bias for the RWE
                    n.t = NULL, equal.weight = T, # size for the RCT (NULL uses the default)
                    overlap = T) {
  if(!is.null())set.seed(seed)
  if (overlap == T) {
    X1 <- rnorm(n, mean.x, 1)
    X2 <- rnorm(n, mean.x, 1)
    X3 <- rnorm(n, mean.x, 1) #+ sin(X1)-sin(X2)
    X1.temp <- rnorm(n, mean.x, 1)
    X2.temp <- rnorm(n, mean.x, 1)
    X3.temp <- rnorm(n, mean.x, 1) #+sin(X1.temp) - sin(X2.temp)
    ## construct functions
    h.X <- cbind(X1, X2, X3)
    h.Xtemp <- cbind(X1.temp, X2.temp, X3.temp)
    f.X <- cbind(X1, X2)
    f.Xtemp <- cbind(X1.temp, X2.temp)

    ## construct Y for RCT and RWE
    ### RCT
    Y0 <- cbind(1, h.X) %*% beta0 + rnorm(n, 0, 1)
    Y1 <- cbind(1, h.X) %*% beta0 + cbind(1, f.X) %*% psi0 + rnorm(n, 0, 1)
    ### RWE
    Y0temp <- cbind(1, h.Xtemp) %*% beta0 + rnorm(n, 0, 1)
    Y1temp <- cbind(1, h.Xtemp) %*% beta0 + cbind(1, f.Xtemp) %*% psi0 + rnorm(n, 0, 1)
  } else {
    # change mean to 0
    X1 <- rnorm(n, 0, 1)
    X2 <- rnorm(n, 0, 1)
    X3 <- rnorm(n, 0, 1) #+ sin(X1)-sin(X2)
    X1.temp <- rnorm(n, 0, 1)
    X2.temp <- rnorm(n, 0, 1)
    X3.temp <- rnorm(n, 0, 1) #+sin(X1.temp) - sin(X2.temp)
    ## construct functions
    h.X <- cbind(X1, X2, X3)
    h.Xtemp <- cbind(X1.temp, X2.temp, X3.temp)
    f.X <- cbind(X1, X2)
    f.Xtemp <- cbind(X1.temp, X2.temp)

    ## construct Y for RCT and RWE
    ### RCT # beta0 = (0, 1, 1, 0)
    Y0 <- cbind(1, h.X) %*% beta0 + rnorm(n, 0, 1)
    Y1 <- cbind(1, h.X) %*% beta0 + cbind(1, X1, abs(X2)) %*% psi0 + rnorm(n, 0, 1)
    ### RWE
    Y0temp <- cbind(1, h.Xtemp) %*% beta0 + rnorm(n, 0, 1)
    Y1temp <- cbind(1, h.Xtemp) %*% beta0 + cbind(1, X1.temp, abs(X2.temp)) %*% psi0 + rnorm(n, 0, 1)
  }
  # true ATE
  psi.pop <- lm((Y1 - Y0) ~ f.X)$coeff
  true.tau <- mean(Y1) - mean(Y0)
  taui <- Y1 - Y0

  ## RCT data
  if (is.null(n.t)) {
    eS <- expoit(-4.5 - 2 * X1 - 2 * X2)
  } else {
    e.opt <- uniroot(function(e) {
      eS <- expoit(e - 2 * X1 - 2 * X2)
      mean(eS) * n - n.t
    }, c(-100, 100))$root

    eS <- expoit(e.opt - 2 * X1 - 2 * X2)
  }
  # eS <- expoit(-5.5 -2 * X1 -2 *  X2)
  mean(eS) * n
  S <- sapply(eS, rbinom, n = 1, size = 1)
  S.ind <- which(S == 1)
  n.t <- length(S.ind)
  ### choose the selected patients
  X.t <- f.X[S.ind, ]
  Y0.t <- Y0[S.ind]
  Y1.t <- Y1[S.ind]
  eS.t <- eS[S.ind]
  if (overlap == F) {
    keep <- which(X.t[, 2] > 0) # change to larger than 0
    X.t <- X.t[keep, ]
    Y0.t <- Y0.t[keep]
    Y1.t <- Y1.t[keep]
    n.t <- length(Y1.t)
  }
  # randomly assign treatment to trial participant
  if (all(equal.weight == T)) {
    ps.t <- rep(0.5, n.t)
    A.t <- rbinom(n.t, 1, ps.t)
  } else {
    # equal.weight <- c(-3,-3)
    # equal.weight <- c(-3, -3)
    ps.t <- expoit(-2 + X.t %*% equal.weight)
    hist(ps.t)
    A.t <- rbinom(n.t, 1, ps.t)
  }
  Y.t <- Y1.t * A.t + Y0.t * (1 - A.t)

  ## OS data
  P.ind <- sample(1:n, size = m)
  X.os <- f.Xtemp[P.ind, ]
  X.confounder.os <- X3.temp[P.ind]
  Y0.os <- Y0temp[P.ind]
  Y1.os <- Y1temp[P.ind]

  # mantain at about .5
  a.opt <- uniroot(function(a) {
    ps <- expoit(a - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)
    mean(ps) - .5
  }, c(-100, 100))$root
  eA <- expoit(a.opt - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)
  # eA <-   expoit(-1 - 1*X.os[,1] - tlocalpar*X.confounder.os )
  mean(eA)
  # hist(eA)
  A.os <- sapply(eA, rbinom, n = 1, size = 1)
  Y.os <- Y1.os * A.os + Y0.os * (1 - A.os)

  # return the RT and RW data-sets
  dat.t <- list(
    Y = Y.t, A = A.t, X = X.t, q = rep(1, n.t),
    ps = ps.t,
    ml.ps = ps.t
  )
  dat.os <- list(Y = Y.os, A = A.os, X = X.os, q = rep(1, m))
  list(
    RT = dat.t,
    RW = dat.os
  )
}

# score functions
ee1 <- function(par, dat) {
  psi <- par
  Y <- dat$Y
  A <- dat$A
  X <- dat$X
  q <- dat$q
  ps <- dat$ps
  mu0 <- dat$mu0
  H <- Y - A * cbind(1, X) %*% psi
  apply(cbind(1, X) * matrix((H - mu0) * (A - ps) * q, length(Y), 1 + dim(X)[2], byrow = FALSE), 2, mean)
}

# score functions divided for sieve approximation
ee1.ml <- function(par, dat) {
  psi <- par
  Y <- dat$Y
  A <- dat$A
  X <- dat$X
  q <- dat$q
  ps <- dat$ml.ps
  mu0 <- dat$ml.mu0
  H <- Y - A * cbind(1, X) %*% psi
  apply(cbind(1, X) * matrix((H - mu0) * (A - ps) * q, length(Y), 1 + dim(X)[2], byrow = FALSE), 2, mean)
}

# score functions divided by sigma for sieve approximation
ee1.ml.new <- function(par, dat) {
  psi <- par
  Y <- dat$Y
  A <- dat$A
  X <- dat$X
  q <- dat$q
  ps <- dat$ml.ps
  mu0 <- dat$ml.mu0
  sigma0 <- dat$ml.sigma0
  sigma1 <- dat$ml.sigma1
  H <- Y - A * cbind(1, X) %*% psi
  apply(sigma0^(-1) * cbind(1, X) * matrix((H - mu0) * (A - ps) * q, length(Y), 1 + dim(X)[2],
    byrow = FALSE
  ), 2, mean)
}

# handy functions
expoit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}
logit <- function(x) {
  return(log(x / (1 - x)))
}
