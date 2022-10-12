#' Elastic integrative analysis for heterogenous treatment effect
#'
#'@description
#' `elasticHTE()` is a test-based dynamic borrowing framework combining
#' a randomized clinical trial (RCT) and a real-world evidence (RWE) study,
#' in which a preliminary test statistic is used to gauge the comparability
#' and reliability of the RWE and decide whether or not to use the RWE in
#' an integrative analysis. The parameter of interest is `psi`, which quantifies
#' how the treatment effect varies over the treatment modifiers.
#'
#'
#' @param dat.t A list contains X, Y and A for RCT.
#' @param dat.os A list contains X, Y and A fro RWE.
#' @param thres.psi A numerical value indicating the threshold for constructing adaptive confidence interval.
#' @param fixed A boolean regarding how to select the tuning parameter `c_gamma`.
#' * `FALSE` (the default): use adaptive selection strategy
#' * `TRUE`: use fixed threshold strategy.
#'   The default fixed threshold is [qchisq(1 -0.05, df = p)],
#'   in which `p` is the dimension of `X` plus one (i.e., add the intercept).
#' @returns A list with components:
#' * est: estimated `psi` associated with the treatment modifiers.
#' * ve: estimated standard error for `psi`.
#' * CIs.inf, CIs.sup: estimated confidence intervals for `psi`.
#' * nuispar: a vector of estimated nuisance parameters
#'
#' @examples
#' m <- 2000 # size of RW dataset
#' # generate the data
#' datasets.all <- GenData(psi0 = c(0,1,1),
#'                              tlocalpar = 2,
#'                              m = m)
#' # datasets.all$RT: the dataset for RT
#' # datasets.all$RW: the dataset for RW
#' # $Y: is a vector of observed outcome
#' # $X: is a matrix of covariates without intercept (n x (p-1))
#' # $A: is a vector of received treatment
#' # $q: is a vector of 1, for permutation-based variance estimation
#' # RT$ps, RT$ml.ps: is a vector of known propensity score (only available in RT)
#'
#' # choice of kappa_n, default is tau * sqrt(log(m)) (similar to BIC criteria)
#' thres.psi  <-  sqrt(log(m)) # threshold for ECI psi
#' # conduct the elastic integrative analysis
#' result.elastic <- elasticHTE(dat.t = datasets.all$RT,
#'                              dat.os = datasets.all$RW,
#'                              fixed = FALSE)
#' result.elastic
#' # return: a list called result.elastic
#' # $est: estimates of the psi in the finite population
#' # $ve: variance estimation of the psi by permutation-based estimation
#' # $CIs.inf, CIs.sup: 95% confidence intervals
#'   ## Wald CIs for regular estimators;
#'   ## ECIs for non-regular estimators, namely elas.v1 and elas.v2
#'   ## elas.v1 performs better in our simulation studies.
#' # $nuispar: estimates for various nuisance parameter
#' @export
elasticHTE <- function(dat.t, # RT
                       dat.os, # RW
                       thres.psi = sqrt(log(m)), # threshold for adaptive CI
                       fixed = FALSE # fixed c_gamma
                       ) {
  n.t <- length(dat.t$Y)
  m <- length(dat.os$Y)
  p <- ncol(dat.t$X) + 1 # add the intercept
  # initialization
  est <- NULL
  ve <- NULL
  nuispar <- NULL
  ##################################################
  ## Step 1: psi estimation
  ## covj.t.: AIPW-adjustement outcome approach
  ## opt.ee.: joint semiparametric efficiency score
  ## ee.rt(ml).: RT only semiparametric efficiency score via sieve models
  ## opt.ee(ml).: joint semiparametric efficiency score via sieve models
  ## S.os1: eta.hat, for detecting assumption violation
  ##################################################

  # for the psi estimation
  psi_est <- function(dat.os = dat.os,
                      dat.t = dat.t) {

    # ---------------------------
    # propensity score estimation
    # OLS for RWE
    glm.out.ols <- glm(A ~ X,
                       family = quasibinomial,
                       weights = q, data = dat.os
    )
    dat.os$ps <- glm.out.ols$fitted.values
    ## construct the ML-x variables
    dat.os$ml.X.os <- cbind(
      dat.os$X[, 1],
      dat.os$X[, 1]^2,
      dat.os$X[, 2],
      dat.os$X[, 2]^2,
      dat.os$X[, 1] * dat.os$X[, 2]
    )

    glm.out.sieve <- glm(A ~ ml.X.os,
                         family = quasibinomial,
                         weights = q, data = dat.os
    )
    dat.os$ml.ps <- glm.out.sieve$fitted.values

    if (any(dat.t$ps != 0.5)) # if not completely RCT
    {
      # OLS for RCT
      glm.out.ols <- glm(A ~ X,
                         family = quasibinomial,
                         weights = q, data = dat.t
      )
      dat.t$ps <- glm.out.ols$fitted.values
      # sieve method
      ## construct the ML-x variables
      dat.t$ml.X.t <- cbind(
        dat.t$X[, 1],
        dat.t$X[, 1]^2,
        dat.t$X[, 2],
        dat.t$X[, 2]^2,
        dat.t$X[, 1] * dat.t$X[, 2]
      )

      glm.out.sieve <- glm(A ~ ml.X.t,
                           family = quasibinomial,
                           weights = q, data = dat.t
      )
      dat.t$ml.ps <- glm.out.sieve$fitted.values
    }
    # ---------------------------

    # ---------------------------
    # outcome model estimation
    # OLS for RCT
    dat.t$ml.X.t <- cbind(
      dat.t$X[, 1],
      dat.t$X[, 1]^2,
      dat.t$X[, 2],
      dat.t$X[, 2]^2,
      dat.t$X[, 1] * dat.t$X[, 2]
    )

    # OLS for mu_0 for RCT and RWE
    mu0.out.t <- glm(Y[which(A == 0)] ~ X[which(A == 0), ],
                     weights = q[which(A == 0)],
                     data = dat.t
    )
    dat.t$mu0 <- cbind(1, dat.t$X) %*% mu0.out.t$coeff

    mu0.out.os <- glm(Y[which(A == 0)] ~ X[which(A == 0), ],
                      weights = q[which(A == 0)],
                      data = dat.os
    )
    dat.os$mu0 <- cbind(1, dat.os$X) %*% mu0.out.os$coeff

    # sieve method for mu for RCT and RWE
    ## RCT
    mu1.out.t <- glm(Y[which(A == 1)] ~ ml.X.t[which(A == 1), ],
                     weights = q[which(A == 1)],
                     data = dat.t
    )
    dat.t$ml.mu1 <- cbind(1, dat.t$ml.X.t) %*%
      mu1.out.t$coeff

    mu0.out.t <- glm(Y[which(A == 0)] ~ ml.X.t[which(A == 0), ],
                     weights = q[which(A == 0)],
                     data = dat.t
    )
    dat.t$ml.mu0 <- cbind(1, dat.t$ml.X.t) %*%
      mu0.out.t$coeff

    ## RWE
    mu1.out.os <- glm(Y[which(A == 1)] ~ ml.X.os[which(A == 1), ],
                      weights = q[which(A == 1)],
                      data = dat.os
    )
    dat.os$ml.mu1 <- cbind(1, dat.os$ml.X.os) %*%
      mu1.out.os$coeff


    mu0.out.os <- glm(Y[which(A == 0)] ~ ml.X.os[which(A == 0), ],
                      weights = q[which(A == 0)],
                      data = dat.os
    )
    dat.os$ml.mu0 <- cbind(1, dat.os$ml.X.os) %*%
      mu0.out.os$coeff

    # obtain the sigma_Y for mu_0 and mu_1 RWE
    dat.os$ml.sigma1 <- summary(mu1.out.os)$dispersion
    dat.os$ml.sigma0 <- summary(mu0.out.os)$dispersion
    # ---------------------------


    # output
    dat.integ <- list(
      Y = c(dat.os$Y, dat.t$Y),
      A = c(dat.os$A, dat.t$A),
      X = rbind(dat.os$X, dat.t$X),
      q = c(dat.os$q, dat.t$q),
      mu0 = c(dat.os$mu0, dat.t$mu0),
      ps = c(dat.os$ps, dat.t$ps),
      ml.mu0 = c(dat.os$ml.mu0, dat.t$ml.mu0),
      ml.ps = c(dat.os$ml.ps, dat.t$ml.ps),
      ml.X = rbind(
        dat.os$ml.X.os,
        dat.t$ml.X.t
      )
    )
    ## covariate adjustment (trial) --> aipw adjusted approach
    Y.tadj <- dat.t$A * (dat.t$Y - dat.t$ml.mu1) / dat.t$ml.ps -
      (1 - dat.t$A) * (dat.t$Y - dat.t$ml.mu0) / (1 - dat.t$ml.ps) + dat.t$ml.mu1 -
      dat.t$ml.mu0
    reg.t <- glm(Y.tadj ~ dat.t$X, weights = dat.t$q)$coeff
    # EFF
    opt.integ <-
      rootSolve::multiroot(f = ee1, start = reg.t, dat = dat.integ)$root
    # RT ml
    opt.ml.t <-
      rootSolve::multiroot(f = ee1.ml, start = reg.t, dat = dat.t)$root
    # RW ml
    opt.ml.integ <-
      rootSolve::multiroot(f = ee1.ml, start = reg.t, dat = dat.integ)$root

    # score function
    S.os1 <- ee1.ml.new(par = opt.ml.t, dat = dat.os)

    # return the psi estimation
    list(
      reg.t = reg.t,
      opt.integ = opt.integ,
      opt.ml.t = opt.ml.t,
      opt.ml.integ = opt.ml.integ,
      S.os1 = S.os1
    )
  }
  psi.list <- psi_est(
    dat.os = dat.os,
    dat.t = dat.t
  )

  est[paste("covj.t.", 1:3, sep = "")] <- psi.list$reg.t
  est[paste("opt.ee.", 1:3, sep = "")] <- psi.list$opt.integ
  est[paste("ee.rt(ml).", 1:3, sep = "")] <- psi.list$opt.ml.t
  est[paste("opt.ee(ml).", 1:3, sep = "")] <- psi.list$opt.ml.integ
  S.os1 <- psi.list$S.os1

  ##################################################
  ## Step 2: asymptotic variance
  ##################################################

  ## permutation-based estimation
  permutation_est <- function(dat.os = dat.os,
                              dat.t = dat.t,
                              nptb = 100) {
    # begin our permutation estimation
    cnames <- c(
      paste("covj.t.", 1:3, sep = ""), paste("opt.ee.", 1:3, sep = ""),
      paste("ee.rt(ml).", 1:3, sep = ""), paste("opt.ee(ml).", 1:3, sep = "")
    )
    ptb.S.os1 <- matrix(0, nptb, p)
    ptb <- matrix(0, nptb, length(cnames))
    colnames(ptb) <- cnames
    for (kkk in 1:nptb) {
      dat.t$q <- rexp(n.t, 1)
      dat.os$q <- rexp(m, 1)
      psi.list.p <- psi_est(
        dat.os = dat.os,
        dat.t = dat.t
      )

      ptb[kkk, paste("covj.t.", 1:3, sep = "")] <- psi.list.p$reg.t
      ptb[kkk, paste("opt.ee.", 1:3, sep = "")] <- psi.list.p$opt.integ
      ptb[kkk, paste("ee.rt(ml).", 1:3, sep = "")] <- psi.list.p$opt.ml.t
      ptb[kkk, paste("opt.ee(ml).", 1:3, sep = "")] <- psi.list.p$opt.ml.integ
      ptb.S.os1[kkk, ] <- psi.list.p$S.os1
    }
    return(list(
      ptb = ptb,
      ptb.S.os1 = ptb.S.os1
    ))
  }

  ## for sanity reason
  flag <- TRUE
  ntime <- 0
  while (flag) {
    ntime <- ntime + 1
    # begin permutation estimation
    list.ptb <- permutation_est(
      dat.os = dat.os,
      dat.t = dat.t,
      nptb = 100
    )
    ptb <- list.ptb$ptb
    ptb.S.os1 <- list.ptb$ptb.S.os1
    # compute Vrt and Veff
    ptboptee1 <- ptb[, paste("ee.rt(ml).", 1:3, sep = "")]
    ptboptee2 <- ptb[, paste("opt.ee(ml).", 1:3, sep = "")]
    Vrt <- var(ptboptee1) * m
    Veff <- var(ptboptee2) * m
    Vrt
    Veff
    Vrteff <- Vrt - Veff
    if (all(eigen(Vrteff)$values > 0)) {
      flag <- FALSE
    }
    if (ntime > 50) {
      break
      cat("wrong \n")
    }
  }

  ## compute the variance of psi
  ve <- apply(ptb, 2, var)

  ##################################################
  ## Step 3: Adaptive selection of gamma (or c_gamma)
  ##################################################

  ## compute the required terms
  rho <- n.t / m
  I.rt <- MASS::ginv(Vrt) / rho
  I.rw <- MASS::ginv(Veff) - MASS::ginv(Vrt)
  Gamma <- MASS::ginv(I.rt) %*% I.rw * rho^(-1 / 2)
  Sigma.S1 <- t(Gamma) %*% I.rt %*% Gamma + I.rw
  iSigSS <- MASS::ginv(Sigma.S1)
  sqrtiSigSS <- expm::sqrtm(iSigSS)
  sqrtVrteff <- Veff %*% expm::sqrtm(Sigma.S1)
  sqrtVeff <- expm::sqrtm(Veff)

  ## test statistics and localparameter for psi
  Tstat1 <- m * t(apply(ptb.S.os1, 2, mean)) %*%
    MASS::ginv(Sigma.S1) %*%
    (apply(ptb.S.os1, 2, mean))
  localpar <- apply(ptb.S.os1, 2, mean) * sqrt(m)
  noncp <- t(localpar) %*% iSigSS %*% (localpar)

  ## simulation for c_gamma selection
  muz1 <- expm::sqrtm(iSigSS) %*% localpar
  muz2 <- expm::sqrtm(Veff) %*% localpar
  ngen <- 1000
  z1.samples <- mvtnorm::rmvnorm(ngen, mean = muz1, diag(3)) # z1
  z2.samples <- mvtnorm::rmvnorm(ngen, mean = muz2, diag(3)) # z2
  z1.samples <- t(z1.samples)
  z2.samples <- t(z2.samples)

  ## construct the asymptotic distribution for RT and RT&RW
  z.rt.samples <- +sqrtVrteff %*% z1.samples - sqrtVeff %*% z2.samples # rej
  z.eff.samples <- -sqrtVeff %*% z2.samples # accept

  ## setup the gamma candidates
  allgamm <- (seq(1 - 1e-10, 1e-10, length.out = 50))
  Lgamma <- length(allgamm)
  critbias <- critvar <- critmse <- matrix(NA, Lgamma, p)
  critmse.tau <- numeric(Lgamma)

  ## begin our selection
  for (kkk in 1:Lgamma) {
    gamm <- allgamm[kkk]
    Icomb <- Tstat1 < qchisq(1 - gamm, df = p)
    id.yes <- apply(z1.samples * z1.samples, 2, sum) < qchisq(1 - gamm, df = p) # integrate or not
    mat.yes <- matrix(id.yes, p, ngen, byrow = TRUE)

    # generate the elastic values
    gen.value <- (z.rt.samples * (1 - mat.yes) +
      z.eff.samples * (mat.yes))
    gen.value <- t(gen.value)

    ## bias
    ## -------------------------------------------
    critbias[kkk, ] <- apply(gen.value, 2, mean, na.rm = TRUE)

    ## variance
    ## -------------------------------------------
    critvar[kkk, ] <- apply(gen.value, 2, var, na.rm = TRUE)

    ## mse
    ## -------------------------------------------
    critmse[kkk, ] <- apply((gen.value - 0)^2, 2, mean, na.rm = TRUE)
  }

  ## chosen by rejective sampling
  kkchosen <- apply(critmse, 2, which.min)
  kkchosen.all <- rep(which.min(apply(critmse, 1, sum)), 3)

  ## select the optimal gamma values (by simulation or analytic form)
  if (fixed) {
    nuispar[paste0("gamma", 1:p)] <- gamm.selected <- rep(0.05, p)
    nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- qchisq(1 - gamm.selected, df = p)
  } else {
    # choose the gamma and c_gamma by minimizing the MSE of the intercept
    nuispar[paste0("gamma", 1:p)] <- gamm.selected <- sapply(
      rep(kkchosen[1], p),
      function(x) allgamm[x])
    nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- sapply(
      rep(kkchosen[1], p),
      function(x) {
        qchisq(1 - allgamm[x], df = p, ncp = 0)})
  }

  ## output the localpar, eta
  nuispar[paste0("eta", 1:p)] <- localpar
  nuispar[paste0("Icomb", 1:p)] <- Icomb <- sapply(cgamma.selected, function(x) x > Tstat1)
  ## combine the SES estimates with the integrated estimator
  est[paste("elas.", 1:p, sep = "")] <- (1 - Icomb) * est[paste("ee.rt(ml).", 1:p, sep = "")] +
    Icomb * est[paste("opt.ee(ml).", 1:p, sep = "")] #- biastemp
  ## expected bias for the elastic estimator
  biastemp <- -pchisq(qchisq(1 - gamm.selected, df = p),
    df = p + 2, ncp = t(muz1) %*% muz1
  ) * Veff %*% localpar
  ## debiased elastic estimator
  # est[paste("elas.",1:p,'.debiased',sep="" )] <-
  #   est[paste("elas.",1:p,sep="" )]  - biastemp/sqrt(m)

  ##################################################
  ## Step 4: elastic integration
  ##################################################

  mat.yes <- sapply(cgamma.selected, function(x) {
    apply(
      z1.samples * z1.samples,
      2, sum
    ) < x
  })
  mat.yes <- t(mat.yes)
  gen.value <- (z.eff.samples * mat.yes + z.rt.samples * (1 - mat.yes))
  gen.value <- t(gen.value)
  ve[paste("elas.", 1:3, sep = "")] <- apply(gen.value, 2, var, na.rm = TRUE) / m

  ##################################################
  ## Step 5: inference
  ##################################################

  ## error allocations for psi
  alpha1.vect <- rep(0.025, p)
  alpha2.vect <- 0.05 - alpha1.vect

  ## upper and lower bounds for mu_1
  UU.muz1 <- sapply(1:p, function(v) {
    alpha2 <- alpha2.vect[v]
    mvtnorm::qmvnorm(1 - alpha2,
      mean = c(muz1), sigma = diag(p),
      tail = "both.tails"
    )$quantile
  })

  LL.muz1 <- -UU.muz1

  ## simulated elastic estimator distribution
  generate_elastic <- function(muz1.new) {
    z0 <- mvtnorm::rmvnorm(ngen, mean = rep(0, 3), diag(3))
    zz0 <- mvtnorm::rmvnorm(ngen, mean = rep(0, 3), diag(3))
    z0 <- t(z0)
    zz0 <- t(zz0)
    z1 <- z0 + c(muz1.new)
    z2 <- zz0 + c(muz2)
    mat.yes <- sapply(cgamma.selected, function(x) apply(z1 * z1, 2, sum) < x)
    mat.yes <- t(mat.yes)
    z.rt <- +sqrtVrteff %*% z1 - sqrtVeff %*% z2 # rej
    z.eff <- -sqrtVeff %*% z2 # accept
    gen.value <- (z.eff * mat.yes + z.rt * (1 - mat.yes))
    gen.value <- t(gen.value)
    gen.value
  }

  ## construct ECI
  nngen <- 100
  qq1 <- qq2 <- matrix(NA, nngen, p)
  Allgen.value <- NULL
  for (nn in 1:nngen) {

    ## combine all the samples
    muz1.new <- mvtnorm::rmvnorm(1, mean = muz1, diag(p))
    Allgen.value <- rbind(
      Allgen.value,
      generate_elastic(muz1.new)
    )

    ## select the 95% CIs for muz1
    while (sum(muz1.new > UU.muz1 | muz1.new < LL.muz1) > 0) {
      muz1.new <- mvtnorm::rmvnorm(1, mean = muz1, diag(p))
    }
    ## re-generate
    gen.value <- generate_elastic(muz1.new)
    qq1[nn, ] <- sapply(1:p, function(v) {
      quantile(gen.value[, v], probs = alpha1.vect[v] / 2, type = 5, na.rm = TRUE)
    })
    qq2[nn, ] <- sapply(1:p, function(v) {
      quantile(gen.value[, v], probs = 1 - alpha1.vect[v] / 2, type = 5, na.rm = TRUE)
    })
  }



  ## construct the Wald CI for the regular estimators
  regular.est.name <- c(
    paste0("covj.t.", 1:3),
    paste0("opt.ee.", 1:3),
    paste0("ee.rt(ml).", 1:3),
    paste0("opt.ee(ml).", 1:3),
    paste0("elas.", 1:3)
  )

  CIs.inf <- CIs.sup <- est[regular.est.name]
  CIs.inf <-
    est[regular.est.name] - qnorm(1 - 0.05 / 2) *
      sqrt(ve[regular.est.name])

  CIs.sup <-
    est[regular.est.name] + qnorm(1 - 0.05 / 2) *
      sqrt(ve[regular.est.name])
  ## begin construct the ECI
  nuispar["conservative"] <- Tstat1 < thres.psi
  nuispar["Tstat.psi"] <- Tstat1

  if (Tstat1 < thres.psi) {

    # simulation
    ## version 1 (CLT)
    CIs.inf[paste("elas.v1.", 1:p, sep = "")] <- est[paste0("elas.", 1:3)] + apply(qq1, 2, min) / sqrt(m)
    CIs.sup[paste("elas.v1.", 1:p, sep = "")] <- est[paste0("elas.", 1:3)] + apply(qq2, 2, max) / sqrt(m)

    ## version 2 (adaptive CI)
    CIs.inf[paste("elas.v2.", 1:p, sep = "")] <-
      est[paste0("elas.", 1:3)] +
      apply(Allgen.value, 2, quantile, probs = 0.05 / 2, type = 5, na.rm = TRUE) / sqrt(m)
    CIs.sup[paste("elas.v2.", 1:p, sep = "")] <-
      est[paste0("elas.", 1:3)] +
      apply(Allgen.value, 2, quantile, probs = (1 - 0.05 / 2), type = 5, na.rm = TRUE) / sqrt(m)
  }
  if (Tstat1 >= thres.psi) {
    gen.value <- generate_elastic(muz1)
    # both two versions are the same
    CIs.inf[paste("elas.v1.", 1:p, sep = "")] <- CIs.inf[paste("elas.v2.", 1:p, sep = "")] <-
      est[paste0("elas.", 1:3)] + apply(gen.value, 2, quantile, probs = 0.025, type = 5, na.rm = TRUE) / sqrt(m)

    CIs.sup[paste("elas.v1.", 1:3, sep = "")] <- CIs.sup[paste("elas.v2.", 1:3, sep = "")] <-
      est[paste0("elas.", 1:3)] + apply(gen.value, 2, quantile, probs = 0.975, type = 5, na.rm = TRUE) / sqrt(m)
  }

  return(list(
    est = est, ve = ve,
    CIs.inf = CIs.inf,
    CIs.sup = CIs.sup,
    nuispar = nuispar
  ))
}
