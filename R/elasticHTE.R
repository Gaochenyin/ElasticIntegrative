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
#' @param dat.os A list contains X, Y and A for RWE.
#' @param thres.psi A numerical value indicating the threshold for constructing adaptive confidence interval.
#' @param fixed A boolean regarding how to select the tuning parameter `c_gamma`.
#' * `FALSE` (the default): use adaptive selection strategy
#' * `TRUE`: use fixed threshold strategy.
#'   The default fixed threshold is [qchisq(1 -0.05, df = p)],
#'   in which `p` is the dimension of `X` plus one (i.e., add the intercept).
#' @param alpha a vector for estimating the linear combination of the coefficients for treatment modifiers (depreciated).
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
#'   ## ACIs for non-regular estimators, namely elas.v1 and elas.v2
#'   ## elas.v1 performs better in our simulation studies.
#' # $nuispar: estimates for various nuisance parameter
#' @export
elasticHTE <- function(dat.t, # RT
                       dat.os, # RW
                       thres.psi = sqrt(log(m)), # threshold for adaptive CI
                       fixed = FALSE, # fixed c_gamma
                       alpha = c(0, 0, 0)
                       )
{

  n.t <- length(dat.t$Y)
  m <- length(dat.os$Y)
  p <- ncol(dat.t$X)+1
  # initialization
  est <-NULL
  ve  <-NULL
  est.tau<-NULL
  ve.tau  <-NULL
  nuispar<-NULL

  # psi estimation
  psi_est <- function(dat.os = dat.os,
                      dat.t = dat.t)
  {
    # for the psi estimation
    # OLS
    glm.out.ols <- glm(A~X,
                       family=quasibinomial,
                       weights=q,data=dat.os)
    dat.os$ps <- glm.out.ols$fitted.values
    # sieve method
    ## construct the ML-x variables
    dat.os$ml.X.os <- cbind(dat.os$X[,1],
                            dat.os$X[,1]^2,
                            dat.os$X[,2],
                            dat.os$X[,2]^2,
                            dat.os$X[,1]*dat.os$X[,2])

    glm.out.sieve <- glm(A~ml.X.os,
                         family=quasibinomial,
                         weights=q,data=dat.os)
    dat.os$ml.ps <- glm.out.sieve$fitted.values

    if(any(dat.t$ps!=0.5)) # if not eqaully weightted
    {
      # OLS
      glm.out.ols <- glm(A~X,
                         family=quasibinomial,
                         weights=q,data=dat.t)
      dat.t$ps <- glm.out.ols$fitted.values
      # sieve method
      ## construct the ML-x variables
      dat.t$ml.X.t <- cbind(dat.t$X[,1],
                            dat.t$X[,1]^2,
                            dat.t$X[,2],
                            dat.t$X[,2]^2,
                            dat.t$X[,1]*dat.t$X[,2])

      glm.out.sieve <- glm(A~ml.X.t,
                           family=quasibinomial,
                           weights=q,data=dat.t)
      dat.t$ml.ps <- glm.out.sieve$fitted.values
    }
    # OLS
    dat.t$ml.X.t <- cbind(dat.t$X[,1],
                          dat.t$X[,1]^2,
                          dat.t$X[,2],
                          dat.t$X[,2]^2,
                          dat.t$X[,1]*dat.t$X[,2])

    # OLS for mu_0 for RCT and RWE
    mu0.out.t <- glm(Y[which(A==0)]~X[which(A==0),],
                     weights=q[which(A==0)],
                     data=dat.t)
    dat.t$mu0 <- cbind(1,dat.t$X)%*%mu0.out.t$coeff

    mu0.out.os <- glm(Y[which(A==0)]~X[which(A==0),],
                      weights=q[which(A==0)],
                      data=dat.os)
    dat.os$mu0 <- cbind(1,dat.os$X)%*%mu0.out.os$coeff

    # sieve method for mu for RCT and RWE
    ## RCT
    mu1.out.t <- glm(Y[which(A==1)]~ml.X.t[which(A==1),],
                     weights=q[which(A==1)],
                     data=dat.t)
    dat.t$ml.mu1 <- cbind(1,dat.t$ml.X.t)%*%
      mu1.out.t$coeff

    mu0.out.t <- glm(Y[which(A==0)]~ml.X.t[which(A==0),],
                     weights=q[which(A==0)],
                     data=dat.t)
    dat.t$ml.mu0 <- cbind(1,dat.t$ml.X.t)%*%
      mu0.out.t$coeff

    ## RWE
    mu1.out.os <- glm(Y[which(A==1)]~ml.X.os[which(A==1),],
                      weights=q[which(A==1)],
                      data=dat.os)
    dat.os$ml.mu1 <- cbind(1,dat.os$ml.X.os)%*%
      mu1.out.os$coeff


    mu0.out.os <- glm(Y[which(A==0)]~ml.X.os[which(A==0),],
                      weights=q[which(A==0)],
                      data=dat.os)
    dat.os$ml.mu0 <- cbind(1,dat.os$ml.X.os)%*%
      mu0.out.os$coeff

    # obtain the sigma_Y for mu_0 and mu_1 RWE
    dat.os$ml.sigma1 <- summary(mu1.out.os)$dispersion
    dat.os$ml.sigma0 <- summary(mu0.out.os)$dispersion
    dat.integ<-list( Y = c(dat.os$Y, dat.t$Y),
                     A = c(dat.os$A, dat.t$A),
                     X = rbind(dat.os$X, dat.t$X),
                     q = c(dat.os$q, dat.t$q),
                     mu0 = c(dat.os$mu0, dat.t$mu0),
                     ps = c(dat.os$ps, dat.t$ps),
                     ml.mu0 = c(dat.os$ml.mu0, dat.t$ml.mu0),
                     ml.ps = c(dat.os$ml.ps, dat.t$ml.ps),
                     ml.X = rbind(dat.os$ml.X.os,
                                  dat.t$ml.X.t))
    ## covariate adjustment (trial) --> aipw adjusted approach
    Y.tadj<- dat.t$A*(dat.t$Y-dat.t$ml.mu1)/dat.t$ml.ps-
      (1-dat.t$A)*(dat.t$Y-dat.t$ml.mu0)/(1-dat.t$ml.ps)+dat.t$ml.mu1-
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
    S.os1 <- ee1.ml.new(par=opt.ml.t, dat=dat.os)

    # return the psi estimation
    list(reg.t = reg.t,
         opt.integ = opt.integ,
         opt.ml.t = opt.ml.t,
         opt.ml.integ = opt.ml.integ,
         S.os1 = S.os1)
  }
  psi.list <- psi_est(dat.os = dat.os,
                      dat.t = dat.t)

  est[paste("covj.t.",1:3,sep="" )] <- psi.list$reg.t
  est[paste("opt.ee.",1:3,sep="" )] <- psi.list$opt.integ
  est[paste("ee.rt(ml).",1:3,sep="" )]<- psi.list$opt.ml.t
  est[paste("opt.ee(ml).",1:3,sep="" )]<- psi.list$opt.ml.integ
  S.os1 <- psi.list$S.os1

  # permutation-based estimation
  permutation_est <- function(dat.os = dat.os,
                              dat.t = dat.t,
                              nptb = 100)
  {
    n.t <- length(dat.t$Y); m <- length(dat.os$Y)
    # begin our permutation estimation
    cnames<-c(paste("covj.t.",1:3,sep="" ),   paste("opt.ee.",1:3,sep="" ),
              paste("ee.rt(ml).",1:3,sep="" ),paste("opt.ee(ml).",1:3,sep="" ))
    ptb.S.os1<-matrix(0, nptb, p)
    ptb <- matrix(0,nptb,length(cnames))
    colnames(ptb)<-cnames
    for(kkk in 1:nptb){
      dat.t$q <-rexp(n.t,1)
      dat.os$q<-rexp(m,1)
      psi.list.p <- psi_est(dat.os = dat.os,
                            dat.t = dat.t)

      ptb[kkk,paste("covj.t.",1:3,sep="" )] <- psi.list.p$reg.t
      ptb[kkk,paste("opt.ee.",1:3,sep="" )] <- psi.list.p$opt.integ
      ptb[kkk,paste("ee.rt(ml).",1:3,sep="" )] <- psi.list.p$opt.ml.t
      ptb[kkk,paste("opt.ee(ml).",1:3,sep="" )] <- psi.list.p$opt.ml.integ
      ptb.S.os1[kkk,] <- psi.list.p$S.os1
    }
    return(list(ptb = ptb,
                ptb.S.os1 = ptb.S.os1))
  }

  flag <- TRUE
  ntime <- 0
  while(flag)
  {
    ntime <- ntime+1
    # begin permutation estimation
    list.ptb <- permutation_est(dat.os = dat.os,
                                dat.t = dat.t,
                                nptb = 100)
    ptb <- list.ptb$ptb
    ptb.S.os1 <- list.ptb$ptb.S.os1
    # compute Vrt and Veff
    ptboptee1<-ptb[,paste("ee.rt(ml).",1:3,sep="" )]
    ptboptee2<-ptb[,paste("opt.ee(ml).",1:3,sep="" )]
    Vrt <- var(ptboptee1)*m; Veff <- var(ptboptee2)*m
    Vrt;Veff
    Vrteff <- Vrt-Veff
    if(all(eigen(Vrteff)$values>0)){flag <- FALSE;}
    if(ntime>50){break; cat('wrong \n')}
  }

  # compute the variance of psi
  ve <- apply(ptb,2,var)

  # compute sigma.S by formula
  rho <- n.t/m
  I.rt <- MASS::ginv(Vrt)/rho
  I.rw <- MASS::ginv(Veff) - MASS::ginv(Vrt)
  Gamma <-  MASS::ginv(I.rt)%*%I.rw*rho^(-1/2)
  Sigma.S1 <- t(Gamma)%*%I.rt%*%Gamma + I.rw
  iSigSS<- MASS::ginv(Sigma.S1)
  sqrtiSigSS<- expm::sqrtm( iSigSS )
  sqrtVrteff <- Veff%*%expm::sqrtm(Sigma.S1)
  sqrtVeff <- expm::sqrtm(Veff)
  # test statistics and localparameter for psi
  Tstat1<- m*t(apply(ptb.S.os1,2,mean))%*%
    MASS::ginv(Sigma.S1)%*%
    (apply(ptb.S.os1,2,mean))
  localpar<-apply(ptb.S.os1,2,mean)*sqrt(m)
  noncp <-t(localpar)%*%iSigSS%*%(localpar)

  # test statistics for tau
  Tstat.tau <- (alpha%*%localpar)*
    MASS::ginv(alpha%*%Sigma.S1%*%alpha)*
    (alpha%*%localpar)


  # simulation for c_gamma selection
  muz1 <- expm::sqrtm(iSigSS)%*%localpar
  muz2 <- expm::sqrtm(Veff)%*%localpar
  ngen <- 1000
  z1.samples <- mvtnorm::rmvnorm(ngen, mean=muz1, diag(3)) # z1
  z2.samples <- mvtnorm::rmvnorm(ngen, mean=muz2, diag(3)) # z2
  z1.samples <-t(z1.samples);z2.samples <- t(z2.samples)

  # a sanity check
  mu.RT <- sqrtVrteff%*%muz1 - sqrtVeff%*%muz2; var.RT <- sqrtVrteff%*%diag(3)%*%t(sqrtVrteff)+
    sqrtVeff%*%diag(3)%*%t(sqrtVeff)
  mu.eff <- -sqrtVeff%*%muz2; var.eff <- sqrtVeff%*%diag(3)%*%t(sqrtVeff)


  z.rt.samples <-  + sqrtVrteff%*%z1.samples- sqrtVeff%*%z2.samples  #rej
  z.eff.samples <-  -sqrtVeff%*%z2.samples # accept

  # setup the gamma candidates
  allgamm<- (seq(1-1e-10 ,1e-10,length.out=50))
  Lgamma<-length(allgamm)
  critbias<-critbias2<-critvar <- critvar2 <-critmse<- critmse2 <- matrix(NA,Lgamma,p)
  critmse.tau <- numeric(Lgamma)

  for(kkk in 1:Lgamma){
    gamm <- allgamm[kkk]
    Icomb<- Tstat1 < qchisq(1-gamm,df=p)
    id.yes <-  apply(z1.samples*z1.samples,2,sum) < qchisq(1-gamm,df=p) # integrate or not
    mat.yes <- matrix(id.yes,p,ngen,byrow=TRUE)


    R <- expm::sqrtm(Sigma.S1)%*%alpha%*%
      MASS::ginv(alpha%*%Sigma.S1%*%alpha)%*%
      alpha%*%expm::sqrtm(Sigma.S1)

    Tstat.tau.samples <- apply((R%*%z1.samples)^2,2, sum)
    # Tstat.tau.samples <-
    #   (alpha%*%expm::sqrtm(Sigma.S1)%*%z1.samples)^2/
    #   c(alpha%*%Sigma.S1%*%alpha)

    id.yes.tau <- Tstat.tau.samples < qchisq(1-gamm,df=1) # degree of freedom should be 1

    gen.value.tau <- alpha%*%z.rt.samples * (1-id.yes.tau)+
      alpha%*%z.eff.samples * id.yes.tau

    # generate the elastic values
    gen.value<- (z.rt.samples *(1-mat.yes)  +
                   z.eff.samples * (mat.yes) )
    gen.value<-t(gen.value)
    # analytic form
    NormalTruncatedFirstMom <- function(mu,
                                        p=1,
                                        a=0,
                                        b=Inf)
    {
      if(a==b) return(rep(0,p))
      else {
        ncp.mu <- t(mu)%*%mu
        dem <- pchisq(b, df=p, ncp = ncp.mu)-
          pchisq(a, df=p, ncp = ncp.mu)
        num <- mu*(
          pchisq(b, df=p+2, ncp = ncp.mu)-
            pchisq(a, df=p+2, ncp = ncp.mu)
        )
        if(dem==0) return(rep(0,p))
        else return(num/dem)

      }
    }
    NormalTruncatedSecondMom <- function(mu,
                                         p=1,
                                         a=0,
                                         b=Inf)
    {
      if(a==b) return(matrix(0, p, p))
      else {
        ncp.mu <- t(mu)%*%mu
        p <- length(mu)
        dem <- pchisq(b, df=p, ncp = ncp.mu)-
          pchisq(a, df=p, ncp = ncp.mu)
        num1 <- diag(p)*(pchisq(b, df=p+2, ncp = ncp.mu)-
                           pchisq(a, df=p+2, ncp = ncp.mu))
        num2 <- mu%*%t(mu)*
          (pchisq(b, df=p+4, ncp = ncp.mu)-
             pchisq(a, df=p+4, ncp = ncp.mu))
        if(dem==0) return(matrix(0, p, p))
        else return((num1+num2)/dem)
        # xdta <- sapply(mu, function(x)rnorm(1000, mean=x))
        # mean(xdta[a<xdta^2&xdta^2<b]^2)
      }
    }

    mu.z1.truncated <- NormalTruncatedFirstMom(muz1, p=p, a = qchisq(1-gamm,df=3),
                                               b= Inf)*
      pchisq(qchisq(1-gamm,df=p),
             df=p, ncp = t(muz1)%*%muz1,
             lower.tail = FALSE)
    mu.tap.analytic <- sqrtVrteff%*%mu.z1.truncated-
      sqrtVeff%*%muz2

    critbias[kkk,]<-apply(gen.value,2,mean,na.rm = TRUE)
    critbias2[kkk,]<-  mu.tap.analytic #biasformula(gamm,noncp,localpar,Veff,Vrteff)

    # variance
    # # -------------------------------------------
    var.z1.truncated <- NormalTruncatedSecondMom(muz1, p=p, a = qchisq(1-gamm,df=p),
                                                 b= Inf)*
      pchisq(qchisq(1-gamm,df=p),
             df=p, ncp = t(muz1)%*%muz1,
             lower.tail = FALSE)-mu.z1.truncated%*%t(mu.z1.truncated)
    var.tap.analytic <- sqrtVrteff%*%var.z1.truncated%*%t(sqrtVrteff) +
      Veff%*%diag(length(muz1))

    critvar[kkk,]<-apply(gen.value,2,var,na.rm = TRUE)
    critvar2[kkk,] <- diag(var.tap.analytic)


    # # -------------------------------------------
    # mse
    mse.tap.analytic <- var.tap.analytic + mu.tap.analytic%*%t(mu.tap.analytic)
    critmse[kkk,]<-apply((gen.value-0)^2,2,mean,na.rm = TRUE)
    critmse2[kkk,] <- diag(mse.tap.analytic)
    critmse.tau[kkk] <- mean((gen.value.tau)^2)
    # ---------------------------------------------------
  }
  mse_ee <- function(cgamma, muz1, muz2,
                     sqrtVrteff, sqrtVeff,
                     alpha = alpha,
                     v = 1)
  {
    # bias
    gamm <- 1-pchisq(cgamma,df=p)
    mu.z1.truncated <- NormalTruncatedFirstMom(muz1, p=p, a = cgamma,
                                               b= Inf)*
      pchisq(cgamma,
             df=p, ncp = t(muz1)%*%muz1,
             lower.tail = FALSE)
    mu.tap.analytic <- sqrtVrteff%*%mu.z1.truncated-
      sqrtVeff%*%muz2
    # variance
    var.z1.truncated <- NormalTruncatedSecondMom(muz1, p=p, a = cgamma,
                                                 b= Inf)*
      pchisq(cgamma,
             df=p, ncp = t(muz1)%*%muz1,
             lower.tail = FALSE)-mu.z1.truncated%*%t(mu.z1.truncated)
    var.tap.analytic <- sqrtVrteff%*%var.z1.truncated%*%t(sqrtVrteff) + Veff%*%diag(length(muz1))
    mse.tap.analytic <- var.tap.analytic + mu.tap.analytic%*%t(mu.tap.analytic)

    if(is.na(v))
      alpha%*%diag(mse.tap.analytic)
    else
      mse.tap.analytic[v,v]
  }
  kkchosen <- apply(critmse, 2, which.min)
  kkchosen.tau <- which.min(critmse.tau) # critmse%*%alpha
  kkchosen.all <- rep(which.min(apply(critmse, 1, sum)), 3)

  # chosen by analytic form
  cgamma.selected.analytic <- sapply(1:3, function(x)optim(1,
                                                           fn = mse_ee,
                                                           method = 'Brent',
                                                           muz1 = muz1, muz2 = muz2,
                                                           sqrtVrteff = sqrtVrteff, sqrtVeff = sqrtVeff,
                                                           v=x,
                                                           lower = 1e-10, upper = 50)$par)

  gamma.selected.analytic <- 1 - pchisq(cgamma.selected.analytic,df=p)

  # chosen by analytic form (in terms of tau)
  cgamma.selected.analytic.all <- optim(1,
                                        fn = mse_ee,
                                        method = 'Brent',
                                        muz1 = muz1, muz2 = muz2,
                                        sqrtVrteff = sqrtVrteff, sqrtVeff = sqrtVeff,
                                        v=NA, alpha = c(1,1,1),
                                        lower =1e-10, upper = 200)$par

  gamma.selected.analytic.all <- 1 - rep(pchisq(cgamma.selected.analytic.all,df=p),p)

  nuispar["gamma.tau"] <- gamm.selected.tau <-  0.05
  nuispar["c_gamma.tau"] <- cgamma.selected.tau <-
    qchisq(1-gamm.selected.tau,df=1)

  # select the optimal gamma values (by simulation or analytic form)
  if(fixed)
  {
    nuispar[paste0("gamma", 1:p)] <- gamm.selected <-  rep(0.05, p)
    nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- qchisq(1-gamm.selected,df=p)


  }else{
    nuispar[paste0("gamma", 1:p)] <- gamm.selected <-  sapply(rep(kkchosen[1],p),
                                                              function(x)allgamm[x])
    nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- sapply(rep(kkchosen[1], p),
                                                                 function(x)
                                                                   qchisq(1-allgamm[x],df=p,ncp=0))

  }

  # output the localpar, eta
  nuispar[paste0("eta",1:p)] <- localpar
  nuispar[paste0("Icomb", 1:p)]<-Icomb<- sapply(cgamma.selected, function(x)x>Tstat1)
  nuispar['Icomb.tau']<-Icomb.tau<- cgamma.selected.tau>Tstat.tau
  # pvalue smoothed
  nuispar[paste0("Icomb.pval", 1:p)] <- Icomb.pval <- 1 - pchisq(Tstat1, df = p)
  # combine the AIPW with the integrated estimator
  est[paste("elas.",1:p,sep="" )]<-(1-Icomb)*est[paste("covj.t.",1:p,sep="" )]+
    Icomb*est[paste("opt.ee(ml).",1:p,sep="" )]#- biastemp
  est['elas.tau']<-(1-Icomb.tau)*alpha%*%est[paste("covj.t.",1:p,sep="" )]+
    Icomb.tau*alpha%*%est[paste("opt.ee(ml).",1:p,sep="" )]#- biastemp
  # expected bias for the elastic estimator
  biastemp <- -pchisq(qchisq(1-gamm.selected,df=p),
                      df=p+2, ncp = t(muz1)%*%muz1)*Veff%*%localpar
  biastemp.tau <- alpha%*%Veff%*%expm::sqrtm(Sigma.S1)%*%{
    (muz1-R%*%muz1)*{1-pchisq(qchisq(1-gamm.selected.tau,df=p),
                              df=p, ncp = t(muz1)%*%R%*%muz1)}+
      R%*%muz1*{1-pchisq(qchisq(1-gamm.selected.tau,df=p),
                         df=p+2, ncp = t(muz1)%*%R%*%muz1)}
  }-alpha%*%expm::sqrtm(Veff)%*%muz2

  est[paste("elas.",1:p,'.debiased',sep="" )] <-
    est[paste("elas.",1:p,sep="" )]  - biastemp/sqrt(m)

  est['elas.tau.debiased'] <-
    est['elas.tau']  - biastemp.tau/sqrt(m)
  # begin our elastic process
  ## indicator combing
  mat.yes <- sapply(cgamma.selected, function(x)apply(z1.samples*z1.samples,
                                                      2,sum)<x)
  mat.yes <- t(mat.yes)
  gen.value<- (z.eff.samples *mat.yes  + z.rt.samples * (1-mat.yes) )
  gen.value<-t(gen.value)
  ve[paste("elas.",1:3,sep="" )] <- apply(gen.value,2,var,na.rm = TRUE)/m

  ve['elas.tau'] <- var((gen.value%*%alpha), na.rm = TRUE)/m



  # inference
  # by default for psi
  alpha1.vect <- rep(0.025, p)
  alpha2.vect <- 0.05 - alpha1.vect

  UU.muz1 <- sapply(1:p, function(v){
    alpha2 <- alpha2.vect[v]
    # mvtnorm::qmvnorm(1-alpha2/2,
    #       mean = apply(ptb.S.os1,2,mean)*sqrt(m), sigma = Sigma.S1,
    #       tail = 'lower.tail')$quantile

    mvtnorm::qmvnorm(1-alpha2,
            mean = c(muz1), sigma = diag(p),
            tail = 'both.tails')$quantile
  })

  LL.muz1 <- -UU.muz1

  # by default for psi.tau
  alpha1.tau <- 0.025
  alpha2.tau <- 0.05 - alpha1.tau
  UU.muz1.tau <- qnorm(1-alpha2.tau/2,
                       mean = alpha%*%(muz1), sd = sqrt(p))

  LL.muz1.tau <- qnorm(alpha2.tau/2,
                       mean = alpha%*%(muz1), sd = sqrt(p))

  UU.localpar <- mvtnorm::qmvnorm(1-alpha2.tau/2,
                         mean = localpar, sigma = Sigma.S1,
                         tail = 'both.tails')$quantile
  LL.localpar <- -UU.localpar

  # by simulation
  ## for each psi
  generate_elastic <- function(muz1.new)
  {
    z0 <-mvtnorm::rmvnorm(ngen, mean=rep(0,3), diag(3))
    zz0<-mvtnorm::rmvnorm(ngen, mean=rep(0,3), diag(3))
    z0<-t(z0);zz0<-t(zz0)
    # muz1<-as.vector(sqrtiSigSS%*%c(localpar.sample))
    # muz2<-as.vector(expm::sqrtm(Veff)%*%c(localpar.sample))
    z1 <- z0+ c(muz1.new)
    # muz2.new <- mvtnorm::rmvnorm(1, mean=muz2, diag(3))
    z2 <- zz0 + c(muz2)
    mat.yes <- sapply(cgamma.selected, function(x)apply(z1*z1,2,sum)<x)
    mat.yes <- t(mat.yes)
    z.rt <-  + sqrtVrteff%*%z1- sqrtVeff%*%z2  #rej
    z.eff <-  -sqrtVeff%*%z2 # accept
    gen.value<- (z.eff *mat.yes  + z.rt * (1-mat.yes) )
    gen.value<-t(gen.value)
    gen.value
  }
  ## for tau
  generate_elastic_tau <- function(localpar.new,
                                   conservative = F)
  {

    z1 <- t(mvtnorm::rmvnorm(ngen, mean=expm::sqrtm(iSigSS)%*%c(localpar.new),
                       diag(p)))
    if(conservative)
    {
      z2<- t(mvtnorm::rmvnorm(ngen, mean=expm::sqrtm(Veff)%*%c(localpar.new),
                       diag(3)))
    }else{
      z2 <- t(mvtnorm::rmvnorm(ngen, mean=muz2, diag(p)))
    }

    # Tstat.tau.samples <-
    #   (alpha%*%expm::sqrtm(Sigma.S1)%*%z1.samples)^2/
    #   c(alpha%*%Sigma.S1%*%alpha)

    # apply((R%*%z1)^2,2, sum)

    Tstat.tau.mat <- c(alpha%*%expm::sqrtm(Sigma.S1)%*%z1)^2/
      c(alpha%*%Sigma.S1%*%alpha)

    gen.value.tau <- alpha%*%Veff%*%expm::sqrtm(Sigma.S1)%*%z1*
      (Tstat.tau.mat>cgamma.selected.tau)-
      alpha%*%sqrtVeff%*%z2
    gen.value.tau
  }


  # simulation for elastic inference
  nngen <- 100
  qq1<-qq2 <-matrix(NA,nngen,p)
  qq.tau.1 <- qq.tau.2 <- numeric(nngen)
  Allgen.value <- NULL
  for(nn in 1:nngen){
    # combine all the samples
    muz1.new <-mvtnorm::rmvnorm(1, mean=muz1, diag(p))
    Allgen.value <- rbind(Allgen.value,
                          generate_elastic(muz1.new))
    # select the 95% for muz1
    while(sum(muz1.new>UU.muz1|muz1.new<LL.muz1)>0){
      muz1.new <-mvtnorm::rmvnorm(1, mean=muz1, diag(p))
    }
    # re-generate
    gen.value <- generate_elastic(muz1.new)
    qq1[nn,]<-sapply(1:p, function(v)quantile(gen.value[,v],
                                              probs=alpha1.vect[v]/2,type=5,na.rm = TRUE))
    qq2[nn,]<-sapply(1:p, function(v)quantile(gen.value[,v],
                                probs=1-alpha1.vect[v]/2,type=5,na.rm = TRUE))

    # adaptive for the tau
    localpar.new <- mvtnorm::rmvnorm(1, mean = localpar,
                            Sigma.S1)
    while(sum(localpar.new>UU.localpar|localpar.new<LL.localpar)>0){
      localpar.new <- mvtnorm::rmvnorm(1, mean = localpar,
                              Sigma.S1)
    }

    qq.tau.1[nn] <- quantile(generate_elastic_tau(localpar.new,
                                    conservative = T),
               probs=alpha1.tau/2,type=5, na.rm = TRUE)
    qq.tau.2[nn] <- quantile(generate_elastic_tau(localpar.new,
                                    conservative = T),
               probs=(1-alpha1.tau/2),type=5, na.rm = TRUE)
  }

  bootq1 <- bootq2 <- est

  ## construct the Wald CI for the regular estimators
  regular.est.name <- c(
    paste0("covj.t.", 1:3),
    paste0("opt.ee.", 1:3),
    paste0("ee.rt(ml).", 1:3),
    paste0("opt.ee(ml).", 1:3),
    paste0("elas.", 1:3)
  )

  bootq1 <-
    est[regular.est.name] - qnorm(1 - 0.05 / 2) *
    sqrt(ve[regular.est.name])

  bootq2 <-
    est[regular.est.name] + qnorm(1 - 0.05 / 2) *
    sqrt(ve[regular.est.name])


  # begin construct the ECI
  nuispar["conservative"]<-Tstat1< thres.psi

  nuispar['Tstat.psi'] <- Tstat1
  nuispar['Tstat.tau'] <- Tstat.tau

  if(Tstat1< thres.psi ){

    # simulation
    ## version 1
    bootq1[paste("elas1.",1:p,sep="" )]<- apply(qq1,2,min)/sqrt(m)
    bootq2[paste("elas1.",1:p,sep="" )]<- apply(qq2,2,max)/sqrt(m)

    ## version 2
    bootq1[paste("elas2.",1:p,  sep="" )] <- #est[paste("elas.",1:p,sep="" )] +
      apply(Allgen.value,2,quantile,probs=0.05/2,type=5,na.rm = TRUE)/sqrt(m)
    bootq2[paste("elas2.",1:p, sep="" )]<- #est[paste("elas.",1:p,sep="" )] +
      apply(Allgen.value,2,quantile,probs=(1-0.05/2),type=5,na.rm = TRUE)/sqrt(m)

  }
  # if  statistics over the threshold
  if(Tstat1>= thres.psi){

    z.rt <- t(mvtnorm::rmvnorm(1000, mean=rep(0,3), Vrt))
    z.rt <- t(z.rt)
    gen.value <- z.rt
    # both two versions are the same
    bootq1[paste("elas1.",1:p,sep="" )] <- bootq1[paste("elas2.",1:p,sep="" )] <- #est[paste("elas.",1:p,sep="" )] +
      apply(gen.value,2,quantile,probs=0.025,type=5,na.rm = TRUE)/sqrt(m)

    bootq2[paste("elas1.",1:3,sep="" )] <- bootq2[paste("elas2.",1:3,sep="" )] <- #est[paste("elas.",1:p,sep="" )] +
      apply(gen.value,2,quantile,probs=0.975,type=5,na.rm = TRUE)/sqrt(m)

  }


  return(list(
    est = est, ve = ve,
    CIs.inf = bootq1,
    CIs.sup = bootq2,
    nuispar = nuispar))
}
