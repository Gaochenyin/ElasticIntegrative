
main<-function(seed){
  Data.list <- GenerateData(beta0 = beta0, # for the mu0 function
                            psi0 = psi, # for constrast function
                            om.id = 0, es.id = 0,
                            n = 1e5, mean.x = 1,  # the whole population
                            n.t = n.t, # use the default for RCT
                            m = m, tlocalpar = tlocalpar, # for the RWE
                            seed = seed, overlap = overlap,
                            equal.weight = equal.weight)
  out<-ElasticEst(dat.t = Data.list$RT,  # RT
                  dat.os = Data.list$RW, # RW
                  es.id = 0, om.id = 0,
                  alpha = c(1,1,1), # for tau = \alpha^T \psi
                  fixed = fixed, cutoff = cutoff,
                  analytic= FALSE, # way to selection c_gamma
                  joint = joint, # jointly select c_gamma
                  .which = 1L, # criteria for c_gamma selection. 1L: Here we choose based on mse(psi_1)
                  epsilon = epsilon, # construct the smoothed version of the estimators
                  conservative = FALSE, # FALSE: project only on the irregular part, i.e. mu_1
                  only.RT = only.RT, sample.split = sample.split
  )
  est<-out$est
  ve<-out$ve
  bootq1<-out$bootq1
  bootq2<-out$bootq2
  nuispar<-out$nuispar
  sample.split.est <- out$split.smooth
  return(list(est=est,ve=ve,
              bootq1=bootq1,bootq2=bootq2,
              nuispar=nuispar,
              sample.split.est = sample.split.est))
}

# handy functions
expoit <- function(x) {return (exp(x)/(1+exp(x)))}
logit <- function(x) {return (log(x/(1-x)))}

GenerateData <- function(beta0 = c(0, 1, 1, 1), # for the mu0 function
                         psi0 = c(0, 1, 1), # for constrast function
                         om.id = 0, es.id = 0,
                         n = 1e5, mean.x = 1,  # the whole population
                         m = 2000, tlocalpar = 0, # size and bias for the RWE
                         n.t = NULL, equal.weight = T, # size for the RCT (NULL uses the default)
                         seed = 1, 
                         overlap = T)
{
  set.seed(seed)
  if(overlap == T)
  {
    # change mean to 1; intercept
    X1 <- rnorm(n, mean.x, 1); X2 <- rnorm(n, mean.x, 1);
    X3 <- rnorm(n, mean.x, 1) #+ sin(X1)-sin(X2)
    X1.temp <- rnorm(n, mean.x, 1); X2.temp <- rnorm(n, mean.x, 1);
    X3.temp <- rnorm(n, mean.x, 1) #+sin(X1.temp) - sin(X2.temp)
    ## construct functions
    h.X <- cbind(X1, X2, X3);h.Xtemp <- cbind(X1.temp, X2.temp, X3.temp)
    f.X <- cbind(X1, X2);f.Xtemp <- cbind(X1.temp, X2.temp)
    
    ## construct Y for RCT and RWE
    ### RCT
    Y0   <-  cbind(1, h.X)%*%beta0+rnorm(n,0,1)
    Y1   <-  cbind(1, h.X)%*%beta0+cbind(1, f.X)%*%psi0+rnorm(n,0,1)
    ### RWE
    Y0temp   <-  cbind(1, h.Xtemp)%*%beta0+rnorm(n,0,1)
    Y1temp   <-  cbind(1, h.Xtemp)%*%beta0+cbind(1, f.Xtemp)%*%psi0+rnorm(n,0,1)
    

  }else
  {
    # change mean to 0; intercept
    X1 <- rnorm(n, 0, 1); X2 <- rnorm(n, 0, 1);
    X3 <- rnorm(n, 0, 1) #+ sin(X1)-sin(X2)
    X1.temp <- rnorm(n, 0, 1); X2.temp <- rnorm(n, 0, 1);
    X3.temp <- rnorm(n, 0, 1) #+sin(X1.temp) - sin(X2.temp)
    ## construct functions
    h.X <- cbind(X1, X2, X3);h.Xtemp <- cbind(X1.temp, X2.temp, X3.temp)
    f.X <- cbind(X1, X2);f.Xtemp <- cbind(X1.temp, X2.temp)
    
    ## construct Y for RCT and RWE
    ### RCT # beta0 = (0, 1, 1, 0)
    Y0   <-  cbind(1, h.X)%*%beta0+rnorm(n,0,1)
    Y1   <-  cbind(1, h.X)%*%beta0+cbind(1, X1, abs(X2))%*%psi0+rnorm(n,0,1)
    ### RWE
    Y0temp   <-  cbind(1, h.Xtemp)%*%beta0+rnorm(n,0,1)
    Y1temp   <-  cbind(1, h.Xtemp)%*%beta0+cbind(1, X1.temp, abs(X2.temp))%*%psi0+rnorm(n,0,1)
    }
  # true ATE
  psi.pop<-lm( (Y1-Y0)~f.X)$coeff
  true.tau<-mean(Y1)-mean(Y0)
  true.tau
  taui<-Y1-Y0
  
  ## RCT data
  if(is.null(n.t))
  {
    eS <- expoit(-4.5 -2 * X1 -2 *  X2)
  }else
  {
    e.opt <- uniroot(function(e)
    {eS <- expoit(e  -2 * X1 -2 *  X2)
    mean(eS)*n-n.t}, c(-100,100))$root
    
    eS <- expoit(e.opt -2 * X1 -2 *  X2)
  }
  # eS <- expoit(-5.5 -2 * X1 -2 *  X2)
  mean(eS)*n
  S <- sapply(eS,rbinom,n = 1, size = 1)
  S.ind <- which(S==1)
  n.t <- length(S.ind)
  ### choose the selected patients
  X.t <- f.X[S.ind,]
  Y0.t <- Y0[S.ind]
  Y1.t <- Y1[S.ind]
  eS.t <- eS[S.ind]
  if(overlap == F)
  {
    keep <-which(X.t[,2]>0) # change to larger than 0
    X.t <- X.t[keep,]
    Y0.t <- Y0.t[keep]
    Y1.t <- Y1.t[keep]
    n.t <- length(Y1.t)
  }
  # randomly assign treatment to trial participant
  if(all(equal.weight == T))
  {
    ps.t <- rep(0.5, n.t)
    A.t  <- rbinom(n.t,1, ps.t) 
  }else{
    # equal.weight <- c(-3,-3)
    # equal.weight <- c(-3, -3)
    ps.t <- expoit(-2+X.t%*%equal.weight)
    hist(ps.t)
    A.t  <- rbinom(n.t,1,ps.t)
  }
  Y.t  <-  Y1.t*A.t+Y0.t*(1-A.t)
  
  ## OS data
  P.ind <- sample(1:n,size = m)
  X.os <- f.Xtemp[P.ind,]
  X.confounder.os <- X3.temp[P.ind]
  Y0.os <- Y0temp[P.ind]
  Y1.os <- Y1temp[P.ind]
  
  # mantain at about .5
  a.opt <- uniroot(function(a)
  {ps <- expoit(a - 1*X.os[,1] + 1*X.os[,2] - tlocalpar*X.confounder.os )
  mean(ps)-.5}, c(-100,100))$root
  eA <-   expoit(a.opt - 1*X.os[,1] + 1*X.os[,2] - tlocalpar*X.confounder.os )
  # eA <-   expoit(-1 - 1*X.os[,1] - tlocalpar*X.confounder.os )
  mean(eA)
  #hist(eA)
  A.os <- sapply(eA,rbinom,n=1,size=1)
  Y.os <- Y1.os*A.os+Y0.os*(1-A.os)
  
  # return the RT and RW data-sets
  dat.t  <- list(Y = Y.t, A = A.t, X = X.t, q = rep(1,n.t),
                 ps = ps.t, 
                 ml.ps = ps.t)
  dat.os  <- list(Y = Y.os, A = A.os, X = X.os, q = rep(1,m))
  list(RT = dat.t,
       RW = dat.os)
}
# # # # a realization
# data.list <- GenerateData(tlocalpar = 0, equal.weight=F)
# dat.t <- data.list$RT;
# dat.os <- data.list$RW
# es.id = 0; om.id = 0;
# alpha = c(1,1,1)
# analytic= FALSE; # way to selection c_gamma
# joint = TRUE; .which = 1# choose only one c_gamma
# fixed = FALSE; # fixed the c_gamma, no selection
# conservative = FALSE

# elastic estimation
ElasticEst <- function(dat.t,  # RT
                       dat.os, # RW
                       es.id = 0, om.id = 0,
                       alpha = c(1,1,1), # for tau(X)
                       analytic= FALSE, # way to selection c_gamma
                       joint = FALSE, .which = 1, # choose only one c_gamma
                       epsilon = 1, # smooth the non-regularity
                       fixed = FALSE, cutoff = 0.05, # fixed the c_gamma, no selection
                       conservative = FALSE, # use conservative CI for psi_sum
                       only.RT = F, sample.split = T
)
{
  
  n.t <- dat.t$Y%>%length()
  m <- dat.os$Y%>%length()
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
    if(es.id == 0)
    {
      # OLS
      glm.out.ols <- glm(A~X,
                         family=quasibinomial,
                         weights=q,data=dat.os)
      dat.os$ps <- glm.out.ols$fitted.values
      # sieve method
      library(dplyr)
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
    }
    if(om.id == 0)
    {
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
      
    }
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
    library(rootSolve)
    opt.integ <- 
      multiroot(f = ee1, start = reg.t, dat = dat.integ)$root
    # RT ml
    opt.ml.t <- 
      multiroot(f = ee1.ml, start = reg.t, dat = dat.t)$root
    # RW ml
    opt.ml.integ <- 
      multiroot(f = ee1.ml, start = reg.t, dat = dat.integ)$root
    
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
  ve['AIPW.tau'] <- (ptb[,paste("covj.t.",1:3, sep="" )]%*%alpha)%>%var
  ve['rt.tau'] <- (ptb[,paste("ee.rt(ml).",1:3, sep="" )]%*%alpha)%>%var
  ve['ee.tau'] <- (ptb[,paste("opt.ee(ml).",1:3, sep="" )]%*%alpha)%>%var
  if(only.RT)return(list(est=est,ve=ve))
  # compute sigma.S by formula
  library(MASS)
  library(expm)
  rho <- n.t/m
  I.rt <- ginv(Vrt)/rho
  I.rw <- ginv(Veff) - ginv(Vrt)
  Gamma <-  ginv(I.rt)%*%I.rw*rho^(-1/2)
  Sigma.S1 <- t(Gamma)%*%I.rt%*%Gamma + I.rw
  iSigSS<- ginv(Sigma.S1)
  sqrtiSigSS<- sqrtm( iSigSS )
  sqrtVrteff <- Veff%*%sqrtm(Sigma.S1)
  sqrtVeff <- sqrtm(Veff)
  # test statistics and localparameter for psi
  Tstat1<- m*t(apply(ptb.S.os1,2,mean))%*%
    ginv(Sigma.S1)%*%
    (apply(ptb.S.os1,2,mean))
  localpar<-apply(ptb.S.os1,2,mean)*sqrt(m)
  noncp <-t(localpar)%*%iSigSS%*%(localpar)
  
  # test statistics for tau
  Tstat.tau <- (alpha%*%localpar)*
    ginv(alpha%*%Sigma.S1%*%alpha)*  
    (alpha%*%localpar)
  
  
  # simulation for c_gamma selection
  muz1 <- sqrtm(iSigSS)%*%localpar
  muz2 <- sqrtm(Veff)%*%localpar
  ngen <- 1000
  library(mvtnorm)
  z1.samples <- rmvnorm(ngen, mean=muz1, diag(3)) # z1
  z2.samples <- rmvnorm(ngen, mean=muz2, diag(3)) # z2
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
    
    
    R <- sqrtm(Sigma.S1)%*%alpha%*%
      ginv(alpha%*%Sigma.S1%*%alpha)%*%
      alpha%*%sqrtm(Sigma.S1)
    
    Tstat.tau.samples <- apply((R%*%z1.samples)^2,2, sum)
    # Tstat.tau.samples <- 
    #   (alpha%*%sqrtm(Sigma.S1)%*%z1.samples)^2/
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
  kkchosen.all <- which.min(apply(critmse, 1, sum))%>%
    rep(3)
  
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
  
  # select the optimal gamma values (by simulation or analytic form)
  if(fixed)
  {
    nuispar[paste0("gamma", 1:p)] <- gamm.selected <-  rep(cutoff, p)
    nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- qchisq(1-gamm.selected,df=p)
    
    nuispar["gamma.tau"] <- gamm.selected.tau <-  cutoff
    nuispar["c_gamma.tau"] <- cgamma.selected.tau <- 
      qchisq(1-gamm.selected.tau,df=1)
  }else{
    nuispar["gamma.tau"] <- gamm.selected.tau <-  allgamm[kkchosen.tau]
    nuispar["c_gamma.tau"] <- cgamma.selected.tau <- 
      qchisq(1-gamm.selected.tau,df=1,ncp=0)
    if(analytic)
    {
      if(joint)
      {
        if(.which == 'sum'){
          nuispar[paste0("gamma", 1:p)] <- gamm.selected <-  gamma.selected.analytic.all
          nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- cgamma.selected.analytic.all
        }else{
          nuispar[paste0("gamma", 1:p)] <- gamm.selected <-  rep(gamma.selected.analytic[.which], p)
          nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- rep(cgamma.selected.analytic[.which], p)
        }
        
      }else{
        nuispar[paste0("gamma", 1:p)] <- gamm.selected <-  gamma.selected.analytic
        nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- cgamma.selected.analytic
      }
      
    }else{
      if(joint)
      {
        if(.which == 'sum'){
          nuispar[paste0("gamma", 1:p)] <- gamm.selected <-  sapply(kkchosen.all, function(x)allgamm[x])
          nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- sapply(kkchosen.all, function(x)
            qchisq(1-allgamm[x],df=p,ncp=0))
        }else{
          nuispar[paste0("gamma", 1:p)] <- gamm.selected <-  sapply(rep(kkchosen[.which],p),
                                                                    function(x)allgamm[x])
          nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- sapply(rep(kkchosen[.which], p),
                                                                       function(x)
                                                                         qchisq(1-allgamm[x],df=p,ncp=0))
        }
      }else{
        nuispar[paste0("gamma", 1:p)] <- gamm.selected <-  sapply(kkchosen, function(x)allgamm[x])
        nuispar[paste0("c_gamma", 1:p)] <- cgamma.selected <- sapply(kkchosen, function(x)
          qchisq(1-allgamm[x],df=p,ncp=0))
      }
    }
  }
  
  # output the localpar, eta
  nuispar[paste0("eta",1:p)] <- localpar
  nuispar[paste0("Icomb", 1:p)]<-Icomb<- sapply(cgamma.selected, function(x)x>Tstat1)
  nuispar['Icomb.tau']<-Icomb.tau<- cgamma.selected.tau>Tstat.tau
  # pvalue smoothed
  nuispar[paste0("Icomb.pval", 1:p)] <- Icomb.pval <- 1 - pchisq(Tstat1, df = p)
  # icdf smoothed
  nuispar[paste0("Icomb.icdf", 1:p)] <- Icomb.icdf <- sapply(cgamma.selected,
                                                             function(x)pnorm(x-Tstat1, sd = epsilon))
  # combine the AIPW with the integrated estimator
  est[paste("elas.",1:p,sep="" )]<-(1-Icomb)*est[paste("covj.t.",1:p,sep="" )]+
    Icomb*est[paste("opt.ee(ml).",1:p,sep="" )]#- biastemp
  est['elas.tau']<-(1-Icomb.tau)*alpha%*%est[paste("covj.t.",1:p,sep="" )]+
    Icomb.tau*alpha%*%est[paste("opt.ee(ml).",1:p,sep="" )]#- biastemp
  # pvalue smoothed
  est[paste("elas.pval",1:p,sep="" )]<-(1-Icomb.pval)*est[paste("covj.t.",1:p,sep="" )]+
    Icomb.pval*est[paste("opt.ee(ml).",1:p,sep="" )]#- biastemp
  # icdf smoothed
  est[paste("elas.icdf",1:p,sep="" )]<-(1-Icomb.icdf)*est[paste("covj.t.",1:p,sep="" )]+
    Icomb.icdf*est[paste("opt.ee(ml).",1:p,sep="" )]#- biastemp
  # expected bias for the elastic estimator
  biastemp <- -pchisq(qchisq(1-gamm.selected,df=p),
                      df=p+2, ncp = t(muz1)%*%muz1)*Veff%*%localpar
  biastemp.tau <- alpha%*%Veff%*%sqrtm(Sigma.S1)%*%{
    (muz1-R%*%muz1)*{1-pchisq(qchisq(1-gamm.selected.tau,df=p),
                              df=p, ncp = t(muz1)%*%R%*%muz1)}+
      R%*%muz1*{1-pchisq(qchisq(1-gamm.selected.tau,df=p),
                         df=p+2, ncp = t(muz1)%*%R%*%muz1)}
  }-alpha%*%sqrtm(Veff)%*%muz2
  
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
  
  ve['elas.tau'] <- (gen.value%*%alpha)%>%var(na.rm = TRUE)/m
  # soft thresholding (smooth) using p-value
  mat.yes.pval <- 1 - pchisq(apply(z1.samples*z1.samples,2,sum), df = p)
  gen.value.pval <- t(z.eff.samples) * mat.yes.pval + t(z.rt.samples) * (1-mat.yes.pval)
  ve[paste("elas.pval.",1:3,sep="" )] <- apply(gen.value.pval,
                                               2,var,na.rm = TRUE)/m
  ## soft thresholding (smooth) using inverse cdf
  ### requires the hyper-parameter epsilon
  mat.yes.icdf <- sapply(cgamma.selected,
                             function(x)
                             {
                               pnorm(x-apply(z1.samples*z1.samples,2,sum), sd = epsilon)
                             })
  mat.yes.icdf <- t(mat.yes.icdf)
  gen.value.icdf <-(z.eff.samples *mat.yes.icdf  + z.rt.samples * (1-mat.yes.icdf) )
  gen.value.icdf <- t(gen.value.icdf)
  ve[paste("elas.icdf.",1:3,sep="" )] <- apply(gen.value.icdf,
                                                 2,var,na.rm = TRUE)/m
  ## sample-splitting + fixed c_gamma + cross-validation
  B <- 50
  if(sample.split == T)
  {
    psi.split.list <- sapply(1:B, function(b)
    {
      set.seed(b)
      # Step1: bootstrap the dat.os and dat.t
      idx.t <- sample(1:n.t, n.t, replace = T)
      idx.os <- sample(1:m, m, replace = T)
      dat.t.boot <- lapply(dat.t, function(x)
      {
        if(is.null(dim(x)))x[idx.t]
        else x[idx.t,]
      }
      )
      dat.os.boot <- lapply(dat.os, function(x)
      {
        if(is.null(dim(x)))x[idx.os]
        else x[idx.os,]
      }
      )
      # Step2: split more folds
      library(caret)
      nfolds <- 10
      folds.t <- createFolds(1:n.t, nfolds, list = T)
      folds.os <- createFolds(1:m, nfolds, list = T)
      
      # ## RCT
      # idx.t.fold1 <- sample(1:n.t, n.t/2, replace = F)
      # idx.t.fold2 <- setdiff(1:n.t, idx.t.fold1)
      # ## RWE
      # idx.os.fold1 <- sample(1:m, m/2, replace = F)
      # idx.os.fold2 <- setdiff(1:m, idx.os.fold1)   
      
      # Step3: conduct tests on each fold
      test.statistc <- function(list.ptb, n.t, m)
      {
        # compute sigma.S by formula
        ptb <- list.ptb$ptb
        ptb.S.os1 <- list.ptb$ptb.S.os1
        # RT and Eff
        ptboptee1<-ptb[,paste("ee.rt(ml).",1:3,sep="" )]
        ptboptee2<-ptb[,paste("opt.ee(ml).",1:3,sep="" )]
        # use the Sigma.S1 
        psi.rt <- apply(ptboptee1, 2, mean)
        psi.ee <- apply(ptboptee2, 2, mean)
        Tstat <- m*t(apply(ptb.S.os1,2,mean))%*%
          ginv(Sigma.S1)%*%(apply(ptb.S.os1,2,mean))
        return(c(psi.rt, psi.ee, test = Tstat))
      }
      folded.estimation <- function(idx.t, idx.os)
      {
        dat.t.fold <- lapply(dat.t.boot, function(x)
        {
          if(is.null(dim(x)))x[idx.t]
          else x[idx.t,]
        })
        dat.os.fold <- lapply(dat.os.boot, function(x)
        {
          if(is.null(dim(x)))x[idx.os]
          else x[idx.os,]
        })
        ptb.fold <- permutation_est(dat.t = dat.t.fold,
                                    dat.os = dat.os.fold, 
                                    nptb = 50)
        test.est.fold <- test.statistc(list.ptb = ptb.fold, 
                                       n.t = length(idx.t),
                                       m = length(idx.os))
        return(test.est.fold)
      }
      ## RCT
      psi.fold.list <- mapply(function(fold.idx.t, fold.idx.os)
      {
        test.est.fold <- folded.estimation(fold.idx.t, fold.idx.os)
        test.est.fold.rest <- folded.estimation(-fold.idx.t, -fold.idx.os)
        # fixed c_gamma
        c_gamma_fixed <- qchisq(1 - 0.05, df = 3)
        if(test.est.fold['test'] < c_gamma_fixed)
        {
          psi.fold <- test.est.fold.rest[paste('opt.ee(ml).', 1:3, sep = '')]
        }else{
          psi.fold <- test.est.fold.rest[paste('ee.rt(ml).', 1:3, sep = '')]
        }
        return(psi.fold)
      }, fold.idx.t = folds.t, fold.idx.os = folds.os)
      apply(psi.fold.list, 1, mean)
    })
  }else{
    psi.split.list <- NA
  }

  
  # inference
  # by default for psi
  alpha1.vect <- rep(0.025, p)
  alpha2.vect <- 0.05 - alpha1.vect
  
  UU.muz1 <- sapply(1:p, function(v){
    alpha2 <- alpha2.vect[v]
    # qmvnorm(1-alpha2/2,
    #       mean = apply(ptb.S.os1,2,mean)*sqrt(m), sigma = Sigma.S1,
    #       tail = 'lower.tail')$quantile
    
    qmvnorm(1-alpha2,
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
  
  UU.localpar <- qmvnorm(1-alpha2.tau/2,
                         mean = localpar, sigma = Sigma.S1,
                         tail = 'both.tails')$quantile
  LL.localpar <- -UU.localpar
  
  # by simulation
  ## for each psi
  generate_elastic <- function(muz1.new)
  {
    z0 <-rmvnorm(ngen, mean=rep(0,3), diag(3))
    zz0<-rmvnorm(ngen, mean=rep(0,3), diag(3))
    z0<-t(z0);zz0<-t(zz0)
    # muz1<-as.vector(sqrtiSigSS%*%c(localpar.sample))
    # muz2<-as.vector(sqrtm(Veff)%*%c(localpar.sample))
    z1 <- z0+ c(muz1.new)
    # muz2.new <- rmvnorm(1, mean=muz2, diag(3))
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
    
    z1 <- rmvnorm(ngen, mean=sqrtm(iSigSS)%*%c(localpar.new), 
                  diag(p))%>%t()
    if(conservative)
    {
      z2<-rmvnorm(ngen, mean=sqrtm(Veff)%*%c(localpar.new),
                  diag(3))%>%t()
    }else{
      z2 <- rmvnorm(ngen, mean=muz2, diag(p))%>%t()
    }
    
    # Tstat.tau.samples <- 
    #   (alpha%*%sqrtm(Sigma.S1)%*%z1.samples)^2/
    #   c(alpha%*%Sigma.S1%*%alpha)
    
    # apply((R%*%z1)^2,2, sum)
    
    Tstat.tau.mat <- c(alpha%*%sqrtm(Sigma.S1)%*%z1)^2/
      c(alpha%*%Sigma.S1%*%alpha)
    
    gen.value.tau <- alpha%*%Veff%*%sqrtm(Sigma.S1)%*%z1*
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
    muz1.new <-rmvnorm(1, mean=muz1, diag(p))
    Allgen.value <- rbind(Allgen.value,
                          generate_elastic(muz1.new))
    # select the 95% for muz1
    while(sum(muz1.new>UU.muz1|muz1.new<LL.muz1)>0){
      muz1.new <-rmvnorm(1, mean=muz1, diag(p))
    }
    # re-generate
    gen.value <- generate_elastic(muz1.new)
    qq1[nn,]<-sapply(1:p, function(v)gen.value[,v]%>%
                       quantile(probs=alpha1.vect[v]/2,type=5,na.rm = TRUE))
    qq2[nn,]<-sapply(1:p, function(v)gen.value[,v]%>%
                       quantile(probs=1-alpha1.vect[v]/2,type=5,na.rm = TRUE))
    
    # adaptive for the tau
    localpar.new <- rmvnorm(1, mean = localpar,
                            Sigma.S1)
    while(sum(localpar.new>UU.localpar|localpar.new<LL.localpar)>0){
      localpar.new <- rmvnorm(1, mean = localpar,
                              Sigma.S1)
    }
    if(conservative)
    {
      qq.tau.1[nn] <- generate_elastic_tau(localpar.new,
                                           conservative = T)%>%
        quantile(probs=alpha1.tau/2,type=5, na.rm = TRUE)
      qq.tau.2[nn] <- generate_elastic_tau(localpar.new,
                                           conservative = T)%>%
        quantile(probs=(1-alpha1.tau/2),type=5, na.rm = TRUE)
      
    }else{
      qq.tau.1[nn] <- generate_elastic_tau(localpar.new)%>%
        quantile(probs=alpha1.tau/2,type=5, na.rm = TRUE)
      qq.tau.2[nn] <- generate_elastic_tau(localpar.new)%>%
        quantile(probs=(1-alpha1.tau/2),type=5, na.rm = TRUE)
    }
  }
  
  bootq1 <- bootq2 <- est
  # begin construct the ECI
  nuispar["conservative"]<-Tstat1< thres.psi
  nuispar["conservative.tau"]<- Tstat.tau < thres.tau
  
  nuispar['Tstat.psi'] <- Tstat1
  nuispar['Tstat.tau'] <- Tstat.tau
  
  if(Tstat1< thres.psi ){
    
    # simulation
    ## version 1
    bootq1[paste("elas1.",1:p,sep="" )]<- apply(qq1,2,min)/sqrt(m)
    bootq2[paste("elas1.",1:p,sep="" )]<- apply(qq2,2,max)/sqrt(m)
    
    ## version 2
    bootq1[paste("elas2.",1:p,  sep="" )] <- 
      apply(Allgen.value,2,quantile,probs=0.05/2,type=5,na.rm = TRUE)/sqrt(m)
    bootq2[paste("elas2.",1:p, sep="" )]<- 
      apply(Allgen.value,2,quantile,probs=(1-0.05/2),type=5,na.rm = TRUE)/sqrt(m)
    
  }
  if(Tstat.tau < thres.tau)
  {
    bootq1['elas1.tau']<- min(qq.tau.1)/sqrt(m)
    bootq2['elas1.tau']<- max(qq.tau.2)/sqrt(m)
    
    bootq1['elas2.tau'] <- (Allgen.value%*%alpha)%>%
      quantile(probs=0.05/2,type=5,na.rm = TRUE)/sqrt(m)
    bootq2['elas2.tau']<- (Allgen.value%*%alpha)%>%
      quantile(probs=1-0.05/2,type=5,na.rm = TRUE)/sqrt(m)
  }
  # if either statistics over the threshold
  
  if(Tstat1>= thres.psi| Tstat.tau>= thres.tau)
  {
    z.rt <- t(rmvnorm(1000, mean=rep(0,3), Vrt))
    z.rt <- t(z.rt)
    gen.value <- z.rt
  }
  
  if(Tstat1>= thres.psi){
    
    # both two versions are the same
    bootq1[paste("elas1.",1:p,sep="" )] <- bootq1[paste("elas2.",1:p,sep="" )] <-
      apply(gen.value,2,quantile,probs=0.025,type=5,na.rm = TRUE)/sqrt(m)
    
    bootq2[paste("elas1.",1:3,sep="" )] <- bootq2[paste("elas2.",1:3,sep="" )] <- 
      apply(gen.value,2,quantile,probs=0.975,type=5,na.rm = TRUE)/sqrt(m)
    
  }
  if(Tstat.tau >= thres.tau)
  {
    # for both verision of tau
    bootq1['elas1.tau'] <- bootq1['elas2.tau'] <-
      (gen.value%*%alpha)%>%quantile(probs=0.025,type=5)/sqrt(m)
    bootq2['elas1.tau'] <- bootq2['elas2.tau'] <- 
      (gen.value%*%alpha)%>%quantile(probs=0.975,type=5)/sqrt(m)
  }
  
  return(list(est=est,ve=ve,
              bootq1=bootq1,bootq2=bootq2,
              nuispar=nuispar,
              split.smooth = psi.split.list))
}


# score functions
ee1<-function(par,dat){
  psi<-par
  Y<-dat$Y
  A<-dat$A
  X<-dat$X
  q<-dat$q
  ps<-dat$ps
  mu0<-dat$mu0
  H<-Y-A*cbind(1,X)%*%psi
  apply(cbind(1,X)* matrix( (H-mu0)*(A-ps)*q,length(Y),1+dim(X)[2], byrow=FALSE),2,mean)
}

ee1.ml.new<-function(par,dat){
  psi<-par
  Y<-dat$Y
  A<-dat$A
  X<-dat$X
  q<-dat$q
  ps<-dat$ml.ps
  mu0<-dat$ml.mu0
  sigma0 <- dat$ml.sigma0
  sigma1 <- dat$ml.sigma1
  H<-Y-A*cbind(1,X)%*%psi
  # apply(((A==0)*sigma0^(-1)+
  #          (A==1)*sigma1^(-1))*cbind(1,X)* matrix( (H-mu0)*(A-ps)*q,length(Y),1+dim(X)[2],
  #                                                  byrow=FALSE),2,mean)
  apply(sigma0^(-1)*cbind(1,X)* matrix( (H-mu0)*(A-ps)*q,length(Y),1+dim(X)[2],
                                        byrow=FALSE),2,mean)
}

ee1.ml<-function(par,dat){
  psi<-par
  Y<-dat$Y
  A<-dat$A
  X<-dat$X
  q<-dat$q
  ps<-dat$ml.ps
  mu0<-dat$ml.mu0
  H<-Y-A*cbind(1,X)%*%psi
  apply(cbind(1,X)* matrix( (H-mu0)*(A-ps)*q,length(Y),1+dim(X)[2], byrow=FALSE),2,mean)
}

UniquePanelCoords <- ggplot2::ggproto(
  "UniquePanelCoords", ggplot2::CoordCartesian,
  
  num_of_panels = 1,
  panel_counter = 1,
  panel_ranges = NULL,
  
  setup_layout = function(self, layout, params) {
    self$num_of_panels <- length(unique(layout$PANEL))
    self$panel_counter <- 1
    layout
  },
  
  setup_panel_params =  function(self, scale_x, scale_y, params = list()) {
    if (!is.null(self$panel_ranges) & length(self$panel_ranges) != self$num_of_panels)
      stop("Number of panel ranges does not equal the number supplied")
    
    train_cartesian <- function(scale, limits, name, given_range = NULL) {
      if (is.null(given_range)) {
        expansion <- ggplot2:::default_expansion(scale, expand = self$expand)
        range <- ggplot2:::expand_limits_scale(scale, expansion,
                                               coord_limits = self$limits[[name]])
      } else {
        range <- given_range
      }
      
      out <- list(
        ggplot2:::view_scale_primary(scale, limits, range),
        sec = ggplot2:::view_scale_secondary(scale, limits, range),
        arrange = scale$axis_order(),
        range = range
      )
      names(out) <- c(name, paste0(name, ".", names(out)[-1]))
      out
    }
    
    cur_panel_ranges <- self$panel_ranges[[self$panel_counter]]
    if (self$panel_counter < self$num_of_panels)
      self$panel_counter <- self$panel_counter + 1
    else
      self$panel_counter <- 1
    
    c(train_cartesian(scale_x, self$limits$x, "x", cur_panel_ranges$x),
      train_cartesian(scale_y, self$limits$y, "y", cur_panel_ranges$y))
  }
)

coord_panel_ranges <- function(panel_ranges, expand = TRUE, default = FALSE, clip = "on") {
  ggplot2::ggproto(NULL, UniquePanelCoords, panel_ranges = panel_ranges, 
                   expand = expand, default = default, clip = clip)
}
# plot
plt.res <- function(bias, variance)
{
  library(latex2exp)
  library(ggplot2)
  bias <- abs(bias)
  mse <- variance + bias^2
  library(tidyr)
  est_bias_m <- reshape2::melt(bias, varnames = c('b', 'param'))%>%
    cbind(name = 'bias')%>%separate(param, c('param', 'dim'), sep = '[.]')
  
  est_var_m <- reshape2::melt(variance, varnames = c('b', 'param'))%>%
    cbind(name = 'variance')%>%separate(param, c('param', 'dim'), sep = '[.]')
  
  est_mse_m <- reshape2::melt(mse, varnames = c('b', 'param'))%>%
    cbind(name = 'mse')%>%separate(param, c('param', 'dim'), sep = '[.]')
  
  est_all_m <- rbind(est_bias_m, est_var_m, est_mse_m)
  xy.labs <- c(c('psi[2]',
                 'psi[3]'),
               c('bias', 'var', 'mse'))
  names(xy.labs) <- c('2', '3',
                      'bias', 'variance', 'mse')
  
  library(ggplot2)
  # delete the elast 0
  # est_all_m <- subset(est_all_m, subset = dim!=0)
  est_all_m$param <- factor(est_all_m$param, levels = c('AIPW', 'RT', 'EE', 'ELAS'))
  # begin our plots
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  est_all_m <- est_all_m[!est_all_m$dim=='tau',]
  p1 <-
    ggplot(est_all_m, aes(x=b, y=value, 
                          color = param,
                          shape = param))+
    geom_line(size = 0.8)+geom_point(stroke = 2)+
    facet_wrap(c('dim','name'), scales = "free",
               labeller =as_labeller(xy.labs,
                                     default = label_parsed))+theme(legend.position = 'top')+
    ylab('')+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 10),
          text = element_text(size = 12),
          legend.text=element_text(size=12))+
    coord_panel_ranges(panel_ranges = list(
      list(y=c(-0.01,0.2)),
      list(y=c(0,0.022)),
      list(y=c(0,0.04)),
      
      list(y=c(-0.01,0.2)),
      list(y=c(0,0.022)),
      list(y=c(0,0.04))))+
    scale_shape_discrete(solid = T,
                         name = '', label = c('aipw',
                                              'rt', 'eff',
                                              'elas'))+
    scale_color_manual(values = c('#67D5B5',
                                  '#EE7785',
                                  '#C89EC4',
                                  '#84B1ED'),
                       name = '', label = c('aipw',
                                            'rt', 'eff',
                                            'elas'))+
    # c('#67D5B5',
    #   '#EE7785',
    #   '#C89EC4',
    #   '#84B1ED')
    
    theme(plot.margin=unit(c(0.3,0.3,0,0),"cm"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  
  p2 <- ggplot(est_all_m, aes(x=b, y=value, 
                              color = param,
                              shape = param))+
    geom_line()+geom_point(stroke = .25)+
    theme(axis.text = element_text(size = 12),
          text = element_text(size = 15),
          legend.text=element_text(size=12))+
    facet_grid(c('dim','name'), scales = "free",
               labeller =as_labeller(xy.labs,
                                     default = label_parsed))+theme(legend.position = 'top')
  
  library(grid)
  library(gtable) 
  gt1 <-  ggplot_gtable(ggplot_build(p1))
  gt2 <-  ggplot_gtable(ggplot_build(p2))
  gt1$grobs[grep('strip-t.+1$', gt1$layout$name)] = gt2$grobs[grep('strip-t', gt2$layout$name)]
  
  gt.side1 = gtable_filter(gt2, 'strip-r-1')
  gt.side2 = gtable_filter(gt2, 'strip-r-2')
  gt.side3 = gtable_filter(gt2, 'strip-r-3')
  
  gt1 = gtable_add_cols(gt1, widths=gt.side1$widths[1], pos = -1)
  gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))
  
  panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
  gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side3, t = panel_id$t[3], l = ncol(gt1))
  
  
  gt1 <- gtable_add_cols(gt1, widths = unit(0.15, 'cm'), pos = -1)
  grid.newpage()
  grid.draw(gt1)
  
}

plt.res.combine <- function(bias.psi0, bias.psi1,
                            variance.psi0, variance.psi1,
                            AIPW = T)
{
  library(latex2exp)
  library(ggplot2)
  bias.psi0 <- abs(bias.psi0);bias.psi1 <- abs(bias.psi1)
  mse.psi0 <- variance.psi0 + bias.psi0^2
  mse.psi1 <- variance.psi1 + bias.psi1^2
  library(tidyr)
  
  est_bias_m0 <- reshape2::melt(bias.psi0, varnames = c('b', 'param'))%>%
    cbind(name = 'bias')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.0'))
  
  est_bias_m1 <- reshape2::melt(bias.psi1, varnames = c('b', 'param'))%>%
    cbind(name = 'bias')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.1'))
  
  est_var_m0 <- reshape2::melt(variance.psi0, varnames = c('b', 'param'))%>%
    cbind(name = 'variance')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.0'))
  
  est_var_m1 <- reshape2::melt(variance.psi1, varnames = c('b', 'param'))%>%
    cbind(name = 'variance')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.1'))
  
  est_mse_m0 <- reshape2::melt(mse.psi0, varnames = c('b', 'param'))%>%
    cbind(name = 'mse')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.0'))
  
  est_mse_m1 <- reshape2::melt(mse.psi1, varnames = c('b', 'param'))%>%
    cbind(name = 'mse')%>%separate(param, c('param', 'dim'), sep = '[.]')%>%
    mutate(dim.name = paste0(dim, '.1'))
  
  est_all_m <- rbind(est_bias_m0, est_bias_m1,
                     est_var_m0, est_var_m1,
                     est_mse_m0, est_mse_m1)
  xy.labs <- c(c('psi[1]==0',
                 'psi[2]==0',
                 'psi[1]==1',
                 'psi[2]== 1'),
               c('bias', 'var', 'MSE'))
  names(xy.labs) <- c('2.0', '3.0',
                      '2.1', '3.1',
                      'bias', 'variance', 'mse')
  
  library(ggplot2)
  # delete the elast 0
  if(AIPW == F)
  {
    line.color <- c('#EE7785',
                    '#C89EC4',
                    '#84B1ED')
    legend.label <- c('RT', 'Eff', 'Elastic')
    est_all_m <- est_all_m[est_all_m$param!='AIPW',]
    # est_all_m <- subset(est_all_m, subset = dim!=0)
    est_all_m$param <- factor(est_all_m$param, levels = c('RT', 'EE','ELAS'),
                              labels = legend.label)
  }else
  {
    line.color <- c('#EE7785',
                    '#67D5B5',
                    '#C89EC4',
                    '#84B1ED')
    legend.label <- c('RT.AIPW','RT.EE', 'Eff','Elas')
    est_all_m$param <- factor(est_all_m$param, levels = c('AIPW','RT', 'EE','ELAS'),
                              labels = legend.label)
  }
  
  # begin our plots
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  # est_all_m <- est_all_m[!est_all_m$dim=='tau',]
  est_all_m$dim.name <- factor(est_all_m$dim.name, levels = c('2.0','3.0',
                                                              '2.1','3.1'))
  p1 <-
    ggplot(est_all_m, aes(x=b, y=value,
                          color = param,
                          shape = param))+
    geom_line(size = 0.8)+geom_point(stroke = 2)+
    facet_wrap(~name+dim.name, scales = "free",
               labeller =as_labeller(xy.labs,
                                     default = label_parsed))+
    theme(legend.position = 'top')+
    scale_shape_discrete(solid = T,
                         name = '', label = legend.label)+
    scale_color_manual(values = line.color,
                       name = '', label = legend.label)+
    ylab('')+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text = element_text(size = 10),
          text = element_text(size = 12),
          legend.text=element_text(size=12))+
    coord_panel_ranges(panel_ranges = list(
      list(y=c(-0.01,0.15)),
      list(y=c(-0.01,0.15)),
      list(y=c(-0.01,0.15)),
      list(y=c(-0.01,0.15)),
      
      list(y=c(0,0.022)),
      list(y=c(0,0.022)),
      list(y=c(0,0.022)),
      list(y=c(0,0.022)),
      
      list(y=c(0,0.03)),
      list(y=c(0,0.03)),
      list(y=c(0,0.03)),
      list(y=c(0,0.03))
    ))+
    # c('#67D5B5',
    #   '#EE7785',
    #   '#C89EC4',
    #   '#84B1ED')
    theme(plot.margin=unit(c(0.3,0.3,0,0),"cm"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))
  
  p2 <-
    ggplot(est_all_m, aes(x=b, y=value, 
                          color = param,
                          shape = param))+
    geom_line()+geom_point(stroke = .25)+
    theme(axis.text = element_text(size = 12),
          text = element_text(size = 15),
          legend.text=element_text(size=12))+
    facet_grid(name~dim.name, scales = "free",
               labeller =as_labeller(xy.labs,
                                     default = label_parsed))+theme(legend.position = 'top')
  
  library(grid)
  library(gtable) 
  gt1 <-  ggplot_gtable(ggplot_build(p1))
  gt2 <-  ggplot_gtable(ggplot_build(p2))
  gt1$grobs[grep('strip-t.+1$', gt1$layout$name)] = gt2$grobs[grep('strip-t', gt2$layout$name)]
  
  gt.side1 = gtable_filter(gt2, 'strip-r-1')
  gt.side2 = gtable_filter(gt2, 'strip-r-2')
  gt.side3 = gtable_filter(gt2, 'strip-r-3')
  
  gt1 = gtable_add_cols(gt1, widths=gt.side1$widths[1], pos = -1)
  gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))
  
  panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
  gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))
  gt1 = gtable_add_grob(gt1, gt.side3, t = panel_id$t[3], l = ncol(gt1))
  
  
  gt1 <- gtable_add_cols(gt1, widths = unit(0.15, 'cm'), pos = -1)
  grid.newpage()
  grid.draw(gt1)
  
  
}
