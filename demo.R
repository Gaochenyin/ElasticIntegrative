rm(list=ls())
run=1
source('R\\functions_elastic.R')

m <- 2000 # size of RW dataset
datasets.all <- GenerateData(psi0 = c(0,1,1),
                             tlocalpar = 2,
                             m = m, seed = 1)

# datasets.all$RT: the dataset for RT
# datasets.all$RW: the dataset for RW
# $Y: is a vector of observed outcome
# $X: is a matrix of covariates without intercept (n x (p-1))
# $A: is a vector of received treatment
# $q: is a vector of 1, for permutation-based variance estimation
# RT$ps, RT$ml.ps: is a vector of known propensity score (only available in RT)

# choice of kappa_n, default is tau * sqrt(log(m)) (similar to BIC criteria)
tau <- 1 # tuning parameter for adaptive confidence interval
thres.psi  <-  tau * sqrt(log(m)) # threshold for ECI psi
result.elastic <- ElasticEst(dat.t = datasets.all$RT,
                             dat.os = datasets.all$RW)

# return: a list called result.elastic
# $est: estimates of the psi in the finite population
# $ve: variance estimation of the psi by permutation-based estimation
# $CIs.inf, CIs.sup: 95% confidence intervals
## Wald CIs for regular estimators;
## ECIs for non-regular estimators, namely elas.v1 and elas.v2
## elas.v1 performs better in our simulation studies.
## $nuispar: estimates for various nuisance parameters
