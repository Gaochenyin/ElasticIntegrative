
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ElasticIntegrative

<!-- badges: start -->
<!-- badges: end -->

## Overview

The goal of *ElasticIntegrative* is to implement integrative analyses
for the heterogenous treatment effect combining a randomized trial and a
real-world evidence study from [Yang et al.,
(2022)](https://arxiv.org/abs/2005.10579).

Two datasets

-   The randomized trial contains observations on (A,X,Y), where the
    treatment assignment A is randomized.

-   The real-world evidence study contains observations on (A,X,Y),
    where the treatment assignment A may be confounded.

## Installation with `devtools`:

``` r
# it takes around 1 minute for the downloads, including the vignettes
devtools::install_github("Gaochenyin/ElasticIntegrative", 
                         build_vignettes = TRUE) 
library(ElasticIntegrative)
```

## Usage

This is an example for illustration

``` r
library(ElasticIntegrative)
## basic example code for data generation
set.seed(2333)
## setups
beta0 <- c(0, 1, 1, 1) # for the mu0 function
psi0 <- c(0, 1, 1) # for contrast function
n <- 1e5; mean.x <- 1 # the whole population
m <- 2000; tlocalpar <- 0 # size and bias for the RWE
# -------------------------------
## begin the data generation
X1.rt <- rnorm(n, mean.x, 1)
X2.rt <- rnorm(n, mean.x, 1)
X3.rt <- rnorm(n, mean.x, 1) 
X1.os <- rnorm(n, mean.x, 1)
X2.os <- rnorm(n, mean.x, 1)
X3.os <- rnorm(n, mean.x, 1) 
## construct the basis functions
h.X.rt <- cbind(X1.rt, X2.rt, X3.rt)
h.X.os <- cbind(X1.os, X2.os, X3.os)
f.X.rt <- cbind(X1.rt, X2.rt)
f.X.os <- cbind(X1.os, X2.os)
# -------------------------------
## construct Y for RCT and RWE populations
### RCT population
Y0.rt <- cbind(1, h.X.rt) %*% beta0 + rnorm(n, 0, 1)
Y1.rt <- cbind(1, h.X.rt) %*% beta0 + cbind(1, f.X.rt) %*% psi0 + rnorm(n, 0, 1)
### RCT sample
eS <- exp(-4.5 - 2 * X1.rt - 2 * X2.rt)/
  (1 + exp(-4.5 - 2 * X1.rt - 2 * X2.rt))
S <- sapply(eS, rbinom, n = 1, size = 1)
S.ind <- which(S == 1)
n.t <- length(S.ind)
### choose the selected patients
X.rt <- f.X.rt[S.ind, ]
Y0.rt <- Y0.rt[S.ind]
Y1.rt <- Y1.rt[S.ind]
eS.rt <- eS[S.ind]
### assign treatment to trial participant with equal probabilities
ps.t <- rep(0.5, n.t)
A.rt <- rbinom(n.t, 1, ps.t)
Y.rt <- Y1.rt * A.rt + Y0.rt * (1 - A.rt)
# -------------------------------
### RWE population
Y0.os <- cbind(1, h.X.os) %*% beta0 + rnorm(n, 0, 1)
Y1.os <- cbind(1, h.X.os) %*% beta0 + cbind(1, f.X.os) %*% psi0 + rnorm(n, 0, 1)
### OS sample
P.ind <- sample(1:n, size = m)
X.os <- f.X.os[P.ind, ]
X.confounder.os <- X3.os[P.ind]
Y0.os <- Y0.os[P.ind]
Y1.os <- Y1.os[P.ind]
# a.opt is chosen to maintain the proportion of the treated at about .5
a.opt <- uniroot(function(a) {
  ps <- expoit(a - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)
  mean(ps) - .5
}, c(-100, 100))$root
eA <- exp(a.opt - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)/
  {1+exp(a.opt - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)}
A.os <- sapply(eA, rbinom, n = 1, size = 1)
Y.os <- Y1.os * A.os + Y0.os * (1 - A.os)
# organize the RT and RW datasets
dat.t <- list(
  Y = Y.rt, A = A.rt, X = X.rt, q = rep(1, n.t),
  ps = ps.t,
  ml.ps = ps.t
)
dat.os <- list(
  Y = Y.os, A = A.os, X = X.os, q = rep(1, m)
  )
# begin the elastic integrative analysis
# choice of kappa_n, default is sqrt(log(m)) (similar to the BIC criteria)
result.elastic <- elasticHTE(dat.t = dat.t,
                             dat.os = dat.os,
                             thres.psi = sqrt(log(m)), # threshold for ACI psi
                             fixed = FALSE)
```

## Value

| **Argument** | **Description**                                                                                                                                                                        |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| est          | the HTE estimators including                                                                                                                                                           |
| i\)          | `covj.t` the covariate adjustment estimator;                                                                                                                                           |
| ii\)         | `ee.rt(ml)` the RCT estimating equation estimator with sieve models;                                                                                                                   |
| iii\)        | `opt.ee(ml)` the RT&RW estimating equation estimator with sieve models;                                                                                                                |
| iv\)         | `elas`the elastic integrative estimator.                                                                                                                                               |
| ve           | the bootstrap variance estimates for est.                                                                                                                                              |
| CI           | 95% percentile bootstrap confidence intervals. Note `covj.t`, `ee.rt(ml)`, and `opt.ee(ml)` use the Wald CI, and `elas` uses the adaptive confidence with threshold `sqrt(log(2000))`. |

    #> est
    #>     covj.t.1     covj.t.2     covj.t.3  ee.rt(ml).1  ee.rt(ml).2  ee.rt(ml).3 
    #>   0.02851533   1.00451768   0.94310969   0.05562283   1.00080762   0.97255215 
    #> opt.ee(ml).1 opt.ee(ml).2 opt.ee(ml).3       elas.1       elas.2       elas.3 
    #>   0.07641376   0.94788773   1.05146075   0.07641376   0.94788773   1.05146075
    #> ve
    #>     covj.t.1     covj.t.2     covj.t.3  ee.rt(ml).1  ee.rt(ml).2  ee.rt(ml).3 
    #>  0.030520919  0.016696668  0.016015331  0.032511804  0.017733931  0.016544565 
    #> opt.ee(ml).1 opt.ee(ml).2 opt.ee(ml).3       elas.1       elas.2       elas.3 
    #>  0.004569076  0.003832585  0.002932012  0.004462569  0.003891839  0.002854052
    #> CI
    #>                    lower     upper
    #> covj.t.1     -0.31389503 0.3709257
    #> covj.t.2      0.75125978 1.2577756
    #> covj.t.3      0.69507293 1.1911465
    #> ee.rt(ml).1  -0.29777886 0.4090245
    #> ee.rt(ml).2   0.73980156 1.2618137
    #> ee.rt(ml).3   0.72045045 1.2246539
    #> opt.ee(ml).1 -0.05606988 0.2088974
    #> opt.ee(ml).2  0.82655051 1.0692249
    #> opt.ee(ml).3  0.94533250 1.1575890
    #> elas.1       -0.12530744 0.2137688
    #> elas.2        0.86336552 1.1751818
    #> elas.3        0.84138533 1.1157781

## More examples

``` r
browseVignettes("ElasticIntegrative")
```

-   [Simulation: adaptive
    selection](https://gaochenyin.github.io/ElasticIntegrative/doc/sim_psi011_111)
-   [Simulation: comparing AIPW and
    SES](https://gaochenyin.github.io/ElasticIntegrative/doc/sim_AIPWvsSES)
-   [Simulation: fixed
    threshold](https://gaochenyin.github.io/ElasticIntegrative/doc/sim_psi011_111_fixed)
