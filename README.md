
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ElasticIntegrative

<!-- badges: start -->
<!-- badges: end -->

## Overview

The goal of *ElasticIntegrative* is toimplement integrative analyses for
the heterogenous treatment effect combining a randomized trial and a
real-world evidence study.

Two datasets

-   The randomized trial contains observations on (A,X,Y), where the
    treatment assignment A is randomized.

-   The real-world evidence study contains observations on (A,X,Y),
    where the treatment assignment A may be confounded.

## Installation with `devtools`:

``` r
devtools::install_github("Gaochenyin/ElasticIntegrative")
library(ElasticIntegrative)
```

## Usage

This is an example for illustration

``` r
library(ElasticIntegrative)
## basic example code for data generation
set.seed(1234)
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
## construct functions
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
# 0.058 is chosen to maintain the proportion of the treated at about .5
eA <- exp(0.058 - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)/
  {1+exp(0.058 - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)}
A.os <- sapply(eA, rbinom, n = 1, size = 1)
Y.os <- Y1.os * A.os + Y0.os * (1 - A.os)
# sort the RT and RW datasets
dat.t <- list(
  Y = Y.rt, A = A.rt, X = X.rt, q = rep(1, n.t),
  ps = ps.t,
  ml.ps = ps.t
)
dat.os <- list(Y = Y.os, A = A.os, X = X.os, q = rep(1, m))
# begin the elastic integrative analysis
# choice of kappa_n, default is tau = 1 (similar to BIC criteria)
tau <- 1
result.elastic <- elasticHTE(dat.t = dat.t,
                             dat.os = dat.os,
                             thres.psi = tau * sqrt(log(m)), # threshold for ACI psi
                             fixed = FALSE)
```

    #> est
    #>  ee.rt(ml).1  ee.rt(ml).2  ee.rt(ml).3 opt.ee(ml).1 opt.ee(ml).2 opt.ee(ml).3 
    #>  -0.06010648   0.92470952   0.96831571   0.01705285   1.03611492   1.00334316 
    #>       elas.1       elas.2       elas.3 
    #>   0.01705285   1.03611492   1.00334316
    #> ve
    #>  ee.rt(ml).1  ee.rt(ml).2  ee.rt(ml).3 opt.ee(ml).1 opt.ee(ml).2 opt.ee(ml).3 
    #>  0.027115933  0.017241212  0.016890441  0.005879353  0.004517951  0.003443317 
    #>       elas.1       elas.2       elas.3 
    #>  0.005813113  0.004243598  0.003365946
    #> CI
    #>                   lower     upper
    #> ee.rt(ml).1  -0.3828521 0.2626391
    #> ee.rt(ml).2   0.6673549 1.1820642
    #> ee.rt(ml).3   0.7135925 1.2230390
    #> opt.ee(ml).1 -0.1332312 0.1673369
    #> opt.ee(ml).2  0.9043746 1.1678553
    #> opt.ee(ml).3  0.8883329 1.1183534
    #> elas.1       -0.1323822 0.1664879
    #> elas.2        0.9084372 1.1637927
    #> elas.3        0.8896324 1.1170539

## More examples

``` r
example(ElasticIntegrative::elasticHTE)
```

``` r
utils::browseVignettes("ElasticIntegrative")
# remember to use build_readme() after every update
```
