
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

This is an example for illustration when
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
is a continuous outcome

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
    #>  0.030100673  0.017169858  0.020132153  0.030455816  0.017169809  0.020936258 
    #> opt.ee(ml).1 opt.ee(ml).2 opt.ee(ml).3       elas.1       elas.2       elas.3 
    #>  0.006598755  0.004010907  0.003080497  0.005862381  0.004043671  0.003188111
    #> CI
    #>                    lower     upper
    #> covj.t.1     -0.31152951 0.3685602
    #> covj.t.2      0.74769614 1.2613392
    #> covj.t.3      0.66501468 1.2212047
    #> ee.rt(ml).1  -0.28642214 0.3976678
    #> ee.rt(ml).2   0.74398645 1.2576288
    #> ee.rt(ml).3   0.68895777 1.2561465
    #> opt.ee(ml).1 -0.08279945 0.2356270
    #> opt.ee(ml).2  0.82375983 1.0720156
    #> opt.ee(ml).3  0.94267838 1.1602431
    #> elas.1       -0.13145638 0.2803858
    #> elas.2        0.86455094 1.1823684
    #> elas.3        0.84800403 1.1233989

## Extension

This is additional example for illustration when
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
is a binary outcome

``` r
library(ElasticIntegrative)
library(SuperLearner)
## basic example code for data generation
set.seed(2333)
## setups
p <- 3; beta0 <- c(-1, 1, 1, 1) # for the mu0 function
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
h.X.rt <- cbind(X1 = X1.rt, X2 = X2.rt, X3 = X3.rt)
h.X.os <- cbind(X1 = X1.os, X2 = X2.os, X3 = X3.os)
f.X.rt <- cbind(X1 = X1.rt, X2 = X2.rt)
f.X.os <- cbind(X1 = X1.os, X2 = X2.os)
# -------------------------------
## construct Y for RCT and RWE populations
### RCT population
prob.A0.rt <- exp(cbind(1, h.X.rt) %*% beta0)/
         (1 + exp(cbind(1, h.X.rt) %*% beta0))
prob.A1.rt <- exp(cbind(1, h.X.rt) %*% beta0 + cbind(1, f.X.rt) %*% psi0)/
         (1 + exp(cbind(1, h.X.rt) %*% beta0 + cbind(1, f.X.rt) %*% psi0))
Y0.rt <- rbinom(n = n, size = 1, 
       prob = prob.A0.rt)
Y1.rt <- rbinom(n = n, size = 1,
       prob = prob.A1.rt)
truth <- mean(prob.A1.rt - prob.A0.rt) # truth is approximated numerically
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
Y0.os <- rbinom(n = n, size = 1,
                prob = exp(cbind(1, h.X.os) %*% beta0)/
                  (1 + exp(cbind(1, h.X.os) %*% beta0)))
Y1.os <- rbinom(n = n, size = 1,
                prob = exp(cbind(1, h.X.os) %*% beta0 + cbind(1, f.X.os) %*% psi0)/
                  (1 + exp(cbind(1, h.X.os) %*% beta0 + cbind(1, f.X.os) %*% psi0)))
### OS sample
P.ind <- sample(1:n, size = m)
X.os <- f.X.os[P.ind, ]
X.confounder.os <- X3.os[P.ind]
Y0.os <- Y0.os[P.ind]
Y1.os <- Y1.os[P.ind]
# a.opt is chosen to maintain the proportion of the treated at about .5
a.opt <- uniroot(function(a) {
  ps <- exp(a - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)/
    {1+exp(a - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)}
  mean(ps) - .5
}, c(-100, 100))$root
eA <- exp(a.opt - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)/
  {1+exp(a.opt - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)}
A.os <- sapply(eA, rbinom, n = 1, size = 1)
Y.os <- Y1.os * A.os + Y0.os * (1 - A.os)
# organize the RT and RW datasets
dat.t <- data.frame(
  Y = Y.rt, A = A.rt, X.rt, q = rep(1, n.t),
  ps = ps.t,
  ml.ps = ps.t
)

dat.os <- data.frame(
  Y = Y.os, A = A.os, X.os, q = rep(1, m)
  )
# begin the elastic integrative analysis
# choice of kappa_n, default is sqrt(log(m)) (similar to the BIC criteria)
# ignore the warnings for using weighted logistic regression
options(warn=-1)
result.elastic <- elasticHTE(mainName = c('X1', 'X2'),
                             contName = c('X1', 'X2'),
                             propenName = c('X1', 'X2'),
                             dat.t = dat.t,
                             dat.os = dat.os,
                             family.Y = binomial(),
                             thres.psi = sqrt(log(m)), # threshold for ACI psi
                             fixed = FALSE,
                             cvControl = list(V = 1L))
```

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
