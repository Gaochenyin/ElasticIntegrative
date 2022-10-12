
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

## Vignettes

``` r
utils::browseVignettes("ElasticIntegrative")
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
n = 1e5; mean.x = 1 # the whole population
m = 2000; tlocalpar = 0 # size and bias for the RWE
equal.weight = T # size for the RCT (NULL uses the default)

## begin the data generation
X1 <- rnorm(n, mean.x, 1)
X2 <- rnorm(n, mean.x, 1)
X3 <- rnorm(n, mean.x, 1) 
X1.temp <- rnorm(n, mean.x, 1)
X2.temp <- rnorm(n, mean.x, 1)
X3.temp <- rnorm(n, mean.x, 1) 
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
# true ATE
psi.pop <- lm((Y1 - Y0) ~ f.X)$coeff
true.tau <- mean(Y1) - mean(Y0)
taui <- Y1 - Y0

## RCT data
eS <- expoit(-4.5 - 2 * X1 - 2 * X2)
S <- sapply(eS, rbinom, n = 1, size = 1)
S.ind <- which(S == 1)
n.t <- length(S.ind)
### choose the selected patients
X.t <- f.X[S.ind, ]
Y0.t <- Y0[S.ind]
Y1.t <- Y1[S.ind]
eS.t <- eS[S.ind]
# assign treatment to trial participant with equal probabilities
ps.t <- rep(0.5, n.t)
A.t <- rbinom(n.t, 1, ps.t)
Y.t <- Y1.t * A.t + Y0.t * (1 - A.t)

## OS data
P.ind <- sample(1:n, size = m)
X.os <- f.Xtemp[P.ind, ]
X.confounder.os <- X3.temp[P.ind]
Y0.os <- Y0temp[P.ind]
Y1.os <- Y1temp[P.ind]
# maintain the proportion of the treated at about .5
a.opt <- uniroot(function(a) {
  ps <- expoit(a - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)
  mean(ps) - .5
}, c(-100, 100))$root
eA <- expoit(a.opt - 1 * X.os[, 1] + 1 * X.os[, 2] - tlocalpar * X.confounder.os)
A.os <- sapply(eA, rbinom, n = 1, size = 1)
Y.os <- Y1.os * A.os + Y0.os * (1 - A.os)
# sort the RT and RW datasets
dat.t <- list(
  Y = Y.t, A = A.t, X = X.t, q = rep(1, n.t),
  ps = ps.t,
  ml.ps = ps.t
)
dat.os <- list(Y = Y.os, A = A.os, X = X.os, q = rep(1, m))
```

## More examples

``` r
example(ElasticIntegrative::elasticHTE)
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
