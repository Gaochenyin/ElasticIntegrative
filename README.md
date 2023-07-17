
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
