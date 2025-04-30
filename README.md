
# PGLMMdredge

<!-- badges: start -->
<!-- badges: end -->

The goal of PGLMMdredge is to adaptat the dredge function from 'MuMIn' to phylogenetic generalized linear mixed models (pglmm from 'phyr'). It generates a model selection table of models based on the dredge function from 'MuMIn' with combinations (subsets) of fixed effect terms in the global model, with optional model inclusion rules. 

## Installation

You can install the development version of PGLMMdredge from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("LucasDLalande/PGLMMdredge")
```

## Usage

Define your random effects, fixed effects, data, model selection criterion, distribution family, and whether to display estimates and/or standard-errors.
``` r
dredge_pglmm(
  formulaRE,
  fixed,
  data,
  rank = c("AICc", "AIC"),
  family = c("gaussian", "binomial", "poisson"),
  cov_ranef,
  estimate = T,
  std.err = F,
  round = 3
)
```
It returns a data frame with ordered models and respective fixed effects structure (and estimates, standard-errors), based on AIC or AICc.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PGLMMdredge)
## basic example code using iris dataset
dredge_pglmm(
   formulaRE = "Sepal.Length ~ 1 + (1 | Species)",
   fixed = c("Sepal.Width", "Petal.Length", "Sepal.Width:Petal.Length"),
   data = iris,
   rank = "AIC",
   family = "gaussian"
 )
```

