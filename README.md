
# PGLMMdredge

<!-- badges: start -->
<!-- badges: end -->

The goal of PGLMMdredge is to adaptat the dredge function from 'MuMIn' to phylogenetic generalized linear mixed models (pglmm from 'phyr'). It generates a model selection table of models as would do the dredge function from 'MuMIn' with combinations (subsets) of fixed effect terms in the global model, with optional model inclusion rules. 

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

This is a basic example using simulated package data (see documentation for the 'senescence' dataset and 'tree_ultra' ultrametric phylogenetic tree)

``` r
# Write the null model using pglmm() from 'phyr'
mod1 <- pglmm(log_onset ~ 1 + (1|Species__),
             data = senescence, family = "gaussian", cov_ranef = list(Species = tree_ultra),
             REML = TRUE, verbose = FALSE, s2.init = .1)

# Returns a model selection table based on the null model and evaluating all possible models based on the fixed effect structure provided
dredge_pglmm(
        formulaRE = "log_onset ~ 1 + (1 | Species)",
        fixed = c("log_afr", "log_mass"),
        data = senescence,
        rank = "AICc",
        family = "gaussian",
        cov_ranef = list(Species = tree_ultra)
)
```
## Citation

To cite the ‘PGLMMdredge’ package in your publications, please use:

  Lalande LD (2025). _PGLMMdredge: A model selection for pglmm_. R package version 0.1.0,
  <https://github.com/LucasDLalande/PGLMMdredge.git>.
