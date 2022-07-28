
# SimInfLME

<!-- badges: start -->
<!-- badges: end -->

## Overview

Package **SimInfLME** implements methods to construct individual and
simultaneous confidence intervals for linear mixed effect using large
sample approximations, parametric and semiparametric bootstrap.

Reference: Reluga and Sperlich (2022). *Simple bootstrap for linear
mixed effects under model misspecification*.

Available at [arXiv.org](https://arxiv.org/abs/2207.12455).

## Installation

You can install the most recent version of **SimInfLME** from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")

remotes::install_github("KatarzynaReluga/postcAIC")
remotes::install_github("KatarzynaReluga/SimInfLME")
```

## Example

This is a basic example which shows you how to construct individual and
simultaneous confidence intervals for linear mixed effects under
nested-error regression model (NERM).

First we fix parameters and generate covariates.

``` r
library(SimInfLME)

# Set seed ----------------------------------------------------------------------
set.seed(11992)

# Set the number of clusters and units in each cluster -------------------------
m = 25
n_j = 5
n = m * n_j

# Set fixed (regression) parameters and variance parameters --------------------
scaling_factor_e = 1
scaling_factor_u = sqrt(0.5)
x_ij = runif(n, 0, 1)
X = cbind(X0 = rep(1, n), X1 = x_ij)
t_X  = t(X)
id_cluster = rep(1:m, each = n_j)
beta = c(1, 1)

# Create a factor vector with cluster labels -----------------------------------
cluster_means <- as.matrix(aggregate(X, list(id_cluster), FUN = mean)[, -1])

# Formula to fit NERM ----------------------------------------------------------
formula_y <-  y ~ -1 + X0 + X1 + (1| id_cluster)
```

Second, we generate outcomes from NERM.

``` r
data_sample <- generate_NERM(generate_u = list(type = "chisquared",
                                               scaling_factor = 1,
                                               dg = 6),
                             generate_e = list(type = "chisquared",
                                               scaling_factor = 1,
                                               dg = 6),
                             beta = beta,
                             X = X,
                             id_cluster = id_cluster,
                             start_seed = 1,
                             no_sim = 1,
                             cluster_means = cluster_means)
```

Third, we compute individual and simultaneous intervals for linear mixed
effects.

``` r
intervals_lme <- construct_intervals(type_method = c("semiparametric"),
                                    type_var_estimator = c("mse_b"),
                                    model_type = c("NERM"),
                                    formula_y, data_sample,
                                    cluster_means = cluster_means,
                                    id_cluster = id_cluster,
                                    alpha = 0.05,
                                    n_boot = 500,
                                    boot_seed = 1)
```

## Simulation study from Reluga and Sperlich (2022)

All
[simulations](https://github.com/KatarzynaReluga/SimInfLME/tree/main/simulations)
in [Reluga and Sperlich (2022)](https://arxiv.org/abs/2207.12455) can be
reproduced using functions in package **SimInfLME**.
