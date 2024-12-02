
# FoCo² <img src="man/figures/logo.svg" title="S'i' fosse foco, arderei 'l mondo (Sonetti, 86) Cecco Angiolieri, Italian poet." alt="logo" align="right" width="150" style="border: none; float: right;"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/FoCo2)](https://CRAN.R-project.org/package=FoCo2)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![devel
version](https://img.shields.io/badge/devel%20version-0.1.0.001-blue.svg)](https://github.com/daniGiro/FoCo2)
[![License:
GPL-3](https://img.shields.io/badge/license-GPL--3-forestgreen.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![R-CMD-check](https://github.com/danigiro/FoCo2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/danigiro/FoCo2/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> *S’i’ fosse foco, arderei ’l mondo* (Sonetti, 86) Cecco Angiolieri,
> Italian poet.

**FoCo²** (**Co**herent **Fo**recast **Co**mbination) is a forecasting
package designed to handle multiple time series forecasts from different
experts, subject to linear constraints. It offers both optimal and
heuristic methods for combining expert forecasts and reconciling them
through a multi-task approach. This process either simultaneously (in
the optimal case) or sequentially (in the heuristic cases) integrates
forecasts from multiple experts while incorporating *a priori*
constraints to produce coherent forecasts.

The package is designed to manage different levels of input complexity:
for example, some experts might produce forecasts for different subsets
of the variables. Whether the input consists of forecasts for all
variables from each expert (balanced case) or partial (unbalanced case)
forecasts for a subset of variables, **FoCo²** efficiently organizes the
data and applies combination and reconciliation methods. It also
includes advanced tools for handling forecast errors. The package
provides functions to estimate the forecast error covariance matrix,
supporting simple, (block) diagonal (assuming uncorrelation) and full
covariance matrices.

## Installation

You can install the **stable** version on [R
CRAN](https://cran.r-project.org/)

``` r
install.packages("FoCo2")
```

You can also install the **development** version from
[Github](https://github.com/daniGiro/FoCo2)

``` r
# install.packages("devtools")
devtools::install_github("danigiro/FoCo2")
```

<!-- ## Other forecast reconciliation package -->
<!-- [**FoReco**](https://github.com/daniGiro/FoReco) <img src="man/figures/foreco.svg" alt="foreco" width="75" height="75" style="border: none; float: left;"/>\ -->
<!-- Classical (bottom-up and top-down), optimal combination and heuristic point ([Di Fonzo and Girolimetto, 2023](https://doi.org/10.1016/j.ijforecast.2021.08.004)) and probabilistic ([Girolimetto et al. 2023](https://doi.org/10.1016/j.ijforecast.2023.10.003)) forecast reconciliation procedures for linearly constrained time series (e.g., hierarchical or grouped time series) in cross-sectional, temporal, or cross-temporal frameworks. -->

## Issues and Contributions

If you encounter any bugs or have suggestions for improvements, please
feel free to report them on [GitHub Issues
page](https://github.com/daniGiro/FoCo2/issues). Contributions are also
welcome!
