# Cross-sectional sequential combination-reconciliation

This function performs a two-step process designed to first combine
forecasts from multiple models or experts and then apply reconciliation
techniques to ensure coherence.

## Usage

``` r
csscr(base, fc = "sa", comb = "ols", res = NULL, mse = TRUE, shrink = TRUE,
      nnw = FALSE, factorized = FALSE, ...)
```

## Arguments

- base:

  A list of \\p\\ numeric (\\h \times n\\) matrix or multivariate time
  series (`mts` class) containing the base forecasts to be reconciled;
  \\h\\ is the forecast horizon, \\n\\ is the total number of time
  series (\\n = n_u + n_b\\) and \\p\\ is the total number of experts.

- fc:

  A string specifying the combination method:

  - "`sa`" - (*default*) simple average (equal weights).

  - "`var`" - (uses `res`) weights derived from the inverse of forecasts
    variances/MSE as proposed by Bates and Granger (1969).

  - "`cov`" - (uses `res`) weights derived using the whole forecast
    error covariance matrix, as proposed by Newbold and Granger (1974).

- comb:

  A string specifying the reconciliation method: `"ols"`, `"wls"`,
  `"shr"`, `"sam"` (see [`FoReco`](https://danigiro.github.io/FoReco/)).
  If `comb = "none"`, no reconciliation is performed and the combined
  forecasts are directly returned.

- res:

  A list of \\p\\ numeric (\\N \times n\\) matrix containing the
  in-sample residuals. This input is used to compute some covariance
  matrices.

- mse:

  If `TRUE` (*default*) the residuals used to compute the covariance
  matrix are not mean-corrected.

- shrink:

  If `TRUE` (*default*), the covariance matrix for `fc = "cov"` is
  shrunk.

- nnw:

  If `TRUE` for `fc = "cov"`, the weights are constrained to be
  non-negative (Conflitti et al., 2015). The *default* is `FALSE`.

- factorized:

  Value to be passed to the
  [`quadprog::solve.QP`](https://rdrr.io/pkg/quadprog/man/solve.QP.html),
  only when `nnw = TRUE`.

- ...:

  Arguments passed on to
  [[`FoReco::csrec`](https://danigiro.github.io/FoReco/reference/csrec.html)](https://danigiro.github.io/FoReco/reference/csrec.html)
  (e.g., `agg_mat` or `cons_mat`).

## Value

A (\\h \times n\\) numeric matrix of cross-sectional combined and
reconciled forecasts.

## References

Bates, J. and Granger, C. W. (1969), The combination of forecasts,
*Operations Research Quarterly*, 20, 451–468.
[doi:10.1057/jors.1969.103](https://doi.org/10.1057/jors.1969.103) .

Conflitti, C., De Mol, C., and Giannone, D. (2015), Optimal combination
of survey forecasts. *International Journal of Forecasting*, 31(4),
1096–1103.
[doi:10.1016/j.ijforecast.2015.03.009](https://doi.org/10.1016/j.ijforecast.2015.03.009)
.

Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination
for linearly constrained multiple time series.
[doi:10.48550/arXiv.2412.03429](https://doi.org/10.48550/arXiv.2412.03429)
.

Newbold, P. and Granger, C. W. (1974), Experience with forecasting
univariate time series and the combination of forecasts, *Journal of the
Royal Statistical Society, A*, 137, 131–146.
[doi:10.2307/2344546](https://doi.org/10.2307/2344546)

## See also

Sequential coherent combination:
[`cssrc()`](https://danigiro.github.io/FoCo2/reference/cssrc.md)

## Examples

``` r
set.seed(123)
# (2 x 3) base forecasts matrix (simulated), expert 1
base1 <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
# (10 x 3) in-sample residuals matrix (simulated), expert 1
res1 <- t(matrix(rnorm(n = 30), nrow = 3))

# (2 x 3) base forecasts matrix (simulated), expert 2
base2 <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
# (10 x 3) in-sample residuals matrix (simulated), expert 2
res2 <- t(matrix(rnorm(n = 30), nrow = 3))

# Base forecasts' and residuals' lists
base <- list(base1, base2)
res <- list(res1, res2)

# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
reco <- csscr(base = base, agg_mat = A, comb = "wls", res = res, fc = "sa")

# Zero constraints matrix for Z - X - Y = 0
C <- t(c(1, -1, -1))
reco <- csscr(base = base, cons_mat = C, comb = "wls", res = res, fc = "sa") # same results

# Incoherent combined forecasts
fc_comb <- csscr(base = base, comb = "none", fc = "sa")
round(C %*% t(fc_comb), 3) # Incoherent forecasts
#>         h-1    h-2
#> [1,] -0.484 -0.626
```
