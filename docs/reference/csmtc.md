# Cross-sectional optimal multi-task forecast combination

This function computes the optimal multi-task linear forecast
combination, as described in Girolimetto and Di Fonzo (2024)

## Usage

``` r
csmtc(base, comb = "ols", res = NULL, approach = "proj",
      nn = NULL, settings = NULL, bounds = NULL, agg_mat = NULL, ...)
```

## Arguments

- base:

  A list of \\p\\ numeric (\\h \times n\\) matrix or multivariate time
  series (`mts` class) containing the base forecasts to be reconciled;
  \\h\\ is the forecast horizon, \\n\\ is the total number of time
  series (\\n = n_u + n_b\\) and \\p\\ is the total number of experts.

- comb:

  A string specifying the reconciliation method. For details, see
  [cscov](https://danigiro.github.io/FoCo2/reference/cscov.md).

- res:

  A list of \\p\\ numeric (\\N \times n\\) matrix containing the
  in-sample residuals. This input is used to compute some covariance
  matrices.

- approach:

  A string specifying the approach used to compute the reconciled
  forecasts. Options include:

  - "`proj`" (*default*): zero-constrained projection approach.

  - "`osqp`": OSQP solver (Stellato et al., 2020).

- nn:

  A string specifying the algorithm to compute non-negative forecasts:

  - "`osqp`": OSQP solver (Stellato et al., 2020).

  - "`sntz`": heuristic "set-negative-to-zero".

- settings:

  An object of class `osqpSettings` specifying settings for the
  [osqp](https://osqp.org/) solver. For details, refer to the [osqp
  documentation](https://osqp.org/) (Stellato et al., 2020).

- bounds:

  A matrix (see `set_bounds` in
  [FoReco](https://danigiro.github.io/FoReco/reference/set_bounds.html))
  with 3 columns (\\i,lower,upper\\), such that

  - Column 1 represents the cross-sectional series (\\i = 1, \dots,
    n\\).

  - Column 2 indicates the *lower* bound.

  - Column 3 indicates the *upper* bound.

- agg_mat:

  A (\\n_u \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix, mapping the \\n_b\\ bottom-level (free) variables
  into the \\n_u\\ upper (constrained) variables.

- ...:

  Arguments passed on to
  [cscov](https://danigiro.github.io/FoCo2/reference/cscov.md).

## Value

A (\\h \times n\\) numeric matrix of cross-sectional multi-task combined
forecasts.

## References

Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination
for linearly constrained multiple time series,
[doi:10.48550/arXiv.2412.03429](https://doi.org/10.48550/arXiv.2412.03429)
.

## See also

Other Optimal combination:
[`cscov()`](https://danigiro.github.io/FoCo2/reference/cscov.md),
[`csocc()`](https://danigiro.github.io/FoCo2/reference/csocc.md),
[`occmat()`](https://danigiro.github.io/FoCo2/reference/occmat.md)

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

## BALANCED PANEL OF FORECASTS
# Base forecasts' and residuals' lists
brc <- list(base1, base2)
erc <- list(res1, res2)

# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
rrc <- csocc(base = brc, agg_mat = A, comb = "shr", res = erc)

yc <- csmtc(base = brc, agg_mat = A, comb = "shr", res = erc)
M <- occmat(base = brc, agg_mat = A, comb = "shr", p = 2, res = erc)$M
M%*%t(yc)-t(rrc)
#> 3 x 2 Matrix of class "dgeMatrix"
#>               h-1 h-2
#> [1,] 0.000000e+00   0
#> [2,] 1.776357e-15   0
#> [3,] 0.000000e+00   0

## UNBALANCED PANEL OF FORECASTS
base2[, 2] <- res2[, 2] <-  NA

# Base forecasts' and residuals' lists
bgc <- list(base1, base2)
egc <- list(res1, res2)
matNA <- matrix(1, 3, 2)
matNA[2,2] <- 0

# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
rgc <- csocc(base = bgc, agg_mat = A, comb = "shr", res = egc)

yc <- csmtc(base = bgc, agg_mat = A, comb = "shr", res = egc)
M <- occmat(base = bgc, agg_mat = A, comb = "shr", p = 2, res = egc, matNA = matNA)$M
M%*%t(yc)-t(rgc)
#> 3 x 2 Matrix of class "dgeMatrix"
#>               h-1 h-2
#> [1,] 0.000000e+00   0
#> [2,] 1.776357e-15   0
#> [3,] 0.000000e+00   0
```
