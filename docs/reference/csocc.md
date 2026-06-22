# Cross-sectional optimal coherent forecast combination

This function implements the optimal linear coherent forecast
combination for a linearly constrained (e.g., hierarchical/grouped)
multiple time series combining forecasts from multiple sources
(experts). The method minimizes forecast errors while ensuring that the
resulting forecasts are coherent, meaning they satisfy the constraints
across the time series. Using a linear constrained optimization
approach, csocc combines forecasts from multiple experts, accounting for
the relations between variables and constraints. In addition, linear
inequality constraints (e.g. non-negative forecasts) can be imposed, if
needed.

## Usage

``` r
csocc(base, agg_mat, cons_mat, comb = "ols", res = NULL,
      approach = "proj", nn = NULL, settings = NULL, bounds = NULL, ...)
```

## Arguments

- base:

  A list of \\p\\ numeric (\\h \times n\\) matrix or multivariate time
  series (`mts` class) containing the base forecasts to be reconciled;
  \\h\\ is the forecast horizon, \\n\\ is the total number of time
  series (\\n = n_u + n_b\\) and \\p\\ is the total number of experts.

- agg_mat:

  A (\\n_u \times n_b\\) numeric matrix representing the cross-sectional
  aggregation matrix, mapping the \\n_b\\ bottom-level (free) variables
  into the \\n_u\\ upper (constrained) variables.

- cons_mat:

  A (\\n_u \times n\\) numeric matrix representing the cross-sectional
  zero constraints: each row represents a constraint equation, and each
  column represents a variable. The matrix can be of full rank, meaning
  the rows are linearly independent, but this is not a strict
  requirement, as the function allows for redundancy in the constraints.

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

  - "`strc`": structural approach.

- nn:

  A string specifying the algorithm to compute non-negative forecasts:

  - "`osqp`": OSQP solver (Stellato et al., 2020).

  - "`sntz`": heuristic "set-negative-to-zero" (Di Fonzo and
    Girolimetto, 2023).

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

- ...:

  Arguments passed on to
  [cscov](https://danigiro.github.io/FoCo2/reference/cscov.md).

## Value

A (\\h \times n\\) numeric matrix of cross-sectional combined and
reconciled forecasts.

## Details

If an expert does not provide a forecast for a given variable, the
missing value can be represented as `NA`. For instance, consider a case
where the number of experts is \\p = 4\\, with \\n = 3\\ variables and
\\h = 2\\. If the second and third experts do not provide forecasts for,
respectively, the second and first variable, the corresponding matrix in
the base list will be as follows: \$\$\texttt{base\[\[1\]\]} =
\begin{bmatrix} \hat{y}\_{1,1}^{1} & \hat{y}\_{2,1}^{1} &
\hat{y}\_{3,1}^{1} \\ \hat{y}\_{1,2}^{1} & \hat{y}\_{2,2}^{1} &
\hat{y}\_{3,2}^{1} \end{bmatrix}, \quad \texttt{base\[\[2\]\]} =
\begin{bmatrix} \hat{y}\_{1,1}^{2} & \texttt{NA} & \hat{y}\_{3,1}^{2} \\
\hat{y}\_{1,2}^{2} & \texttt{NA} & \hat{y}\_{3,2}^{2} \end{bmatrix},\$\$
\$\$\texttt{base\[\[3\]\]} = \begin{bmatrix} \texttt{NA} &
\hat{y}\_{2,1}^{3} & \hat{y}\_{3,1}^{3} \\ \texttt{NA} &
\hat{y}\_{2,2}^{3} & \hat{y}\_{3,2}^{3} \end{bmatrix}, \quad
\texttt{base\[\[4\]\]} = \begin{bmatrix} \hat{y}\_{1,1}^{4} &
\hat{y}\_{2,1}^{4} & \hat{y}\_{3,1}^{4} \\ \hat{y}\_{1,2}^{4} &
\hat{y}\_{2,2}^{4} & \hat{y}\_{3,2}^{4} \end{bmatrix}.\$\$

## References

Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation
of solar forecasts, *Solar Energy*, 251, 13–29.
[doi:10.1016/j.solener.2023.01.003](https://doi.org/10.1016/j.solener.2023.01.003)

Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination
for linearly constrained multiple time series,
[doi:10.48550/arXiv.2412.03429](https://doi.org/10.48550/arXiv.2412.03429)
.

Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020),
OSQP: An Operator Splitting solver for Quadratic Programs, *Mathematical
Programming Computation*, 12, 4, 637-672.
[doi:10.1007/s12532-020-00179-2](https://doi.org/10.1007/s12532-020-00179-2)

## See also

Other Optimal combination:
[`cscov()`](https://danigiro.github.io/FoCo2/reference/cscov.md),
[`csmtc()`](https://danigiro.github.io/FoCo2/reference/csmtc.md),
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
rrc <- csocc(base = brc, agg_mat = A, comb = "wls", res = erc)

# Zero constraints matrix for Z - X - Y = 0
C <- t(c(1, -1, -1))
rrc <- csocc(base = brc, cons_mat = C, comb = "wls", res = erc) # same results

## UNBALANCED PANEL OF FORECASTS
base2[, 2] <- res2[, 2] <-  NA

# Base forecasts' and residuals' lists
bgc <- list(base1, base2)
egc <- list(res1, res2)

# Aggregation matrix for Z = X + Y
A <- t(c(1,1))
rgc <- csocc(base = bgc, agg_mat = A, comb = "shrbe", res = egc)
```
