# Matrices for the optimal coherent forecast combination

This function computes the matrices required for the optimal coherent
forecast combination
[csocc](https://danigiro.github.io/FoCo2/reference/csocc.md), as
described in Girolimetto and Di Fonzo (2024). These matrices serve as
the foundation for building forecasts that effectively combines the
individual information from multiple experts while ensuring coherence
across the variables.

## Usage

``` r
occmat(agg_mat, cons_mat, p = NULL, matNA = NULL,
       comb = "ols", res = NULL, approach = "proj", ...)
```

## Arguments

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

- p:

  Total number of experts, \\p\\.

- matNA:

  A (\\n \times p\\) matrix consisting of 0s and 1s, where each element
  indicates whether expert \\j\\ (column) has provided a forecast for
  variable \\i\\ (row). If expert \\j\\ has provided a forecast for
  variable \\i\\, the corresponding element \\(i,j)\\ is 1; otherwise,
  it is 0.

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

- ...:

  Arguments passed on to
  [cscov](https://danigiro.github.io/FoCo2/reference/cscov.md).

## Value

A list of matrices:

- M:

  Projection matrix.

- Omega:

  Matrix of the combination weights of the optimal linear multi-task
  forecast combination.

- W:

  Forecast error covariance matrix of the base forecasts.

- Wc:

  Forecast error covariance matrix of the combined forecasts.

- Wtilde:

  Forecast error covariance matrix of the reconciled combined forecasts.

- K:

  Matrix that replicates a vector (see Girolimetto and Di Fonzo, 2024).

## References

Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination
for linearly constrained multiple time series,
[doi:10.48550/arXiv.2412.03429](https://doi.org/10.48550/arXiv.2412.03429)
.

## See also

Other Optimal combination:
[`cscov()`](https://danigiro.github.io/FoCo2/reference/cscov.md),
[`csmtc()`](https://danigiro.github.io/FoCo2/reference/csmtc.md),
[`csocc()`](https://danigiro.github.io/FoCo2/reference/csocc.md)
