# Cross-sectional covariance matrix approximation

Extended version of the
[[`FoReco::cscov`](https://danigiro.github.io/FoReco/reference/cscov.html)](https://danigiro.github.io/FoReco/reference/cscov.html)
function, introducing two new approximations for the covariance matrix
(both shrunk and sample versions). Specifically, `shrbe`/`sambe` assume
no correlation between experts, while `shrbv`/`sambv` assume no
correlation between variables.

## Usage

``` r
# S3 method for class 'shrbe'
cscov(comb = "shrbe", ..., n = NULL, p = NULL, matNA = NULL,
      res = NULL, mse = TRUE, shrink_fun = NULL)

# S3 method for class 'sambe'
cscov(comb = "sambe", ..., n = NULL, p = NULL, matNA = NULL,
      res = NULL, mse = TRUE)

# S3 method for class 'shrbv'
cscov(comb = "shrbv", ..., n = NULL, p = NULL, matNA = NULL,
      res = NULL, mse = TRUE, shrink_fun = NULL)

# S3 method for class 'sambv'
cscov(comb = "sambv", ..., n = NULL, p = NULL, matNA = NULL,
      res = NULL, mse = TRUE)
```

## Arguments

- comb:

  A string specifying the reconciliation method.

  - [`FoReco`](https://danigiro.github.io/FoReco/) approaches: `"ols"`,
    `"wls"`, `"shr"`, `"sam"`.

  - "`shrbe`"/"`sambe`" - shrunk/sample block-diagonal covariance by
    experts.

  - "`shrbv`"/"`sambv`" - shrunk/sample block-diagonal covariance by
    variables.

- ...:

  Arguments passed on to
  [[`FoReco::cscov`](https://danigiro.github.io/FoReco/reference/cscov.html)](https://danigiro.github.io/FoReco/reference/cscov.html).

- n:

  Total number of variables, \\n\\.

- p:

  Total number of experts, \\p\\.

- matNA:

  A (\\n \times p\\) matrix consisting of 0s and 1s, where each element
  indicates whether expert \\j\\ (column) has provided a forecast for
  variable \\i\\ (row). If expert \\j\\ has provided a forecast for
  variable \\i\\, the corresponding element \\(i,j)\\ is 1; otherwise,
  it is 0.

- res:

  A list of \\p\\ numeric (\\N \times n\\) matrix containing the
  in-sample residuals. This input is used to compute some covariance
  matrices.

- mse:

  If `TRUE` (*default*) the residuals used to compute the covariance
  matrix are not mean-corrected.

- shrink_fun:

  Shrinkage function of the covariance matrix,
  [[`FoReco::shrink_estim`](https://danigiro.github.io/FoReco/reference/shrink_estim.html)](https://danigiro.github.io/FoReco/reference/shrink_estim.html)
  (*default*).

## Value

A (\\m \times m\\) symmetric positive (semi-)definite matrix, with \\m =
\sum\_{j = 1}^p n_j\\, \\n_j \leq n\\.

## See also

Other Optimal combination:
[`csmtc()`](https://danigiro.github.io/FoCo2/reference/csmtc.md),
[`csocc()`](https://danigiro.github.io/FoCo2/reference/csocc.md),
[`occmat()`](https://danigiro.github.io/FoCo2/reference/occmat.md)
