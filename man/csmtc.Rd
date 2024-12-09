% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimal.R
\name{csmtc}
\alias{csmtc}
\title{Cross-sectional optimal multi-task forecast combination}
\usage{
csmtc(base, agg_mat = NULL, comb = "ols", res = NULL, ...)
}
\arguments{
\item{base}{A list of \eqn{p} numeric (\eqn{h \times n}) matrix or multivariate time series
(\code{mts} class) containing the base forecasts to be reconciled; \eqn{h} is the forecast
horizon, \eqn{n} is the total number of time series (\eqn{n = n_u + n_b}) and
\eqn{p} is the total number of experts.}

\item{agg_mat}{A (\eqn{n_u \times n_b}) numeric matrix representing the cross-sectional
aggregation matrix, mapping the \eqn{n_b} bottom-level (free)
variables into the \eqn{n_u} upper (constrained) variables.}

\item{comb}{A string specifying the reconciliation method. For details, see \link{cscov}.}

\item{res}{A list of \eqn{p} numeric (\eqn{N \times n}) matrix containing the
in-sample residuals. This input is used to compute some covariance matrices.}

\item{...}{Arguments passed on to \link{cscov}.}
}
\value{
A (\eqn{h \times n}) numeric matrix of cross-sectional multi-task combined forecasts.
}
\description{
This function computes the optimal multi-task linear forecast combination, as described
in Girolimetto and Di Fonzo (2024)
}
\examples{
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
M\%*\%t(yc)-t(rrc)

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
M\%*\%t(yc)-t(rgc)

}
\references{
Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination for linearly
constrained multiple time series, \emph{mimeo}.
}
\seealso{
Other Optimal combination: 
\code{\link{cscov}()},
\code{\link{csocc}()},
\code{\link{occmat}()}
}
\concept{Optimal combination}
