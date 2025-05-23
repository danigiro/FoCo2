% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimal.R
\name{csocc}
\alias{csocc}
\title{Cross-sectional optimal coherent forecast combination}
\usage{
csocc(base, agg_mat, cons_mat, comb = "ols", res = NULL,
      approach = "proj", nn = NULL, settings = NULL, bounds = NULL, ...)
}
\arguments{
\item{base}{A list of \eqn{p} numeric (\eqn{h \times n}) matrix or multivariate time series
(\code{mts} class) containing the base forecasts to be reconciled; \eqn{h} is the forecast
horizon, \eqn{n} is the total number of time series (\eqn{n = n_u + n_b}) and
\eqn{p} is the total number of experts.}

\item{agg_mat}{A (\eqn{n_u \times n_b}) numeric matrix representing the cross-sectional
aggregation matrix, mapping the \eqn{n_b} bottom-level (free)
variables into the \eqn{n_u} upper (constrained) variables.}

\item{cons_mat}{A (\eqn{n_u \times n}) numeric matrix representing the cross-sectional
zero constraints: each row represents a constraint equation, and each column represents
a variable. The matrix can be of full rank, meaning the rows are linearly independent,
but this is not a strict requirement, as the function allows for redundancy in the constraints.}

\item{comb}{A string specifying the reconciliation method. For details, see \link{cscov}.}

\item{res}{A list of \eqn{p} numeric (\eqn{N \times n}) matrix containing the
in-sample residuals. This input is used to compute some covariance matrices.}

\item{approach}{A string specifying the approach used to compute the reconciled
forecasts. Options include:
\itemize{
\item "\code{proj}" (\emph{default}): zero-constrained projection approach.
\item "\code{strc}": structural approach.
}}

\item{nn}{A string specifying the algorithm to compute non-negative forecasts:
\itemize{
\item "\code{osqp}": OSQP solver (Stellato et al., 2020).
\item "\code{sntz}": heuristic "set-negative-to-zero" (Di Fonzo and Girolimetto, 2023).
}}

\item{settings}{An object of class \code{osqpSettings} specifying settings
for the \href{https://osqp.org/}{\pkg{osqp}} solver. For details, refer to the
\href{https://osqp.org/}{\pkg{osqp} documentation} (Stellato et al., 2020).}

\item{bounds}{A matrix (see \code{set_bounds} in \href{https://danigiro.github.io/FoReco/reference/set_bounds.html}{\pkg{FoReco}})
with 3 columns (\eqn{i,lower,upper}), such that
\itemize{
\item Column 1 represents the cross-sectional series (\eqn{i = 1, \dots, n}).
\item Column 2 indicates the \emph{lower} bound.
\item Column 3 indicates the \emph{upper} bound.
}}

\item{...}{Arguments passed on to \link{cscov}.}
}
\value{
A (\eqn{h \times n}) numeric matrix of cross-sectional combined and reconciled forecasts.
}
\description{
This function implements the optimal linear coherent forecast combination for a
linearly constrained (e.g., hierarchical/grouped) multiple time series
combining forecasts from multiple sources (experts). The method minimizes forecast
errors while ensuring that the resulting forecasts are coherent, meaning
they satisfy the constraints across the time series. Using a linear constrained
optimization approach, \link{csocc} combines forecasts from multiple experts,
accounting for the relations between variables and constraints. In addition, linear
inequality constraints (e.g. non-negative forecasts) can be imposed, if needed.
}
\details{
If an expert does not provide a forecast for a given variable,
the missing value can be represented as \code{NA}.
For instance, consider a case where the number of experts is
\eqn{p = 4}, with \eqn{n = 3} variables and \eqn{h = 2}. If the
second and third experts do not provide forecasts for, respectively, the second and first variable,
the corresponding matrix in the base list will be as follows:
\deqn{\texttt{base[[1]]} = \begin{bmatrix} \hat{y}_{1,1}^{1} & \hat{y}_{2,1}^{1} & \hat{y}_{3,1}^{1} \\
\hat{y}_{1,2}^{1} & \hat{y}_{2,2}^{1} & \hat{y}_{3,2}^{1} \end{bmatrix}, \quad
\texttt{base[[2]]} = \begin{bmatrix} \hat{y}_{1,1}^{2} & \texttt{NA} & \hat{y}_{3,1}^{2} \\
\hat{y}_{1,2}^{2} & \texttt{NA} & \hat{y}_{3,2}^{2} \end{bmatrix},}
\deqn{\texttt{base[[3]]} = \begin{bmatrix} \texttt{NA} & \hat{y}_{2,1}^{3} & \hat{y}_{3,1}^{3} \\
\texttt{NA} & \hat{y}_{2,2}^{3} & \hat{y}_{3,2}^{3} \end{bmatrix}, \quad
\texttt{base[[4]]} = \begin{bmatrix} \hat{y}_{1,1}^{4} & \hat{y}_{2,1}^{4} & \hat{y}_{3,1}^{4} \\
\hat{y}_{1,2}^{4} & \hat{y}_{2,2}^{4} & \hat{y}_{3,2}^{4} \end{bmatrix}.}
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

}
\references{
Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation of solar
forecasts, \emph{Solar Energy}, 251, 13–29. \doi{10.1016/j.solener.2023.01.003}

Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination for linearly
constrained multiple time series, \emph{mimeo}.

Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020), OSQP:
An Operator Splitting solver for Quadratic Programs,
\emph{Mathematical Programming Computation}, 12, 4, 637-672.
\doi{10.1007/s12532-020-00179-2}
}
\seealso{
Other Optimal combination: 
\code{\link{cscov}()},
\code{\link{csmtc}()},
\code{\link{occmat}()}
}
\concept{Optimal combination}
