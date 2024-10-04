#' Cross-sectional optimal coherent forecast combination
#'
#' This function implements the optimal linear coherent forecast combination for a
#' linearly constrained (e.g., hierarchical/grouped) multiple time series
#' combining forecasts from multiple sources (experts). The method minimizes forecast
#' errors while ensuring that the resulting forecasts are coherent, meaning
#' they satisfy the constraints across the time series. Using a regression-based
#' optimization approach, it combines forecasts from multiple experts,
#' accounting for the relations between variables and constraints.
#'
#' @usage csocc(base, agg_mat, cons_mat, comb = "ols", res = NULL,
#'       approach = "proj", nn = NULL, settings = NULL, bounds = NULL, ...)
#'
#' @param base A list of \eqn{p} numeric (\eqn{h \times n}) matrix or multivariate time series
#' (\code{mts} class) containing the base forecasts to be reconciled; \eqn{h} is the forecast
#' horizon, \eqn{n} is the total number of time series (\eqn{n = n_u + n_b}) and
#' \eqn{p} is the total number of experts.
#' @param agg_mat A (\eqn{n_u \times n_b}) numeric matrix representing the cross-sectional
#' aggregation matrix. It maps the \eqn{n_b} bottom-level (free)
#' variables into the \eqn{n_u} upper (constrained) variables.
#' @param cons_mat A (\eqn{n_u \times n}) numeric matrix representing the cross-sectional
#' zero constraints. It spans the null space for the reconciled forecasts.
#' @param comb A string specifying the reconciliation method. For details, see [cscov].
#' @param res A list of \eqn{p} numeric (\eqn{N \times n}) matrix containing the
#' in-sample residuals. This input is used to compute some covariance matrices.
#' @param approach A string specifying the approach used to compute the reconciled
#' forecasts. Options include:
#'   \itemize{
#'   \item "\code{proj}" (\emph{default}): zero-constrained approach.
#'   \item "\code{strc}": Structural approach.
#'   }
#' @param nn A string specifying the algorithm to compute non-negative reconciled forecasts:
#'   \itemize{
#'   \item "\code{osqp}": OSQP solver (Stellato et al., 2020).
#'   \item "\code{sntz}": heuristic "set-negative-to-zero" (Di Fonzo and Girolimetto, 2023).
#'   }
#' @param settings An object of class \code{osqpSettings} specifying settings
#' for the \href{https://osqp.org/}{\pkg{osqp}} solver. For details, refer to the
#' \href{https://osqp.org/}{\pkg{osqp} documentation} (Stellato et al., 2020).
#' @param bounds A (\eqn{n \times 2}) numeric matrix specifying the cross-sectional bounds.
#' The first column represents the lower bound, and the second column represents the upper bound.
#' @param ... Arguments passed on to [cscov].
#'
#' @details
#' If an expert does not provide a forecast for a given variable,
#' the missing value can be represented as \code{NA}.
#' For instance, consider a case where the number of experts is
#' \eqn{p = 4}, with \eqn{n = 3} variable and \eqn{h = 2}. If the
#' second expert does not provide predictions for the second variable,
#' the corresponding matrix in the base list will be as follows:
#' \deqn{\texttt{base[[2]]} = \begin{bmatrix} \hat{y}_{1,1} & \texttt{NA} & \hat{y}_{3,2} \\
#' \hat{y}_{1,2} & \texttt{NA} & \hat{y}_{3,2} \end{bmatrix}.}
#'
#'
#' @references
#' Di Fonzo, T. and Girolimetto, D. (2023), Spatio-temporal reconciliation of solar
#' forecasts, \emph{Solar Energy}, 251, 13–29. \doi{10.1016/j.solener.2023.01.003}
#'
#' Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination for linearly
#' constrained multiple time series, \emph{mimeo}.
#'
#' Stellato, B., Banjac, G., Goulart, P., Bemporad, A. and Boyd, S. (2020), OSQP:
#' An Operator Splitting solver for Quadratic Programs,
#' \emph{Mathematical Programming Computation}, 12, 4, 637-672.
#' \doi{10.1007/s12532-020-00179-2}
#'
#' @returns A (\eqn{h \times n}) numeric matrix of cross-sectional combined and reconciled forecasts.
#'
#' @family Optimal coherent combination
#'
#' @examples
#' set.seed(123)
#' # (2 x 3) base forecasts matrix (simulated), expert 1
#' base1 <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
#' # (10 x 3) in-sample residuals matrix (simulated), expert 1
#' res1 <- t(matrix(rnorm(n = 30), nrow = 3))
#'
#' # (2 x 3) base forecasts matrix (simulated), expert 2
#' base2 <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
#' # (10 x 3) in-sample residuals matrix (simulated), expert 2
#' res2 <- t(matrix(rnorm(n = 30), nrow = 3))
#'
#' ## RECTANGULAR CASE
#' # Base forecasts' and residuals' lists
#' brc <- list(base1, base2)
#' erc <- list(res1, res2)
#'
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' rrc <- csocc(base = brc, agg_mat = A, comb = "wls", res = erc)
#'
#' # Zero constraints matrix for Z - X - Y = 0
#' C <- t(c(1, -1, -1))
#' rrc <- csocc(base = brc, cons_mat = C, comb = "wls", res = erc) # same results
#'
#' ## NON-RECTANGULAR CASE
#' base2[, 2] <- res2[, 2] <-  NA
#'
#' # Base forecasts' and residuals' lists
#' bgc <- list(base1, base2)
#' egc <- list(res1, res2)
#'
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' rgc <- csocc(base = bgc, agg_mat = A, comb = "shrbe", res = egc)
#'
#' @export
csocc <- function(base, agg_mat, cons_mat,
                   comb = "ols", res = NULL, approach = "proj",
                   nn = NULL, settings = NULL, bounds = NULL, ...){

  # Check if either 'agg_mat' or 'cons_mat' is specified
  if(missing(agg_mat) && missing(cons_mat)){
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.", call = NULL)
  } else if(!missing(agg_mat)){
    tmp <- cstools(agg_mat = agg_mat)
  } else {
    tmp <- cstools(cons_mat = cons_mat)
  }

  n <- tmp$dim[["n"]]
  strc_mat <- tmp$strc_mat
  cons_mat <- tmp$cons_mat
  p <- length(base)

  # Check if 'base' is provided and its dimensions match with the data
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  ina <- sapply(base, function(bmat){
    is.na(colSums(bmat))
  })
  base <- lapply(base, rbind)
  base <- do.call(cbind, base)

  if(NCOL(base) != n*p){
    cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
  }
  base <- base[, !as.vector(ina), drop = FALSE]

  # Compute covariance
  if(!is.null(res)){
    res <- do.call(cbind, res)
    res <- res[, !as.vector(ina), drop = FALSE]
  }

  cov_mat <- cscov(comb = comb, n = NCOL(base), matNA = ina, p = p, nv = n,
                   agg_mat = rbind(do.call(rbind, rep(list(strc_mat), p-1)), agg_mat),
                   res = res, ...)

  if(NROW(cov_mat) != NCOL(base) | NCOL(cov_mat) != NCOL(base)){
    if(any(as.vector(ina))){
      cov_mat <- cov_mat[!as.vector(ina), !as.vector(ina), drop = FALSE]
    }else{
      cli_abort(c("Incorrect covariance dimensions.",
                  "i"="Check {.arg res} dimensions."), call = NULL)
    }
  }

  reco_mat <- resemble(base = base,
                       cov_mat = cov_mat,
                       strc_mat = strc_mat,
                       cons_mat = cons_mat,
                       approach = approach,
                       nn = nn,
                       ina = as.vector(ina),
                       p = p,
                       bounds = bounds,
                       settings = settings)

  rownames(reco_mat) <- paste0("h-", 1:NROW(reco_mat))
  if(is.null(colnames(base))){
    colnames(reco_mat) <- paste0("s-", 1:NCOL(reco_mat))
  } else {
    colnames(reco_mat) <- colnames(base[,1:NCOL(reco_mat), drop = FALSE])
  }

  attr(reco_mat, "FoReco") <- list2env(list(info = attr(reco_mat, "info"),
                                            framework = "Cross-sectional",
                                            forecast_horizon = NROW(reco_mat),
                                            comb = comb,
                                            cs_n = n,
                                            rfun = "csocc"))
  attr(reco_mat, "info") <- NULL
  return(reco_mat)
}
