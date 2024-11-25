#' Cross-sectional optimal coherent forecast combination
#'
#' This function implements the optimal linear coherent forecast combination for a
#' linearly constrained (e.g., hierarchical/grouped) multiple time series
#' combining forecasts from multiple sources (experts). The method minimizes forecast
#' errors while ensuring that the resulting forecasts are coherent, meaning
#' they satisfy the constraints across the time series. Using a linear constrained
#' optimization approach, [csocc] combines forecasts from multiple experts,
#' accounting for the relations between variables and constraints. In addition, linear
#' inequality constraints (e.g. non-negative forecasts) can be imposed, if needed.
#'
#' @usage csocc(base, agg_mat, cons_mat, comb = "ols", res = NULL,
#'       approach = "proj", nn = NULL, settings = NULL, bounds = NULL, ...)
#'
#' @param base A list of \eqn{p} numeric (\eqn{h \times n}) matrix or multivariate time series
#' (\code{mts} class) containing the base forecasts to be reconciled; \eqn{h} is the forecast
#' horizon, \eqn{n} is the total number of time series (\eqn{n = n_u + n_b}) and
#' \eqn{p} is the total number of experts.
#' @param agg_mat A (\eqn{n_u \times n_b}) numeric matrix representing the cross-sectional
#' aggregation matrix, mapping the \eqn{n_b} bottom-level (free)
#' variables into the \eqn{n_u} upper (constrained) variables.
#' @param cons_mat A (\eqn{n_u \times n}) numeric matrix representing the cross-sectional
#' zero constraints: each row represents a constraint equation, and each column represents
#' a variable. The matrix can be of full rank, meaning the rows are linearly independent,
#' but this is not a strict requirement, as the function allows for redundancy in the constraints.
#' @param comb A string specifying the reconciliation method. For details, see [cscov].
#' @param res A list of \eqn{p} numeric (\eqn{N \times n}) matrix containing the
#' in-sample residuals. This input is used to compute some covariance matrices.
#' @param approach A string specifying the approach used to compute the reconciled
#' forecasts. Options include:
#'   \itemize{
#'   \item "\code{proj}" (\emph{default}): zero-constrained projection approach.
#'   \item "\code{strc}": structural approach.
#'   }
#' @param nn A string specifying the algorithm to compute non-negative forecasts:
#'   \itemize{
#'   \item "\code{osqp}": OSQP solver (Stellato et al., 2020).
#'   \item "\code{sntz}": heuristic "set-negative-to-zero" (Di Fonzo and Girolimetto, 2023).
#'   }
#' @param settings An object of class \code{osqpSettings} specifying settings
#' for the \href{https://osqp.org/}{\pkg{osqp}} solver. For details, refer to the
#' \href{https://osqp.org/}{\pkg{osqp} documentation} (Stellato et al., 2020).
#' @param bounds A matrix (see \code{set_bounds} in \href{https://danigiro.github.io/FoReco/reference/set_bounds.html}{\pkg{FoReco}})
#' with 3 columns (\eqn{i,lower,upper}), such that
#' \itemize{
#'   \item Column 1 represents the cross-sectional series (\eqn{i = 1, \dots, n}).
#'   \item Column 2 indicates the \emph{lower} bound.
#'   \item Column 3 indicates the \emph{upper} bound.
#' }
#' @param ... Arguments passed on to [cscov].
#'
#' @details
#' If an expert does not provide a forecast for a given variable,
#' the missing value can be represented as \code{NA}.
#' For instance, consider a case where the number of experts is
#' \eqn{p = 4}, with \eqn{n = 3} variables and \eqn{h = 2}. If the
#' second and third experts do not provide forecasts for, respectively, the second and first variable,
#' the corresponding matrix in the base list will be as follows:
#' \deqn{\texttt{base[[1]]} = \begin{bmatrix} \hat{y}_{1,1}^{1} & \hat{y}_{2,1}^{1} & \hat{y}_{3,1}^{1} \\
#' \hat{y}_{1,2}^{1} & \hat{y}_{2,2}^{1} & \hat{y}_{3,2}^{1} \end{bmatrix}, \quad
#' \texttt{base[[2]]} = \begin{bmatrix} \hat{y}_{1,1}^{2} & \texttt{NA} & \hat{y}_{3,1}^{2} \\
#' \hat{y}_{1,2}^{2} & \texttt{NA} & \hat{y}_{3,2}^{2} \end{bmatrix},}
#' \deqn{\texttt{base[[3]]} = \begin{bmatrix} \texttt{NA} & \hat{y}_{2,1}^{3} & \hat{y}_{3,1}^{3} \\
#' \texttt{NA} & \hat{y}_{2,2}^{3} & \hat{y}_{3,2}^{3} \end{bmatrix}, \quad
#' \texttt{base[[4]]} = \begin{bmatrix} \hat{y}_{1,1}^{4} & \hat{y}_{2,1}^{4} & \hat{y}_{3,1}^{4} \\
#' \hat{y}_{1,2}^{4} & \hat{y}_{2,2}^{4} & \hat{y}_{3,2}^{4} \end{bmatrix}.}
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
#' @family Optimal combination
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
#' ## BALANCED PANEL OF FORECASTS
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
#' ## UNBALANCED PANEL OF FORECASTS
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

  # Check bounds cs
  if(!is.null(bounds)){
    if(is.vector(bounds)){
      bounds <- matrix(bounds, ncol = length(bounds))
    }
    if(NCOL(bounds) != 3){
      cli_abort("{.arg bounds} is not a matrix with 3 columns.", call = NULL)
    }
    bounds_approach <- attr(bounds, "approach")

    bounds <- bounds[bounds[,1] <= n, , drop = FALSE]

    if(is.null(bounds)){
      cli_warn("No valid bounds", call = NULL)
    }else{
      attr(bounds, "approach") <- bounds_approach
    }
  }

  # Compute covariance
  if(!is.null(res)){
    res <- do.call(cbind, res)
    res <- res[, !as.vector(ina), drop = FALSE]
  }

  cov_mat <- cscov(comb = comb, n = ifelse(comb == "ols", NCOL(base), n),
                   matNA = ina, p = p,
                   agg_mat = rbind(do.call(rbind, rep(list(strc_mat), p-1)), agg_mat),
                   res = res, ...)

  if(NROW(cov_mat) != NCOL(base) | NCOL(cov_mat) != NCOL(base)){
    if(any(as.vector(ina))){
      if(NROW(cov_mat) != length(ina)){
        cli_abort(c("Incorrect covariance dimensions.",
                    "i"="Check {.arg res} dimensions."), call = NULL)
      }
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
  }else{
    colnames(reco_mat) <- rownames(ina)
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

#' Cross-sectional optimal multi-task forecast combination
#'
#' This function computes the optimal multi-task linear forecast combination, as described
#' in Girolimetto and Di Fonzo (2024)
#'
#' @usage csmtc(base, agg_mat = NULL, comb = "ols", res = NULL, ...)
#'
#' @inheritParams csocc
#'
#' @returns A (\eqn{h \times n}) numeric matrix of cross-sectional multi-task combined forecasts.
#'
#' @references
#' Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination for linearly
#' constrained multiple time series, \emph{mimeo}.
#'
#' @family Optimal combination
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
#' ## BALANCED PANEL OF FORECASTS
#' # Base forecasts' and residuals' lists
#' brc <- list(base1, base2)
#' erc <- list(res1, res2)
#'
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' rrc <- csocc(base = brc, agg_mat = A, comb = "shr", res = erc)
#'
#' yc <- csmtc(base = brc, agg_mat = A, comb = "shr", res = erc)
#' M <- occmat(base = brc, agg_mat = A, comb = "shr", p = 2, res = erc)$M
#' M%*%t(yc)-t(rrc)
#'
#' ## UNBALANCED PANEL OF FORECASTS
#' base2[, 2] <- res2[, 2] <-  NA
#'
#' # Base forecasts' and residuals' lists
#' bgc <- list(base1, base2)
#' egc <- list(res1, res2)
#' matNA <- matrix(1, 3, 2)
#' matNA[2,2] <- 0
#'
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' rgc <- csocc(base = bgc, agg_mat = A, comb = "shr", res = egc)
#'
#' yc <- csmtc(base = bgc, agg_mat = A, comb = "shr", res = egc)
#' M <- occmat(base = bgc, agg_mat = A, comb = "shr", p = 2, res = egc, matNA = matNA)$M
#' M%*%t(yc)-t(rgc)
#'
#' @export
csmtc <- function(base, agg_mat = NULL,
                  comb = "ols", res = NULL, ...){

  # Check if 'base' is provided and its dimensions match with the data
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  base <- lapply(base, rbind)
  n <- NCOL(base[[1]])
  p <- length(base)

  if(!is.null(agg_mat)){
    tmp <- cstools(agg_mat = agg_mat)
    agg_mat <- tmp$agg_mat
    strc_mat <- tmp$strc_mat
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

  cov_mat <- cscov(comb = comb,  n = ifelse(comb == "ols", NCOL(base), n), matNA = ina, p = p,
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

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(n))[!as.vector(ina), , drop = FALSE]

  if(isDiagonal(cov_mat)){
    cov_inv <- Matrix::.sparseDiagonal(x = Matrix::diag(cov_mat)^(-1))
    cov_k <- Matrix::crossprod(k_mat, cov_inv)%*%k_mat
    cov_k_inv <- Matrix::.sparseDiagonal(x = Matrix::diag(cov_k)^(-1))
    base_comp <- Matrix::tcrossprod(base, cov_inv)%*%Matrix::tcrossprod(k_mat, cov_k_inv)
  }else{
    cov_k <- lin_sys(cov_mat, k_mat)
    k_cov_k <- methods::as(Matrix::crossprod(k_mat, cov_k), "CsparseMatrix")
    ls2 <- lin_sys(k_cov_k, t(cov_k))
    base_comp <- methods::as(Matrix::tcrossprod(base, ls2), "CsparseMatrix")
  }

  rownames(base_comp) <- paste0("h-", 1:NROW(base_comp))
  if(is.null(colnames(base))){
    colnames(base_comp) <- paste0("s-", 1:NCOL(base_comp))
  }else{
    colnames(base_comp) <- rownames(ina)
  }

  return(base_comp)
}
