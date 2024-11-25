#' Matrices for the optimal coherent forecast combination
#'
#' This function computes the matrices required for the optimal coherent
#' forecast combination [csocc], as described in Girolimetto and Di Fonzo (2024).
#' These matrices serve as the foundation for building forecasts that effectively
#' combines the individual information from multiple experts while ensuring
#' coherence across the variables.
#'
#' @usage occmat(agg_mat, cons_mat, p = NULL, matNA = NULL,
#'        comb = "ols", res = NULL, approach = "proj", ...)
#'
#' @inheritParams csocc
#' @inheritParams cscov
#'
#' @return A list of matrices:
#' \item{M}{Projection matrix.}
#' \item{Omega}{Matrix of the combination weights of the optimal linear multi-task forecast combination.}
#' \item{W}{Forecast error covariance matrix of the base forecasts.}
#' \item{Wc}{Forecast error covariance matrix of the combined forecasts.}
#' \item{Wtilde}{Forecast error covariance matrix of the reconciled combined forecasts.}
#' \item{K}{Matrix that replicates a vector (see Girolimetto and Di Fonzo, 2024).}
#'
#' @references
#' Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination for linearly
#' constrained multiple time series, \emph{mimeo}.
#'
#' @family Optimal combination
#'
#' @export
occmat <- function(agg_mat, cons_mat, p = NULL, matNA = NULL,
                   comb = "ols", res = NULL, approach = "proj", ...){

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
  if(is.null(p)){
    if(is.null(res)){
      cli_abort("Argument {.arg p} is missing,
              with no default.", call = NULL)
    }else{
      p <- length(res)
    }
  }else{
    if(!is.null(res) && p != length(res)){
      cli_abort("Argument {.arg p} is not equal to the length of {.arg res}.",
                call = NULL)
    }
  }

  if(!is.null(matNA)){
    ina <- is.na(matNA) | matNA==0
  }else{
    ina <- matrix(FALSE, p, n)
  }

  # Compute covariance
  if(!is.null(res)){
    res <- do.call(cbind, res)
    res <- res[, !as.vector(ina), drop = FALSE]
  }

  cov_mat <- cscov(comb = comb,  n = ifelse(comb == "ols", n*p, n), matNA = ina, p = p,
                   agg_mat = rbind(do.call(rbind, rep(list(strc_mat), p-1)), agg_mat),
                   res = res, ...)

  if(NROW(cov_mat) != sum(!ina) | NCOL(cov_mat) != sum(!ina)){
    if(any(as.vector(ina))){
      if(NROW(cov_mat) != length(ina)){
        cli_abort(c("Incorrect covariance dimensions.",
                    "i"="Check {.arg res} dimensions."), call = NULL)
      }
      cov_mat <- cov_mat[!as.vector(ina), , drop = FALSE][, !as.vector(ina), drop = FALSE]
    }else{
      cli_abort(c("Incorrect covariance dimensions.",
                  "i"="Check {.arg res} dimensions."), call = NULL)
    }
  }

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(NCOL(cons_mat)))[!as.vector(ina), , drop = FALSE]

  if(isDiagonal(cov_mat)){
    cov_k <- lin_sys(cov_mat, k_mat)
    k_cov_k <- methods::as(Matrix::crossprod(k_mat, cov_k), "CsparseMatrix")
    cov_c <- lin_sys(k_cov_k, t(cons_mat))
    c_cov_c <- methods::as(cons_mat%*%cov_c, "CsparseMatrix")
    ls1 <- lin_sys(c_cov_c, cons_mat)
    Omega <- t(lin_sys(k_cov_k, t(cov_k)))
    M <- Diagonal(n) - cov_c%*%ls1
  }else{
    cov_k <- lin_sys(cov_mat, k_mat)
    k_cov_k <- methods::as(Matrix::crossprod(k_mat, cov_k), "CsparseMatrix")
    cov_c <- lin_sys(k_cov_k, t(cons_mat))
    c_cov_c <- methods::as(cons_mat%*%cov_c, "CsparseMatrix")
    ls1 <- lin_sys(c_cov_c, cons_mat)
    Omega <- t(lin_sys(k_cov_k, t(cov_k)))
    M <- Diagonal(n) - cov_c%*%ls1
  }

  Wc <- solve(k_cov_k)

  return(list(M = M, Omega = Omega, W = cov_mat,
              Wc = Wc, K = k_mat, Wtilde = M%*%Wc))
}
