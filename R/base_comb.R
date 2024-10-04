#' Base forecast component of the optimal coherent forecast combination
#'
#' This function computes the base forecast component of the optimal coherent
#' forecast combination [csocc], as described in Girolimetto and Di Fonzo (2024)
#'
#' @usage occbase(base, agg_mat = NULL, comb = "ols", res = NULL, ...)
#'
#' @inheritParams csocc
#'
#' @returns A (\eqn{h \times n}) numeric matrix of cross-sectional combined forecasts.
#'
#' @references
#' Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination for linearly
#' constrained multiple time series, \emph{mimeo}.
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
#' rrc <- csocc(base = brc, agg_mat = A, comb = "shr", res = erc)
#'
#' yc <- occbase(base = brc, agg_mat = A, comb = "shr", res = erc)
#' M <- occmat(base = brc, agg_mat = A, comb = "shr", p = 2, res = erc)$M
#' M%*%t(yc)-t(rrc)
#'
#' ## NON-RECTANGULAR CASE
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
#' yc <- occbase(base = bgc, agg_mat = A, comb = "shr", res = egc)
#' M <- occmat(base = bgc, agg_mat = A, comb = "shr", p = 2, res = egc, matNA = matNA)$M
#' M%*%t(yc)-t(rgc)
#'
#' @export
occbase <- function(base, agg_mat = NULL,
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
  } else {
    colnames(base_comp) <- colnames(base[,1:NCOL(base_comp), drop = FALSE])
  }

  return(base_comp)
}
