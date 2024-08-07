#' Base forecast component
#'
#' @inheritParams csocc
#'
#' @returns A (\eqn{h \times n}) numeric matrix of cross-sectional combined forecasts.
#' @export
csbase <- function(base, agg_mat = NULL, block_diag = "fr", p = NULL,
                   comb = "ols", res = NULL, ...){

  block_diag <- match.arg(block_diag, c("fr", "fc", "none"))

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

  base <- do.call(cbind, base)
  if(NCOL(base) != n*p){
    cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
  }

  # Compute covariance
  if(block_diag %in% c("fr", "fc")){
    if(block_diag == "fc" && !is.null(res)){
      res <- simplify2array(res)
      res <- lapply(1:(dim(res)[2]), function(i) res[,i,])
    }
    cov_mat <- lapply(1:length(res), function(id){
      cscov(comb = comb, n = n, agg_mat = agg_mat,
            res = res[[id]], ...)
    })
    cov_mat <- bdiag(cov_mat)
    if(block_diag == "fc" && !is.null(res)){
      P <- commat(n, p)
      cov_mat <- t(P)%*%cov_mat%*%P
    }
  }else{
    if(!is.null(res)){
      res <- do.call(cbind, res)
    }
    cov_mat <- cscov(comb = comb, n = n*p,
                     agg_mat = rbind(do.call(rbind, rep(list(strc_mat), p-1)), agg_mat),
                     res = res, ...)
  }

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(n))

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
