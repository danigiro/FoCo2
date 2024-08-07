#' Projection matrix
#'
#' @inheritParams csocc
#'
#' @return A list
#' @export
csprojmat <- function(agg_mat, cons_mat, block_diag = "fr", p = NULL,
                      comb = "ols", res = NULL, approach = "proj",
                      ...){

  block_diag <- match.arg(block_diag, c("fr", "fc", "none"))

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

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(NCOL(cons_mat)))

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

  return(list(M = M, Omega = Omega, W = cov_mat, Wc = solve(k_cov_k)))
}
