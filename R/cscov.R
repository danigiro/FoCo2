#' Cross-sectional covariance matrix approximation
#'
#' Extended version of the
#' \href{https://danigiro.github.io/FoReco/reference/cscov.html}{\code{FoReco::cscov}}
#' function, introducing two new approximations for the covariance matrix
#' (both shrunk and sample versions). Specifically, \code{shrbe}/\code{sambe}
#' assume no correlation between experts, while \code{shrbv}/\code{sambv}
#' assume no correlation between variables.
#'
#' @usage
#' \method{cscov}{shrbe}(comb = "shrbe", ..., n = NULL, p = NULL, matNA = NULL,
#'       res = NULL, mse = TRUE, shrink_fun = NULL)
#'
#' @inheritParams csocc
#' @param p Total number of experts, \eqn{p}.
#' @param matNA A (\eqn{n \times p}) matrix consisting of 0s and 1s,
#' where each element indicates whether expert \eqn{j} (column) has
#' provided a forecast for variable \eqn{i} (row). If expert \eqn{j}
#' has provided a forecast for variable \eqn{i}, the corresponding
#' element \eqn{(i,j)} is 1; otherwise, it is 0.
#' @param n Total number of variables, \eqn{n}.
#' @param comb A string specifying the reconciliation method.
#'   \itemize{
#'      \item \href{https://danigiro.github.io/FoReco/}{\code{FoReco}} approaches:
#'            \code{"ols"}, \code{"wls"}, \code{"shr"}, \code{"sam"}.
#'      \item "\code{shrbe}"/"\code{sambe}" - shrunk/sample block-diagonal covariance by experts.
#'      \item "\code{shrbv}"/"\code{sambv}" - shrunk/sample block-diagonal covariance by variables.
#'    }
#' @param mse If \code{TRUE} (\emph{default}) the residuals used to compute the covariance
#' matrix are not mean-corrected.
#' @param shrink_fun Shrinkage function of the covariance matrix,
#' \href{https://danigiro.github.io/FoReco/reference/shrink_estim.html}{\code{FoReco::shrink_estim}}
#' (\emph{default}).
#' @param ... Arguments passed on to
#' \href{https://danigiro.github.io/FoReco/reference/cscov.html}{\code{FoReco::cscov}}.
#'
#' @returns A (\eqn{m \times m}) symmetric positive (semi-)definite matrix, with
#' \eqn{m = \sum_{j = 1}^p n_j}, \eqn{n_j \leq n}.
#'
#' @rdname cscov
#' @name cscov
#'
#' @family Optimal combination
#'
#' @exportS3Method FoReco::cscov
cscov.shrbe <- function(comb = "shrbe", ..., n = NULL, p = NULL, matNA = NULL,
                        res = NULL, mse = TRUE, shrink_fun = NULL){

  if(is.null(n)){
    if(is.list(res)){
      n <- NCOL(res[[1]])
    }else{
      cli_abort("Argument {.arg n} is NULL.", call = NULL)
    }
  }

  if(is.null(p)){
    if(is.list(res)){
      p <- length(res)
    }else{
      cli_abort("Argument {.arg p} is NULL.", call = NULL)
    }
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(is.list(res)){
    res <- do.call(cbind, res)
  }
  id <- rep(1:p, each = n)
  if(!is.null(matNA)){
    if(is.logical(matNA)){
      ina <- as.vector(!matNA)
    }else{
      ina <- as.vector(is.na(matNA) | matNA!=0)
    }
    if(NCOL(res) != sum(ina)){
      res <- res[,ina, drop = FALSE]
    }
    id <- id[ina]
  }

  if(is.null(shrink_fun)){
    if(!is.null(matNA)){
      if(sum(ina) == n*p){
        shrink_fun <- shrink_estim
      }else{
        shrink_fun <- shrink_estim_na
      }
    }else{
      shrink_fun <- shrink_estim_na
    }
  }
  Wlist <- lapply(sort(unique(id)), function(i){
    shrink_fun(res[, id == i, drop = FALSE], mse = mse)
  })
  bdiag(Wlist)
}

#' @usage
#' \method{cscov}{sambe}(comb = "sambe", ..., n = NULL, p = NULL, matNA = NULL,
#'       res = NULL, mse = TRUE)
#' @rdname cscov
#' @exportS3Method FoReco::cscov
cscov.sambe <- function(comb = "sambe", ..., n = NULL, p = NULL, matNA = NULL,
                        res = NULL, mse = TRUE){
  if(is.null(n)){
    if(is.list(res)){
      n <- NCOL(res[[1]])
    }else{
      cli_abort("Argument {.arg n} is NULL.", call = NULL)
    }
  }

  if(is.null(p)){
    if(is.list(res)){
      p <- length(res)
    }else{
      cli_abort("Argument {.arg p} is NULL.", call = NULL)
    }
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(is.list(res)){
    res <- do.call(cbind, res)
  }
  id <- rep(1:p, each = n)
  if(!is.null(matNA)){
    if(is.logical(matNA)){
      ina <- as.vector(!matNA)
    }else{
      ina <- as.vector(is.na(matNA) | matNA!=0)
    }
    if(NCOL(res) != sum(ina)){
      res <- res[,ina, drop = FALSE]
    }
    id <- id[ina]
  }

  Wlist <- lapply(sort(unique(id)), function(i){
    sample_estim(res[, id == i, drop = FALSE], mse = mse)
  })
  bdiag(Wlist)
}

#' @usage
#' \method{cscov}{shrbv}(comb = "shrbv", ..., n = NULL, p = NULL, matNA = NULL,
#'       res = NULL, mse = TRUE, shrink_fun = NULL)
#' @rdname cscov
#' @exportS3Method FoReco::cscov
cscov.shrbv <- function(comb = "shrbv", ..., n = NULL, p = NULL, matNA = NULL,
                        res = NULL, mse = TRUE, shrink_fun = NULL){
  if(is.null(n)){
    if(is.list(res)){
      n <- NCOL(res[[1]])
    }else{
      cli_abort("Argument {.arg n} is NULL.", call = NULL)
    }
  }

  if(is.null(p)){
    if(is.list(res)){
      p <- length(res)
    }else{
      cli_abort("Argument {.arg p} is NULL.", call = NULL)
    }
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(is.list(res)){
    res <- do.call(cbind, res)
  }
  id <- rep(1:n, p)
  if(!is.null(matNA)){
    if(is.logical(matNA)){
      ina <- as.vector(!matNA)
    }else{
      ina <- as.vector(is.na(matNA) | matNA!=0)
    }
    if(NCOL(res) != sum(ina)){
      res <- res[,ina, drop = FALSE]
    }
    id <- id[ina]
  }

  if(is.null(shrink_fun)){
    if(!is.null(matNA)){
      if(sum(ina) == n*p){
        shrink_fun <- shrink_estim
      }else{
        shrink_fun <- shrink_estim_na
      }
    }else{
      shrink_fun <- shrink_estim_na
    }
  }

  Slist <- lapply(sort(unique(id)), function(i){
    shrink_fun(res[, id == i, drop = FALSE], mse = mse)
  })
  cov_mat <- bdiag(Slist)

  P <- commat(n, p)
  if(!is.null(matNA)){
    if(any(!ina)){
      P <- P[, ina, drop = FALSE]
      P <- P[rowSums(P)!=0, , drop = FALSE]
    }
  }
  t(P)%*%cov_mat%*%P

}

#' @usage
#' \method{cscov}{sambv}(comb = "sambv", ..., n = NULL, p = NULL, matNA = NULL,
#'       res = NULL, mse = TRUE)
#' @rdname cscov
#' @exportS3Method FoReco::cscov
cscov.sambv <- function(comb = "sambv", ..., n = NULL, p = NULL, matNA = NULL,
                        res = NULL, mse = TRUE){
  if(is.null(n)){
    if(is.list(res)){
      n <- NCOL(res[[1]])
    }else{
      cli_abort("Argument {.arg n} is NULL.", call = NULL)
    }
  }

  if(is.null(p)){
    if(is.list(res)){
      p <- length(res)
    }else{
      cli_abort("Argument {.arg p} is NULL.", call = NULL)
    }
  }

  if(is.null(res)){
    cli_abort("Argument {.arg res} is NULL.", call = NULL)
  }

  if(is.list(res)){
    res <- do.call(cbind, res)
  }
  id <- rep(1:n, p)
  if(!is.null(matNA)){
    if(is.logical(matNA)){
      ina <- as.vector(!matNA)
    }else{
      ina <- as.vector(is.na(matNA) | matNA!=0)
    }
    if(NCOL(res) != sum(ina)){
      res <- res[,ina, drop = FALSE]
    }
    id <- id[ina]
  }

  Slist <- lapply(sort(unique(id)), function(i){
    sample_estim(res[, id == i, drop = FALSE], mse = mse)
  })
  cov_mat <- bdiag(Slist)

  P <- commat(n, p)
  if(!is.null(matNA)){
    if(any(!ina)){
      P <- P[, ina, drop = FALSE]
      P <- P[rowSums(P)!=0, , drop = FALSE]
    }
  }
  t(P)%*%cov_mat%*%P
}
