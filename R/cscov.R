#' Cross-sectional covariance matrix approximation
#'
#' Extended version of the
#' \href{https://danigiro.github.io/FoReco/reference/cscov.html}{\code{FoReco::cscov}},
#' introducing two new block approximations for the covariance matrix
#' (shrunk and sample version). Specifically, \code{shrbe}/\code{sambe}
#' assume no correlation between experts, while \code{shrbv}/\code{sambv}
#' assume no correlation between variables.
#'
#' @usage
#' \method{cscov}{shrbe}(comb = "shrbe", ..., nv = NULL, p = NULL, matNA = NULL,
#'       res = NULL, mse = TRUE, shrink_fun = shrink_estim)
#'
#' @inheritParams csocc
#' @param p Total number of experts.
#' @param matNA A (\eqn{n \times p}) matrix consisting of 0s and 1s,
#' where each element indicates whether expert \eqn{j} (column) has
#' provided a forecast for variable \eqn{i} (row). If expert \eqn{j}
#' has provided a forecast for variable \eqn{i}, the corresponding
#' element \eqn{(i,j)} is 1; otherwise, it is 0.
#' @param nv Total number of variables.
#' @param comb A string specifying the reconciliation method.
#'   \itemize{
#'      \item \href{https://danigiro.github.io/FoReco/}{\code{FoReco}} approximations:
#'            \code{"ols"}, \code{"wls"}, \code{"shr"}, \code{"sam"}.
#'      \item "\code{shrbe}"/"\code{sambe}" - shrunk/sample block diagonal covariance by experts.
#'      \item "\code{shrbv}"/"\code{sambv}" - shrunk/sample block diagonal covariance by variables.
#'    }
#' @param mse If \code{TRUE} (\emph{default}) the residuals used to compute the covariance
#' matrix are not mean-corrected.
#' @param shrink_fun Shrinkage function of the covariance matrix,
#' \href{https://danigiro.github.io/FoReco/reference/shrink_estim.html}{\code{FoReco::shrink_estim}}
#' (\emph{default}).
#' @param ... Arguments passed on to
#' \href{https://danigiro.github.io/FoReco/reference/cscov.html}{\code{FoReco::cscov}}.
#'
#' @returns A (\eqn{m \times m}) symmetric positive (semi-)definite matrix.
#'
#' @rdname cscov
#' @name cscov
#'
#' @family Optimal coherent combination
#'
#' @exportS3Method FoReco::cscov
cscov.shrbe <- function(comb = "shrbe", ..., nv = NULL, p = NULL, matNA = NULL,
                        res = NULL, mse = TRUE, shrink_fun = shrink_estim){

  if(is.null(nv)){
    if(is.list(res)){
      nv <- NCOL(res[[1]])
    }else{
      cli_abort("Argument {.arg nv} is NULL.", call = NULL)
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

  id <- rep(1:p, each = nv)
  if(!is.null(matNA)){
    if(is.logical(matNA)){
      ina <- as.vector(t(matNA))
      print(dim(matNA))
    }else{
      ina <- t(is.na(matNA) | matNA==0)
    }
    if(NCOL(res) != sum(!ina)){
      res <- res[,!ina, drop = FALSE]
    }
    id <- id[!ina]
  }
  Wlist <- lapply(1:p, function(i){
    shrink_fun(res[, id == i, drop = FALSE], mse = mse)
  })
  bdiag(Wlist)
}

#' @usage
#' \method{cscov}{sambe}(comb = "sambe", ..., nv = NULL, p = NULL, matNA = NULL,
#'       res = NULL, mse = TRUE)
#' @rdname cscov
#' @exportS3Method FoReco::cscov
cscov.sambe <- function(comb = "sambe", ..., nv = NULL, p = NULL, matNA = NULL,
                        res = NULL, mse = TRUE){
  if(is.null(nv)){
    if(is.list(res)){
      nv <- NCOL(res[[1]])
    }else{
      cli_abort("Argument {.arg nv} is NULL.", call = NULL)
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
  id <- rep(1:p, each = nv)
  if(!is.null(matNA)){
    if(is.logical(matNA)){
      ina <- as.vector(t(matNA))
    }else{
      ina <- t(is.na(matNA) | matNA==0)
    }

    if(NCOL(res) != sum(!ina)){
      res <- res[,!ina, drop = FALSE]
    }
    id <- id[!ina]
  }
  Wlist <- lapply(1:p, function(i){
    sample_estim(res[, id == i, drop = FALSE], mse = mse)
  })
  bdiag(Wlist)
}

#' @usage
#' \method{cscov}{shrbv}(comb = "shrbv", ..., nv = NULL, p = NULL, matNA = NULL,
#'       res = NULL, mse = TRUE, shrink_fun = shrink_estim)
#' @rdname cscov
#' @exportS3Method FoReco::cscov
cscov.shrbv <- function(comb = "shrbv", ..., nv = NULL, p = NULL, matNA = NULL,
                        res = NULL, mse = TRUE, shrink_fun = shrink_estim){
  if(is.null(nv)){
    if(is.list(res)){
      nv <- NCOL(res[[1]])
    }else{
      cli_abort("Argument {.arg nv} is NULL.", call = NULL)
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
  id <- rep(1:nv, p)
  if(!is.null(matNA)){
    if(is.logical(matNA)){
      ina <- as.vector(t(matNA))
    }else{
      ina <- t(is.na(matNA) | matNA==0)
    }

    if(NCOL(res) != sum(!ina)){
      res <- res[,!ina, drop = FALSE]
    }
    id <- id[!ina]
  }
  Slist <- lapply(1:nv, function(i){
    shrink_fun(res[, id == i, drop = FALSE], mse = mse)
  })
  cov_mat <- bdiag(Slist)

  P <- commat(nv, p)
  if(any(ina)){
    P <- P[!ina, , drop = FALSE]
    P <- P[, colSums(P)!=0, drop = FALSE]
  }
  t(P)%*%cov_mat%*%P

}

#' @usage
#' \method{cscov}{sambv}(comb = "sambv", ..., nv = NULL, p = NULL, matNA = NULL,
#'       res = NULL, mse = TRUE)
#' @rdname cscov
#' @exportS3Method FoReco::cscov
cscov.sambv <- function(comb = "sambv", ..., nv = NULL, p = NULL, matNA = NULL,
                        res = NULL, mse = TRUE){
  if(is.null(nv)){
    if(is.list(res)){
      nv <- NCOL(res[[1]])
    }else{
      cli_abort("Argument {.arg nv} is NULL.", call = NULL)
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
  id <- rep(1:nv, p)
  if(!is.null(matNA)){
    if(is.logical(matNA)){
      ina <- as.vector(t(matNA))
    }else{
      ina <- t(is.na(matNA) | matNA==0)
    }

    if(NCOL(res) != sum(!ina)){
      res <- res[,!ina, drop = FALSE]
    }
    id <- id[!ina]
  }
  Slist <- lapply(1:nv, function(i){
    sample_estim(res[, id == i, drop = FALSE], mse = mse)
  })
  cov_mat <- bdiag(Slist)

  P <- commat(nv, p)
  if(any(ina)){
    P <- P[!ina, , drop = FALSE]
    P <- P[, colSums(P)!=0, drop = FALSE]
  }
  t(P)%*%cov_mat%*%P
}
