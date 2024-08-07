resemble <- function(approach, base, immutable = NULL, nn = NULL, ...){
  tsp(base) <- NULL # Remove ts

  class_base <- approach

  # Set class of 'base' to include 'approach' and reconcile
  class(approach) <- c(class(approach), class_base)
  rmat <- .resemble(approach = approach, base = base, nn = nn, ...)

  # Check if 'nn' is provided and adjust 'rmat' accordingly
  if(!is.null(nn)){
    if(nn == "osqp"){
      nn <- paste(approach, nn, sep = "_")
    }

    if(!all(rmat >= -sqrt(.Machine$double.eps))){
      class(approach)[length(class(approach))] <- nn
      rmat <- .resemble(approach = approach, base = base, nn = nn, reco = rmat, ...)
    } else if(!all(rmat >= 0)){
      rmat[rmat < 0] <- 0
    }
  }
  return(rmat)
}

.resemble <- function(approach, ...){
  UseMethod("resemble", approach)
}

resemble.proj <- function(base, cons_mat, cov_mat, p, ...){

  # check input
  if(missing(base) | missing(cons_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(NCOL(cons_mat)))

  if(NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }
  if(isDiagonal(cov_mat)){
    cov_inv <- Matrix::.sparseDiagonal(x = Matrix::diag(cov_mat)^(-1))
    cov_k <- Matrix::crossprod(k_mat, cov_inv)%*%k_mat
    cov_k_inv <- Matrix::.sparseDiagonal(x = Matrix::diag(cov_k)^(-1))
    base_comp <- Matrix::tcrossprod(base, cov_inv)%*%Matrix::tcrossprod(k_mat, cov_k_inv)
    lm_sx1 <- methods::as(cons_mat%*%Matrix::tcrossprod(cov_k_inv, cons_mat), "CsparseMatrix")
    lm_dx1 <- methods::as(Matrix::tcrossprod(cons_mat, base_comp), "CsparseMatrix")
    reco <- base_comp - t(Matrix::tcrossprod(cov_k_inv, cons_mat)%*%lin_sys(lm_sx1, lm_dx1))
  }else{
    cov_k <- lin_sys(cov_mat, k_mat)
    k_cov_k <- methods::as(Matrix::crossprod(k_mat, cov_k), "CsparseMatrix")
    cov_c <- lin_sys(k_cov_k, t(cons_mat))
    c_cov_c <- methods::as(cons_mat%*%cov_c, "CsparseMatrix")
    ls1 <- lin_sys(c_cov_c, cons_mat)
    ls2 <- lin_sys(k_cov_k, t(cov_k))
    base_comp <- methods::as(Matrix::tcrossprod(base, ls2), "CsparseMatrix")
    reco <- base_comp - t(cov_c%*%ls1%*%t(base_comp))
  }
  return(as.matrix(reco))
}

resemble.strc <- function(base, strc_mat, cov_mat, p, ...){
  # check input
  if(missing(base) | missing(strc_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(is.null(strc_mat)){
    cli_abort("Please provide a valid {.arg agg_mat} for the structural approach.",
              call = NULL)
  }

  strc_matp <- Matrix::kronecker(rep(1, p), strc_mat)

  if(NROW(strc_matp) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  # Point reconciled forecasts
  if(isDiagonal(cov_mat)){
    cov_mat_inv <- .sparseDiagonal(x = diag(cov_mat)^(-1))
    StWm <- Matrix::crossprod(strc_matp, cov_mat_inv)
    lm_sx1 <- methods::as(StWm %*% strc_matp, "CsparseMatrix")
    lm_dx1 <- methods::as(Matrix::tcrossprod(StWm, base), "CsparseMatrix")
    reco <- t(strc_mat %*% lin_sys(lm_sx1, lm_dx1))
    return(as.matrix(reco))
  } else {
    Q <- lin_sys(cov_mat, strc_matp)
    lm_sx1 <- methods::as(t(strc_matp) %*% Q, "CsparseMatrix")
    lm_dx1 <- methods::as(t(base %*% Q), "CsparseMatrix")
    reco <- t(strc_mat %*% lin_sys(lm_sx1, lm_dx1))
    return(as.matrix(reco))
  }
}

resemble.sntz <- function(base, reco, strc_mat, cov_mat, id_nn = NULL, settings = NULL, ...){
  # Check input
  if(missing(strc_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(missing(reco)){
    reco <- base
  }

  if(is.null(strc_mat)){
    cli_abort(c("Argument {.arg agg_mat} is missing. The {.strong sntz} approach
                is available only for hierarchical/groupped time series."), call = NULL)
  }

  if(is.null(id_nn)){
    bts <- find_bts(strc_mat)
    id_nn <- rep(0, NCOL(reco))
    id_nn[bts] <- 1
  }

  bts <- reco[, id_nn == 1, drop = FALSE]
  if(is.null(settings$type)){
    sntz_type <- "bu"
  }else{
    sntz_type <- settings$type
  }
  tol <- sqrt(.Machine$double.eps)
  switch(sntz_type,
         bu = {
           bts[bts<tol] <- 0
         })

  as.matrix(bts %*% t(strc_mat))
}
