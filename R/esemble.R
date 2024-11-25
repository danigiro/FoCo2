.esemble <- function(y, weights, res = NULL){
  if(NROW(weights) == length(y)){
    tmp <- simplify2array(y)
    tmp[is.na(tmp)] <- 0
    sapply(1:NCOL(weights), function(i) tmp[,i,]%*%weights[,i])
  }else{
    tmp <- do.call(cbind, y)
    tmp%*%weights
  }
}

.weights <- function(y = NULL, fc = "sa", res = NULL, mse = TRUE, shrink = TRUE,
                     factorized = FALSE, nnw = FALSE, ...){
  if(fc == "sa"){
    tmp <- simplify2array(y)
    w <- apply(tmp, 2, wfoco_sa)
  }else if(fc %in% c("bg", "var")){
    tmp <- simplify2array(res)
    w <- apply(tmp, 2, wfoco_var, mse = mse)
  }else if(fc %in% c("nb", "ng", "cov")){
    tmp <- simplify2array(res)
    w <- apply(tmp, 2, wfoco_cov, mse = mse, shrink = shrink, nn = nnw)
  }else{
    w <- extract_omega(fc = fc, p = length(y), n = NCOL(y[[1]]), res = res, ...)
    #cli::cli_abort("Not available fc.")
  }
  return(w)
}

wfoco_sa <- function(fc){
  xna <- !is.na(fc)
  xna <- apply(xna, 2, all)
  xna*(1/sum(xna))
}

wfoco_var <- function(res, mse = TRUE){
  idna <- colSums(!is.na(res))==0
  res <- res[, !idna, drop = FALSE]
  w <- apply(res, 2, function(x) ifelse(mse, sum(x^2, na.rm = TRUE)/sum(!is.na(x)),
                                        var(x, na.rm = TRUE)))
  w <- w^(-1)

  w <- w/sum(w)
  w_out <- rep(0, length(idna))
  w_out[!idna] <- w
  w_out
}

wfoco_cov <- function(res, mse = TRUE, shrink = TRUE, nn = FALSE, factorized = FALSE){
  idna <- colSums(!is.na(res))==0
  res <- res[, !idna, drop = FALSE]
  if(shrink){
    cov_mat <- FoReco::shrink_estim(res, mse = mse)
  }else{
    cov_mat <- sample_estim(res, mse = mse)
  }
  p <- NCOL(cov_mat)
  w <- lin_sys(cov_mat, rep(1, p))

  if(nn & any(w<0)){
    w <-  tryCatch({
      wqp_sys(Dmat = cov_mat, p = p,
              factorized = factorized, nearPD = FALSE, scale = FALSE)
    }, error = function(cond){ # e.g. err: cov_mat is not positive definite!
      tryCatch({
        wqp_sys(Dmat = cov_mat, p = p,
                factorized = factorized, nearPD = TRUE, scale = FALSE)
      }, error = function(cond){ # e.g. err: constraints are inconsistent, no solution!
        wqp_sys(Dmat = cov_mat, p = p,
                factorized = factorized, nearPD = TRUE, scale = TRUE)
      })
    })
    w[w<0] <- 0
  }

  w <- w/sum(w)
  w_out <- rep(0, length(idna))
  w_out[!idna] <- w
  w_out
}

extract_omega <- function(fc = "ols", res = NULL, p = NULL, n = NULL, agg_mat = NULL, comb = NULL, ...){
  comb <- fc

  if(!is.null(agg_mat)){
    tmp <- cstools(agg_mat = agg_mat)
    agg_mat <- tmp$agg_mat
    strc_mat <- tmp$strc_mat
  }


  if(is.null(p) || is.null(n)){
    if(is.null(res)){
      cli_abort("{.arg p} and/or {.arg n} are missing,
              with no default.", call = NULL)
    }else{
      n <- NCOL(res[[1]])
      p <- length(res)
    }
  }

  ina <- matrix(FALSE, n, p)

  # Compute covariance
  if(comb %in% c("wls", "shr", "sam")){
    res <- do.call(cbind, res)
    res <- res[, !as.vector(ina), drop = FALSE]
  }

  cov_mat <- cscov(comb = comb, n = ifelse(comb == "ols", n*p, n), matNA = ina, p = p,
                   agg_mat = rbind(do.call(rbind, rep(list(strc_mat), p-1)), agg_mat),
                   res = res, ...)

  if(NROW(cov_mat) != sum(!ina) | NCOL(cov_mat) != sum(!ina)){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} dimensions."), call = NULL)
  }

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(n))
  cov_k <- lin_sys(cov_mat, k_mat)
  k_cov_k <- methods::as(Matrix::crossprod(k_mat, cov_k), "CsparseMatrix")
  return(t(lin_sys(k_cov_k, t(cov_k))))
}
