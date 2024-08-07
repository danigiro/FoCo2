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

.weights <- function(y = NULL, fc = "sa", res = NULL, mse = TRUE, shrink = TRUE, ...){
  if(fc == "sa"){
    tmp <- simplify2array(y)
    w <- apply(tmp, 2, wfoco_sa)
  }else if(fc == "bg"){
    tmp <- simplify2array(res)
    w <- apply(tmp, 2, wfoco_bg, mse = mse)
  }else if(fc %in% c("nb", "ng")){
    tmp <- simplify2array(res)
    w <- apply(tmp, 2, wfoco_nb, mse = mse, shrink = shrink)
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

wfoco_bg <- function(res, mse = TRUE){
  res <- na.omit(res)
  w <- apply(res, 2, function(x) ifelse(mse, sum(x^2)/length(x), var(x)))
  w <- w^(-1)
  w/sum(w)
}

wfoco_nb <- function(res, mse = TRUE, shrink = TRUE){
  if(shrink){
    cov_mat <- FoReco::shrink_estim(res, mse = mse)
  }else{
    cov_mat <- FoReco:::sample_estim(res, mse = mse)
  }

  w <- lin_sys(cov_mat, rep(1, NCOL(cov_mat)))
  w/sum(w)
}

extract_omega <- function(fc = "ols", res = NULL, p = NULL, n = NULL, agg_mat = NULL, comb = NULL, ...){
  fc <- strsplit(fc, "-")[[1]]
  if(length(fc) == 1){
    block_diag <- "none"
    comb <- fc
  }else{
    block_diag <- fc[1]
    comb <- fc[2]
  }

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
  cov_k <- lin_sys(cov_mat, k_mat)
  k_cov_k <- methods::as(Matrix::crossprod(k_mat, cov_k), "CsparseMatrix")
  return(t(lin_sys(k_cov_k, t(cov_k))))
}
