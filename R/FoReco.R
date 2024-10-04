## Function extracted from FoReco
# Solve a System of Equations: Robust function
lin_sys <- function(msx, mdx){
  if(NCOL(msx)>100){
    if(!is(msx, "symmetricMatrix")){
      msx <- forceSymmetric(msx)
      mdx <- methods::as(mdx, "CsparseMatrix")
    }
  }

  out <- tryCatch(solve(msx, mdx), error = function(cond){
    tryCatch(solve(qr(msx), mdx), error = function(cond){
      backsolve(chol(msx), mdx)
    })
  })
  out[is.na(out)] <- 0
  return(out)
}

# Remove NA values (row and columns)
remove_na <- function(x){
  out <- stats::na.omit(x)
  if(NROW(out) == 0){
    x <- x[, !(colSums(!is.na(x)) == 0)]
    out <- stats::na.omit(x)
  }
  return(out)
}

# Sample covariance matrix
sample_estim <- function(x, mse = TRUE){
  if(mse){
    if(any(is.na(x))){
      x <- remove_na(x)
    }
    crossprod(x) / NROW(x)
  }else{
    stats::var(x, na.rm = TRUE)
  }
}

# find the bottom time series given strc_mat
find_bts <- function(strc_mat){
  strc_mat <- Matrix(strc_mat, sparse = TRUE)
  strc_mat@i[strc_mat@p[-1]] + 1
}


## Function extracted from FoCo
# Solve a Quadratic problem: Robust function
wqp_sys <- function(Dmat, p, dvec = NULL, factorized = FALSE, nearPD = FALSE, scale = FALSE){
  if(nearPD){
    cov_mat <- Matrix::nearPD(Dmat)$mat
  }

  if(factorized){
    Dmat <- solve(chol(Dmat))
  }

  if(is.null(dvec)){
    dvec = rep(0, p)
  }

  if(scale){
    sc <- norm(Dmat,"2")
    Dmat <- Dmat/sc
  }

  tmp <- quadprog::solve.QP(Dmat = Dmat,
                            dvec = dvec,
                            Amat = cbind(rep(1, p),
                                         diag(p)),
                            bvec=c(1, rep(0, p)),
                            meq = 1, factorized = factorized)
  return(tmp$solution)
}
