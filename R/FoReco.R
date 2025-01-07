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
      tryCatch(backsolve(chol(msx), mdx), error = function(cond){
        cli_warn("Nearest positive definite matrix transformation.", call = NULL)
        solve(nearPD(msx)$mat, mdx)
      })
    })
  })
  out[is.na(out)] <- 0
  return(out)
}

# Remove NA values (row and columns)
remove_na <- function(x){
  inax <- is.na(x)
  if(any(inax)){
    out <- stats::na.omit(x)
    if(NROW(out) == 0){
      x <- x[, !(colSums(!inax) == 0)]
      inax <- is.na(x)
      #out <- stats::na.omit(x)
    }
    row_na <- rowSums(inax)
    x <- x[!(rowSums(inax) == NCOL(x)), , drop = FALSE]
  }
  return(x)
}

# Sample covariance matrix
sample_estim <- function(x, mse = TRUE){
  if(mse){
    if(any(is.na(x))){
      x <- remove_na(x)
    }

    if(any(is.na(x))){
      if(is.vector(x)){
        n <- sum(!is.na(x))
      }else{
        n <- colSums(!is.na(x))
      }
      n <- tcrossprod(n, rep(1, length(n)))
      x[is.na(x)] <- 0
      crossprod(x) / pmin(t(n), n)
    }else{
      crossprod(x) / NROW(x)
    }
  }else{
    stats::var(x, na.rm = TRUE, use = "na.or.complete")
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
    Dmat <- Matrix::nearPD(Dmat)$mat
  }

  if(factorized){
    Dmat <- solve(chol(Dmat))
  }

  if(is.null(dvec)){
    dvec <- rep(0, p)
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

## Modify version of FoReco::shrink_estim to deal with NA
shrink_estim_na <- function(x, mse = TRUE){
  if(is.matrix(x) == TRUE && is.numeric(x) == FALSE){
    cli_abort("{.arg x} is not a numeric matrix.", call = NULL)
  }
  x <- remove_na(x)
  n <- colSums(!is.na(x))
  # Full
  covm <- sample_estim(x = x, mse = mse)

  # Target
  tar <- Diagonal(x = diag(covm))

  id <- colSums(!is.na(x))/max(n)

  if(all(n>3) & all(id>0.75)){
    n <- tcrossprod(n, rep(1, length(n)))
    n <- pmin(t(n), n)

    # Lambda
    xs <- scale(x, center = FALSE, scale = sqrt(diag(covm)))
    xs[is.nan(xs)] <- xs[is.na(xs)] <- 0
    #xs <- xs[stats::complete.cases(xs), ]
    vS <- (1/(n*(n-1))) * (crossprod(xs^2) - ((1 / n) * (crossprod(xs))^2))
    diag(vS) <- 0
    corm <- cov2cor(covm)
    corm[is.nan(corm)] <- 0
    diag(corm) <- diag(corm)-1
    corm <- corm^2
    lambda <- sum(vS) / sum(corm)
    if(is.nan(lambda)){
      lambda <- 1
    }
    lambda <- max(min(lambda, 1), 0)
  }else{
    lambda <- 1
  }

  # Shrinkage
  shrink_cov <- lambda * tar + (1 - lambda) * covm
  shrink_cov <- drop0(shrink_cov)
  attr(shrink_cov, "lambda") <- lambda

  return(shrink_cov)
}
