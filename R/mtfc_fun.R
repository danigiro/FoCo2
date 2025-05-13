mtfc <- function(approach, base, nn = NULL, bounds = NULL, ...){
  tsp(base) <- NULL # Remove ts

  class_base <- approach

  # Set class of 'base' to include 'approach' and reconcile
  class(approach) <- c(class(approach), class_base)
  mtfore <- .mtfc(approach = approach, base = base, nn = nn, bounds = bounds, ...)

  # Check if 'nn' is provided and adjust 'mtfore' accordingly
  if(!is.null(nn)){
    if(nn %in% c("osqp", TRUE)){
      nn <- "osqp"
    }

    if(!all(mtfore >= -sqrt(.Machine$double.eps), na.rm = TRUE)){
      class(approach)[length(class(approach))] <- nn
      mtfore <- .mtfc(approach = approach, base = base, nn = nn, mtfore = mtfore,
                          bounds = bounds, ...)
    }else if(!all(mtfore >= 0, na.rm = TRUE)){
      mtfore[mtfore < 0] <- 0
    }
  }

  if(!is.null(bounds)){
    nbid <- bounds[,1,drop = TRUE]

    checkb <- apply(mtfore, 1, function(x){
      idl <- any(x[nbid]<bounds[,2,drop = TRUE] - sqrt(.Machine$double.eps))
      idb <- any(x[nbid]>bounds[, 3, drop = TRUE] + sqrt(.Machine$double.eps))

      idl0 <- any(x[nbid]<bounds[,2,drop = TRUE])
      idb0 <- any(x[nbid]>bounds[, 3, drop = TRUE])
      c(any(c(idl, idb)), any(c(idl0, idb0)))
    })

    if(any(checkb[1,], na.rm = TRUE)){
      if(is.null(attr(bounds, "approach")) || attr(bounds, "approach") == "osqp"){
        attr(bounds, "approach") <- paste(approach, "osqp", sep = "_")
      }

      class(approach)[length(class(approach))] <- attr(bounds, "approach")
      mtfore <- .mtfc(approach = approach, base = base, nn = nn, mtfore = mtfore,
                          bounds = bounds, ...)
    }else if(any(checkb[2,], na.rm = TRUE)){
      mtfore <- t(apply(mtfore, 1, function(x){
        id <- x[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
        x[nbid][id] <- bounds[,2,drop = TRUE][id]

        id <- x[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
        x[nbid][id] <- bounds[, 3, drop = TRUE][id]
        x
      }))
    }
  }
  return(mtfore)
}

.mtfc <- function(approach, ...){
  UseMethod("mtfc", approach)
}

mtfc.proj <- function(base, n, cov_mat, p, ina, ...){
  # check input
  if(missing(base) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base} and {.arg cov_mat}.",
              call = NULL)
  }

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(n))[!ina, , drop = FALSE]

  if(NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  if(any(ina)){
    if(any(colSums(k_mat) == 0)){
      cli_abort("Each variable must have at least one base forecasts.", call = NULL)
    }
  }

  if(isDiagonal(cov_mat)){
    cov_inv <- Matrix::.sparseDiagonal(x = Matrix::diag(cov_mat)^(-1))
    cov_k <- Matrix::crossprod(k_mat, cov_inv)%*%k_mat
    cov_k_inv <- Matrix::.sparseDiagonal(x = Matrix::diag(cov_k)^(-1))
    basec <- Matrix::tcrossprod(base, cov_inv)%*%Matrix::tcrossprod(k_mat, cov_k_inv)
  }else{
    cov_k <- lin_sys(cov_mat, k_mat)
    k_cov_k <- methods::as(Matrix::crossprod(k_mat, cov_k), "CsparseMatrix")
    ls2 <- lin_sys(k_cov_k, t(cov_k))
    basec <- methods::as(Matrix::tcrossprod(base, ls2), "CsparseMatrix")
  }

  return(as.matrix(basec))
}

mtfc.osqp <- function(base, n, cov_mat, p, ina,
                      nn = NULL, id_nn = NULL, bounds = NULL,
                      mtfore = NULL, settings = NULL, ...){

  # check input
  if(missing(base) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base} and {.arg cov_mat}.",
              call = NULL)
  }

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(n))[!ina, , drop = FALSE]

  if(NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  if(!is.null(nn) & !is.null(mtfore)){
    id <- which(rowSums(mtfore < (-sqrt(.Machine$double.eps))) != 0)
    if(!is.null(bounds)){
      id_b <- which(apply(mtfore[, bounds[,1], drop = FALSE], 1,
                          function(x) any(x <= bounds[,2]) | any(x >= bounds[,3])))
      if(length(id_b) > 0){
        id <- sort(unique(c(id, id_b)))
      }
    }

    if(length(id) == 0){
      mtfore[mtfore < 0] <- 0
      return(mtfore)
    }
  } else {
    id <- 1:NROW(base)
  }

  # Linear constrains H = 0
  l <- NULL
  u <- NULL
  A <- NULL

  # P matrix
  if(isDiagonal(cov_mat)){
    cov_inv <- Diagonal(x = diag(cov_mat)^(-1))
    cov_k1 <- cov_inv%*%k_mat
  } else {
    cov_k1 <- lin_sys(cov_mat, k_mat)
  }
  P <- methods::as(Matrix::crossprod(k_mat, cov_k1), "CsparseMatrix")

  # nn constraints (only on the building block variables)
  if(!is.null(nn)){
    if(!(nn %in% c("osqp", TRUE, "proj_osqp"))){
      cli_warn("Non-negative combined forecasts obtained with osqp.", call = NULL)
    }
    A <- rbind(A, .sparseDiagonal(n))
    l <- c(l, rep(0, n))
    u <- c(u, rep(Inf, n))
  }

  # other constraints
  if(!is.null(bounds)){
    A <- rbind(A, Diagonal(n)[bounds[,1,drop = TRUE], ])
    l <- c(l, bounds[,2,drop = TRUE])
    u <- c(u, bounds[,3,drop = TRUE])
  }

  if(is.null(settings)){
    settings <- osqpSettings(
      verbose = FALSE,
      eps_abs = 1e-5,
      eps_rel = 1e-5,
      polish_refine_iter = 100,
      polish = TRUE
    )
  }

  # OSQP
  osqp_step <- apply(base[id, , drop = FALSE], 1, function(x){
    q <- (-1) * t(cov_k1) %*% as.vector(x)
    mtf <- solve_osqp(P, q, A, l, u, settings)

    # Fix a problem of osqp
    if(mtf$info$status_val == -4){
      u[u == Inf] <- max(x)*100
      mtf <- solve_osqp(P, q, A, l, u, settings)
    }

    out <- list()
    out$mtfore <- mtf$x

    if(mtf$info$status_val != 1){
      cli_warn(c("x"="OSQP failed: check the results.",
                 "i"="OSQP flag = {mtf$info$status_val}",
                 "i"="OSQP pri_res = {mtf$info$pri_res}"), call = NULL)
    }

    if(!is.null(bounds)){
      nbid <- bounds[,1,drop = TRUE]
      id <- out$mtfore[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
      out$mtfore[nbid][id] <- bounds[,2,drop = TRUE][id]

      id <- out$mtfore[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
      out$mtfore[nbid][id] <- bounds[, 3, drop = TRUE][id]
    }

    out$info <- c(mtf$info$obj_val, mtf$info$run_time, mtf$info$iter,
                  mtf$info$pri_res, mtf$info$status_val, mtf$info$status_polish)

    return(out)
  })
  osqp_step <- do.call("rbind", osqp_step)

  # Point reconciled forecasts
  if(!is.null(mtfore)){
    mtfore[id, ] <- do.call("rbind", osqp_step[, "mtfore"])
  }else{
    mtfore <- do.call("rbind", osqp_step[, "mtfore"])
  }

  if(!is.null(nn)){
    mtfore[which(mtfore <= sqrt(.Machine$double.eps))] <- 0
  }

  class(mtfore) <- setdiff(class(mtfore), "proj_osqp")

  info <- do.call("rbind", osqp_step[, "info"])
  colnames(info) <- c(
    "obj_val", "run_time", "iter", "pri_res",
    "status", "status_polish"
  )
  rownames(info) <- id
  attr(mtfore, "info") <- info
  return(mtfore)
}

mtfc.sntz <- function(base, reco, strc_mat, ...){

  if(missing(reco)){
    mtfore <- base
  }

  tol <- sqrt(.Machine$double.eps)
  mtfore[mtfore<tol] <- 0

  as.matrix(mtfore)
}

mtfc.sftb <- function(base, mtfore, bounds = NULL, ...){

  if(missing(mtfore)){
    mtfore <- base
  }

  if(is.null(bounds)){
    return(mtfore)
  }

  nbid <- bounds[,1,drop = TRUE]
  mtfore <- t(apply(mtfore, 1, function(x){
    id <- x[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
    x[nbid][id] <- bounds[,2,drop = TRUE][id]

    id <- x[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
    x[nbid][id] <- bounds[, 3, drop = TRUE][id]
    x
  }))

  as.matrix(mtfore)
}
