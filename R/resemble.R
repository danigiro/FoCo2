resemble <- function(approach, base, nn = NULL, bounds = NULL, ...){
  tsp(base) <- NULL # Remove ts

  class_base <- approach

  # Set class of 'base' to include 'approach' and reconcile
  class(approach) <- c(class(approach), class_base)
  rmat <- .resemble(approach = approach, base = base, nn = nn, bounds = bounds, ...)

  # Check if 'nn' is provided and adjust 'rmat' accordingly
  if(!is.null(nn)){
    if(nn == "osqp"){
      nn <- paste(approach, nn, sep = "_")
    }

    if(!all(rmat >= -sqrt(.Machine$double.eps))){
      class(approach)[length(class(approach))] <- nn
      rmat <- .resemble(approach = approach, base = base, nn = nn, reco = rmat,
                        bounds = bounds, ...)
    } else if(!all(rmat >= 0)){
      rmat[rmat < 0] <- 0
    }
  }

  if(!is.null(bounds)){
    nbid <- bounds[,1,drop = TRUE]

    checkb <- apply(rmat, 1, function(x){
      idl <- any(x[nbid]<bounds[,2,drop = TRUE] - sqrt(.Machine$double.eps))
      idb <- any(x[nbid]>bounds[, 3, drop = TRUE] + sqrt(.Machine$double.eps))

      idl0 <- any(x[nbid]<bounds[,2,drop = TRUE])
      idb0 <- any(x[nbid]>bounds[, 3, drop = TRUE])
      c(any(c(idl, idb)), any(c(idl0, idb0)))
    })

    if(any(checkb[1,])){
      if(is.null(attr(bounds, "approach")) || attr(bounds, "approach") == "osqp"){
        attr(bounds, "approach") <- paste(approach, "osqp", sep = "_")
      }

      class(approach)[length(class(approach))] <- attr(bounds, "approach")
      rmat <- .resemble(approach = approach, base = base, nn = nn, reco = rmat,
                        bounds = bounds, ...)
    }else if(any(checkb[2,])){
      rmat <- t(apply(rmat, 1, function(x){
        id <- x[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
        x[nbid][id] <- bounds[,2,drop = TRUE][id]

        id <- x[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
        x[nbid][id] <- bounds[, 3, drop = TRUE][id]
        x
      }))
    }
  }
  return(rmat)
}

.resemble <- function(approach, ...){
  UseMethod("resemble", approach)
}

resemble.proj <- function(base, cons_mat, cov_mat, p, ina, ...){
  # check input
  if(missing(base) | missing(cons_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(NCOL(cons_mat)))[!ina, , drop = FALSE]

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

resemble.strc <- function(base, strc_mat, cov_mat, p, ina, ...){
  # check input
  if(missing(base) | missing(strc_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(is.null(strc_mat)){
    cli_abort("Please provide a valid {.arg agg_mat} for the structural approach.",
              call = NULL)
  }

  strc_matp <- Matrix::kronecker(rep(1, p), strc_mat)[!ina, , drop = FALSE]

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

resemble.proj_osqp <- function(base, cons_mat, cov_mat, p, ina,
                               nn = NULL, id_nn = NULL, bounds = NULL,
                               reco = NULL, settings = NULL, ...){

  # check input
  if(missing(base) | missing(cons_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg cons_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  k_mat <- Matrix::kronecker(rep(1, p), .sparseDiagonal(NCOL(cons_mat)))[!ina, , drop = FALSE]

  if(NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  if(is.null(id_nn)){
    id_nn <- rep(1, NCOL(cons_mat))
  }

  if(!is.null(nn) & !is.null(reco)){
    id <- which(rowSums(reco < (-sqrt(.Machine$double.eps))) != 0)
    if(!is.null(bounds)){
      id_b <- which(apply(reco[, bounds[,1], drop = FALSE], 1,
                          function(x) any(x <= bounds[,2]) | any(x >= bounds[,3])))
      if(length(id_b) > 0){
        id <- sort(unique(c(id, id_b)))
      }
    }

    if(length(id) == 0){
      reco[reco < 0] <- 0
      return(reco)
    }
  } else {
    id <- 1:NROW(base)
  }

  c <- ncol(cons_mat)
  r <- nrow(cons_mat)
  # Linear constrains H = 0
  l <- rep(0, r)
  u <- rep(0, r)
  A <- cons_mat

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
      cli_warn("Non-negative reconciled forecasts obtained with osqp.", call = NULL)
    }
    A <- rbind(A, .sparseDiagonal(c)[id_nn == 1, ])
    l <- c(l, rep(0, sum(id_nn)))
    u <- c(u, rep(Inf, sum(id_nn)))
  }

  # other constraints
  if(!is.null(bounds)){
    A <- rbind(A, Diagonal(c)[bounds[,1,drop = TRUE], ])
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
    rec <- solve_osqp(P, q, A, l, u, settings)

    # Fix a problem of osqp
    if(rec$info$status_val == -4){
      u[u == Inf] <- max(x)*100
      rec <- solve_osqp(P, q, A, l, u, settings)
    }

    out <- list()
    out$reco <- rec$x

    if(rec$info$status_val != 1){
      cli_warn(c("x"="OSQP failed: check the results.",
                 "i"="OSQP flag = {rec$info$status_val}",
                 "i"="OSQP pri_res = {rec$info$pri_res}"), call = NULL)
    }

    if(!is.null(bounds)){
      nbid <- bounds[,1,drop = TRUE]
      id <- out$reco[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
      out$reco[nbid][id] <- bounds[,2,drop = TRUE][id]

      id <- out$reco[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
      out$reco[nbid][id] <- bounds[, 3, drop = TRUE][id]
    }

    out$info <- c(rec$info$obj_val, rec$info$run_time, rec$info$iter,
                  rec$info$pri_res, rec$info$status_val, rec$info$status_polish)

    return(out)
  })
  osqp_step <- do.call("rbind", osqp_step)

  # Point reconciled forecasts
  if(!is.null(reco)){
    reco[id, ] <- do.call("rbind", osqp_step[, "reco"])
  }else{
    reco <- do.call("rbind", osqp_step[, "reco"])
  }

  if(!is.null(nn)){
    reco[which(reco <= sqrt(.Machine$double.eps))] <- 0
  }


  class(reco) <- setdiff(class(reco), "proj_osqp")

  info <- do.call("rbind", osqp_step[, "info"])
  colnames(info) <- c(
    "obj_val", "run_time", "iter", "pri_res",
    "status", "status_polish"
  )
  rownames(info) <- id
  attr(reco, "info") <- info
  return(reco)
}

resemble.strc_osqp <- function(base, strc_mat, cov_mat, p, ina,
                               nn = NULL, id_nn = NULL, bounds = NULL,
                               reco = NULL, settings = NULL, ...){
  # check input
  if(missing(base) | missing(strc_mat) | missing(cov_mat)){
    cli_abort("Mandatory arguments: {.arg base}, {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(is.null(strc_mat)){
    cli_abort("Please provide a valid {.arg agg_mat} for the structural approach.",
              call = NULL)
  }

  strc_matp <- Matrix::kronecker(rep(1, p), strc_mat)[!ina, , drop = FALSE]

  if(NROW(strc_matp) != NROW(cov_mat) | NCOL(base) != NROW(cov_mat)){
    cli_abort("The size of the matrices does not match.", call = NULL)
  }

  if(any(is.na(base))){
    ina <- is.na(base[1,])
    cov_mat <- cov_mat[!ina, !ina, drop = FALSE]
    strc_matp <- strc_matp[!ina, , drop = FALSE]
    base <- base[, !ina, drop = FALSE]

    check_na <- matrix(!ina, NROW(strc_mat), p)
    if(any(rowSums(check_na) == 0)){
      cli_abort("Each variable must have at least one base forecasts.", call = NULL)
    }
  }

  if(is.null(id_nn)){
    bts <- find_bts(strc_mat)
    id_nn <- rep(0, NROW(strc_mat))
    id_nn[bts] <- 1
  }

  if(!is.null(nn) & !is.null(reco)){
    id <- which(rowSums(reco < (-sqrt(.Machine$double.eps))) != 0)
    if(!is.null(bounds)){
      id_b <- which(apply(reco, 1, function(x) all(bounds[,1] <= x) & all(bounds[,2] >= x)))
      if(length(id_b) > 0){
        id <- sort(unique(c(id, id_b)))
      }
    }
    if(length(id) == 0){
      reco[reco < 0] <- 0
      return(reco)
    }
  } else {
    id <- 1:NROW(base)
  }

  r <- NROW(strc_mat)
  c <- NCOL(strc_mat)
  A <- NULL
  l <- NULL
  u <- NULL

  # P matrix and q1 vector
  if(isDiagonal(cov_mat)){
    Q <- Diagonal(x = diag(cov_mat)^(-1))
    P <- t(strc_matp) %*% Q %*% strc_matp
    q1 <- (-1) * t(Q %*% strc_matp)
  } else {
    Q <- lin_sys(cov_mat, strc_matp)
    P <- t(strc_matp) %*% Q
    q1 <- (-1) * t(Q)
  }

  if(isDiagonal(cov_mat)){
    cov_inv <- Diagonal(x = diag(cov_mat)^(-1))
    cov_Sp <- cov_inv%*%strc_matp
  } else {
    cov_Sp <- lin_sys(cov_mat, strc_matp)
  }
  P <- methods::as(Matrix::crossprod(strc_matp, cov_Sp), "CsparseMatrix")

  # nn constraints (only on the building block variables - bottom variables)
  if(!is.null(nn)){
    if(!(nn %in% c("osqp", TRUE, "strc_osqp"))){
      cli_warn("Non-negative reconciled forecasts obtained with osqp.", call = NULL)
    }
    A <- .sparseDiagonal(c)
    l <- rep(0, sum(c))
    u <- rep(Inf, sum(c))
  }

  # other constraints
  if(!is.null(bounds)){
    A <- rbind(A, strc_mat[bounds[,1,drop = TRUE], ,drop = FALSE])
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
    q <- (-1) * t(cov_Sp) %*% as.vector(x)
    rec <- solve_osqp(P, q, A, l, u, settings)

    # Fix a problem of osqp
    if(rec$info$status_val == -4){
      u[u == Inf] <- max(x)*100
      rec <- solve_osqp(P, q, A, l, u, settings)
    }

    out <- list()
    out$reco <- as.numeric(strc_mat %*% rec$x)

    if(rec$info$status_val != 1){
      cli_warn(c("x"="OSQP failed: check the results.",
                 "i"="OSQP flag = {rec$info$status_val}",
                 "i"="OSQP pri_res = {rec$info$pri_res}"), call = NULL)
    }

    if(!is.null(bounds)){
      nbid <- bounds[,1,drop = TRUE]
      id <- out$reco[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
      out$reco[nbid][id] <- bounds[,2,drop = TRUE][id]

      id <- out$reco[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
      out$reco[nbid][id] <- bounds[, 3, drop = TRUE][id]
    }

    out$info <- c(
      rec$info$obj_val, rec$info$run_time, rec$info$iter, rec$info$pri_res,
      rec$info$status_val, rec$info$status_polish
    )

    return(out)
  })
  osqp_step <- do.call("rbind", osqp_step)

  # Point reconciled forecasts
  if(!is.null(reco)){
    reco[id, ] <- do.call("rbind", osqp_step[, "reco"])
  }else{
    reco <- do.call("rbind", osqp_step[, "reco"])
  }
  if(!is.null(nn)){
    reco[which(reco <= sqrt(.Machine$double.eps))] <- 0
  }

  class(reco) <- setdiff(class(reco), "strc_osqp")

  info <- do.call("rbind", osqp_step[, "info"])
  colnames(info) <- c(
    "obj_val", "run_time", "iter", "pri_res",
    "status", "status_polish"
  )
  rownames(info) <- id
  attr(reco, "info") <- info
  return(reco)
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

resemble.sftb <- function(base, reco, strc_mat, id_nn = NULL, bounds = NULL, ...){
  # Check input
  if(missing(strc_mat)){
    cli_abort("Mandatory arguments: {.arg strc_mat} and {.arg cov_mat}.",
              call = NULL)
  }

  if(missing(reco)){
    reco <- base
  }

  if(is.null(bounds)){
    return(reco)
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

  nbid <- bounds[,1,drop = TRUE]
  reco <- t(apply(reco, 1, function(x){
    id <- x[nbid]<=bounds[,2,drop = TRUE]+sqrt(.Machine$double.eps)
    x[nbid][id] <- bounds[,2,drop = TRUE][id]

    id <- x[nbid]>=bounds[, 3, drop = TRUE]-sqrt(.Machine$double.eps)
    x[nbid][id] <- bounds[, 3, drop = TRUE][id]
    x
  }))

  bts <- reco[, id_nn == 1, drop = FALSE]
  as.matrix(bts %*% t(strc_mat))
}
