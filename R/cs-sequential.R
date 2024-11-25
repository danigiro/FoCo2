#' Cross-sectional sequential reconciliation-combination
#'
#' This function applies a sequential method that first reconciles the base forecasts
#' from each expert to satisfy the linear constraints, and then combines the reconciled
#' forecasts obtained so far. [cssrc] may be applied only in 'balanced' cases (e.g.,
#' \eqn{n_j = n} \eqn{\forall j}, see Girolimetto and Di Fonzo, 2024)
#'
#' @usage cssrc(base, fc = "sa", comb = "ols", res = NULL, mse = TRUE, shrink = TRUE,
#'       nnw = FALSE, factorized = FALSE, ...)
#'
#' @inheritParams csocc
#' @param fc  A string specifying the combination method:
#'   \itemize{
#'   \item "\code{sa}" - (\emph{default}) simple average (equal weights).
#'   \item "\code{var}" - (uses \code{res}) weights derived from the inverse
#'   of forecasts variances/MSE as proposed by Bates and Granger (1969).
#'   \item "\code{cov}" - (uses \code{res}) weights derived using the whole
#'   forecast error covariance matrix, as proposed by Newbold and Granger (1974).
#'   }
#' @param comb A string specifying the reconciliation method: \code{"ols"}, \code{"wls"},
#' \code{"shr"}, \code{"sam"} (see \href{https://danigiro.github.io/FoReco/}{\code{FoReco}}).
#' @param mse If \code{TRUE} (\emph{default}) the residuals used to compute the covariance
#' matrix are not mean-corrected.
#' @param shrink  If \code{TRUE} (\emph{default}), the covariance matrix
#' for \code{fc = "cov"} is shrunk.
#' @param nnw If \code{TRUE} for \code{fc = "cov"}, the weights are constrained to be
#' non-negative (Conflitti et al., 2015). The \emph{default} is \code{FALSE}.
#' @param factorized Value to be passed to the \code{\link[quadprog:solve.QP]{quadprog::solve.QP}},
#' only when \code{nnw = TRUE}.
#' @param ... Arguments passed on to
#' \href{https://danigiro.github.io/FoReco/reference/csrec.html}{\code{FoReco::csrec}}
#' (e.g., \code{agg_mat} or \code{cons_mat}).
#'
#' @inherit csocc return
#'
#' @references
#' Bates, J. and Granger, C. W. (1969), The combination of forecasts,
#' \emph{Operations Research Quarterly}, 20, 451–468. \doi{10.1057/jors.1969.103}.
#'
#' Conflitti, C., De Mol, C., and Giannone, D. (2015), Optimal combination of survey
#' forecasts. \emph{International Journal of Forecasting}, 31(4), 1096–1103.
#' \doi{10.1016/j.ijforecast.2015.03.009}.
#'
#' Girolimetto, D. and Di Fonzo, T. (2024), Coherent forecast combination for linearly
#' constrained multiple time series, \emph{mimeo}.
#'
#' Newbold, P. and Granger, C. W. (1974), Experience with forecasting
#' univariate time series and the combination of forecasts,
#' \emph{Journal of the Royal Statistical Society, A}, 137, 131–146.
#' \doi{10.2307/2344546}
#'
#' @family Sequential coherent combination
#'
#' @examples
#' set.seed(123)
#' # (2 x 3) base forecasts matrix (simulated), expert 1
#' base1 <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
#' # (10 x 3) in-sample residuals matrix (simulated), expert 1
#' res1 <- t(matrix(rnorm(n = 30), nrow = 3))
#'
#' # (2 x 3) base forecasts matrix (simulated), expert 2
#' base2 <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
#' # (10 x 3) in-sample residuals matrix (simulated), expert 2
#' res2 <- t(matrix(rnorm(n = 30), nrow = 3))
#'
#' # Base forecasts' and residuals' lists
#' base <- list(base1, base2)
#' res <- list(res1, res2)
#'
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' reco <- cssrc(base = base, agg_mat = A, comb = "wls", res = res, fc = "sa")
#'
#' # Zero constraints matrix for Z - X - Y = 0
#' C <- t(c(1, -1, -1))
#' reco <- cssrc(base = base, cons_mat = C, comb = "wls", res = res, fc = "sa") # same results
#'
#' # WARNING!
#' reco_v <- cssrc(base = base, agg_mat = A, comb = "wls", res = res, fc = "var")
#' round(C %*% t(reco_v), 3) # Incoherent forecasts
#'
#' @export
cssrc <- function(base, fc = "sa", comb = "ols", res = NULL, mse = TRUE, shrink = TRUE,
                  nnw = FALSE, factorized = FALSE, ...){

  ina <- sapply(base, function(bmat){
    is.na(colSums(bmat))
  })

  if(any(ina)){
    cli_abort("{.arg cssrc} does not work with unbalanced panel of forecasts.", call = NULL)
  }

  reco <- lapply(1:length(base), function(j){
    FoReco::csrec(base = base[[j]],
                  res = res[[j]],
                  mse = mse,
                  comb = comb,
                  ...)
  })

  w <- .weights(y = reco, fc = fc, res = res, mse = mse, shrink = shrink, nnw = nnw)
  out <- .esemble(reco, weights = w)
  if(NCOL(out) == 1){
    out <- rbind(out)
  }
  attr_par <- as.list(attr(reco[[1]], "FoReco"))
  attr_par <- attr_par[!(names(attr_par)%in%c("info", "rfun"))]
  attr_par$rfun <- "cssrc"
  attr_par$info <- lapply(reco, function(x) attr(x, "FoReco")$info)
  attr_par$fc <- fc
  attr_par$reco <- reco
  attr_par$weights <- w
  attr(out, "FoReco") <- list2env(attr_par)
  dimnames(out) <- dimnames(reco[[1]])
  return(out)
}


#' Cross-sectional sequential combination-reconciliation
#'
#' This function performs a two-step process designed to first combine
#' forecasts from multiple models or experts and then apply reconciliation
#' techniques to ensure coherence.
#'
#' @usage csscr(base, fc = "sa", comb = "ols", res = NULL, mse = TRUE, shrink = TRUE,
#'       nnw = FALSE, factorized = FALSE, ...)
#'
#' @param comb A string specifying the reconciliation method: \code{"ols"}, \code{"wls"},
#' \code{"shr"}, \code{"sam"} (see \href{https://danigiro.github.io/FoReco/}{\code{FoReco}}).
#' If \code{comb = "none"}, no reconciliation is performed and the combined forecasts are
#' directly returned.
#'
#' @inheritParams cssrc
#'
#' @inherit csocc return
#' @inherit cssrc references
#' @family Sequential coherent combination
#'
#' @examples
#' set.seed(123)
#' # (2 x 3) base forecasts matrix (simulated), expert 1
#' base1 <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
#' # (10 x 3) in-sample residuals matrix (simulated), expert 1
#' res1 <- t(matrix(rnorm(n = 30), nrow = 3))
#'
#' # (2 x 3) base forecasts matrix (simulated), expert 2
#' base2 <- matrix(rnorm(6, mean = c(20, 10, 10)), 2, byrow = TRUE)
#' # (10 x 3) in-sample residuals matrix (simulated), expert 2
#' res2 <- t(matrix(rnorm(n = 30), nrow = 3))
#'
#' # Base forecasts' and residuals' lists
#' base <- list(base1, base2)
#' res <- list(res1, res2)
#'
#' # Aggregation matrix for Z = X + Y
#' A <- t(c(1,1))
#' reco <- csscr(base = base, agg_mat = A, comb = "wls", res = res, fc = "sa")
#'
#' # Zero constraints matrix for Z - X - Y = 0
#' C <- t(c(1, -1, -1))
#' reco <- csscr(base = base, cons_mat = C, comb = "wls", res = res, fc = "sa") # same results
#'
#' # Incoherent combined forecasts
#' fc_comb <- csscr(base = base, comb = "none", fc = "sa")
#' round(C %*% t(fc_comb), 3) # Incoherent forecasts
#' @export
csscr <- function(base, fc = "sa", comb = "ols", res = NULL, mse = TRUE, shrink = TRUE,
                  nnw = FALSE, factorized = FALSE, ...){
  base <- lapply(base, rbind)
  vec_names <- colnames(base[[1]])
  w <- .weights(y = base, fc = fc, res = res, mse = mse, shrink = shrink, nnw = nnw)

  base <- .esemble(base, weights = w)
  if(!is.null(res)){
    res <- .esemble(res, weights = w)
  }

  if(comb == "none"){
    reco <- rbind(base)
    if(is.null(vec_names)){
      colnames(reco) <- paste0("s-", 1:NCOL(reco))
    }else{
      colnames(reco) <- vec_names
    }
    rownames(reco) <- paste0("h-", 1:NROW(reco))
  }else{
    reco <- FoReco::csrec(base = base, comb = comb, res = as.matrix(res), mse = mse, ...)
    attr(reco, "FoReco")$rfun <- "csscr"
    attr(reco, "FoReco")$fc <- "fc"
    attr(reco, "FoReco")$base <- base
    attr(reco, "FoReco")$weights <- w
  }

  return(reco)
}






