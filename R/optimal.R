#' Cross-sectional optimal coherent combination
#'
#' @param base A list of \eqn{p} (\eqn{h \times n}) numeric matrix or multivariate time series
#' (\code{mts} class) containing the base forecasts to be reconciled; \eqn{h} is the forecast
#' horizon, and \eqn{n} is the total number of time series (\eqn{n = n_a + n_b}).
#' @param agg_mat A (\eqn{n_a \times n_b}) numeric matrix representing the cross-sectional
#' aggregation matrix. It maps the \eqn{n_b} bottom-level (free)
#' variables into the \eqn{n_a} upper (constrained) variables.
#' @param cons_mat A (\eqn{n_a \times n}) numeric matrix representing the cross-sectional
#' zero constraints. It spans the null space for the reconciled forecasts.
#' @param block_diag A string specifying the block-diagonal structure of the covariance matrix:
#' \code{fr} reconciliation structure, \code{fc} combination structure, \code{none} full covariance matrix.
#' @param comb A string specifying the reconciliation method. For a complete list, see [cscov].
#' @param res A list of \eqn{p} (\eqn{N \times n}) optional numeric matrix containing the
#' in-sample residuals. This matrix is used to compute some covariance matrices.
#' @param approach A string specifying the approach used to compute the reconciled
#' forecasts. Options include:
#'   \itemize{
#'   \item "\code{proj}" (\emph{default}): Projection approach according to Byron (1978, 1979).
#'   \item "\code{strc}": Structural approach as proposed by Hyndman et al. (2011).
#'   }
#' @param nn A string specifying the algorithm to compute non-negative reconciled forecasts:
#'   \itemize{
#'   \item "\code{sntz}": heuristic "set-negative-to-zero" (Di Fonzo and Girolimetto, 2023).
#'   }
#' @param settings Not yet implemented
#' @param bounds Not yet implemented
#' @param ... Arguments passed on to \code{cscov} (\pkg{FoReco}).
#'
#' @returns A (\eqn{h \times n}) numeric matrix of cross-sectional reconciled forecasts.
#' @export
csocc <- function(base, agg_mat, cons_mat, block_diag = "fr",
                  comb = "ols", res = NULL, approach = "proj",
                  nn = NULL, settings = NULL, bounds = NULL, ...){

  block_diag <- match.arg(block_diag, c("fr", "fc", "none"))

  # Check if either 'agg_mat' or 'cons_mat' is specified
  if(missing(agg_mat) && missing(cons_mat)){
    cli_abort("Argument {.arg agg_mat} (or {.arg cons_mat}) is missing,
              with no default.", call = NULL)
  } else if(!missing(agg_mat)){
    tmp <- cstools(agg_mat = agg_mat)
  } else {
    tmp <- cstools(cons_mat = cons_mat)
  }

  n <- tmp$dim[["n"]]
  strc_mat <- tmp$strc_mat
  cons_mat <- tmp$cons_mat
  p <- length(base)

  # Check if 'base' is provided and its dimensions match with the data
  if(missing(base)){
    cli_abort("Argument {.arg base} is missing, with no default.", call = NULL)
  }

  base <- lapply(base, rbind)
  base <- do.call(cbind, base)

  if(NCOL(base) != n*p){
    cli_abort("Incorrect {.arg base} columns dimension.", call = NULL)
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
    cov_mat <- cscov(comb = comb, n = NCOL(base),
                     agg_mat = rbind(do.call(rbind, rep(list(strc_mat), p-1)), agg_mat),
                     res = res, ...)
  }


  if(NROW(cov_mat) != n*p | NCOL(cov_mat) != n*p){
    cli_abort(c("Incorrect covariance dimensions.",
                "i"="Check {.arg res} dimensions."), call = NULL)
  }

  reco_mat <- resemble(base = base,
                       cov_mat = cov_mat,
                       strc_mat = strc_mat,
                       cons_mat = cons_mat,
                       approach = approach,
                       nn = nn,
                       p = p,
                       bounds = bounds,
                       settings = settings)

  rownames(reco_mat) <- paste0("h-", 1:NROW(reco_mat))
  if(is.null(colnames(base))){
    colnames(reco_mat) <- paste0("s-", 1:NCOL(reco_mat))
  } else {
    colnames(reco_mat) <- colnames(base[,1:NCOL(reco_mat), drop = FALSE])
  }

  attr(reco_mat, "FoReco") <- list2env(list(info = attr(reco_mat, "info"),
                                            framework = "Cross-sectional",
                                            forecast_horizon = NROW(reco_mat),
                                            comb = comb,
                                            cs_n = n,
                                            rfun = "csocc"))
  attr(reco_mat, "info") <- NULL
  return(reco_mat)
}
