#' Cross-sectional sequential reconciliation-combination
#'
#' @inheritParams csocc
#' @param fc  A string specifying the combination method
#' @param mse  todo
#' @param shrink  todo
#' @param ... Arguments passed on to \code{csrec} (\pkg{FoReco}).
#'
#' @inherit csocc return
#' @export
cssrc <- function(base, res = NULL, fc = "sa", mse = TRUE, shrink = TRUE, ...){
  reco <- lapply(1:length(base), function(j){
    FoReco::csrec(base = base[[j]],
                   res = res[[j]],
                   mse = mse,
                   ...)
  })

  w <- .weights(y = reco, fc = fc, res = res, mse = mse, shrink = shrink)
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
  attr(out, "FoReco") <- list2env(attr_par)
  dimnames(out) <- dimnames(reco[[1]])
  return(out)
}


#' Cross-sectional sequential combination-reconciliation
#'
#' @inheritParams cssrc
#'
#' @inherit csocc return
#' @export
csscr <- function(base, res = NULL, fc = "sa", mse = TRUE, shrink = TRUE, ...){
  base <- lapply(base, rbind)
  w <- .weights(y = base, fc = fc, res = res, mse = mse, shrink = shrink, ...)

  base <- .esemble(base, weights = w)
  if(!is.null(res)){
    res <- .esemble(res, weights = w)
  }
  reco <- FoReco::csrec(base = base, res = as.matrix(res), mse = mse, ...)

  attr(reco, "FoReco")$rfun <- "csscr"
  attr(reco, "FoReco")$fc <- "fc"
  attr(reco, "FoReco")$base <- base
  return(reco)
}






