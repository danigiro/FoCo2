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

# find the bottom time series given strc_mat
find_bts <- function(strc_mat){
  strc_mat <- Matrix(strc_mat, sparse = TRUE)
  strc_mat@i[strc_mat@p[-1]] + 1
}
