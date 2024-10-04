# test cross-sectional reconciliation
if(require(testthat)){
  A <- matrix(c(1,1,1,1,
                1,1,0,0,
                0,0,1,1), 3, byrow = TRUE)
  set.seed(123)
  res1 <- matrix(rnorm(100*sum(dim(A))), 100, sum(dim(A)))
  base1 <- t(rnorm(sum(dim(A)), 1))
  res2 <- matrix(rnorm(100*sum(dim(A))), 100, sum(dim(A)))
  base2 <- t(rnorm(sum(dim(A)), 1))
  res3 <- matrix(rnorm(100*sum(dim(A))), 100, sum(dim(A)))
  base3 <- t(rnorm(sum(dim(A)), 1))
  C <- cbind(diag(NROW(A)), -A)
  comb <- "shr"
  base <- list(base1, base2, base3)
  res <- list(res1, res2, res3)
  test_that("Optimal cross-sectional coherent combination", {
    r1 <- csocc(base = base, agg_mat = A, comb = comb,
                res = res, approach = "strc")
    r2 <- csocc(base = base, agg_mat = A, comb = comb,
                res = res, approach = "proj")
    r3 <- csocc(base = base, agg_mat = A, comb = comb,
                res = res, approach = "strc_osqp")
    r4 <- csocc(base = base, agg_mat = A, comb = comb,
                res = res, approach = "proj_osqp")
    r5 <- csocc(base = base, cons_mat = C, comb = comb,
                res = res, approach = "strc")
    r6 <- csocc(base = base, cons_mat = C, comb = comb,
                res = res, approach = "proj")

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r1, r5, ignore_attr = TRUE)
    expect_equal(r1, r6, ignore_attr = TRUE)
    expect_equal(max(abs(C%*%t(r1))), 0)
  })

  test_that("Covariance check", {
    for(i in c("ols", "str", "wls", "shr", "shrbe", "shrbv", "sam", "sambe", "sambv")){
      expect_no_error({
        csocc(base = base, agg_mat = A, comb = i,
              res = res, approach = "proj")
      })
    }
  })

  base[[1]][1,NCOL(base1)] <- -10
  test_that("Optimal nonegative cross-sectional reconciliation", {
    r1 <- csocc(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", nn = "strc_osqp")
    r2 <- csocc(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", nn = "proj_osqp")
    r3 <- csocc(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", nn = "sntz")


    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(max(abs(C%*%t(r1))), 0)
    expect_equal(max(abs(C%*%t(r3))), 0)
  })

  base[[1]][1,1] <- NA
  base_err <- base
  base_err[[2]][1,1] <- NA
  test_that("Optimal cross-sectional coherent combination with NA", {
    r1 <- csocc(base = base, agg_mat = A, comb = comb,
                res = res, approach = "strc")
    r2 <- csocc(base = base, agg_mat = A, comb = comb,
                res = res, approach = "proj")
    r3 <- csocc(base = base, agg_mat = A, comb = comb,
                res = res, approach = "strc_osqp")
    r4 <- csocc(base = base, agg_mat = A, comb = comb,
                res = res, approach = "proj_osqp")
    r5 <- csocc(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", nn = "strc_osqp")
    r6 <- csocc(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", nn = "proj_osqp")
    r7 <- csocc(base = base, agg_mat = A, comb = comb, res = res,
                approach = "proj", nn = "sntz")

    expect_equal(r1, r2, ignore_attr = TRUE)
    expect_equal(r1, r3, ignore_attr = TRUE)
    expect_equal(r1, r4, ignore_attr = TRUE)
    expect_equal(r5, r6, ignore_attr = TRUE)
    expect_equal(max(abs(C%*%t(r1))), 0)
    expect_equal(max(abs(C%*%t(r5))), 0)
    expect_equal(max(abs(C%*%t(r7))), 0)
  })

  test_that("Covariance check with NA", {
    expect_no_error({
      csocc(base = base, agg_mat = A, comb = "ols",
            res = res, approach = "proj")
    })

    expect_no_error({
      csocc(base = base, agg_mat = A, comb = "str",
            res = res, approach = "proj")
    })

    expect_no_error({
      csocc(base = base, agg_mat = A, comb = "wls",
            res = res, approach = "proj")
    })

    expect_no_error({
      csocc(base = base, agg_mat = A, comb = "shr",
            res = res, approach = "proj")
    })

    expect_no_error({
      csocc(base = base, agg_mat = A, comb = "sam",
            res = res, approach = "proj")
    })

    expect_no_error({
      csocc(base = base, agg_mat = A, comb = "shrbe",
            res = res, approach = "proj")
    })

    expect_no_error({
      csocc(base = base, agg_mat = A, comb = "shrbv",
            res = res, approach = "proj")
    })

    expect_no_error({
      csocc(base = base, agg_mat = A, comb = "sambe",
            res = res, approach = "proj")
    })

    expect_no_error({
      csocc(base = base, agg_mat = A, comb = "sambe",
            res = res, approach = "proj")
    })
  })

}
