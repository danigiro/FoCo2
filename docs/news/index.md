# Changelog

## FoCo2 0.1.4

- Every reconciliation function (`...`, `...`) now returns an object of
  the new S3 class `foreco`, defined in FoReco. The objects are built
  through FoReco’s exported `new_foreco_class()` constructor and
  therefore integrate seamlessly with FoReco’s
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  `components()` methods.

## FoCo2 0.1.3

CRAN release: 2026-03-13

- Fixed compatibility issues with
  [`osqp`](https://CRAN.R-project.org/package=osqp) 1.0. The package is
  now fully compatible with any version.
- Fixed bugs and improved stability.
- Fixed documentation.

## FoCo2 0.1.2

CRAN release: 2025-06-14

- Release on CRAN

## FoCo2 0.1.0

- Release on Github
