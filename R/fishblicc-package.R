#' The 'fishblicc' package.
#'
#' @description Procedures to fit catch curves to length frequency data based
#'   using length intervals.
#'
#' @details The package supports fitting a Bayesian length interval catch curve,
#'   with selectivity modelled as an optional mixture of logistic and normal
#'   functions. The package is designed to be used with multiple gears with
#'   different selectivities fishing the same stock.
#'
#' @keywords internal
"_PACKAGE"
#' @name fishblicc-package
#' @aliases fishblicc
#' @useDynLib fishblicc, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.26.13. https://mc-stan.org
#'
NULL
