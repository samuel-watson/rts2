#' The 'rts2' package.
#'
#' @description Fitting of Log Gaussian Cox Process models using maximum likelihood including spectral and SPDE approximations 
#' to the latent field, as well as Bayesian model fitting through Stan.
#'
#' @docType package
#' @name rts2-package
#' @aliases rts2
#' @useDynLib rts2, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import lubridate
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'
NULL
