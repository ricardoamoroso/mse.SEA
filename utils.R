# Utility functions used across the package

#' Null coalescing operator
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' Sample from lognormal distribution using mean and CV
#' @param mean Vector of means
#' @param cv Coefficient of variation
#' @param n Number of samples (defaults to length of mean)
rlnorm_cv <- function(mean, cv, n = NULL) {
  if (is.null(n)) n <- length(mean)
  if (cv <= 0) return(rep_len(mean, n))
  sigma2 <- log(1 + cv^2)
  mu     <- log(mean) - 0.5 * sigma2
  rlnorm(n, meanlog = mu, sdlog = sqrt(sigma2))
}

#' Clamp values to [lo, hi]
#' @param x Numeric vector
#' @param lo Lower bound
#' @param hi Upper bound
clamp <- function(x, lo, hi) pmin(hi, pmax(lo, x))
