# =================================================================
# mmsy.R — Multispecies Maximum Sustainable Yield under shared effort
# =================================================================

#' Compute analytic MMSY under a shared effort fleet
#'
#' Given Schaefer dynamics and a single shared effort E across species:
#'   Equilibrium B_i*(E) = K_i * (1 - q_i * E / r_i)
#'   Yield_i(E) = q_i * E * B_i*
#'   E_mmsy = (sum K_i q_i) / (2 * sum K_i q_i^2 / r_i)
#'
#' @param r  Vector of intrinsic growth rates
#' @param K  Vector of carrying capacities
#' @param q  Vector of catchability coefficients
#' @param enforce_Emax Logical; cap E_mmsy at min(r/q) to avoid negative biomass
#' @return Named list with E_mmsy, yield functions, and reference biomasses
#' @export
mmsy_shared_effort <- function(r, K, q, enforce_Emax = TRUE) {
  num    <- sum(K * q)
  den    <- 2 * sum(K * q^2 / r)
  E_star <- if (den > 0) num / den else 0

  if (enforce_Emax) {
    E_max  <- min(r / q)
    E_star <- min(E_star, E_max)
  }

  Y_components <- function(E) {
    mat <- outer(E, K * q, `*`) - outer(E^2, K * (q^2) / r, `*`)
    mat[mat < 0] <- 0
    if (length(E) == 1) drop(mat) else mat
  }

  Y_total <- function(E) {
    Y <- Y_components(E)
    if (is.matrix(Y)) rowSums(Y) else sum(Y)
  }

  Bmsy_multi_vec        <- K * (1 - q * E_star / r)
  Bmsy_multi_vec[Bmsy_multi_vec < 0] <- 0

  list(
    E_mmsy                = max(E_star, 0),
    Y_components          = Y_components,
    Y_total               = Y_total,
    Bmsy_multi_vec        = Bmsy_multi_vec,
    Btot_msy_multi        = sum(Bmsy_multi_vec),
    Bmsy_individual_vec   = K / 2,
    Emsy_i_if_independent = r / (2 * q),
    MSY_i_if_independent  = r * K / 4
  )
}

#' Estimate shared effort from observed depletions
#'
#' Uses equilibrium relationship d_i = 1 - (q_i / r_i) * E
#' to find the E that best explains all depletions simultaneously (LS).
#'
#' @param r        Vector of growth rates
#' @param q        Vector of catchabilities
#' @param dep      Vector of observed depletions B/K in (0, 1)
#' @param weights  Optional species weights
#' @param dep_cap  Upper cap on depletion (default 0.99)
#' @param dep_floor Lower floor on depletion (default 0.01)
#' @return Scalar estimated shared effort
#' @export
estimate_shared_E_from_depletion <- function(r, q, dep,
                                             weights   = NULL,
                                             dep_cap   = 0.99,
                                             dep_floor = 0.01) {
  dep_use <- pmin(dep_cap, pmax(dep_floor, dep))
  x <- q / r
  y <- 1 - dep_use

  if (is.null(weights)) {
    num <- sum(x * y)
    den <- sum(x^2)
  } else {
    w   <- weights / mean(weights)
    num <- sum(w * x * y)
    den <- sum(w * x^2)
  }
  if (den <= 0) return(0)
  num / den
}
