# =================================================================
# mmsy.R — Multispecies Maximum Sustainable Yield under shared effort
# =================================================================

#' Compute MMSY under a shared effort fleet
#'
#' Given Schaefer dynamics and a single shared effort E across species:
#'   Equilibrium B_i*(E) = K_i * (1 - q_i * E / r_i), floored at 0
#'   Yield_i(E) = q_i * E * max(0, B_i*)
#'
#' E_mmsy is found numerically over [0, min(r/q)] — the interval where the
#' total yield curve is unimodal. The analytic closed-form formula
#' (sum(K*q) / (2*sum(K*q^2/r))) is incorrect when some species collapse
#' before the complex optimum, because negative biomasses cancel positive
#' yields from other species, pulling E_mmsy artificially low.
#'
#' @param r  Vector of intrinsic growth rates
#' @param K  Vector of carrying capacities
#' @param q  Vector of catchability coefficients
#' @param enforce_Emax Logical; kept for backwards compatibility, no longer used
#' @return Named list with E_mmsy, yield functions, and reference biomasses
#' @export
mmsy_shared_effort <- function(r, K, q, enforce_Emax = TRUE) {

  # Total yield at effort E, with collapsed species (B < 0) contributing zero
  Y_tmp <- function(E) {
    Beq_i <- pmax(0, K * (1 - q * E / r))
    sum(q * E * Beq_i)
  }

  # Find E_mmsy numerically; [0, min(r/q)] is unimodal so optimise() is reliable
  E_star <- optimise(Y_tmp, interval = c(0, max(r / q)), maximum = TRUE)$maximum

  Y_components <- function(E) {
    Beq_i <- pmax(0, K * (1 - q * E / r))
    q * E * Beq_i
  }

  Y_total <- function(E) sum(Y_components(E))


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
