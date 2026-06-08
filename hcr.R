# =================================================================
# hcr.R — Harvest Control Rules
#
# Structure:
#   1. Mathematical tools  — reusable formulas
#   2. HCR functions       — the actual rules, use tools above
#   3. Registry            — maps type names to HCR functions
#   4. Dispatcher          — single entry point for the simulator
#
# HCR CONTRACT:
# Every HCR function returns a named list:
#   list(
#     E               = scalar,      # recommended effort
#     prop_accessible = vector(S),   # fraction of biomass accessible per species
#     quota           = vector(S)    # catch quota per species (Inf = no quota)
#   )
# Catch in dynamics.R is always:
#   C = min(q * E * B * prop_accessible, quota)
#
# Implementation error is applied per HCR type to the relevant control lever.
# =================================================================


# -----------------------------------------------------------------
# 1. TOOLS
# -----------------------------------------------------------------

#' Hockey-stick shape: converts a depletion signal into effort
#'
#' Below d_lim  → minimum effort (default = closure)
#' Above d_trig → maximum effort
#' Between      → linear ramp
#'
#' @param d      Depletion signal in [0,1]
#' @param d_lim  Lower breakpoint (default 0.2)
#' @param d_trig Upper breakpoint (default 0.4)
#' @param Emax   Effort at full stock (default 1)
#' @param Emin   Effort at collapsed stock (default 0)
#' @export
hcr_hockeystick <- function(d, d_lim = 0.2, d_trig = 0.4,
                            Emax = 1, Emin = 0) {
  if (d <= d_lim)  return(Emin)
  if (d >= d_trig) return(Emax)
  Emin + (Emax - Emin) * (d - d_lim) / (d_trig - d_lim)
}

#' Aggregate a vector of species depletions into one control signal
#'
#' "min"           → controlled by the most depleted species
#' "weighted_mean" → controlled by the weighted average depletion
#'
#' @param depl_vec  Vector of B/K values, one per species
#' @param weights   Optional weights (e.g. by K or economic value)
#' @param type      Aggregation method
#' @export
agg_depletion <- function(depl_vec, weights = NULL,
                          type = c("min", "weighted_mean")) {
  type <- match.arg(type)
  if (type == "min") return(min(depl_vec))
  if (is.null(weights)) weights <- rep(1 / length(depl_vec), length(depl_vec))
  sum(weights * depl_vec)
}

#' Apply lognormal implementation error to a control value
#'
#' @param x     Recommended control value (scalar or vector)
#' @param bias  Systematic bias on log scale (positive = overshoot)
#' @param cv    Coefficient of variation for random component (0 = no noise)
#' @param n     Length of output (defaults to length of x)
impl_error <- function(x, bias = 0, cv = 0, n = length(x)) {
  if (bias == 0 && cv <= 0) return(x)
  sigma <- if (cv > 0) sqrt(log(1 + cv^2)) else 0
  x * exp(bias + rnorm(n, mean = -0.5 * sigma^2, sd = sigma))
}

#' Build a standardised control list (default: no closure, no quota)
#'
#' @param E               Effort scalar
#' @param S               Number of species
#' @param prop_accessible Per-species accessible fraction (default 1)
#' @param quota           Per-species quota (default Inf)
make_control <- function(E, S, prop_accessible = rep(1, S),
                         quota = rep(Inf, S)) {
  list(
    E               = as.numeric(E)[1],
    prop_accessible = as.numeric(prop_accessible),
    quota           = as.numeric(quota)
  )
}


# -----------------------------------------------------------------
# 2. HCR FUNCTIONS
# Each returns make_control(E, S, ...) — same structure always.
# Implementation error is applied inside each function to the
# relevant control lever.
# -----------------------------------------------------------------

# Ignore the stock signal — always fish at the same effort
hcr_fn_constant_effort <- function(params, state) {
  E <- params$E %||% 0
  E <- impl_error(E, bias = params$impl_bias %||% 0,
                  cv   = params$impl_cv   %||% 0)
  make_control(E, state$S)
}

# Respond to the MOST depleted species (most precautionary)
# Uses BK_obs from state (already has assess_cv noise applied in dynamics.R)
hcr_fn_B_hs_min <- function(params, state) {
  n_smooth <- params$n_smooth %||% 1   # 1 = no smoothing
  if (n_smooth > 1 && !is.null(state$BK_hist)) {
    n     <- nrow(state$BK_hist)
    idx   <- max(1, n - n_smooth + 1):n
    signal <- min(colMeans(state$BK_hist[idx, , drop = FALSE]))
  } else {
    signal <- min(state$BK_obs)
  }
  E <- hcr_hockeystick(signal,
                       d_lim  = params$d_lim  %||% 0.2,
                       d_trig = params$d_trig %||% 0.25,
                       Emax   = params$Emax   %||% 1,
                       Emin   = params$Emin   %||% 0)
  E <- impl_error(E, bias = params$impl_bias %||% 0,
                  cv   = params$impl_cv   %||% 0)
  make_control(E, state$S)
}

# Respond to the AVERAGE depletion across species (less precautionary)
# Uses BK_obs from state (already has assess_cv noise applied in dynamics.R)
hcr_fn_B_hs_mean <- function(params, state) {
  n_smooth <- params$n_smooth %||% 1
  if (n_smooth > 1 && !is.null(state$BK_hist)) {
    n     <- nrow(state$BK_hist)
    idx   <- max(1, n - n_smooth + 1):n
    signal <- mean(colMeans(state$BK_hist[idx, , drop = FALSE]))
  } else {
    signal <- mean(state$BK_obs)
  }
  E <- hcr_hockeystick(signal,
                       d_lim  = params$d_lim  %||% 0.2,
                       d_trig = params$d_trig %||% 0.25,
                       Emax   = params$Emax   %||% 1,
                       Emin   = params$Emin   %||% 0)
  E <- impl_error(E, bias = params$impl_bias %||% 0,
                  cv   = params$impl_cv   %||% 0)
  make_control(E, state$S)
}

# Respond to the MOST depleted species using the CPUE INDEX
# Uses observed index I/I0 instead of true B/K
hcr_fn_index_hs_min <- function(params, state) {
  signal <- min(clamp(state$I_t / state$I0, 0, 1))
  E <- hcr_hockeystick(signal,
                       d_lim  = params$d_lim  %||% 0.2,
                       d_trig = params$d_trig %||% 0.4,
                       Emax   = params$Emax   %||% 1,
                       Emin   = params$Emin   %||% 0)
  E <- impl_error(E, bias = params$impl_bias %||% 0,
                  cv   = params$impl_cv   %||% 0)
  make_control(E, state$S)
}

hcr_fn_index_hs_mean <- function(params, state) {
  signal <- mean(clamp(state$I_t / state$I0, 0, 1))
  E <- hcr_hockeystick(signal,
                       d_lim  = params$d_lim  %||% 0.2,
                       d_trig = params$d_trig %||% 0.4,
                       Emax   = params$Emax   %||% 1,
                       Emin   = params$Emin   %||% 0)
  E <- impl_error(E, bias = params$impl_bias %||% 0,
                  cv   = params$impl_cv   %||% 0)
  make_control(E, state$S)
}

# Spatial closure: protects a fraction of biomass per species
# prop_closed is a vector of length S (or scalar, recycled)
# Implementation error = poaching (fleet accesses more of closed area)
hcr_fn_spatial_closure <- function(params, state) {
  S           <- state$S
  prop_closed <- rep(params$prop_closed %||% 0.3, length.out = S)
  E_base      <- params$E_base %||% 1

  # Implementation error on prop_closed: positive bias = more poaching
  prop_closed_actual <- clamp(
    impl_error(prop_closed,
               bias = params$impl_bias %||% 0,
               cv   = params$impl_cv   %||% 0,
               n    = S),
    lo = 0, hi = 1
  )

  make_control(E  = E_base,
               S  = S,
               prop_accessible = 1 - prop_closed_actual)
}

# Respond to the mean B/K of the n_low least resilient species (lowest r)
# n_low = 1 → least resilient only; n_low = S → same as B_hs_mean
hcr_fn_B_hs_low_r <- function(params, state) {
  n_low  <- min(params$n_low %||% 1, state$S)
  idx    <- order(state$r)[1:n_low]   # indices of n_low lowest-r species
  signal <- mean(state$BK_obs[idx])
  E <- hcr_hockeystick(signal,
                       d_lim  = params$d_lim  %||% 0.2,
                       d_trig = params$d_trig %||% 0.25,
                       Emax   = params$Emax   %||% 1,
                       Emin   = params$Emin   %||% 0)
  E <- impl_error(E, bias = params$impl_bias %||% 0,
                  cv   = params$impl_cv   %||% 0)
  make_control(E, state$S)
}

# Implementation error on effort (fleet fishes during closed season)
hcr_fn_seasonal_closure <- function(params, state) {
  prop_season <- params$prop_season %||% 0.25
  E_base      <- params$E_base      %||% 1
  E <- E_base * (1 - prop_season)
  E <- impl_error(E, bias = params$impl_bias %||% 0,
                  cv   = params$impl_cv   %||% 0)
  make_control(E, state$S)
}

# 2-over-3 rule: adjust effort based on ratio of mean CPUE
# last 2 years vs preceding 3 years
hcr_fn_2over3 <- function(params, state) {
  max_change <- params$max_change %||% 0.20
  agg_type   <- params$agg_type   %||% "min"
  E_prev     <- if (is.null(state$E_prev) || all(is.na(state$E_prev))) 1 else state$E_prev[1]
  I_hist     <- state$I_hist
  n          <- if (!is.null(I_hist)) nrow(I_hist) else 0

  if (n < 5) {
    E <- if (is.numeric(E_prev)) E_prev[1] else 1
    return(make_control(E, state$S))
  }

  n_low  <- params$n_low %||% 1
  # Normalise by first year so min/mean pick worst trend, not smallest species
  I_ref  <- pmax(I_hist[1, ], 1e-9)
  I_norm <- sweep(I_hist, 2, I_ref, "/")
  I_agg  <- apply(I_norm, 1, function(row) {
    if (agg_type == "min")        min(row)
    else if (agg_type == "mean")  mean(row)
    else if (agg_type == "low_r") mean(row[order(state$r)[1:min(n_low, state$S)]])
    else mean(row)
  })
  recent   <- mean(I_agg[(n-1):n])
  previous <- mean(I_agg[(n-4):(n-2)])

  if (previous <= 0) return(make_control(E_prev[1], state$S))

  ratio  <- recent / previous
  change <- clamp(ratio, 1 - max_change, 1 + max_change)
  E      <- E_prev[1] * change
  E      <- impl_error(E, bias = params$impl_bias %||% 0,
                       cv   = params$impl_cv   %||% 0)
  make_control(E, state$S)
}

# Slope rule: adjust effort based on slope of log(CPUE) over last n_years
hcr_fn_slope <- function(params, state) {
  n_years    <- params$n_years    %||% 5
  lambda     <- params$lambda     %||% 1
  max_change <- params$max_change %||% 0.20
  agg_type   <- params$agg_type   %||% "min"
  E_prev     <- if (is.null(state$E_prev) || all(is.na(state$E_prev))) 1 else state$E_prev[1]
  I_hist     <- state$I_hist
  n          <- if (!is.null(I_hist)) nrow(I_hist) else 0

  if (n < n_years) {
    E <- if (is.numeric(E_prev)) E_prev[1] else 1
    return(make_control(E, state$S))
  }

  n_low    <- params$n_low %||% 1
  I_recent <- I_hist[(n - n_years + 1):n, , drop = FALSE]
  # Normalise by first year of window so min/mean pick worst trend
  I_ref    <- pmax(I_recent[1, ], 1e-9)
  I_norm   <- sweep(I_recent, 2, I_ref, "/")
  I_agg    <- apply(I_norm, 1, function(row) {
    if (agg_type == "min")        min(row)
    else if (agg_type == "mean")  mean(row)
    else if (agg_type == "low_r") mean(row[order(state$r)[1:min(n_low, state$S)]])
    else mean(row)
  })

  years <- seq_len(n_years)
  slope <- tryCatch(
    coef(lm(log(pmax(I_agg, 1e-6)) ~ years))[2],
    error = function(e) 0
  )

  change <- clamp(1 + lambda * slope, 1 - max_change, 1 + max_change)
  E      <- E_prev[1] * change
  E      <- impl_error(E, bias = params$impl_bias %||% 0,
                       cv   = params$impl_cv   %||% 0)
  make_control(E, state$S)
}


# -----------------------------------------------------------------
# 3. REGISTRY
# -----------------------------------------------------------------

hcr_registry <- list(
  constant_effort  = hcr_fn_constant_effort,
  B_hs_min         = hcr_fn_B_hs_min,
  B_hs_mean        = hcr_fn_B_hs_mean,
  B_hs_low_r       = hcr_fn_B_hs_low_r,
  index_hs_min     = hcr_fn_index_hs_min,
  index_hs_mean = hcr_fn_index_hs_mean,
  spatial_closure  = hcr_fn_spatial_closure,
  seasonal_closure = hcr_fn_seasonal_closure,
  two_over_three   = hcr_fn_2over3,
  slope_rule       = hcr_fn_slope
)


# -----------------------------------------------------------------
# 4. DISPATCHER
# -----------------------------------------------------------------

#' Apply the HCR for a given MP and system state
#'
#' Accepts either an mp_spec object (from mp_custom()) or a plain list
#' with elements \code{type} and \code{params}.
#'
#' @param mp    An mp_spec object OR list(type, params)
#' @param state Named list: B_t, I_t, I0, K_t, E_prev, shared_fleet, S, I_hist
#' @return Named list: E, prop_accessible (length S), quota (length S)
#' @export
hcr_apply <- function(mp, state) {
  hcr_type <- if (!is.null(mp$hcr_type)) mp$hcr_type else mp$type
  if (is.null(hcr_type)) stop("mp must have 'hcr_type' or 'type'")
  fn <- hcr_registry[[hcr_type]]
  if (is.null(fn)) stop("Unknown HCR type: ", hcr_type)
  fn(mp$params, state)
}
