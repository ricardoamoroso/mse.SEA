# =================================================================
# dynamics.R — Core MSE simulator
# =================================================================

#' @import data.table
NULL

#' Run multispecies MSE simulation
#'
#' Simulates Schaefer population dynamics for S species under a
#' shared (or per-species) fleet, with a pluggable HCR.
#' The simulation is split into two periods:
#' \enumerate{
#'   \item \strong{Conditioning} (\code{nyears_cond}): population is driven by
#'         the conditioning specification (e.g. B/K sampling). No HCR.
#'   \item \strong{Projection} (\code{nyears_proj}): HCR is active.
#' }
#' Only projection years are returned in the output.
#'
#' @param nyears_proj  Number of projection years (HCR active)
#' @param nyears_cond  Number of conditioning years (default 0 for bk_range)
#' @param r            Vector of intrinsic growth rates (length S)
#' @param K            Vector of carrying capacities (length S)
#' @param q            Vector of catchability coefficients (length S)
#' @param qI           Vector of index catchability (NULL = q * 0.2)
#' @param species_names Optional character vector of species names
#' @param shared_fleet Logical; single effort for all species?
#' @param hcr          HCR specification: list(type, params)
#' @param conditioning List describing initial state. type = "bk_range":
#'   dep_lo and dep_hi vectors (length S); each sim draws
#'   Dep_ini ~ Uniform(dep_lo, dep_hi) per species.
#' @param cv_r         CV for parameter uncertainty on r
#' @param cv_K         CV for parameter uncertainty on K
#' @param cv_q         CV for parameter uncertainty on q
#' @param cv_qI        CV for parameter uncertainty on qI
#' @param process_cv   CV for process error on growth
#' @param rho_proc     AR1 autocorrelation in process error
#' @param obs_type     Observation model type (currently only "per_species")
#' @param obs_cv       CV for index observation error
#' @param assess_cv    CV for the "assessment". Shortcut to create error in B/K
#' @param Blim_frac    Fraction of K used as Blim threshold
#' @param Btrig_frac   Fraction of K used as Btrig threshold
#' @param nsims        Number of simulations
#' @param seed         Random seed
#' @param burn         Years to exclude from metrics (stabilisation)
#' @param om_id        Optional OM preset id (overrides uncertainty params)
#' @param om_presets   data.table of OM presets (required if om_id is given)
#' @return A data.table with columns sim, year, spp, B, C, I, E, r, K, q, qI.
#'   year is indexed from 1 = first projection year.
#' @export
simulate_multispp <- function(
    nyears_proj   = 60,
    nyears_cond   = 0,
    r, K, q,
    qI            = NULL,
    species_names = NULL,
    shared_fleet  = TRUE,
    hcr           = list(type = "constant_effort", params = list(E = 0.2)),
    conditioning  = list(
      type   = "bk_range",
      dep_lo = NULL,   # vector length S, defaults to rep(0.3, S)
      dep_hi = NULL    # vector length S, defaults to rep(0.8, S)
    ),
    cv_r          = 0.2,
    cv_K          = 0.2,
    cv_q          = 0.2,
    cv_qI         = 0.0,
    process_cv    = 0.1,
    rho_proc      = 0.0,
    obs_type      = c("per_species", "fleet_cpue", "cpue_plus_comp"),
    obs_cv        = 0.2,
    assess_cv     = 0.1,
    Blim_frac     = 0.2,
    Btrig_frac    = 0.4,
    nsims         = 100,
    seed          = 42,
    burn          = 10,
    E_init        = NULL,   # starting effort; NULL = E_mmsy from sim params
    om_spec       = NULL,   # om_spec object from om_presets.R
    om_id         = NULL,
    om_presets    = NULL
) {
  set.seed(seed)
  S      <- length(r)
  nyears <- nyears_cond + nyears_proj
  stopifnot(length(K) == S, length(q) == S)

  if (is.null(qI)) qI <- q
  stopifnot(length(qI) == S)

  if (is.null(species_names))
    species_names <- paste0("sp", seq_len(S))

  obs_type <- match.arg(obs_type)

  # --- Resolve conditioning ---
  cond_type <- conditioning$type %||% "bk_range"
  if (cond_type == "bk_range") {
    dep_lo <- conditioning$dep_lo %||% rep(0.3, S)
    dep_hi <- conditioning$dep_hi %||% rep(0.8, S)
    if (length(dep_lo) == 1) dep_lo <- rep(dep_lo, S)
    if (length(dep_hi) == 1) dep_hi <- rep(dep_hi, S)
    stopifnot(length(dep_lo) == S, length(dep_hi) == S)
    stopifnot(all(dep_lo >= 0), all(dep_hi <= 1), all(dep_lo <= dep_hi))
  } else {
    stop("conditioning type '", cond_type, "' not yet implemented.")
  }

  # --- Override uncertainty params from OM preset ---
  # Priority: om_spec object > om_id lookup > individual CV arguments
  if (!is.null(om_spec)) {
    # om_spec is an om_spec object from om_presets.R
    cv_r       <- om_spec$r_cv
    cv_K       <- om_spec$K_cv
    cv_q       <- om_spec$q_cv
    cv_qI      <- om_spec$qI_cv
    process_cv <- om_spec$proc_cv
    rho_proc   <- om_spec$proc_rho
    obs_cv     <- om_spec$obs_cv
    assess_cv  <- 0   # om_spec controls all uncertainty; no extra assess noise
  } else if (!is.null(om_id)) {
    # Legacy: look up om_id in a data.frame of presets
    if (is.null(om_presets))
      stop("om_id was provided but om_presets is NULL")
    om_row <- om_presets[om_presets$om_id == om_id, ]
    if (nrow(om_row) == 0L)
      stop("Unknown om_id: ", om_id)
    cv_r       <- om_row$r_cv
    cv_K       <- om_row$K_cv
    cv_q       <- om_row$q_cv
    if (!is.null(om_row$qI_cv)) cv_qI <- om_row$qI_cv
    process_cv <- om_row$proc_cv
    rho_proc   <- om_row$proc_rho
    obs_cv     <- om_row$obs_cv
  }

  # --- Collapse dep range to point for fully deterministic OM ---
  # When all CVs are 0, spreading sims across dep_lo..dep_hi creates
  # spurious uncertainty in the plot. Collapse to midpoint.
  if (cv_r == 0 && cv_K == 0 && cv_q == 0 && process_cv == 0 && obs_cv == 0) {
    dep_lo <- (dep_lo + dep_hi) / 2
    dep_hi <- dep_lo
  }

  out_list <- vector("list", nsims)

  draw_params <- function() {
    list(
      r  = rlnorm_cv(r,  cv_r),
      K  = rlnorm_cv(K,  cv_K),
      q  = rlnorm_cv(q,  cv_q),
      qI = rlnorm_cv(qI, cv_qI)
    )
  }

  gen_proc_eps <- function() {
    z   <- matrix(rnorm(nyears * S), nrow = nyears, ncol = S)
    eps <- z
    if (abs(rho_proc) > 0) {
      for (i in 2:nyears)
        eps[i, ] <- rho_proc * eps[i-1, ] + sqrt(1 - rho_proc^2) * z[i, ]
    }
    eps
  }

  for (sim in seq_len(nsims)) {
    par  <- draw_params()
    r_s  <- par$r; K_s <- par$K
    q_s  <- par$q; qI_s <- par$qI

    B  <- matrix(NA_real_, nrow = nyears, ncol = S)
    C  <- matrix(NA_real_, nrow = nyears, ncol = S)
    I  <- matrix(NA_real_, nrow = nyears, ncol = S)
    BK <- matrix(NA_real_, nrow = nyears, ncol = S)  # observed B/K history
    E  <- if (shared_fleet) rep(NA_real_, nyears) else matrix(NA_real_, nyears, S)

    # --- Initial state from conditioning ---
    B[1, ] <- runif(S, min = dep_lo, max = dep_hi) * K_s
    I0     <- qI_s * B[1, ]

    # --- Starting effort for E_prev at first projection year ---
    E_start <- if (!is.null(E_init)) E_init else mmsy_shared_effort(r_s, K_s, q_s)$E_mmsy


    # All types: fish at MMSY (or supplied E_nom) to generate I_hist
    # The B/K range anchors the initial state; process error spreads it
    E_nom <- conditioning$E_nom %||% {
      if (shared_fleet)
        rep(mmsy_shared_effort(r_s, K_s, q_s)$E_mmsy, S)
      else
        r_s / (2 * q_s)
    }
    E_nom <- rep(E_nom, length.out = S)

    proc_eps <- gen_proc_eps()

    for (t in 1:nyears) {
      in_projection <- t > nyears_cond

      # Observation model
      if (obs_type == "per_species") {
        if (obs_cv <= 0) {
          I[t, ] <- qI_s * B[t, ]
        } else {
          tau    <- sqrt(log(1 + obs_cv^2))
          I[t, ] <- qI_s * B[t, ] * exp(rnorm(S, mean = -0.5 * tau^2, sd = tau))
        }
      } else {
        stop("obs_type '", obs_type, "' not yet implemented.")
      }

      # Assessment shortcut: noisy depletion estimate seen by HCR
      if (assess_cv > 0) {
        sigma  <- sqrt(log(1 + assess_cv^2))
        BK_obs <- clamp((B[t, ] / K_s) * exp(rnorm(S, -0.5 * sigma^2, sigma)), 0, 1)
        #BK_obs <- clamp((B[t,] / K_s) * exp(rnorm(S, 0, sigma)), 0, 1)

      } else {
        BK_obs <- B[t, ] / K_s
      }
      BK[t, ] <- BK_obs

      # Effort: HCR only during projection
      if (in_projection) {
        state_t <- list(
          BK_obs       = BK_obs,
          BK_hist      = BK[1:t, , drop = FALSE],
          B_t          = B[t, ],
          I_t          = I[t, ],
          I0           = I0,
          K_t          = K_s,
          r            = r_s,
          E_prev       = if (t == nyears_cond + 1) E_start else if (shared_fleet) E[t-1] else E[t-1, ],
          shared_fleet = shared_fleet,
          S            = S,
          I_hist       = I[1:t, , drop = FALSE]
        )
        ctrl <- hcr_apply(hcr, state_t)
      } else {
        # conditioning period: nominal effort, full access, no quota
        ctrl <- list(
          E               = E_nom[1],
          prop_accessible = rep(1, S),
          quota           = rep(Inf, S)
        )
      }

      # Unpack control action
      Et              <- rep(as.numeric(ctrl$E), length.out = S)
      prop_accessible <- rep(as.numeric(ctrl$prop_accessible), length.out = S)
      quota           <- rep(as.numeric(ctrl$quota), length.out = S)

      if (shared_fleet) {
        E[t] <- ctrl$E
      } else {
        E[t, ] <- Et
      }

      # Population dynamics
      growth_det <- r_s * B[t, ] * (1 - B[t, ] / K_s)
      if (process_cv > 0) {
        sigma  <- sqrt(log(1 + process_cv^2))
        noise  <- proc_eps[t, , drop = FALSE]
        growth <- growth_det * exp(sigma * as.vector(noise) - 0.5 * sigma^2)
      } else {
        growth <- growth_det
      }

      C[t, ] <- pmin(pmax(0, q_s * Et * B[t, ] * prop_accessible), quota)
      B_next <- clamp(B[t, ] + growth - C[t, ], 1e-8, 3 * K_s)
      if (t < nyears) B[t+1, ] <- B_next
    }

    # --- Return all years with phase flag ---
    vec_by_row <- function(M) as.vector(t(M))
    E_all <- if (shared_fleet) E else E

    dt <- data.table::data.table(
      sim   = sim,
      year  = rep(seq_len(nyears), each = S),
      phase = rep(ifelse(seq_len(nyears) <= nyears_cond, "cond", "proj"), each = S),
      spp   = rep(species_names, times = nyears),
      B     = vec_by_row(B),
      C     = vec_by_row(C),
      I     = vec_by_row(I),
      r     = rep(r_s,  times = nyears),
      K     = rep(K_s,  times = nyears),
      q     = rep(q_s,  times = nyears),
      qI    = rep(qI_s, times = nyears)
    )
    if (shared_fleet) {
      dt[, E := rep(E_all, each = S)]
    } else {
      dt[, E := as.vector(t(E_all))]
    }
    out_list[[sim]] <- dt
  }

  res <- data.table::rbindlist(out_list)
  data.table::setattr(res, "config", list(
    nyears_proj  = nyears_proj,
    nyears_cond  = nyears_cond,
    Blim_frac    = Blim_frac,
    Btrig_frac   = Btrig_frac,
    shared_fleet = shared_fleet,
    hcr          = hcr,
    conditioning = conditioning,
    cv_r         = cv_r,
    cv_K         = cv_K,
    cv_q         = cv_q,
    cv_qI        = cv_qI,
    process_cv   = process_cv,
    rho_proc     = rho_proc,
    obs_cv       = obs_cv,
    assess_cv    = assess_cv,
    burn         = burn,
    om_id        = om_id
  ))
  res
}
