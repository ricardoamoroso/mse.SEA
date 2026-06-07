# =================================================================
# conditioning.R — Initialization diagnostics for MSE app
#
# Purpose:
#   Given species parameters and a candidate E_init, show:
#   1. E_mmsy   — analytic multispecies MSY effort
#   2. E_implied — effort back-calculated from stated depletions (LS)
#   3. E_init   — user-supplied historical effort (editable)
#
#   For each of these three effort levels, compute the implied
#   equilibrium depletion per species and compare to stated Dep_ini.
#   The residuals reveal internal consistency of the inputs.
#
# All functions are pure (no Shiny dependencies) so they can be
# tested outside the app.
# =================================================================


# -----------------------------------------------------------------
# 1. CORE DIAGNOSTICS FUNCTION
# -----------------------------------------------------------------

#' Compute initialization diagnostics for a given E_init
#'
#' @param r        Vector of mean growth rates (length S)
#' @param K        Vector of mean carrying capacities (length S)
#' @param q        Vector of catchabilities (length S)
#' @param dep_ini  Vector of stated depletions B/K (length S)
#' @param E_init   Candidate historical effort (scalar)
#' @param species_names Character vector of species names
#' @return A list with:
#'   \item{E_mmsy}{Analytic MMSY effort}
#'   \item{E_implied}{Back-calculated effort from depletions}
#'   \item{E_init}{The supplied E_init}
#'   \item{table}{data.frame with per-species diagnostics}
#'   \item{consistency_flag}{Logical: are residuals all < 0.10?}
conditioning_diagnostics <- function(r, K, q, dep_ini,
                                     E_init, species_names = NULL) {
  S <- length(r)
  if (is.null(species_names)) species_names <- paste0("sp", seq_len(S))

  # --- Reference efforts ---
  mmsy_res  <- mmsy_shared_effort(r, K, q)
  E_mmsy    <- mmsy_res$E_mmsy
  E_implied <- estimate_shared_E_from_depletion(r, q, dep_ini)

  # --- Implied depletion at each effort level ---
  # Schaefer equilibrium: d = 1 - q*E/r  (clamped to [0,1])
  dep_at_E <- function(E) clamp(1 - q * E / r, 0, 1)

  dep_at_Emmsy   <- dep_at_E(E_mmsy)
  dep_at_Eimpl   <- dep_at_E(E_implied)
  dep_at_Einit   <- dep_at_E(E_init)

  # --- Residuals: stated - implied at E_init ---
  residuals <- dep_ini - dep_at_Einit

  # --- Consistency flag: all residuals within ±0.10 ---
  consistency_flag <- all(abs(residuals) <= 0.10)

  # --- Per-species table ---
  tbl <- data.frame(
    Species          = species_names,
    Dep_stated       = round(dep_ini,        3),
    Dep_at_Emmsy     = round(dep_at_Emmsy,   3),
    Dep_at_Eimplied  = round(dep_at_Eimpl,   3),
    Dep_at_Einit     = round(dep_at_Einit,   3),
    Residual         = round(residuals,       3),
    stringsAsFactors = FALSE
  )

  list(
    E_mmsy           = E_mmsy,
    E_implied        = E_implied,
    E_init           = E_init,
    table            = tbl,
    consistency_flag = consistency_flag
  )
}


# -----------------------------------------------------------------
# 2. CONSISTENCY MESSAGE
# -----------------------------------------------------------------

#' Generate a human-readable consistency message
#'
#' @param diag  Output from conditioning_diagnostics()
#' @return HTML string suitable for renderUI()
conditioning_message <- function(diag) {
  tbl  <- diag$table
  bad  <- tbl[abs(tbl$Residual) > 0.10, "Species"]
  n_bad <- length(bad)

  if (diag$consistency_flag) {
    HTML("<span style='color:#2ca02c; font-weight:bold;'>
         &#10003; Good consistency — all species depletions are within
         &plusmn;0.10 of the implied depletion at E_init.</span>")
  } else {
    spp_list <- paste(bad, collapse = ", ")
    HTML(sprintf(
      "<span style='color:#d62728; font-weight:bold;'>
       &#9888; Inconsistency detected for: %s.<br>
       </span>
       <span style='color:#555;'>
       The stated depletion differs by more than 0.10 from what
       E_init would imply at equilibrium. Consider adjusting E_init
       or reviewing the depletion estimates for these species.
       </span>",
      spp_list
    ))
  }
}


# -----------------------------------------------------------------
# 3. EFFORT REFERENCE TABLE
# -----------------------------------------------------------------

#' Summary table of the three reference effort levels
#'
#' @param diag  Output from conditioning_diagnostics()
#' @return data.frame with one row per effort level
effort_reference_table <- function(diag) {
  data.frame(
    Effort_level = c("E_mmsy (analytic)",
                     "E_implied (from depletions)",
                     "E_init (your input)"),
    Value        = round(c(diag$E_mmsy,
                           diag$E_implied,
                           diag$E_init), 6),
    Description  = c(
      "Effort that maximises multispecies yield at equilibrium",
      "Effort back-calculated from stated depletions (least squares)",
      "Historical effort you supply — used to start the simulation"
    ),
    stringsAsFactors = FALSE
  )
}
