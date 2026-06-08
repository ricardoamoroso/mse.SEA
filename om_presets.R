# Operating Model specifications
# OMs define HOW uncertainty is applied in the simulation.
# They do NOT define where central parameter values come from
# (that is the job of conditioning or user CSV input).

#' Create a custom Operating Model specification
#'
#' Base constructor — all other OM functions are wrappers around this.
#'
#' @param r_cv        CV for intrinsic growth rate uncertainty (drawn once per sim)
#' @param K_cv        CV for carrying capacity uncertainty (drawn once per sim)
#' @param q_cv        CV for fishing catchability uncertainty (drawn once per sim)
#' @param qI_cv       CV for index catchability uncertainty (drawn once per sim)
#' @param dep_cv      CV for initial depletion uncertainty (drawn once per sim)
#' @param proc_cv     CV for process error magnitude (varies each year)
#' @param proc_rho    AR1 autocorrelation in process error (within species)
#' @param proc_rho_sp Cross-species correlation in process error (shared environment)
#' @param obs_cv      CV for index observation error
#' @param catch_cv    CV for catch observation error (proxy for IUU/misreporting)
#' @param beta_hyper  Hyperstability coefficient: 0 = linear CPUE, 1 = fully hyperstable
#' @param label       Display name for plots and tables
#' @return An object of class \code{om_spec}
#' @export
om_custom <- function(
  # Parameter uncertainty (drawn once per simulation)
  r_cv        = 0.00,
  K_cv        = 0.00,
  q_cv        = 0.00,
  qI_cv       = 0.00,
  # Process error
  proc_cv     = 0.00,
  proc_rho    = 0.00,
  proc_rho_sp = 0.00,
  # Observation error
  obs_cv      = 0.00,
  catch_cv    = 0.00,
  # Structural assumptions
  beta_hyper  = 0.00,
  assess_cv = 0.00,
  # Label
  label       = "Custom OM"
) {
  # Basic input validation
  stopifnot(
    r_cv        >= 0, K_cv    >= 0,
    q_cv        >= 0, qI_cv   >= 0,
    proc_cv     >= 0,
    proc_rho    >= 0, proc_rho    <= 1,
    proc_rho_sp >= 0, proc_rho_sp <= 1,
    obs_cv      >= 0, catch_cv    >= 0,
    beta_hyper  >= 0, beta_hyper  <= 1
  )

  structure(
    list(
      r_cv        = r_cv,
      K_cv        = K_cv,
      q_cv        = q_cv,
      qI_cv       = qI_cv,
      proc_cv     = proc_cv,
      proc_rho    = proc_rho,
      proc_rho_sp = proc_rho_sp,
      obs_cv      = obs_cv,
      catch_cv    = catch_cv,
      beta_hyper  = beta_hyper,
     assess_cv = assess_cv,
     label       = label
    ),
    class = "om_spec"
  )
}

om_custom <- function(r_cv = 0, K_cv = 0, q_cv = 0, qI_cv = 0,
                      proc_cv = 0, proc_rho = 0, obs_cv = 0,
                      assess_cv = 0, label = "Custom") {
  structure(list(r_cv = r_cv, K_cv = K_cv, q_cv = q_cv,
                 qI_cv = qI_cv, proc_cv = proc_cv, proc_rho = proc_rho,
                 obs_cv = obs_cv, assess_cv = assess_cv, label = label),
            class = "om_spec")
}

#' Print method for om_spec
#' @export
print.om_spec <- function(x, ...) {
  cat("Operating Model:", x$label, "\n")
  cat("  --- Parameter uncertainty ---\n")
  cat("  r_cv:        ", x$r_cv,        "\n")
  cat("  K_cv:        ", x$K_cv,        "\n")
  cat("  q_cv:        ", x$q_cv,        "\n")
  cat("  qI_cv:       ", x$qI_cv,       "\n")
  cat("  --- Process error ---\n")
  cat("  proc_cv:     ", x$proc_cv,     "\n")
  cat("  proc_rho:    ", x$proc_rho,    "\n")
  cat("  proc_rho_sp: ", x$proc_rho_sp, "\n")
  cat("  --- Observation error ---\n")
  cat("  obs_cv:      ", x$obs_cv,      "\n")
  cat("  catch_cv:    ", x$catch_cv,    "\n")
  cat("  --- Structural ---\n")
  cat("  beta_hyper:  ", x$beta_hyper,  "\n")
  invisible(x)
}

# ------------------------------------------------------------------
# Built-in OM presets
# ------------------------------------------------------------------

om_deterministic <- om_custom(
  label = "Deterministic — no uncertainty"
)

om_baseline <- om_custom(
  proc_cv = 0.20,
  label   = "Baseline — process error only"
)

om_param_uncertainty <- om_custom(
  r_cv    = 0.10,
  K_cv    = 0.20,
  proc_cv = 0.20,
  label   = "Parameter + process uncertainty"
)

om_autocorr <- om_custom(
  proc_cv  = 0.30,
  proc_rho = 0.50,
  obs_cv   = 0.10,
  label    = "Autocorrelated observation error"
)

# om_hyperstability <- om_custom(
#   beta_hyper = 0.50,
#   proc_cv    = 0.20,
#   obs_cv     = 0.10,
#   label      = "Hyperstability"
# )



# om_environmental <- om_custom(
#   proc_cv     = 0.30,
#   proc_rho    = 0.40,
#   proc_rho_sp = 0.60,
#   label       = "Correlated environmental variability"
# )
#' List all available built-in OMs
#'
#' Returns a summary table of all preset OMs with their
#' key parameter values. Useful for workshop overview.
#' @return A data.frame summarising built-in OMs
#' @export
# om_list <- function() {
#   oms <- list(
#     om_deterministic,
#     om_baseline,
#     om_param_uncertainty,
#     om_autocorr
#     #om_hyperstability,
#     #om_environmental
#   )
#   data.frame(
#     label       = sapply(oms, `[[`, "label"),
#     proc_cv     = sapply(oms, `[[`, "proc_cv"),
#     proc_rho    = sapply(oms, `[[`, "proc_rho"),
#     #proc_rho_sp = sapply(oms, `[[`, "proc_rho_sp"),
#     obs_cv      = sapply(oms, `[[`, "obs_cv"),
#     #catch_cv    = sapply(oms, `[[`, "catch_cv"),
#     #beta_hyper  = sapply(oms, `[[`, "beta_hyper"),
#     stringsAsFactors = FALSE
#   )
# }

om_list <- list(
  om_deterministic     = om_deterministic,
  om_baseline          = om_baseline,
  om_param_uncertainty = om_param_uncertainty,
  om_autocorr          = om_autocorr,
  om_custom_ui         = list(label = "Custom OM (set parameters below)")
)
