# =================================================================
# mp_presets.R — Management Procedure definitions
#
# Structure:
#   1. mp_custom()  — base constructor
#   2. Reactive MPs — depend on E_mmsy computed from user data
#                     these are defined in the Shiny module, not here
#   3. Fixed MPs    — do not depend on user data, defined here
#
# Note: constant effort MPs (Emsy, 0.8*Emsy, 0.6*Emsy) cannot be
# defined here because E_mmsy depends on the uploaded species
# parameters. They are created at runtime in the Shiny server
# using mmsy_shared_effort().
# =================================================================


# -----------------------------------------------------------------
# 1. BASE CONSTRUCTOR
# -----------------------------------------------------------------

#' Create a custom MP specification
#'
#' @param mp_id    Short identifier string
#' @param name     Display name for plots and tables
#' @param hcr_type HCR type — must match a type in hcr_registry
#' @param ...      HCR parameters (d_lim, d_trig, Emax, Emin, E)
#' @return An object of class mp_spec
#' @export
mp_custom <- function(mp_id, name, hcr_type, ...) {
  structure(
    list(mp_id    = mp_id,
         name     = name,
         hcr_type = hcr_type,
         params   = list(...)),
    class = "mp_spec"
  )
}

#' Print method for mp_spec
#' @export
print.mp_spec <- function(x, ...) {
  cat("MP:", x$name, "\n")
  cat("  id:       ", x$mp_id,    "\n")
  cat("  hcr_type: ", x$hcr_type, "\n")
  if (length(x$params) > 0) {
    cat("  params:\n")
    for (nm in names(x$params))
      cat("    ", nm, ":", x$params[[nm]], "\n")
  }
  invisible(x)
}


# -----------------------------------------------------------------
# 2. FIXED MPs
# These do not depend on user data so can be defined here.
# Biomass and index based HCRs with standard reference points.
# -----------------------------------------------------------------

mp_no_fishing <- mp_custom(
  mp_id    = "MP_no_fishing",
  name     = "No fishing (E = 0)",
  hcr_type = "constant_effort",
  E        = 0
)


# ── Biomass hockey-stick MPs ─────────────────────────────────────


mp_B_hs_mean <- mp_custom(
  mp_id    = "MP_B_hs_mean",
  name     = "Biomass HS (mean depletion)",
  hcr_type = "B_hs_mean",
  d_lim    = 0.2,
  d_trig   = 0.6,
  Emax     = ".E_init",
  Emin     = 0.0
)

mp_index_hs_min <- mp_custom(
  mp_id    = "MP_index_hs_min",
  name     = "Index HS (most depleted species)",
  hcr_type = "index_hs_min",
  d_lim    = 0.2,
  d_trig   = 0.4,
  Emax     = ".E_init",
  Emin     = 0.0
)

#
# ── Seasonal closure ─────────────────────────────────────────────
mp_seasonal_closure <- mp_custom(
  mp_id       = "MP_seasonal_closure",
  name        = "Seasonal closure (X months)",
  hcr_type    = "seasonal_closure",
  prop_season = ".closure_months_prop",
  E_base      = ".E_init"
)
#
# ── Spatial closure ──────────────────────────────────────────────
mp_spatial_closure <- mp_custom(
  mp_id       = "MP_spatial_closure",
  name        = "Spatial closure",
  hcr_type    = "spatial_closure",
  prop_closed = ".prop_closed_csv",
  E_base      = ".E_init"
)

#
# ── 2-over-3 rule ────────────────────────────────────────────────
mp_2over3_min <- mp_custom(
  mp_id      = "MP_2over3_min",
  name       = "2-over-3 rule (min species)",
  hcr_type   = "two_over_three",
  agg_type   = "min",
  max_change = 0.20
)

mp_2over3_mean <- mp_custom(
  mp_id      = "MP_2over3_mean",
  name       = "2-over-3 rule (mean species)",
  hcr_type   = "two_over_three",
  agg_type   = "mean",
  max_change = 0.20
)

#
# # Slope rule responding to most depleted species
mp_slope_min <- mp_custom(
  mp_id      = "MP_slope_min",
  name       = "Slope rule (min species)",
  hcr_type   = "slope_rule",
  agg_type   = "min",
  n_years    = 5,
  lambda     = 1,
  max_change = 0.20
)

mp_slope_mean <- mp_custom(
  mp_id      = "MP_slope_mean",
  name       = "Slope rule (mean species)",
  hcr_type   = "slope_rule",
  agg_type   = "mean",
  n_years    = 5,
  lambda     = 1,
  max_change = 0.20
)

# -----------------------------------------------------------------
# 3. RUNTIME MPs — created in Shiny server, documented here
# These require E_mmsy computed from mmsy_shared_effort()
#
mp_const_emsy <- mp_custom(
   mp_id    = "MP_constEmsy",
   name     = "Constant E = Emsy",
   hcr_type = "constant_effort",
   E        = ".E_mmsy"
 )


# mp_const_08emsy <- mp_custom("MP_const08Emsy", "Constant E = 0.8 Emsy",
#                               "constant_effort", E = mmsy$E_mmsy * 0.8)
#
# mp_const_06emsy <- mp_custom("MP_const06Emsy", "Constant E = 0.6 Emsy",
#                               "constant_effort", E = mmsy$E_mmsy * 0.6)
# -----------------------------------------------------------------


# Bottom of mp_presets.R
mp_all <- Filter(
  function(x) inherits(x, "mp_spec"),
  mget(ls(), envir = environment())
)
