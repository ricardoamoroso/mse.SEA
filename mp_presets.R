# =================================================================
# mp_presets.R — Management Procedure definitions
#
# TOKENS (resolved by resolve_mp_params in app.R):
#   ".E_mmsy"              — multispecies MSY effort
#   ".E_mmsy_08"           — E_mmsy * 0.8
#   ".E_mmsy_06"           — E_mmsy * 0.6
#   ".E_init"              — historical effort from user
#   ".closure_months_prop" — input$closure_months / 12
#   ".prop_closed_csv"     — PropB_area column from CSV
# =================================================================

# -----------------------------------------------------------------
# 1. BASE CONSTRUCTOR
# -----------------------------------------------------------------

mp_custom <- function(mp_id, name, hcr_type, ...) {
  structure(
    list(mp_id    = mp_id,
         name     = name,
         hcr_type = hcr_type,
         params   = list(...)),
    class = "mp_spec"
  )
}

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
# -----------------------------------------------------------------

# No fishing
mp_no_fishing <- mp_custom(
  mp_id    = "MP_no_fishing",
  name     = "No fishing (E = 0)",
  hcr_type = "constant_effort",
  E        = 0
)

# Constant effort at Emsy fractions
mp_const_emsy <- mp_custom(
  mp_id    = "MP_constEmsy",
  name     = "Constant E = Emsy",
  hcr_type = "constant_effort",
  E        = ".E_mmsy"
)

mp_const_08emsy <- mp_custom(
  mp_id    = "MP_const08Emsy",
  name     = "Constant E = 0.8 Emsy",
  hcr_type = "constant_effort",
  E        = ".E_mmsy_08"
)

mp_const_06emsy <- mp_custom(
  mp_id    = "MP_const06Emsy",
  name     = "Constant E = 0.6 Emsy",
  hcr_type = "constant_effort",
  E        = ".E_mmsy_06"
)

# Biomass hockey-stick MPs
mp_B_hs_min <- mp_custom(
  mp_id    = "MP_B_hs_min",
  name     = "Biomass HS (most depleted species)",
  hcr_type = "B_hs_min",
  d_lim    = 0.2,
  d_trig   = 0.4,
  Emax     = ".E_init",
  Emin     = 0.0
)

mp_B_hs_mean <- mp_custom(
  mp_id    = "MP_B_hs_mean",
  name     = "Biomass HS (mean depletion)",
  hcr_type = "B_hs_mean",
  d_lim    = 0.2,
  d_trig   = 0.6,
  Emax     = ".E_init",
  Emin     = 0.0
)

mp_B_hs_low_r <- mp_custom(
  mp_id    = "MP_B_hs_low_r",
  name     = "Biomass HS (n least resilient species)",
  hcr_type = "B_hs_low_r",
  d_lim    = 0.2,
  d_trig   = 0.4,
  Emax     = ".E_init",
  Emin     = 0.0,
  n_low    = 2
)

# Index hockey-stick MPs
mp_index_hs_min <- mp_custom(
  mp_id    = "MP_index_hs_min",
  name     = "Index HS (most depleted species)",
  hcr_type = "index_hs_min",
  d_lim    = 0.2,
  d_trig   = 0.4,
  Emax     = ".E_init",
  Emin     = 0.0
)

mp_index_hs_mean <- mp_custom(
  mp_id    = "MP_index_hs_mean",
  name     = "Index HS (mean species)",
  hcr_type = "index_hs_mean",
  d_lim    = 0.2,
  d_trig   = 0.4,
  Emax     = ".E_init",
  Emin     = 0.0
)

# Seasonal closure
mp_seasonal_closure <- mp_custom(
  mp_id       = "MP_seasonal_closure",
  name        = "Seasonal closure (X months)",
  hcr_type    = "seasonal_closure",
  prop_season = ".closure_months_prop",
  E_base      = ".E_init"
)

# Spatial closure
mp_spatial_closure <- mp_custom(
  mp_id       = "MP_spatial_closure",
  name        = "Spatial closure",
  hcr_type    = "spatial_closure",
  prop_closed = ".prop_closed_csv",
  E_base      = ".E_init"
)

# 2-over-3 rule
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

# Slope rule
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
# 3. AUTO-DISCOVER ALL MPs
# Any mp_custom() object defined above is picked up automatically.
# No manual list maintenance needed.
# -----------------------------------------------------------------
mp_all <- Filter(
  function(x) inherits(x, "mp_spec"),
  mget(ls(), envir = environment())
)
