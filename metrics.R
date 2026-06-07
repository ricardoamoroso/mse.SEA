# =================================================================
# metrics.R — Performance metrics for multispecies MSE
# =================================================================

#' Compute performance metrics over last 10 projection years
#'
#' @param dt A data.table returned by simulate_multispp()
#'   (projection years only, year re-indexed to 1..nyears_proj)
#' @return A list with two elements:
#'   \describe{
#'     \item{by_sim}{Per-simulation metrics (one row per sim)}
#'     \item{summary}{Medians and probabilities across all sims}
#'   }
#' @export
metrics_multispp_10y <- function(dt) {
  cfg <- attr(dt, "config")
  ny  <- cfg$nyears_proj %||% cfg$nyears %||% max(dt$year)

  dt[, depl       := B / K]
  dt[, Blim       := 0.2 * K]
  dt[, Bmsy       := 0.5 * K]
  dt[, below_Blim := as.integer(B < Blim)]
  dt[, below_Bmsy := as.integer(B < Bmsy)]

  dt10 <- dt[year > ny - 10]

  # ---- Catch metrics (total across species per sim-year) ----
  tot10    <- dt10[, .(C_total = sum(C)), by = .(sim, year)]
  per_sim  <- tot10[, .(
    sum_catch_10yr  = sum(C_total),
    mean_catch_10yr = mean(C_total),
    AAV = mean(
      abs(C_total - data.table::shift(C_total)) /
        pmax(1e-9, data.table::shift(C_total)),
      na.rm = TRUE
    )
  ), by = sim]

  # ---- Risk metrics ----
  # For each sim-year: fraction of species below reference point
  # Then average that fraction across years (per sim)
  risk_yr <- dt10[, .(
    frac_below_Blim = mean(below_Blim),
    frac_below_Bmsy = mean(below_Bmsy)
  ), by = .(sim, year)]

  risk10 <- risk_yr[, .(
    # Average depletion risk across all species-years
    prop_below_Blim = mean(frac_below_Blim),
    prop_below_Bmsy = mean(frac_below_Bmsy),
    # Proportion of YEARS in which >= 50% of species are below reference
    prop_years_half_spp_Blim = mean(frac_below_Blim >= 0.5),
    prop_years_half_spp_Bmsy = mean(frac_below_Bmsy >= 0.5)
  ), by = sim]

  # ---- Mean depletion ----
  depl10 <- dt10[, .(mean_depl_10yr = mean(depl)), by = sim]

  # ---- 2/3 rule ----
  # Proportion of years where AT MOST 1/3 of species are below Blim
  # (i.e. at least 2/3 are above Blim)
  rule_2_3 <- risk_yr[, .(
    prop_years_2_3_above_Blim = mean(frac_below_Blim <= 1/3)
  ), by = sim]
  rule_2_3[, meets_2_3_rule := prop_years_2_3_above_Blim >= 0.75]

  # ---- Merge all per-sim metrics ----
  out <- Reduce(
    function(x, y) merge(x, y, by = "sim"),
    list(per_sim, risk10, depl10, rule_2_3)
  )

  # ---- Summary across sims ----
  med <- out[, .(
    med_sum_catch_10yr       = median(sum_catch_10yr),
    med_mean_catch_10yr      = median(mean_catch_10yr),
    med_AAV                  = median(AAV),
    prob_AAV_le_20pct        = mean(AAV <= 0.20),
    med_prop_below_Blim      = median(prop_below_Blim),
    med_prop_below_Bmsy      = median(prop_below_Bmsy),
    prob_half_spp_below_Blim = mean(prop_years_half_spp_Blim >= 0.5),
    prob_half_spp_below_Bmsy = mean(prop_years_half_spp_Bmsy >= 0.5),
    med_mean_depl_10yr       = median(mean_depl_10yr),
    prob_meets_2_3_rule      = mean(meets_2_3_rule)
  )]

  list(by_sim = out, summary = med)
}
