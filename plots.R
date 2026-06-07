# =================================================================
# plots.R — Standard simulation plots
# =================================================================

#' Plot depletion (B/K) — median + 90% CI ribbon per species
#' @param dt data.table from simulate_multispp() — projection years only
#' @export
plot_B_fun <- function(dt) {
  dt <- copy(dt)
  dt[, depl    := B / K]
  dt[, proj_yr := year - min(year) + 1]   # re-index to 1..nyears_proj
  dt[, spp     := factor(spp)]

  ggplot2::ggplot(dt, ggplot2::aes(proj_yr, depl, colour = spp, fill = spp)) +
    ggplot2::stat_summary(
      fun.min  = function(x) quantile(x, 0.05),
      fun.max  = function(x) quantile(x, 0.95),
      geom     = "ribbon", alpha = 0.15, colour = NA
    ) +
    ggplot2::stat_summary(fun = median, geom = "line", linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0.2, linetype = "dashed",
                        colour = "red",    linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dotted",
                        colour = "grey40", linewidth = 0.5) +
    ggplot2::annotate("text", x = -Inf, y = 0.21,
                      label = "Blim",  hjust = -0.1, size = 3, colour = "red") +
    ggplot2::annotate("text", x = -Inf, y = 0.51,
                      label = "Bmsy",  hjust = -0.1, size = 3, colour = "grey40") +
    ggplot2::scale_colour_viridis_d(option = "turbo") +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(x = "Year", y = "Depletion (B/K)",
                  colour = "Species", fill = "Species",
                  title    = "Depletion trajectories",
                  subtitle = "Median \u00b1 90% CI across simulations") +
    ggplot2::scale_y_continuous(limits = c(0, NA))
}

#' Plot catch trajectories — median + 90% CI ribbon per species
#' @param dt data.table from simulate_multispp() — projection years only
#' @export
plot_C_fun <- function(dt) {
  dt <- copy(dt)
  dt[, proj_yr := year - min(year) + 1]
  dt[, spp     := factor(spp)]

  ggplot2::ggplot(dt, ggplot2::aes(proj_yr, C, colour = spp, fill = spp)) +
    ggplot2::stat_summary(
      fun.min  = function(x) quantile(x, 0.05),
      fun.max  = function(x) quantile(x, 0.95),
      geom     = "ribbon", alpha = 0.15, colour = NA
    ) +
    ggplot2::stat_summary(fun = median, geom = "line", linewidth = 1) +
    ggplot2::scale_colour_viridis_d(option = "turbo") +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(x = "Year", y = "Catch",
                  colour = "Species", fill = "Species",
                  title    = "Catch trajectories",
                  subtitle = "Median \u00b1 90% CI across simulations") +
    ggplot2::scale_y_continuous(limits = c(0, NA))
}

#' Plot effort trajectories — median + 90% CI ribbon
#' @param dt data.table from simulate_multispp() — projection years only
#' @export
plot_E_fun <- function(dt) {
  dt      <- copy(dt)
  dt_E    <- unique(dt[, .(sim, year, E)])
  dt_E[, proj_yr := year - min(year) + 1]

  ggplot2::ggplot(dt_E, ggplot2::aes(proj_yr, E)) +
    ggplot2::stat_summary(
      fun.min  = function(x) quantile(x, 0.05),
      fun.max  = function(x) quantile(x, 0.95),
      geom     = "ribbon", alpha = 0.2, fill = "#2166ac", colour = NA
    ) +
    ggplot2::stat_summary(fun = median, geom = "line",
                          linewidth = 1, colour = "#2166ac") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Year", y = "Effort",
                  title    = "Effort trajectories",
                  subtitle = "Median \u00b1 90% CI across simulations") +
    ggplot2::scale_y_continuous(limits = c(0, NA))
}

#' Plot catch composition — proportions summing to 1 per year
#'
#' Computes median catch per species-year first, then converts to
#' proportions of the total median, so the stacked area always sums to 1.
#' Species are ordered by total median catch (largest at bottom) to avoid
#' the visual gaps that occur with character-ordered stacking.
#'
#' @param dt data.table from simulate_multispp() — projection years only
#' @export
plot_Cprop_fun <- function(dt) {
  df         <- as.data.frame(dt)
  df$proj_yr <- df$year - min(df$year) + 1

  # Median catch per species-year across sims
  med_C <- aggregate(C ~ proj_yr + spp, data = df,
                     FUN = function(x) median(x, na.rm = TRUE))
  names(med_C)[names(med_C) == "C"] <- "med_C"

  # Total per year -> proportion
  tot           <- aggregate(med_C ~ proj_yr, data = med_C, FUN = sum)
  names(tot)[2] <- "C_total"
  med_C         <- merge(med_C, tot, by = "proj_yr")
  med_C$C_prop  <- med_C$med_C / pmax(med_C$C_total, 1e-9)

  # Order species: largest total at bottom of stack
  spp_totals <- aggregate(med_C ~ spp, data = med_C, FUN = sum)
  spp_order  <- spp_totals$spp[order(-spp_totals$med_C)]
  med_C$spp  <- factor(med_C$spp, levels = spp_order)
  med_C      <- med_C[order(med_C$proj_yr, med_C$spp), ]

  # Pre-compute cumulative ymin/ymax per year — geom_ribbon has no gap issues
  ribbon_list <- lapply(split(med_C, med_C$proj_yr), function(yr_df) {
    yr_df      <- yr_df[order(yr_df$spp), ]
    yr_df$ymax <- cumsum(yr_df$C_prop)
    yr_df$ymin <- c(0, head(yr_df$ymax, -1))
    yr_df
  })
  ribbon_df <- do.call(rbind, ribbon_list)

  ggplot2::ggplot(ribbon_df,
                  ggplot2::aes(x = proj_yr, ymin = ymin, ymax = ymax,
                               fill = spp)) +
    ggplot2::geom_ribbon() +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0),
                                breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(x    = "Year",
                  y    = "Proportion of total catch",
                  fill = "Species",
                  title    = "Catch composition",
                  subtitle = "Based on median catch across simulations") +
    ggplot2::scale_y_continuous(limits = c(0, NA))
}

#' Plot the HCR shape
#' @param hcr list(type, params) — from mp_spec via selected_mp_r()
#' @export
plot_HCR_fun <- function(hcr) {
  type  <- hcr$type
  p     <- hcr$params
  d_seq <- seq(0, 1.2, length.out = 300)

  if (type == "constant_effort") {
    E_val <- p$E %||% 0
    df    <- data.frame(d = d_seq, E = E_val)
    ggplot2::ggplot(df, ggplot2::aes(d, E)) +
      ggplot2::geom_line(linewidth = 1.5, colour = "#2166ac") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(x = "Depletion signal (B/K)", y = "Effort") +
      ggplot2::ylim(0, max(E_val * 1.2, 0.1))

  } else if (type %in% c("spatial_closure", "seasonal_closure")) {
    prop   <- if (type == "spatial_closure") (p$prop_closed %||% 0.3)
    else                            (p$prop_season %||% 0.25)
    E_base <- p$E_base %||% 1
    E_eff  <- E_base * (1 - prop)
    df     <- data.frame(d = d_seq, E = E_eff)
    ggplot2::ggplot(df, ggplot2::aes(d, E)) +
      ggplot2::geom_line(linewidth = 1.5, colour = "#d6604d") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(x = "Depletion signal (B/K)", y = "Effective effort",
                    subtitle = sprintf("Fixed effort reduction: %.0f%%", prop * 100))

  } else if (type %in% c("B_hs_min", "B_hs_mean", "index_hs_min", "B_hs_low_r")) {
    dlim <- p$d_lim  %||% 0.2
    dtrg <- p$d_trig %||% 0.4
    Emax <- p$Emax   %||% 1
    Emin <- p$Emin   %||% 0
    E_vals <- sapply(d_seq, hcr_hockeystick,
                     d_lim = dlim, d_trig = dtrg, Emax = Emax, Emin = Emin)
    df <- data.frame(d = d_seq, E = E_vals)
    ggplot2::ggplot(df, ggplot2::aes(d, E)) +
      ggplot2::geom_vline(xintercept = c(dlim, dtrg),
                          linetype = "dashed", colour = "grey60") +
      ggplot2::geom_line(linewidth = 1.5, colour = "#2166ac") +
      ggplot2::annotate("text", x = dlim, y = Emax * 0.05,
                        label = "Blim",  hjust = -0.1, size = 3.5) +
      ggplot2::annotate("text", x = dtrg, y = Emax * 0.05,
                        label = "Btrig", hjust = -0.1, size = 3.5) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(x = "Depletion signal (B/K)", y = "Effort") +
      ggplot2::scale_y_continuous(limits = c(0, NA))

  } else {
    ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::labs(
        title    = paste(type, "\u2014 dynamic rule"),
        subtitle = "Effort adjusts each year; no static curve to display."
      )
  }
}
