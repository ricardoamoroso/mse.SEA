# =================================================================
# app.R — Multispecies MSE Shiny Application
# Sources: dynamics.R, hcr.R, metrics.R, mmsy.R,
#          mp_presets.R, om_presets.R, plots.R, utils.R
# =================================================================

# ---- Libraries ----
# These are loaded automatically when the package is attached.
# Listed here so they are attached when running the app standalone.
library(shiny)
library(shinyjs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggpubr)
library(zip)
library(scales)
library(DT)

# ---- Source helper files ----
.app_dir <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile)),
  error = function(e) getwd()
)
for (.f in c("utils.R", "mmsy.R", "hcr.R",
             "mp_presets.R", "om_presets.R",
             "dynamics.R", "metrics.R", "plots.R",
             "conditioning.R")) {
  source(file.path(.app_dir, .f), local = FALSE)
}
rm(.app_dir, .f)
cat("conditioning_diagnostics exists:", exists("conditioning_diagnostics"), "\n")

# Wrap metrics (nyears_proj + nyears_cond)
# and restrict to projection phase only for the app
# Filter simulation output to projection phase only (for plots and metrics table)
# dynamics.R tags every row with phase = "cond" (conditioning years) or phase = "proj" (projection years).
proj_only <- function(dt) {
  if ("phase" %in% names(dt)) dt[phase == "proj"] else dt
}

# connect the new dynamics.R with the old code metrics.R so I do not have to change the way we calculated the metrics before
# Config format compatibility — dynamics.R stores nyears_proj + nyears_cond separately, but metrics.R was written expecting a single nyears. The wrapper resolves whichever format is present.
metrics_multispp_10y_safe <- function(dt) {
  cfg  <- attr(dt, "config")
  # Support both old (nyears) and new (nyears_proj + nyears_cond) config formats
  ny   <- cfg$nyears_proj %||% cfg$nyears %||% max(dt$year)
  nc   <- cfg$nyears_cond %||% 0

  # Work only on projection years (phase == "proj" if column exists, else all)
  if ("phase" %in% names(dt)) {
    dt_proj <- dt[phase == "proj"]
    # Re-index year within projection (1 = first proj year)
    dt_proj[, proj_year := year - nc]
  } else {
    dt_proj <- copy(dt)
    dt_proj[, proj_year := year]
  }

  # Temporarily rename proj_year so metrics.R code works unchanged
  dt_work <- copy(dt_proj)
  dt_work[, year := proj_year]
  data.table::setattr(dt_work, "config",
                      modifyList(as.list(cfg), list(nyears = ny, nyears_proj = ny, nyears_cond = 0)))

  metrics_multispp_10y(dt_work)
}

# ==============================================================
# OM and MP lookup tables — built after helper files are sourced
# ==============================================================
# here is we update the table of the OMs it should update this part.
om_spec_map     <- om_list   # om_list defined in om_presets.R
om_choices_list <- setNames(names(om_list),
                            sapply(om_list, `[[`, "label"))

# Fixed MP objects from mp_presets.R


fixed_mp_list <- c(mp_all, list(
  mp_custom(
    mp_id    = "MP_custom",
    name     = "Custom MP (set parameters below)",
    hcr_type = "constant_effort",
    E        = 1
  )
))

fixed_mp_choices <- setNames(
  sapply(fixed_mp_list, `[[`, "mp_id"),
  sapply(fixed_mp_list, `[[`, "name")
)


# ==============================================================
# Helper: descriptions for sidebar
# ==============================================================
om_description <- HTML(
  "<b>Operating Models (OMs)</b> define how uncertainty is represented:<br>
  <b>Deterministic:</b> no noise at all — pure model behaviour.<br>
  <b>Baseline:</b> year-to-year process error (CV = 0.20).<br>
  <b>Param + process:</b> uncertainty in r and K + process error.<br>
  <b>Autocorr obs:</b> autocorrelated + noisy survey indices."
)

mp_description <- HTML(
  "<b>Management Procedures (MPs):</b><br>
  <b>Const E = Emsy/0.8/0.6 Emsy:</b> Fixed effort at multispecies MSY or a fraction of it.<br>
  <b>Biomass HS (min/mean):</b> Hockey-stick on most-depleted or mean B/K.<br>
  <b>Index HS (min):</b> Hockey-stick on CPUE index of most-depleted species.<br>
  <b>Spatial closure:</b> Part of the fishing ground closed.<br>
  <b>Seasonal closure:</b> X months closed per year (3 months = 25% effort reduction).<br>
  <b>2-over-3 (min/mean):</b> Adjust effort from recent vs past CPUE ratio.<br>
  <b>Slope rule (min/mean):</b> Adjust effort based on CPUE trend slope."
)

# ==============================================================
# CSV READER: handles r/r_low/r_high, K/K_low/K_high, Dep_ini/range
# Reads the uploaded CSV and normalises all parameter columns into a consistent internal format regardless of what columns the user provided.
# For each of r, K, and Dep_ini it handles three cases
# ==============================================================
read_param_csv <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  nms <- names(df)

  # --- r ---
  # For r:
  # r_low + r_high present → range, compute r_mean as midpoint (or read it if present), set r_is_range = TRUE
  # Only r present → point value, r_lo = r_hi = r_mean, r_is_range = FALSE
  if (all(c("r_low", "r_hi") %in% nms)) {
    df$r_lo   <- as.numeric(df$r_low)
    df$r_hi   <- as.numeric(df$r_hi)
    df$r_mean <- if ("r_mean" %in% nms) as.numeric(df$r_mean) else (df$r_lo + df$r_hi) / 2
    df$r_is_range <- TRUE
  } else if ("r" %in% nms) {
    df$r_mean <- as.numeric(df$r)       # ← read "r", not "r_mean"
    df$r_lo   <- df$r_mean
    df$r_hi   <- df$r_mean
    df$r_is_range <- FALSE
  } else {
    stop("CSV must have column 'r_mean' OR columns 'r_low' and 'r_hi'.")
  }

  # --- K --- For K: same logic, also accepts K_mean alone without a range
  if (all(c("K_low", "K_hi") %in% nms)) {
    df$K_lo   <- as.numeric(df$K_low)
    df$K_hi   <- as.numeric(df$K_hi)
    df$K_mean <- if ("K_mean" %in% nms) as.numeric(df$K_mean) else (df$K_lo + df$K_hi) / 2
    df$K_is_range <- TRUE
  } else if ("K_mean" %in% nms) {      # ← accept K_mean alone (no range)
    df$K_mean <- as.numeric(df$K_mean)
    df$K_lo   <- df$K_mean
    df$K_hi   <- df$K_mean
    df$K_is_range <- FALSE
  } else if ("K" %in% nms) {           # ← also accept plain "K"
    df$K_mean <- as.numeric(df$K)
    df$K_lo   <- df$K_mean
    df$K_hi   <- df$K_mean
    df$K_is_range <- FALSE
  } else {
    stop("CSV must have column 'K_mean' OR columns 'K_low' and 'K_hi'.")
  }

  # --- Dep_ini ---
  # dep_explicit_bounds = TRUE means the user gave Dep_ini_low/Dep_ini_high
  # explicitly, so assess_cv should be 0 (no extra noise on depletion signal).
  # dep_explicit_bounds = FALSE means only a point Dep_ini was given, so the
  # app adds ±0.1 range AND assess_cv = 0.1 (dynamics.R default) for HCR noise.

  # --- Dep_ini ---
  if (all(c("Dep_low", "Dep_hi") %in% nms)) {
    df$dep_lo  <- as.numeric(df$Dep_low)
    df$dep_hi  <- as.numeric(df$Dep_hi)
    df$dep_mean <- if ("Dep_ini" %in% nms) as.numeric(df$Dep_ini) else (df$dep_lo + df$dep_hi) / 2
    df$dep_is_range        <- TRUE
    df$dep_explicit_bounds <- TRUE
  } else if ("Dep_ini" %in% nms) {
    dep         <- as.numeric(df$Dep_ini)
    df$dep_mean <- dep
    df$dep_lo   <- pmax(dep - 0.1, 0.05)
    df$dep_hi   <- pmin(dep + 0.1, 0.99)
    df$dep_is_range        <- TRUE
    df$dep_explicit_bounds <- FALSE
  } else {
    stop("CSV must have column 'Dep_ini' OR columns 'Dep_low' and 'Dep_hi'.")
  }

  # --- q (required) ---
  if (!"q" %in% nms) stop("CSV must have column 'q'.")
  df$q <- as.numeric(df$q)

  # --- Species name ---
  if ("Species" %in% nms) {
    df$species_name <- as.character(df$Species)
  } else {
    df$species_name <- paste0("sp", seq_len(nrow(df)))
  }

  # --- Group (optional) ---
  if ("Group" %in% nms) df$Group <- as.character(df$Group)

  # --- PropB_area (optional) ---
  # Vector of species-specific proportions of biomass affected by spatial closure.
  # e.g. c(0.1, 1.0, 0.3) means closing x% of area affects 10%, 100%, 30%
  # of each species' biomass respectively.
  # If not in CSV, defaults to NULL (hcr.R will use scalar 0.3 for all species).
  if ("PropB_area" %in% nms) {
    df$prop_B_area <- as.numeric(df$PropB_area)
  } else {
    df$prop_B_area <- NA  # signals: use scalar default
  }

  df
}

# ==============================================================
# UI
# ==============================================================
ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML("
      .sidebar-section { margin-bottom: 12px; }
      .param-note { font-size: 11px; color: #666; }
      .run-info { background: #f7f7f7; border-radius: 6px;
                  padding: 8px; margin-top: 6px; font-size: 12px; }
    "))
  ),
  titlePanel("Multispecies MSE — Management Strategy Evaluation"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # --- CSV upload ---
      div(class = "sidebar-section",
          fileInput("param_file", "Upload species parameters (CSV)",
                    accept = ".csv"),
          p(class = "param-note",
            "Required columns: Species, q, and either r or r_low/r_high,
           either K or K_low/K_high, and either Dep_ini or Dep_ini_low/Dep_ini_high.")
      ),
      hr(),

      # --- E_init ---
      div(class = "sidebar-section",
          h5("Historical effort (E_init)"),
          numericInput("E_init_input",
                       label    = NULL,
                       value    = NA,
                       min      = 0,
                       step     = 0.0001),
          p(class = "param-note",
            "Leave blank to use E_mmsy as starting effort.
     See the Initialisation tab for guidance on choosing this value."),
          uiOutput("einit_hint")
      ),
      hr(),

      # --- OM selector ---
      div(class = "sidebar-section",
          selectInput("om_id", "Operating Model (OM)",
                      choices = om_choices_list,
                      selected = "om_baseline"),
          om_description
      ),
      hr(),

      # --- Custom OM sliders ---
      conditionalPanel(
        condition = "input.om_id == 'om_custom_ui'",
        sliderInput("cv_r",     "CV on r (growth rate)",      0, 0.5, 0,    step = 0.05),
        sliderInput("cv_K",     "CV on K (carrying capacity)", 0, 0.5, 0,   step = 0.05),
        sliderInput("proc_cv",  "Process error CV",            0, 0.5, 0.2, step = 0.05),
        sliderInput("obs_cv",   "Observation error CV",        0, 0.5, 0.2, step = 0.05),
        sliderInput("proc_rho", "Autocorrelation (process)",   0, 0.9, 0,   step = 0.1),
        sliderInput("assess_cv", "Assessment error CV (depletion signal for HS rules)",
                    0, 0.5, 0.1, step = 0.05)
      ),
      hr(),

      # --- MP selector ---
      div(class = "sidebar-section",
          selectInput("mp_id", "Management Procedure (MP)",
                      choices = fixed_mp_choices,
                      selected = "MP_constEmsy"),
          mp_description,
          # Seasonal closure — months slider
          conditionalPanel(
            condition = "input.mp_id == 'MP_seasonal_closure'",
            br(),
            sliderInput("closure_months",
                        label = "Months closed per year",
                        min = 1, max = 7, value = 3, step = 1,
                        post = " month(s)"),
            uiOutput("closure_prop_label")
          ),

          # Spatial closure: prop_closed is read from PropB_area column in CSV
      ),
      hr(),

      conditionalPanel(
        condition = "input.mp_id == 'MP_custom'",
        br(),
        selectInput("custom_mp_type", "MP type",
                    choices = c(
                      "Constant effort"                  = "constant_effort",
                      "Hockey-stick (min species)"       = "B_hs_min",
                      "Hockey-stick (mean species)"      = "B_hs_mean",
                      "Hockey-stick (n least resilient)" = "B_hs_low_r",
                      "Index HS (min species)"           = "index_hs_min",
                      "Index HS (mean species)"          = "index_hs_mean",
                      "2-over-3"                         = "two_over_three",
                      "Slope rule"                       = "slope_rule"
                    )),

        # Constant effort
        conditionalPanel(
          condition = "input.custom_mp_type == 'constant_effort'",
          sliderInput("custom_E_fraction", "Effort fraction of Emsy",
                      0, 1.5, 1, step = 0.1)
        ),

        # All hockey-stick types share d_lim, d_trig, Emax
        conditionalPanel(
          condition = "input.custom_mp_type == 'B_hs_min' ||
                       input.custom_mp_type == 'B_hs_mean' ||
                       input.custom_mp_type == 'B_hs_low_r' ||
                       input.custom_mp_type == 'index_hs_min' ||
                       input.custom_mp_type == 'index_hs_mean'",
          sliderInput("custom_d_lim",     "Limit depletion (d_lim)",       0,   0.5, 0.2, step = 0.05),
          sliderInput("custom_d_trig",    "Trigger depletion (d_trig)",    0,   1.0, 0.4, step = 0.05),
          sliderInput("custom_Emax_frac", "Max effort (fraction of Emsy)", 0,   1.5, 1.0, step = 0.1)
        ),

        # n_low only for B_hs_low_r
        conditionalPanel(
          condition = "input.custom_mp_type == 'B_hs_low_r'",
          sliderInput("custom_n_low", "Number of most vulnerable species",
                      1, 5, 2, step = 1)
        ),

        # 2-over-3
        conditionalPanel(
          condition = "input.custom_mp_type == 'two_over_three'",
          selectInput("custom_agg_type_2o3", "Apply rule to",
                      choices = c("Most depleted species" = "min",
                                  "Mean across species"   = "mean")),
          sliderInput("custom_max_change", "Max annual change",
                      0, 0.5, 0.2, step = 0.05)
        ),

        # Slope rule
        conditionalPanel(
          condition = "input.custom_mp_type == 'slope_rule'",
          selectInput("custom_agg_type_slope", "Apply rule to",
                      choices = c("Most depleted species" = "min",
                                  "Mean across species"   = "mean")),
          sliderInput("custom_max_change_slope", "Max annual change",
                      0, 0.5, 0.2, step = 0.05),
          sliderInput("custom_n_years",  "Years for slope calculation", 2, 10, 5,   step = 1),
          sliderInput("custom_lambda",   "Sensitivity (lambda)",        0.1, 3, 1,  step = 0.1)
        )
      ),
      hr(),

      # --- Simulation settings ---
      div(class = "sidebar-section",
          h5("Simulation settings"),
          numericInput("nyears_proj", "Projection years", value = 60, min = 10, max = 300),
          numericInput("nyears_cond", "Conditioning years", value = 0, min = 0, max = 50),
          numericInput("nsims", "Number of simulations", value = 100, min = 1, max = 2000),
          checkboxInput("shared_fleet", "Shared fleet (single effort for all species)", TRUE)
      ),
      hr(),

      #Add the run name
      textInput("run_label", "Run label (optional)",
                placeholder = "e.g. Seasonal 3 months, low effort"),

      # --- Run & clear ---
      actionButton("run_sim", "Run simulation", class = "btn-primary btn-block"),
      br(), br(),
      actionButton("clear_runs", "Clear trade-off runs", class = "btn-warning btn-block"),

      # --- Run info box (shown after run) ---
      uiOutput("run_info_box")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",

        # ---- Tab 0: Species table ----
        tab_initialisation <- tabPanel(
          "Initialisation",
          br(),

          # --- Top: three effort reference values ---
          h4("Step 1 — Review reference effort levels"),
          p(style = "color:#555; font-size:13px;",
            "Upload your species CSV and set E_init in the sidebar.
     The table below shows three reference effort levels to help
     you decide what historical effort to use."),
          DTOutput("effort_ref_table"),
          br(),

          # --- Middle: consistency diagnostics ---
          h4("Step 2 — Check consistency with stated depletions"),
          p(style = "color:#555; font-size:13px;",
            "For each effort level the table shows the implied equilibrium
     depletion. Compare to your stated Dep_ini. Large residuals
     (> 0.10) suggest the effort estimate and depletion estimate
     are inconsistent with each other."),
          uiOutput("conditioning_message_ui"),
          br(),
          DTOutput("conditioning_table"),
          br(),

          # --- Bottom: confirm ---
          h4("Step 3 — Confirm and proceed"),
          p(style = "color:#555; font-size:13px;",
            "Once you are satisfied with E_init, select an OM and MP
     in the sidebar and press Run."),
          uiOutput("conditioning_status_ui")
        ),

        # ---- Tab 1: Species table ----


        tabPanel("Species data",
                 br(),
                 uiOutput("species_note_ui"),
                 br(),
                 h4("Species parameters (from CSV)"),
                 DTOutput("species_table"),
                 br(),
                 uiOutput("mmsy_header_ui"),
                 DTOutput("mmsy_table")
        ),

        # ---- Tab 1: HCR ----
        tabPanel("Harvest control rule",
                 br(),
                 downloadButton("download_HCR", "Download plot"),
                 br(), br(),
                 plotOutput("plot_HCR", height = "420px")
        ),

        # ---- Tab 2: Simulation results ----
        tabPanel("Simulation results",
                 br(),
                 fluidRow(
                   column(12,
                          downloadButton("download_B",       "Depletion"),
                          downloadButton("download_C",       "Catch"),
                          downloadButton("download_E",       "Effort"),
                          downloadButton("download_Cprop",   "Catch composition"),
                          downloadButton("download_sim_zip", "All plots (ZIP)")
                   )
                 ),
                 br(),
                 fluidRow(
                   column(6, plotOutput("plot_B",     height = "300px")),
                   column(6, plotOutput("plot_C",     height = "300px"))
                 ),
                 fluidRow(
                   column(6, plotOutput("plot_E",     height = "300px")),
                   column(6, plotOutput("plot_Cprop", height = "300px"))
                 )
        ),

        # ---- Tab 3: Metrics ----
        tabPanel("Metrics",
                 br(),
                 downloadButton("download_metrics_sim",  "Download by-sim table"),
                 downloadButton("download_metrics_summ", "Download summary table"),
                 br(), br(),
                 h4("Summary (medians across simulations)"),
                 DTOutput("metrics_table_summary"),
                 br(),
                 h4("Per-simulation metrics"),
                 DTOutput("metrics_table_sim")
        ),

        # ---- Tab 4: Trade-off plots ----
        tabPanel("Trade-off plots",
                 br(),
                 fluidRow(
                   column(12,
                          downloadButton("download_tradeoff1",   "Catch vs Blim"),
                          downloadButton("download_tradeoff2",   "Catch vs AAV"),
                          downloadButton("download_tradeoff3",   "Catch vs Bmsy"),
                          downloadButton("download_tradeoff4",   "Catch vs 2/3 rule"),
                          downloadButton("download_tradeoff_zip","All trade-offs (ZIP)")
                   )
                 ),
                 br(),
                 h4("Average of last 10 projection years (each point = one run)", style = "font-size:14px;"),
                 fluidRow(
                   column(6, plotOutput("tradeoffPlot",      height = "300px")),
                   column(6, plotOutput("tradeoffPlot_AAV",  height = "300px"))
                 ),
                 fluidRow(
                   column(6, plotOutput("tradeoffPlot_Bmsy", height = "300px")),
                   column(6, plotOutput("tradeoffPlot_2_3",  height = "300px"))
                 )
        )
      )
    )
  )
)

# ==============================================================
# SERVER
# ==============================================================
server <- function(input, output, session) {

  # ---- 1. Load CSV ----
  param_data <- reactive({
    req(input$param_file)
    tryCatch(
      read_param_csv(input$param_file$datapath),
      error = function(e) {
        showNotification(paste("CSV error:", e$message), type = "error", duration = 8)
        NULL
      }
    )
  })

  # ---- 2. Species names ----
  species_names_r <- reactive({
    df <- param_data(); req(df)
    df$species_name
  })

  # ---- 3. MMSY (computed from CSV means) ---- using mmsy_shared_effort() from mmsy.R. Returns E_mmsy, Bmsy per species, etc.
  mmsy_r <- reactive({
    df <- param_data(); req(df)
    mmsy_shared_effort(df$r_mean, df$K_mean, df$q)
  })

  # ---- 3b. E_init resolution ----
  # Priority: user input → E_implied fallback → E_mmsy fallback
  E_init_r <- reactive({
    mmsy <- mmsy_r()
    df   <- param_data()

    # Fallback chain
    E_mmsy_val   <- if (!is.null(mmsy)) mmsy$E_mmsy else NA
    E_implied_val <- if (!is.null(df)) {
      tryCatch(
        estimate_shared_E_from_depletion(df$r_mean, df$q, df$dep_mean),
        error = function(e) NA
      )
    } else NA

    user_val <- input$E_init_input
    if (!is.null(user_val) && !is.na(user_val) && user_val > 0) {
      user_val
    } else if (!is.na(E_implied_val)) {
      E_implied_val
    } else {
      E_mmsy_val
    }
  })

  # Hint shown below the E_init input in the sidebar
  output$einit_hint <- renderUI({
    mmsy <- mmsy_r()
    df   <- param_data()
    if (is.null(mmsy) || is.null(df)) return(NULL)

    E_mmsy_val <- mmsy$E_mmsy
    E_impl_val <- tryCatch(
      estimate_shared_E_from_depletion(df$r_mean, df$q, df$dep_mean),
      error = function(e) NA
    )

    tags$p(style = "color:#555; font-size:11px; margin-top:4px;",
           sprintf("E_mmsy = %.5f", E_mmsy_val),
           br(),
           if (!is.na(E_impl_val))
             sprintf("E_implied = %.5f", E_impl_val)
           else
             "E_implied: unavailable"
    )
  })


  # ---- 4
  # Generic runtime token resolver.
  # Replaces string tokens in mp$params with actual runtime values.
  # To add a new token: define it in the preset as ".token_name"
  # and add it to the token_map here. Nothing else needs to change.
  resolve_mp_params <- function(mp, runtime) {
    token_map <- list(
      .E_mmsy              = runtime$E_mmsy,
      .E_mmsy_08           = runtime$E_mmsy * 0.8,
      .E_mmsy_06           = runtime$E_mmsy * 0.6,
      .E_init              = runtime$E_init,
      .closure_months_prop = runtime$closure_months / 12,
      .prop_closed_csv     = runtime$prop_closed_csv
    )
    mp$params <- lapply(mp$params, function(v) {
      if (is.character(v) && v %in% names(token_map)) token_map[[v]] else v
    })
    mp
  }


  # ---- 5. Get selected MP spec ----
  # Injects E_mmsy into params that need it (E_base, Emax).
  # Also reads the closure sliders and converts to proportions.
  #  gets the chosen MP and injects E_mmsy into parameters that need it:
  # Closure MPs: sets E_base = E_mmsy and reads input$closure_months to set prop_season or prop_closed
  # Hockey-stick MPs: sets Emax = E_mmsy
  # This is necessary because mp_presets.R stores NULL for these fields as a placeholder

  selected_mp_r <- reactive({
    req(input$mp_id)
    id   <- input$mp_id
    mmsy <- mmsy_r(); req(mmsy)
    df   <- param_data()


    # Handle custom MP built from UI sliders
    if (id == "MP_custom") {
      mp <- switch(input$custom_mp_type,
                   constant_effort = mp_custom(
                     mp_id    = "MP_custom",
                     name     = "Custom constant effort",
                     hcr_type = "constant_effort",
                     E        = mmsy$E_mmsy * (input$custom_E_fraction %||% 1)
                   ),
                   B_hs_min = mp_custom(
                     mp_id    = "MP_custom",
                     name     = "Custom HS (min species)",
                     hcr_type = "B_hs_min",
                     d_lim    = input$custom_d_lim   %||% 0.2,
                     d_trig   = input$custom_d_trig  %||% 0.4,
                     Emax     = mmsy$E_mmsy * (input$custom_Emax_frac %||% 1),
                     Emin     = 0
                   ),
                   B_hs_mean = mp_custom(
                     mp_id    = "MP_custom",
                     name     = "Custom HS (mean species)",
                     hcr_type = "B_hs_mean",
                     d_lim    = input$custom_d_lim   %||% 0.2,
                     d_trig   = input$custom_d_trig  %||% 0.4,
                     Emax     = mmsy$E_mmsy * (input$custom_Emax_frac %||% 1),
                     Emin     = 0
                   ),
                   B_hs_low_r = mp_custom(
                     mp_id    = "MP_custom",
                     name     = "Custom HS (n least resilient)",
                     hcr_type = "B_hs_low_r",
                     d_lim    = input$custom_d_lim   %||% 0.2,
                     d_trig   = input$custom_d_trig  %||% 0.4,
                     Emax     = mmsy$E_mmsy * (input$custom_Emax_frac %||% 1),
                     Emin     = 0,
                     n_low    = input$custom_n_low   %||% 2
                   ),
                   index_hs_min = mp_custom(
                     mp_id    = "MP_custom",
                     name     = "Custom index HS (min species)",
                     hcr_type = "index_hs_min",
                     d_lim    = input$custom_d_lim   %||% 0.2,
                     d_trig   = input$custom_d_trig  %||% 0.4,
                     Emax     = mmsy$E_mmsy * (input$custom_Emax_frac %||% 1),
                     Emin     = 0
                   ),
                   index_hs_mean = mp_custom(
                     mp_id    = "MP_custom",
                     name     = "Custom index HS (mean species)",
                     hcr_type = "index_hs_mean",
                     d_lim    = input$custom_d_lim   %||% 0.2,
                     d_trig   = input$custom_d_trig  %||% 0.4,
                     Emax     = mmsy$E_mmsy * (input$custom_Emax_frac %||% 1),
                     Emin     = 0
                   ),
                   two_over_three = mp_custom(
                     mp_id      = "MP_custom",
                     name       = "Custom 2-over-3",
                     hcr_type   = "two_over_three",
                     agg_type   = input$custom_agg_type_2o3   %||% "min",
                     max_change = input$custom_max_change %||% 0.2
                   ),
                   slope_rule = mp_custom(
                     mp_id      = "MP_custom",
                     name       = "Custom slope rule",
                     hcr_type   = "slope_rule",
                     agg_type   = input$custom_agg_type_slope  %||% "min",
                     n_years    = input$custom_n_years         %||% 5,
                     lambda     = input$custom_lambda          %||% 1,
                     max_change = input$custom_max_change_slope %||% 0.2
                   )
      )
      return(mp)
    }


    runtime <- list(
      E_mmsy          = mmsy$E_mmsy,
      E_init          = E_init_r(),
      closure_months  = input$closure_months %||% 3,
      prop_closed_csv = if (!is.null(df) && !any(is.na(df$prop_B_area)))
        df$prop_B_area else NULL
    )

    mp_found <- Filter(function(m) m$mp_id == id, fixed_mp_list)
    if (length(mp_found) == 0) stop("Unknown MP id: ", id)
    resolve_mp_params(mp_found[[1]], runtime)
  })


  # ---- 6. Get selected OM spec ----
  # simply looks up the chosen OM id in om_spec_map and returns the om_spec object.
  selected_om_r <- reactive({
    req(input$om_id)
    if (input$om_id == "om_custom_ui") {
      om_custom(
        r_cv     = input$cv_r,
        K_cv     = input$cv_K,
        proc_cv  = input$proc_cv,
        obs_cv   = input$obs_cv,
        proc_rho = input$proc_rho,
        assess_cv = input$assess_cv %||% 0,
        label    = "Custom OM"
      )
    } else {
      om_spec_map[[input$om_id]]
    }
  })


  # ---- Closure labels ----
  output$closure_prop_label <- renderUI({
    months <- input$closure_months %||% 3
    pct    <- round(months / 12 * 100, 1)
    tags$p(style = "color:#555; font-size:12px; margin-top:-8px;",
           sprintf("= %.1f%% effort reduction", pct))
  })



  # ---- 7. Species data tab ----

  # Small note shown above the table (or prompt to upload)
  output$species_note_ui <- renderUI({
    df <- param_data()
    if (is.null(df)) {
      helpText("Upload a CSV file to see species parameters.")
    } else {
      tags$p(style = "color:#555; font-size:13px;",
             sprintf("%d species loaded.", nrow(df)))
    }
  })

  output$species_table <- renderDT({
    df <- param_data()
    req(df)
    disp <- data.frame(
      Species = df$species_name,
      K       = round(df$K_mean),
      r       = round(df$r_mean, 4),
      Dep_ini = round((df$dep_lo + df$dep_hi) / 2, 3),
      q       = formatC(df$q, format = "e", digits = 3),
      stringsAsFactors = FALSE
    )
    if ("Group" %in% names(df))   disp$Group     <- df$Group
    if (any(df$r_is_range))        disp$r_range   <- paste0("[", round(df$r_lo, 4), ", ", round(df$r_hi, 4), "]")
    if (any(df$K_is_range))        disp$K_range   <- paste0("[", round(df$K_lo),    ", ", round(df$K_hi),    "]")
    if (any(df$dep_is_range))      disp$dep_range <- paste0("[", round(df$dep_lo, 3), ", ", round(df$dep_hi, 3), "]")
    datatable(disp, rownames = FALSE, options = list(
      scrollX = TRUE, pageLength = 25, dom = "tip"
    ))
  })

  # Dynamic header showing E_mmsy
  output$mmsy_header_ui <- renderUI({
    mmsy <- mmsy_r()
    req(mmsy)
    h4(sprintf("Multispecies MSY reference points  (E_mmsy = %.5f)", mmsy$E_mmsy))
  })

  output$mmsy_table <- renderDT({
    mmsy <- mmsy_r()
    req(mmsy)
    df <- param_data(); req(df)
    tbl <- data.frame(
      Species    = df$species_name,
      Bmsy_multi = round(mmsy$Bmsy_multi_vec),
      Bmsy_indep = round(mmsy$Bmsy_individual_vec),
      Emsy_indep = round(mmsy$Emsy_i_if_independent, 4),
      MSY_indep  = round(mmsy$MSY_i_if_independent, 2)
    )
    datatable(tbl, rownames = FALSE, options = list(
      scrollX = TRUE, pageLength = 25, dom = "tip"
    ))
  })

  # ---- 7b. Initialisation tab outputs ----
  diag_r <- reactive({
    df   <- param_data(); req(df)
    mmsy <- mmsy_r();     req(mmsy)
    E_in <- E_init_r();   req(!is.na(E_in))

    # dep_mean: midpoint of stated range
    dep_mean <- (df$dep_lo + df$dep_hi) / 2

    conditioning_diagnostics(
      r            = df$r_mean,
      K            = df$K_mean,
      q            = df$q,
      dep_ini      = dep_mean,
      E_init       = E_in,
      species_names = df$species_name
    )
  })

  output$effort_ref_table <- renderDT({
    diag <- diag_r(); req(diag)
    tbl  <- effort_reference_table(diag)
    datatable(tbl, rownames = FALSE,
              options = list(dom = "t", ordering = FALSE)) |>
      formatRound("Value", digits = 6)
  })

  output$conditioning_table <- renderDT({
    diag <- diag_r(); req(diag)
    tbl  <- diag$table

    # Colour residuals: green if |resid| <= 0.10, red otherwise
    dt <- datatable(tbl, rownames = FALSE,
                    options = list(dom = "t", ordering = FALSE,
                                   pageLength = 25)) |>
      formatRound(c("Dep_stated", "Dep_at_Emmsy",
                    "Dep_at_Eimplied", "Dep_at_Einit",
                    "Residual"), digits = 3) |>
      formatStyle(
        "Residual",
        backgroundColor = styleInterval(
          c(-0.10, 0.10),
          c("#f4cccc", "#ffffff", "#d9ead3")  # red / white / green
        )
      )
    dt
  })

  output$conditioning_message_ui <- renderUI({
    diag <- diag_r(); req(diag)
    conditioning_message(diag)
  })

  output$conditioning_status_ui <- renderUI({
    diag <- diag_r()
    if (is.null(diag)) {
      tags$p(style = "color:#aaa;",
             "Upload a CSV and set E_init to see the diagnostics.")
    } else if (diag$consistency_flag) {
      tags$p(style = "color:#2ca02c; font-weight:bold;",
             sprintf("Ready — using E_init = %.5f", diag$E_init))
    } else {
      tags$p(style = "color:#e6a817; font-weight:bold;",
             sprintf("Proceed with caution — E_init = %.5f
                    but some species are inconsistent.
                    Review the table above.", diag$E_init))
    }
  })

  # ---- 8. HCR plot ----
  plot_HCR_reactive <- reactive({
    mp <- selected_mp_r(); req(mp)
    # Use the hcr_hockeystick from hcr.R
    type  <- mp$hcr_type
    p     <- mp$params
    d_seq <- seq(0, 1.2, length.out = 300)

    if (type == "constant_effort") {
      E_val <- p$E %||% 0
      df_hcr <- data.frame(d = d_seq, E = E_val)
      ggplot(df_hcr, aes(d, E)) +
        geom_line(linewidth = 1.5, colour = "#2166ac") +
        theme_minimal(base_size = 14) +
        labs(x = "Depletion signal (B/K)", y = "Effort",
             title = paste("HCR:", mp$name)) +
        ylim(0, max(E_val * 1.2, 0.1))
    } else if (type %in% c("B_hs_min","B_hs_mean","index_hs_min","index_hs_mean","B_hs_low_r")) {
      dlim <- p$d_lim %||% 0.2; dtrg <- p$d_trig %||% 0.4
      Emax <- p$Emax  %||% 1;   Emin <- p$Emin   %||% 0
      E_vals <- sapply(d_seq, hcr_hockeystick,
                       d_lim = dlim, d_trig = dtrg, Emax = Emax, Emin = Emin)
      df_hcr <- data.frame(d = d_seq, E = E_vals)
      ggplot(df_hcr, aes(d, E)) +
        geom_vline(xintercept = c(dlim, dtrg), linetype = "dashed", colour = "grey60") +
        geom_line(linewidth = 1.5, colour = "#2166ac") +
        annotate("text", x = dlim, y = Emax * 0.05, label = "Blim", hjust = -0.1, size = 3.5) +
        annotate("text", x = dtrg, y = Emax * 0.05, label = "Btrig", hjust = -0.1, size = 3.5) +
        theme_minimal(base_size = 14) +
        labs(x = "Depletion signal (B/K)", y = "Effort",
             title = paste("HCR:", mp$name))
    } else if (type == "seasonal_closure") {
      prop   <- p$prop_season %||% 0.25
      E_base <- p$E_base %||% 1
      E_eff  <- E_base * (1 - prop)
      df_hcr <- data.frame(d = d_seq, E = E_eff)
      ggplot(df_hcr, aes(d, E)) +
        geom_line(linewidth = 1.5, colour = "#d6604d") +
        theme_minimal(base_size = 14) +
        labs(x = "Depletion signal (B/K)", y = "Effective effort",
             title = paste("HCR:", mp$name),
             subtitle = sprintf("Fixed effort reduction: %.0f%%", prop * 100))
    } else if (type == "spatial_closure") {
      df        <- param_data()
      prop_vec  <- if (!is.null(df) && !any(is.na(df$prop_B_area))) df$prop_B_area
      else rep(p$prop_closed %||% 0.3, length.out = 1)
      spp_names <- if (!is.null(df)) df$species_name else paste0("sp", seq_along(prop_vec))
      tbl_df    <- data.frame(Species    = spp_names,
                              PropB_area = prop_vec,
                              stringsAsFactors = FALSE)
      ggplot(tbl_df, aes(x = reorder(Species, -PropB_area), y = PropB_area)) +
        geom_col(fill = "#d6604d", width = 0.6) +
        geom_text(aes(label = scales::percent(PropB_area, accuracy = 1)),
                  vjust = -0.4, size = 4) +
        scale_y_continuous(labels = scales::percent_format(),
                           limits = c(0, 1.1), expand = c(0, 0)) +
        theme_minimal(base_size = 14) +
        labs(x = NULL, y = "Proportion of biomass affected",
             title = "Spatial closure — species-specific biomass impact",
             subtitle = "Proportion of each species' biomass in the closed area (from CSV PropB_area)")
    } else {
      ggplot() + theme_void() +
        labs(title = paste(mp$name, "— shape not displayable (feedback rule)"),
             subtitle = "This MP adjusts effort dynamically; no static HCR curve to show.")
    }
  })

  output$plot_HCR <- renderPlot({ plot_HCR_reactive() })

  output$download_HCR <- downloadHandler(
    filename = function() paste0("HCR_", input$mp_id, ".png"),
    content  = function(file) ggsave(file, plot = plot_HCR_reactive(), width = 7, height = 5)
  )

  # ---- 9. Run simulation ---- Only fires when the Run button is pressed
  sim_res <- eventReactive(input$run_sim, {
    df  <- param_data();  req(df)
    mp  <- selected_mp_r()
    om  <- selected_om_r()

    withProgress(message = "Running MSE simulation...", value = 0, {
      incProgress(0.1, detail = "Setting up parameters")

      # --- Conditioning ---
      cond <- list(
        type   = "bk_range",
        dep_lo = df$dep_lo,
        dep_hi = df$dep_hi
      )

      # --- cv_r / cv_K from CSV ranges (used only when om_spec doesn't override) ---
      # Range in CSV  -> derive CV from the range (uniform distribution formula)
      # No range      -> 0.2 (dynamics.R default)
      # om_spec overrides these inside dynamics.R for all OM-controlled CVs,
      # including setting them to 0 for om_deterministic.
      cv_r_csv <- if (any(df$r_is_range))
        mean((df$r_hi - df$r_lo) / (2 * df$r_mean * sqrt(3)), na.rm = TRUE)
      else 0.2

      cv_K_csv <- if (any(df$K_is_range))
        mean((df$K_hi - df$K_lo) / (2 * df$K_mean * sqrt(3)), na.rm = TRUE)
      else 0.2

      # assess_cv: 0 when explicit dep bounds, 0.1 otherwise
      # (dynamics.R sets it to 0 when om_spec is deterministic)
      assess_cv_use <- if (isTRUE(df$dep_explicit_bounds[1])) 0 else 0.1

      mmsy_vals <- mmsy_r(); req(mmsy_vals)

      incProgress(0.3, detail = "Simulating")

      res <- simulate_multispp(
        nyears_proj   = input$nyears_proj,
        nyears_cond   = input$nyears_cond,
        r             = df$r_mean,
        K             = df$K_mean,
        q             = df$q,
        species_names = df$species_name,
        shared_fleet  = input$shared_fleet,
        hcr           = mp,
        conditioning  = cond,
        om_spec       = om,           # full om_spec object — overrides all CV args
        cv_r          = cv_r_csv,     # fallback if om_spec not provided
        cv_K          = cv_K_csv,     # fallback if om_spec not provided
        assess_cv     = assess_cv_use,
        E_init        = E_init_r(),
        nsims         = input$nsims,
        seed          = 42
      )
      incProgress(1.0, detail = "Done")
      res
    })
  })

  # ---- Run info ----
  output$run_info_box <- renderUI({
    dt <- sim_res()
    if (is.null(dt)) return(NULL)
    cfg <- attr(dt, "config")
    div(class = "run-info",
        tags$b("Last run:"), br(),
        sprintf("MP: %s", input$mp_id), br(),
        sprintf("OM: %s", input$om_id), br(),
        sprintf("Sims: %d | Proj years: %d | Cond years: %d",
                input$nsims, input$nyears_proj, input$nyears_cond)
    )
  })

  # ---- 10. Simulation plots ----
  output$plot_B     <- renderPlot({ dt <- sim_res(); req(dt); plot_B_fun(proj_only(dt)) })
  output$plot_C     <- renderPlot({ dt <- sim_res(); req(dt); plot_C_fun(proj_only(dt)) })
  output$plot_E     <- renderPlot({ dt <- sim_res(); req(dt); plot_E_fun(proj_only(dt)) })
  output$plot_Cprop <- renderPlot({ dt <- sim_res(); req(dt); plot_Cprop_fun(proj_only(dt)) })

  output$download_B <- downloadHandler(
    filename = function() "depletion_plot.png",
    content  = function(file) { dt <- sim_res(); req(dt)
    ggsave(file, plot = plot_B_fun(proj_only(dt)), width = 7, height = 5) }
  )
  output$download_C <- downloadHandler(
    filename = function() "catch_plot.png",
    content  = function(file) { dt <- sim_res(); req(dt)
    ggsave(file, plot = plot_C_fun(proj_only(dt)), width = 7, height = 5) }
  )
  output$download_E <- downloadHandler(
    filename = function() "effort_plot.png",
    content  = function(file) { dt <- sim_res(); req(dt)
    ggsave(file, plot = plot_E_fun(proj_only(dt)), width = 7, height = 5) }
  )
  output$download_Cprop <- downloadHandler(
    filename = function() "catch_composition_plot.png",
    content  = function(file) { dt <- sim_res(); req(dt)
    ggsave(file, plot = plot_Cprop_fun(proj_only(dt)), width = 7, height = 5) }
  )
  output$download_sim_zip <- downloadHandler(
    filename = function() paste0("simulation_plots_", input$mp_id, "_", input$om_id, ".zip"),
    contentType = "application/zip",
    content  = function(file) {
      dt <- sim_res(); req(dt)
      tmp <- tempdir()
      fs <- c(file.path(tmp, "depletion.png"),
              file.path(tmp, "catch.png"),
              file.path(tmp, "effort.png"),
              file.path(tmp, "catch_composition.png"))
      ggsave(fs[1], plot = plot_B_fun(proj_only(dt)),     width = 7, height = 5)
      ggsave(fs[2], plot = plot_C_fun(proj_only(dt)),     width = 7, height = 5)
      ggsave(fs[3], plot = plot_E_fun(proj_only(dt)),     width = 7, height = 5)
      ggsave(fs[4], plot = plot_Cprop_fun(proj_only(dt)), width = 7, height = 5)
      utils::zip(zipfile = file, files = fs, flags = "-j")
    }
  )

  # ---- 11. Metrics ----
  metrics_res <- reactive({
    dt <- sim_res(); req(dt)
    metrics_multispp_10y_safe(proj_only(dt))
  })

  output$metrics_table_summary <- renderDT({
    m <- metrics_res(); req(m)
    df <- as.data.frame(m$summary)
    # Round all numeric columns to 4 decimal places for display
    df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, 4) else x)
    datatable(df, rownames = FALSE,
              options = list(scrollX = TRUE, dom = "tip", pageLength = 5)) |>
      formatRound(columns = names(df)[sapply(df, is.numeric)], digits = 4)
  })

  output$metrics_table_sim <- renderDT({
    m <- metrics_res(); req(m)
    df <- as.data.frame(m$by_sim)
    datatable(df, rownames = FALSE,
              filter = "top",
              options = list(
                scrollX    = TRUE,
                pageLength = 20,
                dom        = "lfrtip"
              )) |>
      formatRound(columns = names(df)[sapply(df, is.numeric)], digits = 4)
  })

  output$download_metrics_sim <- downloadHandler(
    filename = function() paste0("metrics_by_sim_", input$mp_id, "_", input$om_id, ".csv"),
    content  = function(file) { m <- metrics_res(); req(m); write.csv(m$by_sim,  file, row.names = FALSE) }
  )
  output$download_metrics_summ <- downloadHandler(
    filename = function() paste0("metrics_summary_", input$mp_id, "_", input$om_id, ".csv"),
    content  = function(file) { m <- metrics_res(); req(m); write.csv(m$summary, file, row.names = FALSE) }
  )

  # ---- 12. Trade-off accumulator ----
  # Every time you press Run, a new row is appended to tradeoff_store with the median and 10th/90th percentiles of the key metrics for that run.
  tradeoff_store <- reactiveVal(data.frame())

  observeEvent(input$run_sim, {
    dt <- sim_res(); req(dt)
    m  <- metrics_multispp_10y_safe(proj_only(dt))
    bs <- m$by_sim
    if (nrow(bs) == 0) return()

    newrow <- data.frame(
      run_id = paste0("Run_", format(Sys.time(), "%Y%m%d_%H%M%S")),
      MP     = input$mp_id,
      OM     = input$om_id,
      run_label = if (nchar(trimws(input$run_label)) > 0)
        trimws(input$run_label)
      else
        input$mp_id,
      med_AAV          = median(bs$AAV,                    na.rm = TRUE),
      p10_mean_AAV     = quantile(bs$AAV,             0.10, na.rm = TRUE),
      p90_mean_AAV     = quantile(bs$AAV,             0.90, na.rm = TRUE),

      med_mean_catch   = median(bs$mean_catch_10yr,        na.rm = TRUE),
      p10_mean_catch   = quantile(bs$mean_catch_10yr, 0.10, na.rm = TRUE),
      p90_mean_catch   = quantile(bs$mean_catch_10yr, 0.90, na.rm = TRUE),

      med_prop_below_Blim  = median(bs$prop_below_Blim,        na.rm = TRUE),
      p10_prop_below_Blim  = quantile(bs$prop_below_Blim, 0.10, na.rm = TRUE),
      p90_prop_below_Blim  = quantile(bs$prop_below_Blim, 0.90, na.rm = TRUE),

      med_prop_below_Bmsy  = median(bs$prop_below_Bmsy,        na.rm = TRUE),
      p10_prop_below_Bmsy  = quantile(bs$prop_below_Bmsy, 0.10, na.rm = TRUE),
      p90_prop_below_Bmsy  = quantile(bs$prop_below_Bmsy, 0.90, na.rm = TRUE),

      med_prop_years_2_3   = median(bs$prop_years_2_3_above_Blim,        na.rm = TRUE),
      p10_prop_years_2_3   = quantile(bs$prop_years_2_3_above_Blim, 0.10, na.rm = TRUE),
      p90_prop_years_2_3   = quantile(bs$prop_years_2_3_above_Blim, 0.90, na.rm = TRUE),

      stringsAsFactors = FALSE
    )
    tradeoff_store(dplyr::bind_rows(tradeoff_store(), newrow))
  })

  observeEvent(input$clear_runs, { tradeoff_store(data.frame()) })

  # ---- 13. Trade-off plot helpers ----
  tradeoff_theme <- function() {
    theme_minimal(base_size = 13) +
      theme(legend.position = "bottom",
            legend.title = element_text(size = 10))
  }

  plot_tradeoff_base <- function(df, y_col, y_lo, y_hi, y_lab) {
    if (nrow(df) == 0) return(ggplot() + theme_void() + ggtitle("No runs yet — press Run"))
    ggplot(df, aes(x = med_mean_catch, y = .data[[y_col]],
                   colour = run_label, shape = OM, label = run_label)) +
      geom_errorbar( aes(ymin = .data[[y_lo]], ymax = .data[[y_hi]]),
                     width = 0, linewidth = 0.8) +
      geom_errorbarh(aes(xmin = p10_mean_catch, xmax = p90_mean_catch),
                     height = 0, linewidth = 0.8) +
      geom_point(size = 4) +
      tradeoff_theme() +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Median mean catch (last 10 yr)", y = y_lab)
  }

  plot_tradeoff1_fun <- function(df)
    plot_tradeoff_base(df, "med_prop_below_Blim",
                       "p10_prop_below_Blim", "p90_prop_below_Blim",
                       "Proportion below Blim")

  plot_tradeoff2_fun <- function(df)
    plot_tradeoff_base(df, "med_AAV",
                       "p10_mean_AAV", "p90_mean_AAV",
                       "Average Annual Variability (AAV)")

  plot_tradeoff3_fun <- function(df)
    plot_tradeoff_base(df, "med_prop_below_Bmsy",
                       "p10_prop_below_Bmsy", "p90_prop_below_Bmsy",
                       "Proportion below Bmsy")

  plot_tradeoff4_fun <- function(df)
    plot_tradeoff_base(df, "med_prop_years_2_3",
                       "p10_prop_years_2_3", "p90_prop_years_2_3",
                       "Prop. years ≥ 2/3 spp above Blim")

  output$tradeoffPlot      <- renderPlot({ plot_tradeoff1_fun(tradeoff_store()) })
  output$tradeoffPlot_AAV  <- renderPlot({ plot_tradeoff2_fun(tradeoff_store()) })
  output$tradeoffPlot_Bmsy <- renderPlot({ plot_tradeoff3_fun(tradeoff_store()) })
  output$tradeoffPlot_2_3  <- renderPlot({ plot_tradeoff4_fun(tradeoff_store()) })

  make_tradeoff_dl <- function(plot_fn, fname) {
    downloadHandler(
      filename = function() fname,
      content  = function(file) {
        df <- tradeoff_store(); req(nrow(df) > 0)
        ggsave(file, plot = plot_fn(df), width = 7, height = 5)
      }
    )
  }

  output$download_tradeoff1 <- make_tradeoff_dl(plot_tradeoff1_fun, "tradeoff_catch_vs_blim.png")
  output$download_tradeoff2 <- make_tradeoff_dl(plot_tradeoff2_fun, "tradeoff_catch_vs_AAV.png")
  output$download_tradeoff3 <- make_tradeoff_dl(plot_tradeoff3_fun, "tradeoff_catch_vs_Bmsy.png")
  output$download_tradeoff4 <- make_tradeoff_dl(plot_tradeoff4_fun, "tradeoff_catch_vs_2_3rule.png")

  output$download_tradeoff_zip <- downloadHandler(
    filename    = function() "tradeoff_plots.zip",
    contentType = "application/zip",
    content     = function(file) {
      df <- tradeoff_store(); req(nrow(df) > 0)
      tmp <- tempdir()
      fs  <- c(file.path(tmp, "tradeoff_blim.png"),
               file.path(tmp, "tradeoff_AAV.png"),
               file.path(tmp, "tradeoff_Bmsy.png"),
               file.path(tmp, "tradeoff_2_3.png"))
      ggsave(fs[1], plot = plot_tradeoff1_fun(df), width = 7, height = 5)
      ggsave(fs[2], plot = plot_tradeoff2_fun(df), width = 7, height = 5)
      ggsave(fs[3], plot = plot_tradeoff3_fun(df), width = 7, height = 5)
      ggsave(fs[4], plot = plot_tradeoff4_fun(df), width = 7, height = 5)
      utils::zip(zipfile = file, files = fs, flags = "-j")
    }
  )

} # server

shinyApp(ui, server)
