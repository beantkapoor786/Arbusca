# Central bslib theme + reusable UI helpers. Visual rules live in
# www/styles.css; this file defines the palette/typography tokens and the
# small helper functions modules use to build consistent markup.

# ---- Palette (5 colours + neutrals) -----------------------------------------
AMF_COLORS <- list(
  primary = "#0f5f73",   # deep scientific teal-blue: actions, active state
  success = "#2f7d4f",   # muted green: completed steps
  warning = "#a86a00",   # amber (AA on white)
  danger  = "#b3372b",   # muted red
  info    = "#5b6b78",   # slate: neutral/secondary emphasis
  bg      = "#f5f6f7",   # page background
  fg      = "#1f2933"    # body text
)

amf_theme <- bslib::bs_theme(
  version = 5,
  bg = AMF_COLORS$bg,
  fg = AMF_COLORS$fg,
  primary = AMF_COLORS$primary,
  secondary = AMF_COLORS$info,
  success = AMF_COLORS$success,
  warning = AMF_COLORS$warning,
  danger = AMF_COLORS$danger,
  info = AMF_COLORS$info,
  base_font = bslib::font_collection(bslib::font_google("Inter"), "system-ui", "-apple-system", "Segoe UI", "sans-serif"),
  code_font = bslib::font_collection(bslib::font_google("JetBrains Mono"), "SFMono-Regular", "Menlo", "Consolas", "monospace"),
  "font-size-base" = "0.9rem",
  "border-radius" = "0.4rem",
  "border-radius-sm" = "0.3rem",
  "border-radius-lg" = "0.5rem",
  "card-border-color" = "#dde1e5",
  "card-cap-bg" = "#ffffff",
  "card-bg" = "#ffffff",
  "input-border-color" = "#c4cbd1",
  "focus-ring-width" = "0.2rem",
  "focus-ring-opacity" = "0.35"
)

# ---- Helpers -----------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# Renders a status badge span for a given FSM state (QUEUED/RUNNING/SUCCESS/
# ERROR) -- NULL for IDLE, since "hasn't run yet" isn't worth a badge and
# showing one there reads as if the step already exists and just hasn't run.
status_badge <- function(status) {
  status <- toupper(status %||% "IDLE")
  if (identical(status, "IDLE")) return(NULL)
  cls <- switch(status,
    "QUEUED"  = "amf-status-queued",
    "RUNNING" = "amf-status-running",
    "SUCCESS" = "amf-status-success",
    "ERROR"   = "amf-status-error",
    "amf-status-idle"
  )
  shiny::span(class = paste("amf-status-badge", cls), status)
}

# A small "?" circle; hovering/focusing shows help_text in a bslib tooltip.
help_icon <- function(help_text) {
  bslib::tooltip(
    shiny::tags$span(class = "amf-help-icon", tabindex = "0", `aria-label` = "Help", "?"),
    help_text
  )
}

# A form-field label with a help icon. Use in place of a plain string label,
# e.g. numericInput(ns("x"), help_label("X", "what X does"), ...).
help_label <- function(label, help_text) {
  shiny::tags$label(label, help_icon(help_text))
}

# DT's wrapper div contains floated internal elements that collapse to zero
# height under Bootstrap 5 (via bslib), letting whatever comes next in the
# page render up into the space the table visually occupies. `overflow:
# auto` establishes a new block-formatting context that contains the floats,
# which fixes it. Wrap every DT::dataTableOutput() in the app with this.
dt_output <- function(outputId) {
  shiny::div(class = "amf-dt", DT::dataTableOutput(outputId))
}

# ---- Step layout helpers -----------------------------------------------------

# Standard step panel: title + plain-language description on top, then a
# "Parameters" card (left) and stacked result cards (right); stacked = TRUE
# gives a single full-width column (parameters above results). `params` and
# `results` are UI (tagLists); modules keep their own input/output IDs.
step_card <- function(title, description, params = NULL, results, stacked = FALSE) {
  shiny::tagList(
    shiny::div(
      class = "amf-step-head",
      shiny::h4(title, class = "amf-step-title"),
      shiny::p(description, class = "amf-step-desc")
    ),
    if (is.null(params)) return_results_only(results) else bslib::layout_columns(
      col_widths = if (stacked) 12 else c(5, 7),
      fill = FALSE,
      bslib::card(
        class = "amf-params-card",
        bslib::card_header(bsicons::bs_icon("sliders"), "Parameters"),
        bslib::card_body(params)
      ),
      shiny::div(class = "amf-results-col", results)
    )
  )
}

# Standalone parameters card, for steps that place it after another card.
params_card <- function(...) {
  bslib::card(
    class = "amf-params-card",
    bslib::card_header(bsicons::bs_icon("sliders"), "Parameters"),
    bslib::card_body(...)
  )
}

# Results-only body for steps without a shared parameters card.
return_results_only <- function(results) shiny::div(class = "amf-results-col", results)

# A titled card for outputs (plots, tables, log). `actions` (e.g. download
# buttons) is shown quietly in the header, right-aligned.
result_card <- function(title, ..., actions = NULL, icon = NULL) {
  bslib::card(
    class = "amf-result-card",
    bslib::card_header(
      shiny::div(
        class = "d-flex justify-content-between align-items-center gap-2",
        shiny::span(icon, title),
        actions
      )
    ),
    bslib::card_body(...)
  )
}

# A labelled sub-section within the parameters card.
param_group <- function(title, ...) {
  shiny::div(class = "amf-param-group", shiny::div(class = "amf-param-group-title", title), ...)
}

# Run / Cancel / status row. Run is the one prominent action; pass the
# module's own label (e.g. "Run Cutadapt") via run_label.
run_row <- function(ns, run_label = "Run", extra = NULL) {
  shiny::div(
    class = "amf-run-row",
    shiny::actionButton(ns("run"), run_label, icon = bsicons::bs_icon("play-fill"), class = "btn-primary"),
    shiny::actionButton(ns("cancel"), "Cancel", class = "btn-outline-secondary btn-sm"),
    shiny::uiOutput(ns("badge"), inline = TRUE),
    extra
  )
}

# Friendly placeholder for outputs that don't exist yet.
empty_state <- function(text, icon = "inbox") {
  shiny::div(class = "amf-empty", bsicons::bs_icon(icon), shiny::span(text))
}

# Compact summary tiles. `items` is a named list: name = title, value = text.
stat_boxes <- function(items) {
  boxes <- lapply(names(items), function(nm) {
    bslib::value_box(title = nm, value = items[[nm]], class = "amf-vb", height = "auto")
  })
  do.call(bslib::layout_columns, c(
    list(fill = FALSE, col_widths = bslib::breakpoints(sm = 12, md = rep(floor(12 / length(items)), length(items)))),
    boxes
  ))
}

# "About this step": input, method and output in plain language. Purely static.
step_io <- function(input, method, output) {
  param_group("About this step",
    shiny::tags$dl(
      class = "amf-io",
      shiny::tags$dt("Reads"), shiny::tags$dd(input),
      shiny::tags$dt("Method"), shiny::tags$dd(method),
      shiny::tags$dt("Writes"), shiny::tags$dd(shiny::tags$code(output))
    )
  )
}

# Placeholder shown (via CSS) only while the preceding output is still empty.
results_placeholder <- function(text, icon = "hourglass") {
  shiny::div(class = "amf-placeholder", empty_state(text, icon))
}
