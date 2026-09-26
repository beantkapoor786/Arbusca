# Arbusca — AMF 18S/SSU DADA2 pipeline GUI. Raw reads -> ASVs -> MaarjAM
# taxonomy -> phyloseq, then figures and statistics (PERMANOVA, ANCOM-BC2).

library(shiny)
library(bslib)
library(bsicons)

source("R/theme.R")
source("R/app_config.R")
source("R/state_engine.R")
source("R/proc_runner.R")
source("R/paths.R")
source("R/mod_logpanel.R")
source("R/mod_setup.R")
source("R/mod_qcheck.R")
source("R/mod_primer.R")
source("R/mod_qc.R")
source("R/mod_denoise.R")
source("R/mod_merge.R")
source("R/mod_chimera.R")
source("R/mod_taxonomy.R")
source("R/mod_output.R")
source("R/mod_figures.R")
source("R/mod_permanova.R")
source("R/mod_ancombc.R")

STEPS <- c(
  qcheck   = "2. Quality Check",
  primer   = "3. Primer Removal",
  qc       = "4. Filter and Trim",
  denoise  = "5. Learn Error Rate and Denoise",
  merge    = "6. Merge Pairs",
  chimera  = "7. Remove Chimeras",
  taxonomy = "8. Taxonomy Assignment using BLAST",
  output   = "9. Phyloseq Object",
  figures  = "10. Figures",
  permanova = "11. PERMANOVA",
  ancombc  = "12. Differential Abundance (ANCOM-BC2)"
)

ALL_STEPS <- c(setup = "1. Setup & Raw Reads", STEPS)

# Short labels for the pipeline stepper; panel titles keep the full names.
STEP_SHORT <- c(
  setup = "Raw reads", qcheck = "Quality check", primer = "Primer removal",
  qc = "Filtering", denoise = "Denoising", merge = "Merge pairs",
  chimera = "Chimera removal", taxonomy = "Taxonomy", output = "Phyloseq",
  figures = "Figures", permanova = "PERMANOVA",
  ancombc = "ANCOM-BC2"
)

next_step_id <- function(step_id) {
  idx <- match(step_id, names(ALL_STEPS))
  if (is.na(idx) || idx == length(ALL_STEPS)) return(NULL)
  names(ALL_STEPS)[idx + 1]
}

ui <- page_fluid(
  title = "Arbusca",
  theme = amf_theme,
  tags$head(tags$link(rel = "stylesheet", href = paste0("styles.css?v=", as.integer(file.mtime("www/styles.css"))))),
  tags$head(tags$script(HTML(
    "Shiny.addCustomMessageHandler('amf-scroll-log', function(id) {
       var el = document.getElementById(id);
       if (el) el.scrollTop = el.scrollHeight;
     });"
  ))),
  tags$header(
    class = "amf-navbar",
    div(class = "amf-nav-left", span(class = "amf-tagline d-none d-lg-inline", "AMF 18S/SSU DADA2 pipeline")),
    div(class = "amf-brand", span(class = "amf-brand-name", "Arbusca"))
  ),
  uiOutput("stepper"),
  div(class = "amf-main", uiOutput("main_panel")),
  uiOutput("proceed_bar")
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    project_dir = NULL,
    active_step = "setup",   # which step's panel is shown (navigation)
    running_step = NULL,     # which step has a live process handle, if any
    handles = list(),        # step_id -> processx/callr handle
    status = list(
      primer = "IDLE", qc = "IDLE", denoise = "IDLE",
      merge = "IDLE", chimera = "IDLE", taxonomy = "IDLE", output = "IDLE"
    ),
    log = list(),             # step_id -> character vector of log lines
    artifacts = list(),       # step_id -> R return value of a callr step
    phyloseq_rev = 0          # bumped whenever 05_output/ changes (see invalidate_analyses)
  )

  setup <- mod_setup_server("setup", rv)

  setup_complete <- reactive({
    !is.null(rv$project_dir) && isTRUE(setup$all_ok()) && nrow(setup$sample_table()) > 0
  })

  qcheck <- mod_qcheck_server("qcheck", rv, setup$sample_table)
  mod_primer_server("primer", rv, setup$sample_table)
  mod_qc_server("qc", rv, setup$sample_table)
  mod_denoise_server("denoise", rv, setup$sample_table)
  mod_merge_server("merge", rv, setup$sample_table)
  mod_chimera_server("chimera", rv, setup$sample_table)
  mod_taxonomy_server("taxonomy", rv, setup$sample_table)
  mod_output_server("output", rv, setup$sample_table)
  figures <- mod_figures_server("figures", rv, setup$sample_table)
  mod_permanova_server("permanova", rv, setup$sample_table)
  mod_ancombc_server("ancombc", rv, setup$sample_table)

  # Resume support (DESIGN.md 4.4): whenever a project directory is picked
  # (fresh selection, "Last used" shortcut, or Browse to a different one),
  # check which steps' complete output sets are already on disk and mark
  # them SUCCESS. Only jump the view to the resume point if something was
  # actually already done -- a brand-new empty project just stays on Setup
  # as usual, since there's nothing to resume.
  # Keyed off sample_table() (not just project_dir) because the sample
  # delimiter/index inputs that shape it live in a renderUI block and aren't
  # available client-side the instant project_dir changes -- scanning only
  # on project_dir would run against a stale/default-guessed table and never
  # retry once the real inputs sync in.
  observeEvent(setup$sample_table(), {
    req(rv$project_dir)
    st <- setup$sample_table()
    req(nrow(st) > 0)

    already_done <- any(isolate(unlist(rv$status)) == "SUCCESS")
    scan_checkpoints(rv, st)
    now_done <- any(isolate(unlist(rv$status)) == "SUCCESS")

    if (now_done) {
      resume_step <- first_incomplete_step(rv)
      rv$active_step <- resume_step
      if (!already_done) {
        showNotification(paste("Resuming from", STEPS[[resume_step]]), type = "message", duration = 5)
      }
    }
  }, ignoreInit = TRUE)

  # Single polling loop drives progress + live logs for whichever step is
  # RUNNING. invalidateLater only fires while something is active, so an
  # idle app costs nothing.
  observe({
    if (any_running(rv)) {
      invalidateLater(750, session)
      poll_step(rv, rv$running_step)
    }
  })

  # Presentation-only marker: Quality Check is optional/exploratory, so once
  # the user has moved on to a later step for the current project it counts as
  # done in the stepper even if profiles weren't generated this session.
  qcheck_passed <- reactiveVal(NULL)   # project_dir it was passed for
  observe({
    idx <- match(rv$active_step, names(ALL_STEPS))
    if (!is.null(rv$project_dir) && !is.na(idx) && idx > match("qcheck", names(ALL_STEPS))) {
      qcheck_passed(rv$project_dir)
    }
  })

  # Horizontal pipeline tracker across the top of the page. Presentation only:
  # states are derived from the same reactives the old badges used.
  output$stepper <- renderUI({
    step_ids <- names(ALL_STEPS)
    items <- lapply(seq_along(step_ids), function(i) {
      step_id <- step_ids[[i]]
      is_active <- identical(rv$active_step, step_id)
      state <- if (identical(step_id, "setup")) {
        if (!isTRUE(setup$all_ok())) "ERROR"
        else if (isTRUE(setup_complete())) "SUCCESS"
        else "IDLE"
      } else if (is.null(rv$project_dir)) {
        "IDLE"
      } else if (identical(step_id, "qcheck")) {
        if (isTRUE(qcheck$done()) || identical(qcheck_passed(), rv$project_dir)) "SUCCESS" else "IDLE"
      } else if (identical(step_id, "figures")) {
        if (isTRUE(figures$done())) "SUCCESS" else "IDLE"
      } else {
        rv$status[[step_id]] %||% "IDLE"
      }
      marker <- switch(state,
        SUCCESS = bs_icon("check-lg"),
        ERROR   = bs_icon("x-lg"),
        RUNNING = bs_icon("arrow-repeat", class = "amf-spin"),
        QUEUED  = bs_icon("hourglass-split"),
        as.character(i)
      )
      tagList(
        if (i > 1) span(class = paste("amf-step-line", if (state == "SUCCESS") "is-done" else ""), `aria-hidden` = "true"),
        actionLink(
          inputId = paste0("nav_", step_id),
          label = tagList(
            span(class = "amf-step-marker", marker),
            span(class = "amf-step-label", STEP_SHORT[[step_id]])
          ),
          class = paste("amf-step", paste0("is-", tolower(state)), if (is_active) "is-active" else ""),
          `aria-current` = if (is_active) "step" else NULL,
          title = ALL_STEPS[[step_id]]
        )
      )
    })
    div(class = "amf-stepper", role = "navigation", `aria-label` = "Pipeline steps", items)
  })

  lapply(names(ALL_STEPS), function(step_id) {
    observeEvent(input[[paste0("nav_", step_id)]], {
      rv$active_step <- step_id
    })
  })

  output$main_panel <- renderUI({
    ui_fn <- get(paste0("mod_", rv$active_step, "_ui"))
    ui_fn(rv$active_step)
  })

  # Current step's completion state, used to gate the Proceed button.
  current_step_complete <- reactive({
    if (identical(rv$active_step, "setup")) {
      setup_complete()
    } else if (identical(rv$active_step, "qcheck")) {
      TRUE
    } else if (identical(rv$active_step, "figures")) {
      TRUE
    } else {
      identical(rv$status[[rv$active_step]], "SUCCESS")
    }
  })

  output$proceed_bar <- renderUI({
    nxt <- next_step_id(rv$active_step)
    if (is.null(nxt)) return(NULL)
    complete <- isTRUE(current_step_complete())
    div(
      class = "amf-proceed",
      actionButton(
        "proceed_btn",
        label = tagList(paste("Proceed to", ALL_STEPS[[nxt]]), bs_icon("arrow-right")),
        class = paste(if (complete) "btn-primary" else "btn-secondary disabled"),
        disabled = if (!complete) "disabled" else NULL
      )
    )
  })

  observeEvent(input$proceed_btn, {
    req(isTRUE(current_step_complete()))
    nxt <- next_step_id(rv$active_step)
    req(!is.null(nxt))
    rv$active_step <- nxt
  })
}

shinyApp(ui, server)
