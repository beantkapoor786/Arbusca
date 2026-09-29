# Step 2 -- Quality Check. Pure visualization of the raw (untrimmed) reads'
# quality profiles, so the user can eyeball run quality before anything is
# processed -- distinct from Step 5's "Quality Control" (mod_qc.R), which
# filters/trims and happens to preview the *trimmed* reads. Not a pipeline
# stage with artifacts, so it carries no rv$status entry and is always
# treated as complete (see app.R's setup/qcheck special-casing).

mod_qcheck_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "2. Quality Check",
    description = "Plots the base-quality profile of your raw reads so you can judge run quality before anything is trimmed. Nothing is written to disk in this step.",
    params = shiny::tagList(
      shiny::uiOutput(ns("input_status")),
      param_group("Profile options",
        shiny::fluidRow(
          shiny::column(4, shiny::numericInput(ns("n_samples"), "Number of samples to profile", value = 2, min = 1, width = "100%")),
          shiny::column(8,
            shiny::radioButtons(ns("view_mode"), "View", choices = c("Per sample" = "per_sample", "Aggregated" = "aggregate"), selected = "per_sample", inline = TRUE),
            shiny::div(class = "amf-param-hint", "Per sample facets each sample separately; Aggregated pools them into one combined profile.")
          )
        )
      ),
      shiny::div(
        class = "amf-run-row",
        shiny::actionButton(ns("generate"), "Generate quality profiles", icon = bsicons::bs_icon("play-fill"), class = "btn-primary"),
        shiny::uiOutput(ns("badge"), inline = TRUE)
      )
    ),
    results = shiny::tagList(
      result_card("Raw read quality profiles",
        shiny::uiOutput(ns("profile_error_alert")),
        shiny::uiOutput(ns("plots")),
        actions = shiny::uiOutput(ns("downloads"))
      )
    )
  )
}

# sample_table: reactive() -> data.frame(sample, fwd, rev) of basenames from
# mod_setup. This step reads the raw files directly in project_dir, before
# any primer trimming or filtering has happened.
mod_qcheck_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw_files <- shiny::reactive({
      st <- sample_table()
      shiny::req(rv$project_dir)
      raw_files_df(rv$project_dir, st)
    })

    output$input_status <- shiny::renderUI({
      files <- raw_files()
      if (nrow(files) == 0) {
        return(shiny::div(class = "alert alert-warning small", "No raw read pairs found yet -- finish Setup & Raw Reads first."))
      }
      shiny::div(class = "text-muted small mb-2", sprintf("%d raw sample pair(s) available.", nrow(files)))
    })

    # Rendered as a background callr job (like Step 6's "Learn error rate")
    # so the RUNNING badge actually shows while plotQualityProfile() -- which
    # can take a while over many raw reads -- is still computing, instead of
    # freezing the UI with no feedback until it's done.
    profile_status <- shiny::reactiveVal("IDLE")
    profile_handle <- shiny::reactiveVal(NULL)
    profile_result <- shiny::reactiveVal(NULL)  # list(fwd=, rev=) ggplot objects
    profile_error <- shiny::reactiveVal(NULL)
    profile_n <- shiny::reactiveVal(1)  # facets per plot (1 when aggregated)
    output$profile_error_alert <- shiny::renderUI({
      shiny::req(identical(profile_status(), "ERROR"), profile_error())
      shiny::div(class = "alert alert-danger small mb-2", profile_error())
    })

    shiny::observeEvent(input$generate, {
      shiny::req(!identical(profile_status(), "RUNNING"))
      files <- raw_files()
      shiny::req(nrow(files) > 0)
      n <- max(1, min(input$n_samples, nrow(files)))
      aggregate <- identical(input$view_mode, "aggregate")

      profile_job <- function(fwd_files, rev_files, aggregate) {
        library(dada2)
        list(
          fwd = plotQualityProfile(fwd_files, aggregate = aggregate),
          rev = plotQualityProfile(rev_files, aggregate = aggregate)
        )
      }

      h <- callr::r_bg(
        func = profile_job,
        args = list(fwd_files = files$fwd[seq_len(n)], rev_files = files$rev[seq_len(n)], aggregate = aggregate),
        stdout = "|", stderr = "|", supervise = TRUE
      )
      profile_handle(h)
      profile_n(if (aggregate) 1 else n)
      profile_status("RUNNING")
      profile_result(NULL)
    })

    shiny::observe({
      h <- profile_handle()
      shiny::req(!is.null(h), identical(profile_status(), "RUNNING"))
      shiny::invalidateLater(500, session)
      if (!h$is_alive()) {
        result <- tryCatch(h$get_result(), error = function(e) e)
        if (inherits(result, "error")) {
          profile_error(conditionMessage(result))
          profile_status("ERROR")
        } else {
          profile_error(NULL)
          profile_result(result)
          profile_status("SUCCESS")
        }
      }
    })

    output$badge <- shiny::renderUI(status_badge(profile_status()))

    # plotQualityProfile facets with facet_wrap's default grid
    # (ceiling(sqrt(n)) columns), so size the plot to its facet rows. Beyond
    # 2 samples the facets get cramped in half-width columns, so Forward and
    # Reverse are stacked full-width instead of side by side.
    plot_height <- shiny::reactive({
      n <- profile_n()
      max(360, ceiling(n / ceiling(sqrt(n))) * 300)
    })
    output$plots <- shiny::renderUI({
      width <- if (profile_n() > 2) 12 else 6
      shiny::fluidRow(
        shiny::column(width, shiny::h6("Forward"), shiny::plotOutput(ns("quality_plot_fwd"), height = "auto")),
        shiny::column(width, shiny::h6("Reverse"), shiny::plotOutput(ns("quality_plot_rev"), height = "auto"))
      )
    })
    output$quality_plot_fwd <- shiny::renderPlot({
      shiny::req(profile_result())
      profile_result()$fwd
    }, height = function() plot_height())
    output$quality_plot_rev <- shiny::renderPlot({
      shiny::req(profile_result())
      profile_result()$rev
    }, height = function() plot_height())

    # Only offer downloads once profiles actually exist for this session.
    output$downloads <- shiny::renderUI({
      shiny::req(profile_result())
      shiny::div(
        class = "d-flex gap-2 mt-2",
        shiny::downloadButton(ns("dl_fwd"), "Download forward plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm"),
        shiny::downloadButton(ns("dl_rev"), "Download reverse plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm")
      )
    })

    output$dl_fwd <- shiny::downloadHandler(
      filename = function() "quality_profile_forward.png",
      content = function(file) ggplot2::ggsave(file, plot = profile_result()$fwd, dpi = 300, width = 8, height = 6, units = "in")
    )
    output$dl_rev <- shiny::downloadHandler(
      filename = function() "quality_profile_reverse.png",
      content = function(file) ggplot2::ggsave(file, plot = profile_result()$rev, dpi = 300, width = 8, height = 6, units = "in")
    )

    # Exposed so the stepper badge can show SUCCESS only once profiles have
    # actually been generated -- not on launch, before anything has run.
    done <- shiny::reactive(identical(profile_status(), "SUCCESS"))

    list(done = done)
  })
}
