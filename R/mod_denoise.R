# Step 5 -- Learn Error Rate and Denoise (DADA2 core), per DESIGN.md section 5.
# Two sub-stages on one page:
#   5a. Learn error rate -- learnErrors(F/R), shown via plotErrors. Runs as
#       its own self-contained job (not part of the pipeline FSM/progress
#       count), since it's a prerequisite the user reviews before denoising,
#       not the pipeline step itself.
#   5b. Denoise -- dada(F/R) only, gated on 5a's error models being
#       available (this session, or loaded from disk if resuming). This is
#       the real pipeline step ("denoise"). Merging pairs and removing
#       chimeras are their own separate pipeline steps (mod_merge.R,
#       mod_chimera.R), each reading this step's dadaF.rds/dadaR.rds.

mod_denoise_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "5. Learn Error Rate and Denoise",
    description = "DADA2 first learns the error rates of your sequencing run (5a), then uses them to separate true biological sequences from sequencing errors (5b), producing amplicon sequence variants (ASVs).",
    results = shiny::tagList(
      shiny::uiOutput(ns("input_status")),
      result_card("5a. Learn error rate",
        shiny::div(
          class = "amf-run-row mb-3",
          shiny::actionButton(ns("learn"), "Learn error rate", icon = bsicons::bs_icon("play-fill"), class = "btn-primary"),
          shiny::actionButton(ns("learn_cancel"), "Cancel", class = "btn-outline-secondary btn-sm"),
          shiny::uiOutput(ns("learn_badge"), inline = TRUE)
        ),
        shiny::uiOutput(ns("learn_results")),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("learn_log"))
      ),
      result_card("5b. Denoise",
        shiny::uiOutput(ns("denoise_gate_msg")),
        shiny::fluidRow(
          shiny::column(5,
            shiny::tags$label("Pooling mode"),
            shiny::selectInput(
              ns("pool"), NULL,
              choices = c("Independent (default, fastest)" = "independent", "Pseudo-pooling" = "pseudo", "Pool all samples" = "pool"),
              width = "100%"
            )
          ),
          shiny::column(7, shiny::div(class = "amf-param-hint mt-4", "Independent is fastest. Pseudo-pooling recovers rare ASVs shared across samples at modest cost. Pool all samples is the most sensitive but slowest and most memory-intensive."))
        ),
        run_row(ns, "Denoise"),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("log"))
      ),
      shiny::uiOutput(ns("results_section")),
      results_placeholder("Run 5b to see denoised sequence counts per sample.")
    )
  )
}

# sample_table: reactive() -> data.frame(sample, fwd, rev) of basenames from
# mod_setup. Denoising reads Step 4's filtered output (02_filtered/).
mod_denoise_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "denoise"

    input_files <- shiny::reactive({
      st <- sample_table()
      shiny::req(rv$project_dir)
      # Depending on rv$status$qc (even though unused) makes this recompute
      # the moment Quality Control finishes, instead of caching the empty
      # result from before 02_filtered/ had anything in it.
      rv$status$qc
      filtered_files_df(rv$project_dir, st)
    })

    output$input_status <- shiny::renderUI({
      files <- input_files()
      if (nrow(files) == 0) {
        return(shiny::div(class = "alert alert-warning small", "No filtered reads found yet -- finish Filter and Trim first."))
      }
      shiny::div(class = "text-muted small mb-2", sprintf("%d filtered sample(s) ready from 02_filtered/.", nrow(files)))
    })

    # --- 5a: learn error rate -- self-contained, not part of the pipeline FSM ---

    learn_status <- shiny::reactiveVal("IDLE")
    learn_handle <- shiny::reactiveVal(NULL)
    learn_log <- shiny::reactiveVal(character(0))
    learn_result <- shiny::reactiveVal(NULL)  # list(errF=, errR=) once available

    # Resuming mid-step: if error models are already on disk (learned in a
    # previous session), pick them up without re-running.
    checked_disk <- shiny::reactiveVal(FALSE)
    shiny::observe({
      shiny::req(!checked_disk())
      files <- input_files()
      shiny::req(rv$project_dir, nrow(files) > 0)
      out_dir <- denoise_dir_path(rv$project_dir)
      errF_path <- file.path(out_dir, "errF.rds")
      errR_path <- file.path(out_dir, "errR.rds")
      if (fs::file_exists(errF_path) && fs::file_exists(errR_path)) {
        learn_result(list(errF = readRDS(errF_path), errR = readRDS(errR_path)))
        learn_status("SUCCESS")
      }
      checked_disk(TRUE)
    })

    # 5a lives outside rv$status (see comment above), so reset_downstream()
    # -- which reset_downstream fires on every (re-)launch of an upstream
    # step, per state_engine.R -- has no way to reach it directly. Watch
    # rv$status$qc instead: the moment it drops off SUCCESS (Filter and Trim
    # re-run after e.g. a primer change), whatever error model is currently
    # held/showing is stale, so clear it and require "Learn error rate" to
    # be run again rather than silently no-op'ing against an old SUCCESS.
    shiny::observeEvent(rv$status$qc, {
      shiny::req(!identical(rv$status$qc, "SUCCESS"))
      learn_status("IDLE")
      learn_handle(NULL)
      learn_log(character(0))
      learn_result(NULL)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$learn, {
      shiny::req(!any_running(rv))
      shiny::req(identical(learn_status(), "IDLE") || identical(learn_status(), "ERROR"))
      files <- input_files()
      shiny::req(nrow(files) > 0)

      out_dir <- denoise_dir_path(rv$project_dir)
      fs::dir_create(out_dir)
      multithread <- .Platform$OS.type != "windows"

      learn_job <- function(filtFs, filtRs, out_dir, multithread) {
        library(dada2)
        errF <- learnErrors(filtFs, multithread = multithread, verbose = TRUE)
        errR <- learnErrors(filtRs, multithread = multithread, verbose = TRUE)
        saveRDS(errF, file.path(out_dir, "errF.rds"))
        saveRDS(errR, file.path(out_dir, "errR.rds"))
        list(errF = errF, errR = errR)
      }

      h <- callr::r_bg(
        func = learn_job,
        args = list(filtFs = files$fwd, filtRs = files$rev, out_dir = out_dir, multithread = multithread),
        stdout = "|", stderr = "|", supervise = TRUE
      )
      learn_handle(h)
      learn_status("RUNNING")
      learn_log(character(0))
      learn_result(NULL)
    })

    shiny::observeEvent(input$learn_cancel, {
      h <- learn_handle()
      if (!is.null(h) && h$is_alive()) h$kill()
      learn_status("IDLE")
      learn_log(character(0))
      learn_handle(NULL)
    })

    shiny::observe({
      h <- learn_handle()
      shiny::req(!is.null(h), identical(learn_status(), "RUNNING"))
      shiny::invalidateLater(750, session)

      new_out <- tryCatch(h$read_output_lines(), error = function(e) character(0))
      new_err <- tryCatch(h$read_error_lines(), error = function(e) character(0))
      if (length(c(new_out, new_err)) > 0) learn_log(c(learn_log(), new_out, new_err))

      if (!h$is_alive()) {
        result <- tryCatch(h$get_result(), error = function(e) e)
        if (inherits(result, "error")) {
          learn_status("ERROR")
          learn_log(c(learn_log(), paste("ERROR:", conditionMessage(result))))
        } else {
          learn_result(result)
          learn_status("SUCCESS")
        }
      }
    })

    output$learn_badge <- shiny::renderUI(status_badge(learn_status()))
    mod_logpanel_server("learn_log", shiny::reactive(learn_log()))

    output$learn_results <- shiny::renderUI({
      res <- learn_result()
      shiny::req(!is.null(res))
      shiny::fluidRow(
        shiny::column(6, shiny::h6("Forward"), shiny::plotOutput(ns("err_plot_fwd"), height = "420px")),
        shiny::column(6, shiny::h6("Reverse"), shiny::plotOutput(ns("err_plot_rev"), height = "420px"))
      )
    })
    output$err_plot_fwd <- shiny::renderPlot({
      res <- learn_result()
      shiny::req(res)
      dada2::plotErrors(res$errF, nominalQ = TRUE)
    }, res = 96)
    output$err_plot_rev <- shiny::renderPlot({
      res <- learn_result()
      shiny::req(res)
      dada2::plotErrors(res$errR, nominalQ = TRUE)
    }, res = 96)

    # --- 5b: denoise -- the real pipeline step, tied to rv$status$denoise ---

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    output$denoise_gate_msg <- shiny::renderUI({
      if (is.null(learn_result())) {
        shiny::div(class = "text-muted small mb-2", "Run 5a first -- denoising needs a learned error model.")
      } else {
        NULL
      }
    })

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv))
      shiny::req(!is.null(learn_result()))
      files <- input_files()
      shiny::req(nrow(files) > 0)

      out_dir <- denoise_dir_path(rv$project_dir)
      fs::dir_create(out_dir)
      pool_arg <- switch(input$pool, independent = FALSE, pseudo = "pseudo", pool = TRUE, FALSE)
      multithread <- .Platform$OS.type != "windows"
      errF <- learn_result()$errF
      errR <- learn_result()$errR

      denoise_job <- function(samples, filtFs, filtRs, errF, errR, out_dir, pool, multithread) {
        library(dada2)
        names(filtFs) <- samples
        names(filtRs) <- samples

        # callr::r_bg runs this in a bare background R session -- it does
        # NOT inherit helper functions from the app's global environment,
        # even ones defined in this same file, so this has to be self-contained.
        normalize_dada_list <- function(x, sample_names) {
          if (inherits(x, "dada")) stats::setNames(list(x), sample_names) else x
        }

        dadaFs <- dada(filtFs, err = errF, multithread = multithread, pool = pool, verbose = TRUE)
        dadaRs <- dada(filtRs, err = errR, multithread = multithread, pool = pool, verbose = TRUE)
        dadaFs <- normalize_dada_list(dadaFs, samples)
        dadaRs <- normalize_dada_list(dadaRs, samples)

        saveRDS(dadaFs, file.path(out_dir, "dadaF.rds"))
        saveRDS(dadaRs, file.path(out_dir, "dadaR.rds"))

        getN <- function(x) sum(dada2::getUniques(x))
        summary_df <- data.frame(
          sample = samples,
          denoisedF = as.numeric(sapply(dadaFs, getN)[samples]),
          denoisedR = as.numeric(sapply(dadaRs, getN)[samples]),
          row.names = NULL
        )

        list(summary = summary_df)
      }

      launch_step_callr(rv, step_id, denoise_job, list(
        samples = files$sample, filtFs = files$fwd, filtRs = files$rev,
        errF = errF, errR = errR, out_dir = out_dir, pool = pool_arg, multithread = multithread
      ))
    })

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    output$results_section <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res))
      result_card("Denoised sequence counts",
        dt_output(ns("summary_table"))
      )
    })

    output$summary_table <- DT::renderDataTable({
      res <- rv$artifacts[[step_id]]
      if (is.null(res)) {
        empty <- data.frame(Sample = character(0), `Denoised (F)` = integer(0), `Denoised (R)` = integer(0), check.names = FALSE)
        return(DT::datatable(empty, rownames = FALSE, options = list(dom = "t")))
      }
      df <- data.frame(
        Sample = res$summary$sample, `Denoised (F)` = res$summary$denoisedF, `Denoised (R)` = res$summary$denoisedR,
        check.names = FALSE
      )
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = nrow(df) + 1))
    })

    mod_logpanel_server("log", shiny::reactive(rv$log[[step_id]] %||% character(0)))
  })
}
