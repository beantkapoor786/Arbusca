# Step 4 -- Quality Control (dada2::filterAndTrim). Per DESIGN.md section 5:
# quality-profile plots before committing truncLen, sliders for
# truncLen/maxEE, maxN=0 and rm.phix=TRUE fixed (DADA2 requirements, not
# user-editable), maxEE relaxed vs typical 16S defaults since AMF SSU
# amplicons + degenerate primers tend to be noisier.

mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "4. Filter and Trim",
    description = "Removes low-quality reads and trims read ends using DADA2 filterAndTrim, so only reliable sequence is passed on to denoising. Choose truncation lengths from the quality profiles, then run.",
    results = shiny::tagList(
      shiny::uiOutput(ns("input_status")),
      result_card("Quality profiles",
        shiny::fluidRow(
          shiny::column(7,
            shiny::actionButton(ns("preview_quality"), "Preview quality profiles", class = "btn-outline-secondary btn-sm mb-2"),
            shiny::actionButton(ns("choose_samples"), "Select samples...", class = "btn-outline-secondary btn-sm mb-2"),
            shiny::uiOutput(ns("preview_badge"), inline = TRUE)
          ),
          shiny::column(5, shiny::radioButtons(ns("view_mode"), "View", choices = c("Per sample" = "per_sample", "Aggregated" = "aggregate"), selected = "per_sample", inline = TRUE))
        ),
        shiny::uiOutput(ns("selected_samples_caption")),
        shiny::uiOutput(ns("preview_error_alert")),
        shiny::fluidRow(
          shiny::column(6, shiny::h6("Forward"), shiny::plotOutput(ns("quality_plot_fwd"), height = "320px")),
          shiny::column(6, shiny::h6("Reverse"), shiny::plotOutput(ns("quality_plot_rev"), height = "320px"))
        ),
        actions = shiny::uiOutput(ns("downloads"))
      ),
      params_card(
      shiny::fluidRow(
        shiny::column(4,
          param_group("Truncation",
            shiny::fluidRow(
              shiny::column(6, shiny::tagList(
                help_label("truncLen (forward)", "Truncates every forward read to this many bases (cutting off the low-quality tail), before filtering by maxEE. 0 = no truncation. Pick this from where the forward quality profile above starts to drop."),
                shiny::numericInput(ns("truncLen_fwd"), NULL, value = 0, min = 0)
              )),
              shiny::column(6, shiny::tagList(
                help_label("truncLen (reverse)", "Same as truncLen (forward), but for the reverse read. Reverse reads are typically lower quality, so this is often shorter than the forward value."),
                shiny::numericInput(ns("truncLen_rev"), NULL, value = 0, min = 0)
              ))
            ),
            shiny::div(class = "amf-param-hint", "0 means no truncation.")
          )
        ),
        shiny::column(5,
          param_group("Expected errors",
            shiny::fluidRow(
              shiny::column(6, shiny::tagList(
                help_label("maxEE (forward)", "Maximum expected errors allowed in a forward read (summed per-base error probability from quality scores) before it's discarded. Lower = stricter filtering, fewer reads kept."),
                shiny::numericInput(ns("maxEE_fwd"), NULL, value = 2, min = 0, step = 0.5)
              )),
              shiny::column(6, shiny::tagList(
                help_label("maxEE (reverse)", "Same as maxEE (forward), but for the reverse read. Set higher here since reverse reads are noisier and a strict value would discard too many pairs."),
                shiny::numericInput(ns("maxEE_rev"), NULL, value = 4, min = 0, step = 0.5)
              ))
            ),
            shiny::div(class = "amf-param-hint", "Relaxed on the reverse read: AMF SSU with degenerate primers is noisier than typical 16S.")
          )
        ),
        shiny::column(3,
          param_group("Quality trimming",
            help_label("truncQ", "Trims each read at the first base with quality score at or below this value, before truncLen/maxEE are applied. Lower = keeps more of the read (trims less aggressively)."),
            shiny::numericInput(ns("truncQ"), NULL, value = 2, min = 0),
            help_label("minLen", "Discards reads shorter than this many bases after truncation and trimming. Applied after primer removal, so this is the length of the remaining read. 0 = no minimum."),
            shiny::numericInput(ns("minLen"), NULL, value = 50, min = 0),
            shiny::div(class = "amf-param-hint", "Fixed by DADA2 requirements: maxN = 0, rm.phix = TRUE.")
          )
        )
      )
      ),
      result_card("Filter and trim reads",
        run_row(ns, "Filter and trim")
      ),
      shiny::uiOutput(ns("vb_reads")),
      result_card("Read tracking", dt_output(ns("track_table")), results_placeholder("Run this step to see reads in and out for each sample."))
    )
  )
}

# sample_table: reactive() -> data.frame(sample, fwd, rev) of basenames from
# mod_setup (same one mod_primer uses -- QC reads Step 3's trimmed output).
mod_qc_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "qc"

    input_files <- shiny::reactive({
      st <- sample_table()
      shiny::req(rv$project_dir)
      # Depending on rv$status$primer (even though unused) makes this
      # recompute the moment Primer Removal finishes, instead of caching the
      # empty result from before those files existed on disk.
      rv$status$primer
      trimmed_files_df(rv$project_dir, st)
    })

    output$input_status <- shiny::renderUI({
      files <- input_files()
      if (nrow(files) == 0) {
        return(shiny::div(class = "alert alert-warning small", "No trimmed reads found yet -- finish Primer Removal first."))
      }
      shiny::div(class = "text-muted small mb-2", sprintf("%d trimmed sample(s) ready from 01_trimmed/.", nrow(files)))
    })

    # Which samples to plot -- defaults to the first 2 (as before) until the
    # user opens "Select samples..." and picks their own set, of any size.
    selected_samples <- shiny::reactiveVal(character(0))

    shiny::observeEvent(input_files(), {
      files <- input_files()
      if (length(selected_samples()) == 0 && nrow(files) > 0) {
        selected_samples(files$sample[seq_len(min(2, nrow(files)))])
      }
    })

    shiny::observeEvent(input$choose_samples, {
      files <- input_files()
      shiny::req(nrow(files) > 0)
      shiny::showModal(shiny::modalDialog(
        title = "Select samples to plot",
        shiny::checkboxGroupInput(
          ns("sample_picker"), NULL,
          choices = files$sample,
          selected = intersect(selected_samples(), files$sample)
        ),
        footer = shiny::tagList(
          shiny::actionButton(ns("select_all_samples"), "Select all", class = "btn-outline-secondary btn-sm"),
          shiny::actionButton(ns("select_none_samples"), "Select none", class = "btn-outline-secondary btn-sm"),
          shiny::modalButton("Cancel"),
          shiny::actionButton(ns("apply_sample_picker"), "Apply", class = "btn-primary")
        )
      ))
    })

    shiny::observeEvent(input$select_all_samples, {
      shiny::updateCheckboxGroupInput(session, "sample_picker", selected = input_files()$sample)
    })
    shiny::observeEvent(input$select_none_samples, {
      shiny::updateCheckboxGroupInput(session, "sample_picker", selected = character(0))
    })

    shiny::observeEvent(input$apply_sample_picker, {
      selected_samples(input$sample_picker %||% character(0))
      shiny::removeModal()
    })

    output$selected_samples_caption <- shiny::renderUI({
      n <- length(selected_samples())
      shiny::div(
        class = "text-muted small mb-2",
        if (n == 0) "No samples selected -- click \"Select samples...\" to choose which to plot."
        else sprintf("Plotting %d sample(s): %s. Inspect these before choosing truncLen below.", n, paste(selected_samples(), collapse = ", "))
      )
    })

    # Rendered as a background callr job (like Step 5's "Learn error rate")
    # so the RUNNING badge actually shows while plotQualityProfile() is still
    # computing, instead of freezing the UI with no feedback until it's done.
    preview_status <- shiny::reactiveVal("IDLE")
    preview_handle <- shiny::reactiveVal(NULL)
    preview_result <- shiny::reactiveVal(NULL)  # list(fwd=, rev=) ggplot objects
    preview_error <- shiny::reactiveVal(NULL)
    output$preview_error_alert <- shiny::renderUI({
      shiny::req(identical(preview_status(), "ERROR"), preview_error())
      shiny::div(class = "alert alert-danger small mb-2", preview_error())
    })

    shiny::observeEvent(input$preview_quality, {
      shiny::req(!identical(preview_status(), "RUNNING"))
      files <- input_files()
      picked <- files[files$sample %in% selected_samples(), ]
      shiny::req(nrow(picked) > 0)
      aggregate <- identical(input$view_mode, "aggregate")

      preview_job <- function(fwd_files, rev_files, aggregate) {
        library(dada2)
        list(
          fwd = plotQualityProfile(fwd_files, aggregate = aggregate),
          rev = plotQualityProfile(rev_files, aggregate = aggregate)
        )
      }

      h <- callr::r_bg(
        func = preview_job,
        args = list(fwd_files = picked$fwd, rev_files = picked$rev, aggregate = aggregate),
        stdout = "|", stderr = "|", supervise = TRUE
      )
      preview_handle(h)
      preview_status("RUNNING")
      preview_result(NULL)
    })

    shiny::observe({
      h <- preview_handle()
      shiny::req(!is.null(h), identical(preview_status(), "RUNNING"))
      shiny::invalidateLater(500, session)
      if (!h$is_alive()) {
        result <- tryCatch(h$get_result(), error = function(e) e)
        if (inherits(result, "error")) {
          preview_error(conditionMessage(result))
          preview_status("ERROR")
        } else {
          preview_error(NULL)
          preview_result(result)
          preview_status("SUCCESS")
        }
      }
    })

    output$preview_badge <- shiny::renderUI(status_badge(preview_status()))
    output$quality_plot_fwd <- shiny::renderPlot({
      shiny::req(preview_result())
      preview_result()$fwd
    })
    output$quality_plot_rev <- shiny::renderPlot({
      shiny::req(preview_result())
      preview_result()$rev
    })

    # Only offer downloads once profiles actually exist for this session.
    output$downloads <- shiny::renderUI({
      shiny::req(preview_result())
      shiny::div(
        class = "d-flex gap-2 mt-2",
        shiny::downloadButton(ns("dl_fwd"), "Download forward plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm"),
        shiny::downloadButton(ns("dl_rev"), "Download reverse plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm")
      )
    })

    output$dl_fwd <- shiny::downloadHandler(
      filename = function() "quality_profile_forward.png",
      content = function(file) ggplot2::ggsave(file, plot = preview_result()$fwd, dpi = 300, width = 8, height = 6, units = "in")
    )
    output$dl_rev <- shiny::downloadHandler(
      filename = function() "quality_profile_reverse.png",
      content = function(file) ggplot2::ggsave(file, plot = preview_result()$rev, dpi = 300, width = 8, height = 6, units = "in")
    )

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv))
      files <- input_files()
      shiny::req(nrow(files) > 0)

      filt_dir <- filtered_dir_path(rv$project_dir)
      fs::dir_create(filt_dir)
      filt_paths <- filtered_fastq_paths(rv$project_dir, files$sample)

      truncLen <- c(input$truncLen_fwd, input$truncLen_rev)
      maxEE <- c(input$maxEE_fwd, input$maxEE_rev)
      truncQ <- input$truncQ
      minLen <- input$minLen %||% 50
      multithread <- .Platform$OS.type != "windows"

      filter_job <- function(samples, fnFs, fnRs, filtFs, filtRs, truncLen, maxEE, truncQ, minLen, multithread) {
        library(dada2)
        # NB: filterAndTrim's signature is (fwd, filt, rev, filt.rev, ...) --
        # NOT (fwd, rev, filt.fwd, filt.rev). Positional args here previously
        # put fnRs where filt was expected, silently checking existence of
        # the wrong (output) path. Named args avoid that trap for good.
        out <- filterAndTrim(
          fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs,
          truncLen = truncLen, maxEE = maxEE, truncQ = truncQ, minLen = minLen,
          maxN = 0, rm.phix = TRUE, compress = TRUE,
          multithread = multithread, verbose = TRUE
        )
        data.frame(sample = samples, reads_in = out[, 1], reads_out = out[, 2], row.names = NULL)
      }

      launch_step_callr(rv, step_id, filter_job, list(
        samples = files$sample, fnFs = files$fwd, fnRs = files$rev,
        filtFs = filt_paths$fwd, filtRs = filt_paths$rev,
        truncLen = truncLen, maxEE = maxEE, truncQ = truncQ, minLen = minLen, multithread = multithread
      ))
    })

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    output$track_table <- DT::renderDataTable({
      res <- rv$artifacts[[step_id]]
      if (is.null(res)) {
        empty <- data.frame(Sample = character(0), `Reads in` = integer(0), `Reads out` = integer(0),
                             `% kept` = character(0), check.names = FALSE)
        return(DT::datatable(empty, rownames = FALSE, options = list(dom = "t")))
      }
      pct <- ifelse(res$reads_in > 0, sprintf("%.1f%%", 100 * res$reads_out / res$reads_in), "NA")
      df <- data.frame(
        Sample = res$sample, `Reads in` = res$reads_in, `Reads out` = res$reads_out, `% kept` = pct,
        check.names = FALSE
      )
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = nrow(df) + 1))
    })

    # Read-only summary tiles (presentation only; reads rv$artifacts).
    output$vb_reads <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res), is.data.frame(res), all(c("reads_in", "reads_out") %in% names(res)))
      rin <- sum(res$reads_in, na.rm = TRUE)
      rout <- sum(res$reads_out, na.rm = TRUE)
      stat_boxes(list(
        "Reads in" = format(rin, big.mark = ","),
        "Reads kept" = format(rout, big.mark = ","),
        "Retained" = if (rin > 0) sprintf("%.1f%%", 100 * rout / rin) else "NA"
      ))
    })
  })
}
