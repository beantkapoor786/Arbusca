# Step 3 -- Remove Ambiguous Bases (dada2::filterAndTrim with maxN = 0).
# Drops every read pair containing an N before primer removal, so the
# primer counts in Step 4 (vcountPattern with IUPAC matching, where an N in
# a read would match any primer base) aren't inflated. Writes 01_filtN/.

mod_filtn_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "3. Remove Ambiguous Bases",
    description = "Removes every read pair that contains an ambiguous base (N), so primers can be detected and removed reliably in the next step.",
    params = shiny::tagList(
      shiny::uiOutput(ns("input_status")),
      step_io("Raw forward and reverse reads.", "DADA2 filterAndTrim with maxN = 0 (other settings at their defaults).", "01_filtN/")
    ),
    results = shiny::tagList(
      result_card("Remove reads with Ns",
        run_row(ns, "Remove Ns"),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("log"))
      ),
      shiny::uiOutput(ns("vb_reads")),
      result_card("Read counts", dt_output(ns("track_table")), results_placeholder("Run this step to see raw and filtered read counts for each sample."))
    )
  )
}

mod_filtn_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "filtn"

    input_files <- shiny::reactive({
      shiny::req(rv$project_dir)
      raw_files_df(rv$project_dir, sample_table())
    })

    output$input_status <- shiny::renderUI({
      files <- input_files()
      if (nrow(files) == 0) {
        return(shiny::div(class = "alert alert-warning small", "No paired samples detected yet -- finish Setup first."))
      }
      shiny::div(class = "text-muted small mb-2", sprintf("%d raw sample(s) ready.", nrow(files)))
    })

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv))
      files <- input_files()
      shiny::req(nrow(files) > 0)

      fs::dir_create(filtn_dir_path(rv$project_dir))
      out_paths <- filtn_fastq_paths(rv$project_dir, files$sample)

      filtn_job <- function(samples, fnFs, fnRs, filtFs, filtRs, multithread) {
        library(dada2)
        out <- filterAndTrim(
          fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs,
          maxN = 0, multithread = multithread, verbose = TRUE
        )
        data.frame(sample = samples, reads_in = out[, 1], reads_out = out[, 2], row.names = NULL)
      }

      launch_step_callr(rv, step_id, filtn_job, list(
        samples = files$sample, fnFs = files$fwd, fnRs = files$rev,
        filtFs = out_paths$fwd, filtRs = out_paths$rev,
        multithread = .Platform$OS.type != "windows"
      ))
    })

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    output$track_table <- DT::renderDataTable({
      res <- rv$artifacts[[step_id]]
      if (is.null(res)) {
        empty <- data.frame(Sample = character(0), `Raw reads` = integer(0), `Filtered reads` = integer(0),
                             `% kept` = character(0), check.names = FALSE)
        return(DT::datatable(empty, rownames = FALSE, options = list(dom = "t")))
      }
      pct <- ifelse(res$reads_in > 0, sprintf("%.1f%%", 100 * res$reads_out / res$reads_in), "NA")
      df <- data.frame(
        Sample = res$sample, `Raw reads` = res$reads_in, `Filtered reads` = res$reads_out, `% kept` = pct,
        check.names = FALSE
      )
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = nrow(df) + 1))
    })

    output$vb_reads <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res), is.data.frame(res))
      rin <- sum(res$reads_in, na.rm = TRUE)
      rout <- sum(res$reads_out, na.rm = TRUE)
      stat_boxes(list(
        "Raw reads" = format(rin, big.mark = ","),
        "Reads without Ns" = format(rout, big.mark = ","),
        "Retained" = if (rin > 0) sprintf("%.1f%%", 100 * rout / rin) else "NA"
      ))
    })

    mod_logpanel_server("log", shiny::reactive(rv$log[[step_id]] %||% character(0)))
  })
}
