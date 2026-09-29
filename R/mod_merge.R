# Step 7 -- Merge Pairs (DADA2 mergePairs -> makeSequenceTable), per
# DESIGN.md section 5's "Merge mode (AML1/AML2-aware)": overlap merge
# (default, short amplicons), forward-only (drop R2 -- amplicon too long to
# overlap, e.g. nested AML1/AML2), or concatenate (join without overlap).
# Reads Step 6's dadaF.rds/dadaR.rds; writes seqtab.rds. The ASV length
# histogram here is the first QC look at the assembled (pre-chimera) ASVs.

mod_merge_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "7. Merge Pairs",
    description = "Joins each forward read with its reverse partner into one full-length sequence, then builds the table of sequences per sample. Choose the mode that matches whether your reads overlap.",
    params = shiny::tagList(
      shiny::uiOutput(ns("input_status")),
      shiny::fluidRow(
        shiny::column(6,
          param_group("Merge mode",
            help_label("Merge mode", "Overlap merge requires forward and reverse reads to overlap and agree in the overlapping region. Forward-only discards the reverse read. Concatenate joins the reads end to end with no overlap required."),
            shiny::selectInput(
              ns("merge_mode"), NULL,
              choices = c(
                "Overlap merge (default, short amplicons)" = "overlap",
                "Forward-only (amplicon too long to overlap)" = "forward_only",
                "Concatenate (join without overlap)" = "concatenate"
              ),
              width = "100%"
            ),
            shiny::div(
              class = "amf-param-hint",
              "Forward-only or Concatenate is safer for long/nested amplicons (e.g. AML1/AML2) that can't overlap on short reads."
            )
          )
        ),
        shiny::column(6, step_io("Denoised forward and reverse reads from step 6.", "DADA2 mergePairs, then makeSequenceTable.", "seqtab.rds"))
      )
    ),
    results = shiny::tagList(
      result_card("Merge pairs",
        run_row(ns, "Merge pairs"),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("log"))
      ),
      shiny::uiOutput(ns("vb_asvs")),
      shiny::uiOutput(ns("results_section")),
      results_placeholder("Run this step to see the ASV length distribution and merged read counts.")
    )
  )
}

mod_merge_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "merge"

    input_files <- shiny::reactive({
      st <- sample_table()
      shiny::req(rv$project_dir)
      rv$status$qc
      filtered_files_df(rv$project_dir, st)
    })

    dada_ready <- shiny::reactive({
      shiny::req(rv$project_dir)
      # Depending on rv$status$denoise (even though unused) makes this
      # recompute the moment Denoise finishes, instead of caching a FALSE
      # result from before dadaF.rds/dadaR.rds existed on disk.
      rv$status$denoise
      out_dir <- denoise_dir_path(rv$project_dir)
      fs::file_exists(file.path(out_dir, "dadaF.rds")) && fs::file_exists(file.path(out_dir, "dadaR.rds"))
    })

    output$input_status <- shiny::renderUI({
      if (!isTRUE(dada_ready())) {
        return(shiny::div(class = "alert alert-warning small", "No denoised reads found yet -- finish Learn Error Rate and Denoise first."))
      }
      files <- input_files()
      shiny::div(class = "text-muted small mb-2", sprintf("%d sample(s) ready to merge.", nrow(files)))
    })

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv))
      shiny::req(isTRUE(dada_ready()))
      files <- input_files()
      shiny::req(nrow(files) > 0)

      out_dir <- denoise_dir_path(rv$project_dir)

      merge_job <- function(samples, filtFs, filtRs, out_dir, merge_mode) {
        library(dada2)
        names(filtFs) <- samples
        names(filtRs) <- samples

        dadaFs <- readRDS(file.path(out_dir, "dadaF.rds"))
        dadaRs <- readRDS(file.path(out_dir, "dadaR.rds"))

        if (identical(merge_mode, "forward_only")) {
          mergers <- NULL
          seqtab <- makeSequenceTable(dadaFs)
        } else if (identical(merge_mode, "concatenate")) {
          mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, justConcatenate = TRUE, verbose = TRUE)
          seqtab <- makeSequenceTable(mergers)
        } else {
          mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
          seqtab <- makeSequenceTable(mergers)
        }

        saveRDS(seqtab, file.path(out_dir, "seqtab.rds"))

        getN <- function(x) sum(dada2::getUniques(x))
        merged <- if (is.null(mergers)) {
          sapply(dadaFs, getN)
        } else if (is.data.frame(mergers)) {
          getN(mergers)
        } else {
          sapply(mergers, getN)
        }
        merged_df <- data.frame(sample = samples, merged = as.numeric(merged[samples]), row.names = NULL)

        list(merged = merged_df, asv_lengths = nchar(colnames(seqtab)), n_asvs = ncol(seqtab))
      }

      launch_step_callr(rv, step_id, merge_job, list(
        samples = files$sample, filtFs = files$fwd, filtRs = files$rev,
        out_dir = out_dir, merge_mode = input$merge_mode
      ))
    })

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    output$results_section <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res))
      result_card("Merge results",
        shiny::h6(sprintf("ASV length distribution (%d ASVs, pre-chimera-removal)", res$n_asvs)),
        shiny::div(class = "text-muted small mb-1", "Expect a tight band around the target SSU region; a scattered/multi-modal spread signals primer, merge-mode, or truncLen problems."),
        shiny::plotOutput(ns("length_hist"), height = "260px"),
        shiny::h6("Merged read counts"),
        dt_output(ns("merged_table"))
      )
    })

    output$length_hist <- shiny::renderPlot({
      res <- rv$artifacts[[step_id]]
      shiny::req(res)
      counts <- table(factor(res$asv_lengths, levels = seq(min(res$asv_lengths) - 2, max(res$asv_lengths) + 2)))
      graphics::barplot(counts, space = 0.1, main = NULL, xlab = "ASV length (bp)", ylab = "Number of ASVs",
                        col = "#0f766e", border = "white")
    })

    output$merged_table <- DT::renderDataTable({
      res <- rv$artifacts[[step_id]]
      if (is.null(res)) {
        empty <- data.frame(Sample = character(0), Merged = integer(0), check.names = FALSE)
        return(DT::datatable(empty, rownames = FALSE, options = list(dom = "t")))
      }
      df <- data.frame(Sample = res$merged$sample, Merged = res$merged$merged, check.names = FALSE)
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = nrow(df) + 1))
    })

    # Read-only summary tile (presentation only; reads rv$artifacts).
    output$vb_asvs <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res), !is.null(res$n_asvs))
      stat_boxes(list(
        "ASVs after merging" = format(res$n_asvs, big.mark = ","),
        "Samples" = format(nrow(res$merged), big.mark = ",")
      ))
    })

    mod_logpanel_server("log", shiny::reactive(rv$log[[step_id]] %||% character(0)))
  })
}
