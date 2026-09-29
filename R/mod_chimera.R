# Step 8 -- Remove Chimeras (DADA2 removeBimeraDenovo), per DESIGN.md
# section 5. Reads Step 7's seqtab.rds; writes seqtab_nochim.rds and the
# full read-tracking table (raw -> non-chimeric), counting the reads each
# earlier step left on disk and pulling the denoised counts back from Step
# 6's dadaF.rds/dadaR.rds. The final ASV
# length histogram here is the last look before Taxonomy -- a tight band
# confirms merge/chimera-removal didn't distort the target region.

mod_chimera_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "8. Remove Chimeras",
    description = "Removes chimeric sequences, artefacts formed when two different templates join during PCR, using DADA2's de novo method. What remains is your final set of ASVs.",
    params = shiny::tagList(
      shiny::uiOutput(ns("input_status")),
      shiny::fluidRow(
        shiny::column(6,
          param_group("Method",
            help_label("Chimera detection method", "Consensus (default): each sequence is tested in every sample separately, and removed only if it is flagged as chimeric across the samples it occurs in. Pooled: all samples are pooled and each sequence is tested once. Per-sample: each sample is tested separately and a sequence is removed from that sample only where it is flagged."),
            shiny::selectInput(
              ns("method"), NULL,
              choices = c("Consensus (default)" = "consensus", "Pooled" = "pooled", "Per-sample" = "per-sample"),
              width = "100%"
            )
          )
        ),
        shiny::column(6, step_io("Merged sequence table from step 7.", "DADA2 removeBimeraDenovo (de novo, no reference database).", "seqtab_nochim.rds, track.rds"))
      )
    ),
    results = shiny::tagList(
      result_card("Remove chimeras",
        run_row(ns, "Remove chimeras"),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("log"))
      ),
      shiny::uiOutput(ns("vb_asvs")),
      shiny::uiOutput(ns("results_section")),
      results_placeholder("Run this step to see the final ASV count, length distribution and read tracking.")
    )
  )
}

mod_chimera_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "chimera"

    seqtab_ready <- shiny::reactive({
      shiny::req(rv$project_dir)
      # Depending on rv$status$merge (even though unused) makes this
      # recompute the moment Merge Pairs finishes, instead of caching a
      # FALSE result from before seqtab.rds existed on disk.
      rv$status$merge
      fs::file_exists(file.path(denoise_dir_path(rv$project_dir), "seqtab.rds"))
    })

    output$input_status <- shiny::renderUI({
      if (!isTRUE(seqtab_ready())) {
        return(shiny::div(class = "alert alert-warning small", "No merged sequence table found yet -- finish Merge Pairs first."))
      }
      shiny::div(class = "text-muted small mb-2", "Sequence table ready -- checking for chimeras compares each ASV against more abundant ones from the same samples.")
    })

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv))
      shiny::req(isTRUE(seqtab_ready()))
      st <- sample_table()
      shiny::req(nrow(st) > 0)

      out_dir <- denoise_dir_path(rv$project_dir)
      multithread <- .Platform$OS.type != "windows"

      # Forward-read files of each earlier step; pairs are filtered together,
      # so the R1 count is the pair count.
      stage_files <- list(
        raw = file.path(rv$project_dir, st$fwd),
        filtN = filtn_fastq_paths(rv$project_dir, st$sample)$fwd,
        trimmed = trimmed_fastq_paths(rv$project_dir, st$sample)$fwd,
        filtered = filtered_fastq_paths(rv$project_dir, st$sample)$fwd
      )

      chimera_job <- function(samples, out_dir, method, multithread, stage_files) {
        library(dada2)
        seqtab <- readRDS(file.path(out_dir, "seqtab.rds"))
        dadaFs <- readRDS(file.path(out_dir, "dadaF.rds"))
        dadaRs <- readRDS(file.path(out_dir, "dadaR.rds"))

        seqtab_nochim <- removeBimeraDenovo(seqtab, method = method, multithread = multithread, verbose = TRUE)
        saveRDS(seqtab_nochim, file.path(out_dir, "seqtab_nochim.rds"))

        getN <- function(x) sum(dada2::getUniques(x))
        denoisedF <- sapply(dadaFs, getN)
        denoisedR <- sapply(dadaRs, getN)
        merged <- rowSums(seqtab)
        nonchim <- rowSums(seqtab_nochim)
        # NA where a step wrote no file (filterAndTrim skips samples with no reads left).
        count_reads <- function(paths) {
          n <- rep(NA_real_, length(paths))
          ok <- file.exists(paths)
          if (any(ok)) n[ok] <- ShortRead::countFastq(paths[ok])$records
          n
        }

        track <- data.frame(
          sample = samples,
          raw = count_reads(stage_files$raw),
          filtN = count_reads(stage_files$filtN),
          trimmed = count_reads(stage_files$trimmed),
          filtered = count_reads(stage_files$filtered),
          denoisedF = as.numeric(denoisedF[samples]),
          denoisedR = as.numeric(denoisedR[samples]),
          merged = as.numeric(merged[samples]),
          nonchim = as.numeric(nonchim[samples]),
          row.names = NULL
        )
        saveRDS(track, file.path(out_dir, "track.rds"))

        list(
          track = track, asv_lengths = nchar(colnames(seqtab_nochim)),
          n_asvs = ncol(seqtab_nochim), n_asvs_before = ncol(seqtab)
        )
      }

      launch_step_callr(rv, step_id, chimera_job, list(samples = st$sample, out_dir = out_dir, method = input$method %||% "consensus", multithread = multithread, stage_files = stage_files))
    })

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    output$results_section <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res))
      result_card("Chimera removal results",
        shiny::div(class = "mb-2", sprintf("%d ASVs before chimera removal -> %d after (%d removed).",
                                            res$n_asvs_before, res$n_asvs, res$n_asvs_before - res$n_asvs)),
        shiny::h6(sprintf("ASV length distribution (%d non-chimeric ASVs)", res$n_asvs)),
        shiny::div(class = "text-muted small mb-1", "Expect a tight band around the target SSU region; a scattered/multi-modal spread signals primer, merge-mode, or truncLen problems."),
        shiny::plotOutput(ns("length_hist"), height = "260px"),
        shiny::h6("Reads tracked (raw -> non-chimeric)"),
        dt_output(ns("track_table")),
        shiny::downloadButton(ns("download_track"), "Download as CSV", class = "btn-outline-secondary btn-sm mb-2")
      )
    })

    output$length_hist <- shiny::renderPlot({
      res <- rv$artifacts[[step_id]]
      shiny::req(res)
      counts <- table(factor(res$asv_lengths, levels = seq(min(res$asv_lengths) - 2, max(res$asv_lengths) + 2)))
      graphics::barplot(counts, space = 0.1, main = NULL, xlab = "ASV length (bp)", ylab = "Number of ASVs",
                        col = "#0f766e", border = "white")
    })

    output$track_table <- DT::renderDataTable({
      res <- rv$artifacts[[step_id]]
      if (is.null(res)) {
        empty <- data.frame(Sample = character(0), Raw = integer(0), `N-filtered` = integer(0), `Primer-trimmed` = integer(0),
                             Filtered = integer(0), `Denoised (F)` = integer(0), `Denoised (R)` = integer(0),
                             Merged = integer(0), `Non-chimeric` = integer(0), check.names = FALSE)
        return(DT::datatable(empty, rownames = FALSE, options = list(dom = "t")))
      }
      t <- res$track
      df <- data.frame(
        Sample = t$sample, Raw = t$raw, `N-filtered` = t$filtN, `Primer-trimmed` = t$trimmed,
        Filtered = t$filtered, `Denoised (F)` = t$denoisedF, `Denoised (R)` = t$denoisedR,
        Merged = t$merged, `Non-chimeric` = t$nonchim, check.names = FALSE
      )
      DT::datatable(df, rownames = FALSE, options = list(dom = "t", pageLength = nrow(df) + 1))
    })

    output$download_track <- shiny::downloadHandler(
      filename = function() sprintf("chimera_read_tracking_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        res <- rv$artifacts[[step_id]]
        shiny::req(res)
        utils::write.csv(res$track, file, row.names = FALSE)
      }
    )

    # Read-only summary tiles (presentation only; reads rv$artifacts).
    output$vb_asvs <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res), !is.null(res$n_asvs), !is.null(res$n_asvs_before))
      removed <- res$n_asvs_before - res$n_asvs
      stat_boxes(list(
        "ASVs before" = format(res$n_asvs_before, big.mark = ","),
        "Chimeras removed" = format(removed, big.mark = ","),
        "Final ASVs" = format(res$n_asvs, big.mark = ",")
      ))
    })

    mod_logpanel_server("log", shiny::reactive(rv$log[[step_id]] %||% character(0)))
  })
}
