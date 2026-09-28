# Step 9 -- Phyloseq Object. Combines Step 8's ASV table (seqtab_nochim.rds)
# and taxonomy assignments (vt_assignments.rds) with user-supplied sample
# metadata into a single phyloseq object, saved to 05_output/phyloseq.rds.
# Also offers a data transformation on top of that object (rarefaction,
# relative abundance, CLR, or none), saved to 05_output/phyloseq_transformed.rds.

# method: "none" (raw counts), "relabund" (proportions), "clr" (centered
# log-ratio, with a pseudocount since raw ASV counts always contain zeros),
# or "rarefy" (subsample every sample to rarefy_depth reads, dropping any
# sample below it -- phyloseq::rarefy_even_depth()'s own behavior).
apply_transformation <- function(ps, method, rarefy_depth = NULL, rngseed = 42, clr_pseudocount = 1) {
  switch(method,
    none = ps,
    relabund = phyloseq::transform_sample_counts(ps, function(x) x / sum(x)),
    clr = {
      otu <- phyloseq::otu_table(ps)
      # Normalize to a samples-x-taxa matrix regardless of the object's own
      # orientation, so the row-wise (per-sample) centering below is correct.
      mat <- if (phyloseq::taxa_are_rows(ps)) t(as(otu, "matrix")) else as(otu, "matrix")
      logm <- log(mat + clr_pseudocount)
      clr_mat <- logm - rowMeans(logm)
      phyloseq::otu_table(ps) <- phyloseq::otu_table(clr_mat, taxa_are_rows = FALSE)
      ps
    },
    rarefy = phyloseq::rarefy_even_depth(
      ps, sample.size = rarefy_depth, rngseed = rngseed,
      replace = FALSE, trimOTUs = TRUE, verbose = FALSE
    ),
    stop("Unknown transformation method: ", method)
  )
}

TRANSFORM_METHOD_LABELS <- c(
  rarefy   = "Rarefaction (subsample to even depth)",
  clr      = "CLR (centered log ratio)",
  relabund = "Relative abundance (proportion)",
  none     = "No transformation (raw counts)"
)

TRANSFORM_METHOD_HELP <- c(
  rarefy   = "Randomly subsamples every sample down to the same number of reads (dropping samples below that depth), so differences in sequencing depth can't drive apparent differences in diversity. Involves discarding data and is sensitive to the random seed -- CLR is often preferred for compositional comparisons.",
  clr      = "Centered log-ratio: log-transforms each sample's proportions relative to that sample's own geometric mean, correcting the compositional (values-sum-to-a-constant) nature of ASV count data. A small pseudocount is added first since raw counts contain zeros, which log() can't handle.",
  relabund = "Divides each sample's counts by that sample's total, so every sample sums to 1. Removes the effect of uneven sequencing depth, but the values are compositional (not independent) and comparisons of absolute abundance are lost.",
  none     = "Keeps the raw ASV read counts as-is -- no transformation applied."
)

mod_output_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "9. Phyloseq Object",
    description = "Combines the ASV table, taxonomy and your sample metadata into a single phyloseq object, the standard R format for downstream community analysis.",
    params = shiny::tagList(
      shiny::uiOutput(ns("phyloseq_gate_alert")),
      shiny::uiOutput(ns("input_status")),
      param_group("Sample metadata",
        shiny::div(
          class = "amf-param-hint mb-2",
          "Upload a CSV with one row per sample. One of its columns must contain the same sample names used throughout this pipeline (Setup's sample sheet) -- pick which column below."
        ),
        shiny::fluidRow(
          shiny::column(6, shiny::fileInput(ns("metadata_file"), "Metadata CSV", accept = ".csv", width = "100%")),
          shiny::column(6, shiny::uiOutput(ns("sample_col_picker")))
        )
      ),
      shiny::div(
        class = "amf-run-row",
        shiny::actionButton(ns("run"), "Build phyloseq object", icon = bsicons::bs_icon("play-fill"), class = "btn-primary"),
        shiny::uiOutput(ns("badge"), inline = TRUE)
      )
    ),
    results = shiny::tagList(
      shiny::uiOutput(ns("metadata_preview_section")),
      shiny::uiOutput(ns("build_error_alert")),
      shiny::uiOutput(ns("results_section")),
      shiny::uiOutput(ns("transform_section")),
      results_placeholder("Upload metadata and build the phyloseq object to see a summary here.")
    )
  )
}

# sample_table: reactive() -> data.frame(sample, fwd, rev) of basenames from
# mod_setup -- used here only for its list of pipeline sample names, to
# validate the uploaded metadata against.
mod_output_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "output"

    seqtab_path <- shiny::reactive(file.path(denoise_dir_path(rv$project_dir), "seqtab_nochim.rds"))
    assignments_path <- shiny::reactive(file.path(taxonomy_dir_path(rv$project_dir), "vt_assignments.rds"))

    inputs_ready <- shiny::reactive({
      shiny::req(rv$project_dir)
      # Depending on rv$status$taxonomy (even though unused) makes this
      # recompute the moment Taxonomy Assignment finishes, instead of caching
      # a FALSE result from before these files existed on disk.
      rv$status$taxonomy
      fs::file_exists(seqtab_path()) && fs::file_exists(assignments_path())
    })

    output$input_status <- shiny::renderUI({
      if (!isTRUE(inputs_ready())) {
        return(shiny::div(class = "alert alert-warning small", "No taxonomy assignments found yet -- finish Taxonomy Assignment first."))
      }
      NULL
    })

    phyloseq_ready <- shiny::reactive(requireNamespace("phyloseq", quietly = TRUE))
    output$phyloseq_gate_alert <- shiny::renderUI({
      if (isTRUE(phyloseq_ready())) return(NULL)
      shiny::div(class = "alert alert-warning small mb-3",
                 "phyloseq package not found -- install via BiocManager::install(\"phyloseq\") to run this step.")
    })

    metadata_df <- shiny::reactive({
      shiny::req(input$metadata_file)
      utils::read.csv(input$metadata_file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    })

    # Best-guess which column holds sample names: whichever column's values
    # overlap the pipeline's sample list the most.
    output$sample_col_picker <- shiny::renderUI({
      df <- metadata_df()
      st <- sample_table()
      shiny::req(nrow(df) > 0, nrow(st) > 0)
      cols <- names(df)
      overlaps <- vapply(cols, function(cl) sum(as.character(df[[cl]]) %in% st$sample), integer(1))
      best <- cols[which.max(overlaps)]
      shiny::selectInput(ns("sample_col"), "Column containing sample names", choices = cols, selected = best, width = "100%")
    })

    output$metadata_preview_section <- shiny::renderUI({
      df <- metadata_df()
      st <- sample_table()
      shiny::req(nrow(df) > 0, input$sample_col)
      meta_samples <- as.character(df[[input$sample_col]])
      matched <- sum(meta_samples %in% st$sample)
      missing <- setdiff(st$sample, meta_samples)
      result_card("Sample metadata",
        shiny::div(
          class = "text-muted small mb-1",
          sprintf("%d / %d metadata rows match a pipeline sample.", matched, length(meta_samples))
        ),
        if (length(missing) > 0) {
          shiny::div(class = "alert alert-warning small mb-2",
                     sprintf("%d pipeline sample(s) have no metadata row and will be dropped from the phyloseq object: %s.",
                             length(missing), paste(missing, collapse = ", ")))
        },
        dt_output(ns("metadata_table"))
      )
    })

    output$metadata_table <- DT::renderDataTable({
      df <- metadata_df()
      shiny::req(nrow(df) > 0)
      DT::datatable(df, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 10))
    })

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    build_error <- shiny::reactiveVal(NULL)
    output$build_error_alert <- shiny::renderUI({
      shiny::req(identical(rv$status[[step_id]], "ERROR"), build_error())
      shiny::div(class = "alert alert-danger small mb-3", build_error())
    })

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv))
      msg <- if (!isTRUE(phyloseq_ready())) "The phyloseq package is not installed."
        else if (!isTRUE(inputs_ready())) "No taxonomy assignments found yet -- finish Taxonomy Assignment first."
        else if (is.null(input$metadata_file)) "Upload a sample metadata CSV first."
      if (!is.null(msg)) return(shiny::showNotification(msg, type = "error"))
      shiny::req(metadata_df(), input$sample_col)
      st <- sample_table()

      rv$status[[step_id]] <- "RUNNING"

      # The build runs synchronously here, so doing it in this observer would
      # set RUNNING and SUCCESS within one flush and the badge would never show
      # RUNNING. Run it once this flush has reached the browser instead.
      session$onFlushed(function() shiny::isolate({
        result <- tryCatch({
          seqtab_nochim <- readRDS(seqtab_path())
          assignments <- readRDS(assignments_path())

          # assignments$asv_id (ASV1, ASV2, ...) was assigned in colnames(seqtab_nochim)
          # order by write_asv_fasta() back in Step 8, so this rename is positional,
          # not a text match -- it just gives both tables the same taxon names.
          colnames(seqtab_nochim) <- assignments$asv_id

          meta <- metadata_df()
          rownames(meta) <- as.character(meta[[input$sample_col]])
          common <- intersect(rownames(seqtab_nochim), rownames(meta))
          if (length(common) == 0) {
            stop("No sample names in the metadata match the pipeline's sample names -- check the selected column.")
          }

          otu <- phyloseq::otu_table(seqtab_nochim[common, , drop = FALSE], taxa_are_rows = FALSE)

          tax_cols <- intersect(c("family", "genus", "vt", "lineage", "classification"), names(assignments))
          tax_mat <- as.matrix(assignments[, tax_cols, drop = FALSE])
          rownames(tax_mat) <- assignments$asv_id
          tax <- phyloseq::tax_table(tax_mat)

          samp <- phyloseq::sample_data(meta[common, , drop = FALSE])

          ps <- phyloseq::phyloseq(otu, tax, samp)

          out_dir <- file.path(rv$project_dir, "05_output")
          fs::dir_create(out_dir)
          saveRDS(ps, file.path(out_dir, "phyloseq.rds"))

          list(
            phyloseq = ps,
            n_samples = phyloseq::nsamples(ps),
            n_taxa = phyloseq::ntaxa(ps),
            unmatched_samples = setdiff(st$sample, common)
          )
        }, error = function(e) e)

        if (inherits(result, "error")) {
          rv$status[[step_id]] <- "ERROR"
          build_error(conditionMessage(result))
        } else {
          build_error(NULL)
          rv$artifacts[[step_id]] <- result
          rv$status[[step_id]] <- "SUCCESS"
          invalidate_analyses(rv)
        }
      }), once = TRUE)
    })

    output$results_section <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res))
      result_card("Phyloseq object",
        shiny::div(class = "mb-2", sprintf("Phyloseq object built: %d samples x %d taxa.", res$n_samples, res$n_taxa)),
        if (length(res$unmatched_samples) > 0) {
          shiny::div(class = "alert alert-warning small mb-2",
                     sprintf("%d pipeline sample(s) excluded (no metadata row): %s.",
                             length(res$unmatched_samples), paste(res$unmatched_samples, collapse = ", ")))
        },
        shiny::div(class = "text-muted small mb-1", sprintf("Saved to %s.", file.path(rv$project_dir, "05_output", "phyloseq.rds"))),
        shiny::verbatimTextOutput(ns("phyloseq_summary"))
      )
    })

    output$phyloseq_summary <- shiny::renderPrint({
      res <- rv$artifacts[[step_id]]
      shiny::req(res)
      print(res$phyloseq)
    })

    # --- Data transformation -- own local status, not tied to rv$status$output
    # (which just means "the base phyloseq object exists"); re-applying a
    # transformation doesn't need to invalidate anything downstream since
    # this is the last step in the pipeline. ---

    transform_status <- shiny::reactiveVal("IDLE")
    transform_result <- shiny::reactiveVal(NULL)  # list(phyloseq=, method=, table=)

    # Reset whenever the base phyloseq object is rebuilt -- a transformation
    # of the old object no longer applies.
    shiny::observeEvent(rv$artifacts[[step_id]], {
      transform_status("IDLE")
      transform_result(NULL)
    }, ignoreInit = TRUE)

    output$transform_section <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res))
      result_card("Data transformation",
        shiny::div(class = "text-muted small mb-2", "Transform the ASV counts in the phyloseq object above for downstream analysis (diversity, ordination, etc.)."),
        shiny::fluidRow(
          shiny::column(4, shiny::selectInput(
            ns("transform_method"), "Method",
            choices = stats::setNames(names(TRANSFORM_METHOD_LABELS), TRANSFORM_METHOD_LABELS),
            width = "100%"
          )),
          shiny::column(4, shiny::uiOutput(ns("transform_param_ui")))
        ),
        shiny::uiOutput(ns("transform_method_caption")),
        shiny::div(
          class = "d-flex align-items-center gap-2 mb-3",
          shiny::actionButton(ns("apply_transform"), "Apply transformation", icon = bsicons::bs_icon("play-fill"), class = "btn-primary"),
          shiny::uiOutput(ns("transform_badge"), inline = TRUE)
        ),
        shiny::uiOutput(ns("transform_error_alert")),
        shiny::uiOutput(ns("transform_results_section"))
      )
    })

    output$transform_method_caption <- shiny::renderUI({
      shiny::req(input$transform_method)
      shiny::div(class = "text-muted small mb-2", TRANSFORM_METHOD_HELP[[input$transform_method]])
    })

    output$transform_param_ui <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(res, input$transform_method)
      if (identical(input$transform_method, "rarefy")) {
        depths <- phyloseq::sample_sums(res$phyloseq)
        shiny::tagList(
          shiny::numericInput(ns("rarefy_depth"), "Sample size (reads/sample)", value = as.integer(min(depths)), min = 1, width = "100%"),
          shiny::div(class = "text-muted small", sprintf("Sample depths currently range %d - %d reads.", as.integer(min(depths)), as.integer(max(depths))))
        )
      } else if (identical(input$transform_method, "clr")) {
        shiny::numericInput(ns("clr_pseudocount"), "Pseudocount", value = 1, min = 0, step = 0.5, width = "100%")
      } else {
        NULL
      }
    })

    output$transform_badge <- shiny::renderUI(status_badge(transform_status()))

    transform_error <- shiny::reactiveVal(NULL)
    output$transform_error_alert <- shiny::renderUI({
      shiny::req(identical(transform_status(), "ERROR"), transform_error())
      shiny::div(class = "alert alert-danger small mb-3", transform_error())
    })

    shiny::observeEvent(input$apply_transform, {
      shiny::req(!any_running(rv))
      res <- rv$artifacts[[step_id]]
      shiny::req(res, input$transform_method)

      transform_status("RUNNING")

      out <- tryCatch({
        ps2 <- apply_transformation(
          res$phyloseq, input$transform_method,
          rarefy_depth = input$rarefy_depth, clr_pseudocount = input$clr_pseudocount %||% 1
        )
        # Stashed as an attribute (not just implied by the filename) so the
        # Figures step can label which transform is behind
        # phyloseq_transformed.rds without having to re-derive or guess it.
        attr(ps2, "amf_transform_method") <- input$transform_method

        out_dir <- file.path(rv$project_dir, "05_output")
        fs::dir_create(out_dir)
        saveRDS(ps2, file.path(out_dir, "phyloseq_transformed.rds"))

        otu <- phyloseq::otu_table(ps2)
        mat <- if (phyloseq::taxa_are_rows(ps2)) t(as(otu, "matrix")) else as(otu, "matrix")
        table_df <- data.frame(Sample = rownames(mat), mat, check.names = FALSE, row.names = NULL)
        utils::write.csv(table_df, file.path(out_dir, "asv_table_transformed.csv"), row.names = FALSE)

        list(phyloseq = ps2, method = input$transform_method, table = table_df)
      }, error = function(e) e)

      if (inherits(out, "error")) {
        transform_status("ERROR")
        transform_error(conditionMessage(out))
      } else {
        transform_error(NULL)
        transform_result(out)
        transform_status("SUCCESS")
        invalidate_analyses(rv)
      }
    })

    output$transform_results_section <- shiny::renderUI({
      res <- transform_result()
      shiny::req(!is.null(res))
      shiny::tagList(
        shiny::div(class = "mb-2", sprintf(
          "%s applied: %d samples x %d taxa.", TRANSFORM_METHOD_LABELS[[res$method]],
          phyloseq::nsamples(res$phyloseq), phyloseq::ntaxa(res$phyloseq)
        )),
        shiny::div(class = "text-muted small mb-1", sprintf("Saved to %s.", file.path(rv$project_dir, "05_output", "phyloseq_transformed.rds"))),
        shiny::h6("Transformed ASV table"),
        dt_output(ns("transform_table")),
        shiny::downloadButton(ns("download_transform_table"), "Download as CSV", class = "btn-outline-secondary btn-sm mb-2")
      )
    })

    output$transform_table <- DT::renderDataTable({
      res <- transform_result()
      shiny::req(res)
      num_cols <- which(vapply(res$table, is.numeric, logical(1)))
      DT::formatSignif(DT::datatable(res$table, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15)),
                       columns = num_cols, digits = 4)
    })

    output$download_transform_table <- shiny::downloadHandler(
      filename = function() sprintf("asv_table_%s_%s.csv", transform_result()$method, format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        utils::write.csv(transform_result()$table, file, row.names = FALSE)
      }
    )
  })
}
