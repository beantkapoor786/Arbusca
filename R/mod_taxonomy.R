# Step 9 -- Taxonomy Assignment (BLAST + MaarjAM), per DESIGN.md section 5.
# Defaults match the user's stated methodology: blastn vs a local MaarjAM
# BLAST+ database, ASVs with >=97% identity, >=95% query coverage, and
# e-value < 1e-50 are called AM fungal ASVs; everything else is unclassified.
#
# The design doc leaves the canonical MaarjAM download source as an open
# question (no verified official URL to hardcode), so this does not attempt
# to auto-download the reference -- the user either points at an existing
# BLAST+ database or supplies a reference FASTA that the app indexes with
# makeblastdb. Both choices persist via app_config.R like the cutadapt path.

taxonomy_defaults <- list(identity = 97, coverage = 95, evalue = "1e-50", max_target_seqs = 5)

validate_blast_db <- function(db_path) {
  tryCatch({
    out <- suppressWarnings(system2("blastdbcmd", c("-db", db_path, "-info"), stdout = TRUE, stderr = TRUE))
    status <- attr(out, "status")
    if (!is.null(status) && status != 0) return(NULL)
    paste(out[nzchar(out)], collapse = " | ")
  }, error = function(e) NULL)
}

# Writes one FASTA record per ASV (sequence itself, named ASV1, ASV2, ...)
# from a seqtab's column names -- the representative sequences DADA2 uses as
# column names throughout.
write_asv_fasta <- function(seqtab_nochim, out_path) {
  seqs <- colnames(seqtab_nochim)
  ids <- paste0("ASV", seq_along(seqs))
  lines <- as.vector(rbind(paste0(">", ids), seqs))
  writeLines(lines, out_path)
  data.frame(asv_id = ids, sequence = seqs, stringsAsFactors = FALSE)
}

# Reads blastn's headerless outfmt-6 TSV (qseqid sseqid pident length qcovhsp
# evalue bitscore), applies the coverage/e-value thresholds (identity is
# already enforced by blastn's own -perc_identity), and keeps each ASV's
# single best-passing hit (lowest e-value, ties broken by higher bitscore).
best_hits <- function(blast_raw_path, identity_thr, coverage_thr, evalue_thr) {
  cols <- c("qseqid", "sseqid", "pident", "length", "qcovhsp", "evalue", "bitscore")
  if (!file.exists(blast_raw_path) || file.info(blast_raw_path)$size == 0) {
    return(stats::setNames(data.frame(matrix(nrow = 0, ncol = length(cols))), cols))
  }
  raw <- utils::read.delim(blast_raw_path, header = FALSE, col.names = cols, stringsAsFactors = FALSE)
  passing <- raw[raw$pident >= identity_thr & raw$qcovhsp >= coverage_thr & raw$evalue < evalue_thr, , drop = FALSE]
  if (nrow(passing) == 0) return(passing)
  passing <- passing[order(passing$qseqid, passing$evalue, -passing$bitscore), ]
  passing[!duplicated(passing$qseqid), ]
}

# Derives a VT/lineage mapping straight from the BLAST database's own
# sequence titles -- no separate MaarjAM taxonomy download exists, but its
# reference FASTA headers already carry Family, Genus, and Virtual Taxon ID,
# e.g. ">gi|ESA00914|gb|AJ854100| Paraglomeraceae Paraglomus sp. 18S rDNA
# TYPE VTX00001". blastdbcmd's "%t" title output starts with the exact same
# id blastn reports as sseqid (both are just the raw defline up to the first
# space, since the db is built without -parse_seqids), so splitting each
# title on its first whitespace reconstructs the sseqid -> description join
# key with no dependency on the original FASTA file being available.
derive_maarjam_mapping <- function(db_path) {
  titles <- tryCatch(
    suppressWarnings(system2("blastdbcmd", c("-db", db_path, "-entry", "all", "-outfmt", "%t"), stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0)
  )
  titles <- titles[nzchar(titles)]
  if (length(titles) == 0) return(NULL)

  split_at_first_space <- function(x) {
    pos <- regexpr("\\s", x)
    if (pos < 0) return(c(x, ""))
    c(substr(x, 1, pos - 1), trimws(substr(x, pos + 1, nchar(x))))
  }
  parts <- vapply(titles, split_at_first_space, character(2), USE.NAMES = FALSE)
  sseqid <- parts[1, ]
  desc <- sub("^\\|\\s*", "", parts[2, ])  # a few upstream headers have a stray leading '|'

  vt <- regmatches(desc, regexpr("VTX[0-9]+", desc))
  vt[lengths(regmatches(desc, gregexpr("VTX[0-9]+", desc))) == 0] <- NA_character_
  tokens <- strsplit(desc, "\\s+")
  family <- vapply(tokens, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
  genus <- vapply(tokens, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))

  lineage <- mapply(function(f, g, v) paste(stats::na.omit(c(f, g, v)), collapse = ";"),
                     family, genus, vt)

  data.frame(sseqid = sseqid, family = family, genus = genus, vt = vt, lineage = lineage,
             stringsAsFactors = FALSE)
}

# Left-joins best-hit rows onto every ASV (so no-hit and filtered-out ASVs
# still appear, correctly labeled "unclassified"), and onto an optional
# mapping table (any columns, joined on "sseqid") for VT/lineage if supplied.
build_assignments <- function(asv_df, hits_df, mapping_df) {
  merged <- merge(asv_df, hits_df, by.x = "asv_id", by.y = "qseqid", all.x = TRUE)
  merged <- merged[match(asv_df$asv_id, merged$asv_id), ]  # preserve ASV1, ASV2, ... order
  if (!is.null(mapping_df) && "sseqid" %in% names(mapping_df)) {
    merged <- merge(merged, mapping_df, by = "sseqid", all.x = TRUE, sort = FALSE)
    merged <- merged[match(asv_df$asv_id, merged$asv_id), ]
  }
  merged$classification <- ifelse(is.na(merged$sseqid), "unclassified", "AM fungal ASV")
  merged
}

mod_taxonomy_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "9. Taxonomy Assignment using BLAST",
    description = "Compares each ASV against the MaarjAM reference database with BLAST+ and assigns a virtual taxon (VT) when the match passes your identity, coverage and e-value thresholds.",
    params = shiny::tagList(
      shiny::uiOutput(ns("blast_gate_alert")),
      param_group("Reference database", shiny::uiOutput(ns("db_setup"))),
      shiny::uiOutput(ns("input_status")),
      param_group("Match thresholds",
        shiny::fluidRow(
          shiny::column(3, shiny::numericInput(ns("identity"), "Min. % identity", value = taxonomy_defaults$identity, min = 0, max = 100, width = "100%")),
          shiny::column(3, shiny::numericInput(ns("coverage"), "Min. % query coverage", value = taxonomy_defaults$coverage, min = 0, max = 100, width = "100%")),
          shiny::column(3, shiny::textInput(ns("evalue"), "Max. e-value", value = taxonomy_defaults$evalue, width = "100%")),
          shiny::column(3, shiny::numericInput(ns("max_target_seqs"), "Max target seqs", value = taxonomy_defaults$max_target_seqs, min = 1, width = "100%"))
        ),
        shiny::div(class = "amf-param-hint", "Defaults per Davison et al. 2012 / Lekberg et al. 2018: \u226597% identity, \u226595% query coverage, e-value < 1e-50.")
      )
    ),
    results = shiny::tagList(
      result_card("Run BLAST",
        run_row(ns, "Run BLAST"),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("log"))
      ),
      shiny::uiOutput(ns("results_section")),
      results_placeholder("Run this step to see the ASV and taxonomy tables.")
    )
  )
}

mod_taxonomy_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "taxonomy"

    seqtab_ready <- shiny::reactive({
      shiny::req(rv$project_dir)
      # Depending on rv$status$chimera (even though unused) makes this
      # recompute the moment Remove Chimeras finishes, instead of caching a
      # FALSE result from before seqtab_nochim.rds existed on disk.
      rv$status$chimera
      fs::file_exists(file.path(denoise_dir_path(rv$project_dir), "seqtab_nochim.rds"))
    })

    output$input_status <- shiny::renderUI({
      if (!isTRUE(seqtab_ready())) {
        return(shiny::div(class = "alert alert-warning small", "No non-chimeric ASVs found yet -- finish Remove Chimeras first."))
      }
      NULL
    })

    # --- BLAST+ readiness: lazy check, surfaced only here (per earlier
    # decision to not front-load binary checks at Setup) ---
    blast_ready <- shiny::reactive({
      nzchar(Sys.which("blastn")) && nzchar(Sys.which("makeblastdb"))
    })
    output$blast_gate_alert <- shiny::renderUI({
      if (isTRUE(blast_ready())) return(NULL)
      shiny::div(class = "alert alert-warning small mb-3",
                 "blastn/makeblastdb not found -- install via 'conda install -c bioconda blast' to run this step.")
    })

    # --- MaarjAM database: point at an existing BLAST+ db, or index a
    # reference FASTA. Persists like the cutadapt path. No auto-download --
    # DESIGN.md leaves the canonical MaarjAM source as an open question, and
    # a URL wasn't confirmed, so this never guesses one. ---

    db_path <- shiny::reactiveVal(NULL)
    db_info <- shiny::reactiveVal(NULL)

    shiny::isolate({
      cfg <- read_app_config()
      if (!is.null(cfg$maarjam_db_path)) {
        info <- validate_blast_db(cfg$maarjam_db_path)
        if (!is.null(info)) {
          db_path(cfg$maarjam_db_path)
          db_info(info)
        }
      }
    })

    build_status <- shiny::reactiveVal("IDLE")
    build_handle <- shiny::reactiveVal(NULL)
    build_log <- shiny::reactiveVal(character(0))

    shiny::observeEvent(input$validate_db_path, {
      shiny::req(nzchar(input$db_path_input))
      if (grepl(" ", input$db_path_input, fixed = TRUE)) {
        output$db_validation <- shiny::renderUI(shiny::div(class = "text-danger small mt-1",
          "BLAST+ cannot handle paths containing spaces (-db treats spaces as database-name separators). Move the database to a path with no spaces."))
        return()
      }
      info <- validate_blast_db(input$db_path_input)
      if (is.null(info)) {
        output$db_validation <- shiny::renderUI(shiny::div(class = "text-danger small mt-1", "Not a valid BLAST+ database (blastdbcmd -info failed)."))
      } else {
        write_app_config(utils::modifyList(read_app_config(), list(maarjam_db_path = input$db_path_input)))
        db_path(input$db_path_input)
        db_info(info)
        output$db_validation <- shiny::renderUI(NULL)
      }
    })

    shiny::observeEvent(input$build_db, {
      shiny::req(identical(build_status(), "IDLE") || identical(build_status(), "ERROR"))
      fasta_path <- trimws(input$fasta_path_input %||% "")
      if (!nzchar(fasta_path)) {
        build_log("No reference FASTA path entered.")
        build_status("ERROR")
        return()
      }
      if (!file.exists(fasta_path)) {
        build_log(paste("File not found:", fasta_path))
        build_status("ERROR")
        return()
      }
      if (fs::is_dir(fasta_path)) {
        build_log("Path is a directory -- point 'reference FASTA' at the .fasta file itself, not its folder.")
        build_status("ERROR")
        return()
      }
      out_prefix <- sub("/+$", "", tools::file_path_sans_ext(fasta_path))
      if (!nzchar(basename(out_prefix))) {
        build_log("Could not derive a database name/prefix from that path.")
        build_status("ERROR")
        return()
      }
      # makeblastdb's -in/-out both accept space-separated lists of multiple
      # filenames, so it splits on any unescaped space in the path -- even
      # passed as a single, correctly quoted argv element with no shell
      # involved. There is no escaping workaround; the path must be moved.
      if (grepl(" ", fasta_path, fixed = TRUE) || grepl(" ", out_prefix, fixed = TRUE)) {
        build_log("BLAST+ cannot handle paths containing spaces (it treats spaces as filename separators in -in/-out). Move the FASTA to a path with no spaces and try again.")
        build_status("ERROR")
        return()
      }
      h <- processx::process$new("makeblastdb", c("-in", fasta_path, "-dbtype", "nucl", "-out", out_prefix),
                                  stdout = "|", stderr = "|")
      build_handle(h)
      build_status("RUNNING")
      build_log(character(0))
      session$userData[[paste0(id, "_build_out_prefix")]] <- out_prefix
    })

    shiny::observe({
      h <- build_handle()
      shiny::req(!is.null(h), identical(build_status(), "RUNNING"))
      shiny::invalidateLater(750, session)
      new_lines <- c(tryCatch(h$read_output_lines(), error = function(e) character(0)),
                     tryCatch(h$read_error_lines(), error = function(e) character(0)))
      if (length(new_lines) > 0) build_log(c(build_log(), new_lines))
      if (!h$is_alive()) {
        if (identical(h$get_exit_status(), 0L)) {
          out_prefix <- session$userData[[paste0(id, "_build_out_prefix")]]
          info <- validate_blast_db(out_prefix)
          write_app_config(utils::modifyList(read_app_config(), list(maarjam_db_path = out_prefix)))
          db_path(out_prefix)
          db_info(info %||% "built")
          build_status("SUCCESS")
        } else {
          build_status("ERROR")
        }
      }
    })

    shiny::observeEvent(input$reset_db, {
      db_path(NULL)
      db_info(NULL)
    })

    output$db_setup <- shiny::renderUI({
      if (!is.null(db_path())) {
        return(shiny::tagList(
          shiny::div(
            class = "alert alert-success small mb-2 d-flex justify-content-between align-items-center",
            shiny::span(paste("MaarjAM DB ready:", db_path())),
            shiny::actionLink(ns("reset_db"), "Change")
          ),
          shiny::div(class = "text-muted small mb-2 font-monospace", db_info())
        ))
      }

      shiny::div(
        class = "border rounded p-3 mb-3",
        shiny::h6("MaarjAM BLAST+ database"),
        shiny::div(class = "text-muted small mb-2", "Point at an existing BLAST+ database, or supply a reference FASTA to index."),
        shiny::fluidRow(
          shiny::column(9, shiny::textInput(ns("db_path_input"), "Path to existing BLAST+ database (prefix, no extension)", width = "100%")),
          shiny::column(3, shiny::actionButton(ns("validate_db_path"), "Validate & save", class = "btn-outline-primary w-100 mt-4"))
        ),
        shiny::uiOutput(ns("db_validation")),
        shiny::div(class = "text-muted small my-2", "-- or --"),
        shiny::fluidRow(
          shiny::column(9, shiny::textInput(ns("fasta_path_input"), "Path to MaarjAM reference FASTA (will be indexed with makeblastdb)", width = "100%")),
          shiny::column(3, shiny::actionButton(ns("build_db"), "Build index", class = "btn-primary w-100 mt-4"))
        ),
        shiny::div(class = "mb-1", status_badge(build_status())),
        if (length(build_log()) > 0) {
          shiny::tags$div(class = "amf-log-console", shiny::HTML(paste(htmltools::htmlEscape(build_log()), collapse = "<br/>")))
        }
      )
    })

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv))
      evalue_num <- suppressWarnings(as.numeric(input$evalue))
      msg <- if (!isTRUE(blast_ready())) "BLAST+ (blastn/makeblastdb) is not installed."
        else if (is.null(db_path())) "Set up the MaarjAM database first (Reference database, above)."
        else if (!isTRUE(seqtab_ready())) "No non-chimeric ASVs found yet -- finish Remove Chimeras first."
        else if (is.na(evalue_num)) "Max. e-value must be a number, e.g. 1e-50."
      if (!is.null(msg)) return(shiny::showNotification(msg, type = "error"))

      out_dir <- taxonomy_dir_path(rv$project_dir)
      fs::dir_create(out_dir)
      seqtab_nochim <- readRDS(file.path(denoise_dir_path(rv$project_dir), "seqtab_nochim.rds"))
      asv_fasta <- file.path(out_dir, "asv.fasta")
      asv_df <- write_asv_fasta(seqtab_nochim, asv_fasta)

      raw_tsv <- file.path(out_dir, "blast_raw.tsv")

      script <- paste(c(
        "set -e",
        sprintf(
          "blastn -db %s -query %s -outfmt '6 qseqid sseqid pident length qcovhsp evalue bitscore' -perc_identity %s -max_target_seqs %s -out %s",
          shQuote(db_path()), shQuote(asv_fasta), input$identity, input$max_target_seqs, shQuote(raw_tsv)
        ),
        sprintf("echo '%d ASVs searched against %s.'", nrow(asv_df), shQuote(db_path()))
      ), collapse = "\n")

      session$userData[[paste0(id, "_run_ctx")]] <- list(
        asv_df = asv_df, raw_tsv = raw_tsv, identity = input$identity,
        coverage = input$coverage, evalue = evalue_num, db_path = db_path()
      )
      launch_step(rv, step_id, "bash", c("-c", script))
    })

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    # Once the blastn job finishes, do the (fast) R-side assignment logic:
    # threshold filtering, best-hit selection, optional mapping join.
    step_status <- shiny::reactive(rv$status[[step_id]])
    shiny::observeEvent(step_status(), {
      shiny::req(identical(step_status(), "SUCCESS"))
      ctx <- session$userData[[paste0(id, "_run_ctx")]]
      shiny::req(!is.null(ctx))

      hits <- best_hits(ctx$raw_tsv, ctx$identity, ctx$coverage, ctx$evalue)
      utils::write.table(hits, file.path(taxonomy_dir_path(rv$project_dir), "blast_besthit.tsv"),
                          sep = "\t", row.names = FALSE, quote = FALSE)

      # A manually supplied override always wins; otherwise derive the
      # VT/lineage mapping straight from the database's own headers, since
      # no separate MaarjAM taxonomy file exists to require of the user.
      mapping_df <- if (!is.null(ctx$db_path)) derive_maarjam_mapping(ctx$db_path) else NULL
      assignments <- build_assignments(ctx$asv_df, hits, mapping_df)
      saveRDS(assignments, file.path(taxonomy_dir_path(rv$project_dir), "vt_assignments.rds"))
      rv$artifacts[[step_id]] <- assignments
    }, ignoreInit = TRUE)

    # ASV (abundance) table: rows = ASV, columns = its sequence plus one count
    # column per sample -- read straight off seqtab_nochim.rds, whose column
    # order write_asv_fasta() used to assign ASV1, ASV2, ... in the first
    # place, so it lines up with res$asv_id without needing a join.
    asv_count_table <- shiny::reactive({
      res <- rv$artifacts[[step_id]]
      shiny::req(res)
      path <- file.path(denoise_dir_path(rv$project_dir), "seqtab_nochim.rds")
      shiny::req(fs::file_exists(path))
      seqtab_nochim <- readRDS(path)
      counts <- as.data.frame(t(seqtab_nochim), check.names = FALSE)
      cbind(data.frame(ASV = res$asv_id, Sequence = res$sequence, stringsAsFactors = FALSE), counts)
    })

    output$results_section <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res))
      n_assigned <- sum(res$classification == "AM fungal ASV")
      lineage_col <- grep("lineage", names(res), ignore.case = TRUE, value = TRUE)
      result_card("Taxonomy results",
        shiny::div(class = "mb-2", sprintf("%d / %d ASVs classified as AM fungal.", n_assigned, nrow(res))),
        if (length(lineage_col) > 0) {
          shiny::tagList(shiny::h6("Taxonomic distribution"), shiny::plotOutput(ns("lineage_barplot"), height = "260px"))
        } else {
          shiny::div(class = "text-muted small", "No lineage column in the mapping file -- showing BLAST results only.")
        },
        shiny::h6("ASV table"),
        shiny::div(class = "text-muted small mb-1", "Per-ASV read counts by sample."),
        dt_output(ns("asv_table")),
        shiny::downloadButton(ns("download_asv_table"), "Download as CSV", class = "btn-outline-secondary btn-sm mb-3"),
        shiny::h6("Taxonomy table"),
        shiny::div(class = "text-muted small mb-1", "Best BLAST hit and classification per ASV."),
        dt_output(ns("assignments_table")),
        shiny::downloadButton(ns("download_assignments_table"), "Download as CSV", class = "btn-outline-secondary btn-sm mb-2")
      )
    })

    output$lineage_barplot <- shiny::renderPlot({
      res <- rv$artifacts[[step_id]]
      shiny::req(res)
      lineage_col <- grep("lineage", names(res), ignore.case = TRUE, value = TRUE)
      shiny::req(length(lineage_col) > 0)
      assigned <- res[res$classification == "AM fungal ASV", , drop = FALSE]
      shiny::req(nrow(assigned) > 0)
      last_rank <- vapply(strsplit(assigned[[lineage_col[1]]], ";"), function(x) if (length(x) > 0) trimws(utils::tail(x, 1)) else NA_character_, character(1))
      counts <- sort(table(last_rank), decreasing = TRUE)
      graphics::barplot(counts, las = 2, col = "#0f766e", border = "white", ylab = "ASVs")
    })

    # On screen the (400+ bp) sequence is shortened, with the full sequence
    # on hover; the CSV download keeps it whole.
    output$asv_table <- DT::renderDataTable({
      df <- asv_count_table()
      DT::datatable(df, rownames = FALSE, options = list(
        scrollX = TRUE, pageLength = 15,
        columnDefs = list(list(targets = 1, render = DT::JS(
          "function(d, type) { return type === 'display' && d.length > 24 ?",
          "'<span title=\"' + d + '\">' + d.substr(0, 24) + '&hellip;</span>' : d; }")))
      ))
    })

    output$download_asv_table <- shiny::downloadHandler(
      filename = function() sprintf("asv_abundance_table_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        utils::write.csv(asv_count_table(), file, row.names = FALSE)
      }
    )

    taxonomy_table_df <- shiny::reactive({
      res <- rv$artifacts[[step_id]]
      if (is.null(res)) {
        return(data.frame(ASV = character(0), `Best hit` = character(0), `% id` = numeric(0),
                           `% coverage` = numeric(0), `e-value` = numeric(0), Classification = character(0), check.names = FALSE))
      }
      keep_cols <- intersect(c("asv_id", "sseqid", "pident", "qcovhsp", "evalue", "classification",
                                names(res)[!(names(res) %in% c("asv_id", "sequence", "sseqid", "pident", "length", "qcovhsp", "evalue", "bitscore", "classification"))]),
                              names(res))
      df <- res[, keep_cols, drop = FALSE]
      names(df)[names(df) == "asv_id"] <- "ASV"
      names(df)[names(df) == "sseqid"] <- "Best hit"
      names(df)[names(df) == "pident"] <- "% id"
      names(df)[names(df) == "qcovhsp"] <- "% coverage"
      names(df)[names(df) == "evalue"] <- "e-value"
      names(df)[names(df) == "classification"] <- "Classification"
      names(df)[names(df) == "family"] <- "Family"
      names(df)[names(df) == "genus"] <- "Genus"
      names(df)[names(df) == "vt"] <- "Virtual taxon"
      names(df)[names(df) == "lineage"] <- "Lineage"
      if ("e-value" %in% names(df)) df[["e-value"]] <- signif(df[["e-value"]], 3)
      df
    })

    output$assignments_table <- DT::renderDataTable({
      DT::datatable(taxonomy_table_df(), rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 15))
    })

    output$download_assignments_table <- shiny::downloadHandler(
      filename = function() sprintf("taxonomy_table_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        utils::write.csv(taxonomy_table_df(), file, row.names = FALSE)
      }
    )

    mod_logpanel_server("log", shiny::reactive(rv$log[[step_id]] %||% character(0)))
  })
}
