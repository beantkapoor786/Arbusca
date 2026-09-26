# Step 3 -- Primer Removal (cutadapt). Presets come from inst/primer_presets.yaml
# (each with its literature citation); the dropdown auto-selects the pair
# that matches the detected read length, per DESIGN.md section 5.

load_primer_presets <- function() {
  presets <- yaml::read_yaml("inst/primer_presets.yaml")
  stats::setNames(presets, vapply(presets, function(p) p$id, character(1)))
}

# --- custom primer sets remembered across sessions (~/.arbusca/config.yaml) ---
# Saved automatically whenever a custom (non-preset) fwd/rev pair is used for
# a primer test or a run, most-recently-used first, so returning to a project
# with the same custom primers doesn't require retyping them.

MAX_SAVED_PRIMER_SETS <- 10

get_saved_primer_sets <- function() {
  cfg <- read_app_config()
  sets <- cfg$custom_primer_sets
  if (is.null(sets)) return(list())
  Filter(function(s) !is.null(s$fwd) && !is.null(s$rev) && nzchar(s$fwd) && nzchar(s$rev), sets)
}

# Adds/moves (fwd, rev) to the front of the saved list, deduped by exact text
# match, capped at MAX_SAVED_PRIMER_SETS.
save_primer_set <- function(fwd, rev) {
  if (!nzchar(fwd) || !nzchar(rev)) return(invisible())
  sets <- get_saved_primer_sets()
  sets <- Filter(function(s) !(identical(s$fwd, fwd) && identical(s$rev, rev)), sets)
  sets <- c(list(list(fwd = fwd, rev = rev)), sets)
  sets <- utils::head(sets, MAX_SAVED_PRIMER_SETS)
  write_app_config(utils::modifyList(read_app_config(), list(custom_primer_sets = sets)))
  invisible(sets)
}

saved_primer_set_label <- function(s) {
  shorten <- function(x) if (nchar(x) > 24) paste0(substr(x, 1, 24), "…") else x
  sprintf("%s / %s", shorten(s$fwd), shorten(s$rev))
}

# Named-list choices for the preset selectInput: built-in presets grouped
# under "Presets", saved custom pairs (if any) grouped under "My saved
# primers", each keyed "saved_<i>" so apply_preset() can tell them apart.
build_preset_choices <- function(presets, saved_sets) {
  builtin <- c(stats::setNames(names(presets), vapply(presets, function(p) p$label, character(1))), "Custom" = "custom")
  choices <- list(Presets = builtin)
  if (length(saved_sets) > 0) {
    ids <- paste0("saved_", seq_along(saved_sets))
    labels <- vapply(saved_sets, saved_primer_set_label, character(1))
    choices[["My saved primers"]] <- stats::setNames(ids, labels)
  }
  choices
}

# IUPAC-aware reverse complement, for building the "read-through" adapter
# (the other primer's revcomp can appear at the 3' end when the amplicon is
# shorter than the read).
revcomp <- function(seq) {
  comp <- c(
    A = "T", C = "G", G = "C", T = "A", U = "A",
    R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
    B = "V", D = "H", H = "D", V = "B", N = "N"
  )
  chars <- strsplit(toupper(seq), "")[[1]]
  paste(rev(unname(comp[chars])), collapse = "")
}

# Reads the first record of a (optionally gzipped) FASTQ file and returns
# its sequence length, or NA if the file can't be read.
detect_read_length <- function(path) {
  tryCatch({
    con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con))
    readLines(con, n = 1)              # @header
    seq_line <- readLines(con, n = 1)  # sequence
    nchar(seq_line)
  }, error = function(e) NA_integer_)
}

default_preset_id <- function(read_len, presets) {
  if (is.na(read_len)) return(names(presets)[1])
  if (abs(read_len - 250) <= 20) return("nu_ssu_0595_0948")
  if (abs(read_len - 300) <= 20) return("nu_ssu_0450_0899")
  names(presets)[1]
}

# Splits a comma-separated primer field into a trimmed, non-empty character
# vector -- lets a field hold multiple candidate primer variants (e.g. when
# the user isn't sure which of two primers was actually used).
parse_primer_list <- function(text) {
  parts <- trimws(strsplit(text %||% "", ",", fixed = TRUE)[[1]])
  toupper(parts[nzchar(parts)])
}

# NULL if every primer is a plain IUPAC DNA sequence, else a message naming
# the offenders (e.g. a primer *name* like "NS31" typed instead of its bases).
invalid_primer_msg <- function(fwd, rev) {
  bad <- c(fwd, rev)[!grepl("^[ACGTURYSWKMBDHVN]+$", c(fwd, rev))]
  if (length(fwd) == 0 || length(rev) == 0) return("Enter at least one forward and one reverse primer sequence.")
  if (length(bad) == 0) return(NULL)
  sprintf("Not a DNA sequence: %s. Enter primer sequences (IUPAC bases), not primer names.", paste(bad, collapse = ", "))
}

# All four orientations of a primer (IUPAC-aware, via Biostrings) -- Forward
# (as given), Complement, Reverse, and RevComp. Checking all four against the
# raw reads catches a primer that's present but misoriented in the data,
# which checking only the as-entered form would miss.
allOrients <- function(primer) {
  dna <- Biostrings::DNAString(primer)
  orients <- list(
    Forward = dna,
    Complement = Biostrings::complement(dna),
    Reverse = Biostrings::reverse(dna),
    RevComp = Biostrings::reverseComplement(dna)
  )
  vapply(orients, as.character, character(1))
}

# Expands each primer in a role (Forward/Reverse) into its 4 orientations,
# as a long data.frame(role, primer, orientation, sequence).
expand_primer_orientations <- function(primers, role) {
  rows <- lapply(primers, function(p) {
    or <- allOrients(p)
    data.frame(role = role, primer = p, orientation = names(or), sequence = unname(or), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Builds a bash script that checks every candidate primer, in all 4
# orientations, against R1 and R2 *independently* (single-end cutadapt,
# count-only, output discarded), emitting a
# "RESULT|sample|role|primer|orientation|sequence|r1_count|r1_total|r2_count|r2_total"
# line per (sample, primer, orientation) -- not assuming a fixed F/R pairing
# or a fixed orientation is what lets this surface "this primer actually
# shows up reverse-complemented in R2" when the user isn't sure which primer
# was used, or how.
build_primer_test_script <- function(cutadapt_cmd, files_df, fwd_primers, rev_primers, err_rate) {
  variants <- rbind(
    expand_primer_orientations(fwd_primers, "Forward"),
    expand_primer_orientations(rev_primers, "Reverse")
  )
  lines <- c("set -e")
  count_block <- function(var_prefix, seq, read_path) {
    c(
      sprintf("%sOUT=$(%s -g %s -e %s -o /dev/null %s 2>&1)", var_prefix, cutadapt_cmd, shQuote(seq), err_rate, shQuote(read_path)),
      sprintf("%sTOTAL=$(echo \"$%sOUT\" | grep -m1 'Total reads processed' | awk -F':' '{print $NF}' | awk '{gsub(\",\",\"\"); print $1}')", var_prefix, var_prefix),
      sprintf("%sCNT=$(echo \"$%sOUT\" | grep -m1 'Reads with adapters' | awk -F':' '{print $NF}' | awk '{gsub(\",\",\"\"); print $1}')", var_prefix, var_prefix)
    )
  }
  for (i in seq_len(nrow(files_df))) {
    s <- files_df$sample[i]
    in1 <- files_df$fwd[i]
    in2 <- files_df$rev[i]
    lines <- c(lines, sprintf("echo '--- %s ---'", s))
    for (j in seq_len(nrow(variants))) {
      role <- variants$role[j]
      primer <- variants$primer[j]
      orientation <- variants$orientation[j]
      seq <- variants$sequence[j]
      lines <- c(
        lines,
        count_block("R1", seq, in1),
        count_block("R2", seq, in2),
        sprintf("echo \"RESULT|%s|%s|%s|%s|%s|$R1CNT|$R1TOTAL|$R2CNT|$R2TOTAL\"", s, role, primer, orientation, seq)
      )
    }
  }
  paste(lines, collapse = "\n")
}

# Parses RESULT lines out of a primer-test job's log into a results data.frame.
parse_primer_results <- function(log_lines) {
  result_lines <- grep("^RESULT\\|", log_lines, value = TRUE)
  empty <- data.frame(sample = character(0), role = character(0), primer = character(0), orientation = character(0),
                       sequence = character(0), r1_count = numeric(0), r1_total = numeric(0),
                       r2_count = numeric(0), r2_total = numeric(0))
  if (length(result_lines) == 0) return(empty)
  parts <- strsplit(result_lines, "\\|")
  rows <- lapply(parts, function(p) {
    data.frame(
      sample = p[2], role = p[3], primer = p[4], orientation = p[5], sequence = p[6],
      r1_count = suppressWarnings(as.numeric(p[7])), r1_total = suppressWarnings(as.numeric(p[8])),
      r2_count = suppressWarnings(as.numeric(p[9])), r2_total = suppressWarnings(as.numeric(p[10]))
    )
  })
  do.call(rbind, rows)
}

PRIMER_TABLE_PREVIEW_ROWS <- 10

# show_all = FALSE truncates to the first PRIMER_TABLE_PREVIEW_ROWS rows, so
# a long sample list doesn't turn the page into a scroll marathon; callers
# pair this with a "Show full table" button (see primer_table_toggle_ui()).
primer_results_table <- function(df, show_all = FALSE) {
  empty_cols <- data.frame(`Sample` = character(0), `Role` = character(0), `Primer` = character(0),
                            `Orientation` = character(0), `Sequence tested` = character(0),
                            `Found in R1` = character(0), `Found in R2` = character(0), check.names = FALSE)
  if (is.null(df) || nrow(df) == 0) {
    return(DT::datatable(empty_cols, rownames = FALSE, options = list(dom = "t")))
  }
  fmt <- function(count, total) {
    pct <- if (total > 0) round(100 * count / total, 1) else NA
    sprintf("%d (%s%%)", count, if (is.na(pct)) "NA" else pct)
  }
  out <- data.frame(
    Sample = df$sample, Role = df$role, Primer = df$primer,
    Orientation = df$orientation, `Sequence tested` = df$sequence,
    `Found in R1` = mapply(fmt, df$r1_count, df$r1_total),
    `Found in R2` = mapply(fmt, df$r2_count, df$r2_total),
    check.names = FALSE
  )
  if (!show_all && nrow(out) > PRIMER_TABLE_PREVIEW_ROWS) {
    out <- out[seq_len(PRIMER_TABLE_PREVIEW_ROWS), , drop = FALSE]
  }
  DT::datatable(out, rownames = FALSE, options = list(dom = "t", pageLength = nrow(out) + 1))
}

# Renders the "Show full table (N rows)" / "Show first 10 rows" toggle for a
# primer results table, or NULL if there's nothing to truncate.
primer_table_toggle_ui <- function(ns, input_id, df, show_all) {
  n <- if (is.null(df)) 0 else nrow(df)
  if (n <= PRIMER_TABLE_PREVIEW_ROWS) return(NULL)
  label <- if (show_all) "Show first 10 rows" else sprintf("Show full table (%d rows)", n)
  shiny::actionLink(ns(input_id), label, class = "small")
}

mod_primer_ui <- function(id) {
  ns <- shiny::NS(id)
  presets <- load_primer_presets()
  choices <- build_preset_choices(presets, get_saved_primer_sets())

  step_card(
    stacked = TRUE,
    title = "3. Primer Removal",
    description = "Trims the PCR primer sequences off the ends of every read with cutadapt, so primer bases are not mistaken for biological variation. Test your primers on the raw reads first if you are unsure which were used.",
    params = shiny::tagList(
      shiny::uiOutput(ns("cutadapt_setup")),
      shiny::uiOutput(ns("read_length_caption")),
      param_group("Primers",
        shiny::fluidRow(
          shiny::column(4, shiny::selectInput(ns("preset"), "Primer preset", choices = choices, width = "100%")),
          shiny::column(8, shiny::uiOutput(ns("citation_caption")))
        ),
        shiny::fluidRow(
          shiny::column(6, shiny::textInput(ns("fwd_primer"), "Forward primer(s) (5'->3', comma-separated)", placeholder = "e.g. TTGGAGGGCAAGTCTGGTGCC", width = "100%")),
          shiny::column(6, shiny::textInput(ns("rev_primer"), "Reverse primer(s) (5'->3', comma-separated)", placeholder = "e.g. GAACCCAAACACTTTGGTTTCC", width = "100%"))
        ),
        shiny::uiOutput(ns("primer_feedback")),
        shiny::div(class = "amf-param-hint", "Not sure which primer was used? List multiple candidates, separated by commas -- cutadapt tries each.")
      ),
      param_group("Matching",
        shiny::fluidRow(
          shiny::column(4, shiny::tagList(
            help_label("Error rate", "Max fraction of mismatches cutadapt allows when matching the primer against a read (substitutions, insertions, and deletions all count). E.g. with a 20bp primer, an error rate of 0.1 allows up to 2 mismatched bases; 0 requires an exact match."),
            shiny::numericInput(ns("err_rate"), NULL, value = 0, min = 0, max = 1, step = 0.01, width = "100%")
          ))
        )
      )
    ),
    results = shiny::tagList(
      result_card("Test primers",
        shiny::div(
          class = "d-flex align-items-center gap-2 mb-2",
          shiny::actionButton(ns("test_primers"), "Test these primers", class = "btn-primary btn-sm"),
          shiny::uiOutput(ns("test_badge"), inline = TRUE),
          shiny::span(class = "amf-param-hint m-0", "It's always a good idea to check if these primers were used in the library prep.")
        ),
        shiny::div(class = "amf-param-hint mb-2", "Checks whether these primers appear in the raw reads, without writing any output."),
        dt_output(ns("test_results_table")),
        shiny::uiOutput(ns("test_results_toggle_ui")),
        shiny::uiOutput(ns("test_results_download_ui"))
      ),
      result_card("Remove primers with cutadapt",
        run_row(ns, "Run Cutadapt"),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("log"))
      ),
      shiny::uiOutput(ns("after_results_section"))
    )
  )
}

# sample_table: reactive() -> data.frame(sample, fwd, rev) of basenames,
# from mod_setup. rv$project_dir resolves them to full paths.
mod_primer_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "primer"
    presets <- load_primer_presets()
    saved_sets <- shiny::reactiveVal(get_saved_primer_sets())

    apply_preset <- function(preset_id) {
      if (identical(preset_id, "custom")) {
        shiny::updateTextInput(session, "fwd_primer", value = "")
        shiny::updateTextInput(session, "rev_primer", value = "")
        return(invisible())
      }
      if (startsWith(preset_id, "saved_")) {
        idx <- as.integer(sub("^saved_", "", preset_id))
        s <- saved_sets()[[idx]]
        if (is.null(s)) return(invisible())
        shiny::updateTextInput(session, "fwd_primer", value = s$fwd)
        shiny::updateTextInput(session, "rev_primer", value = s$rev)
        return(invisible())
      }
      p <- presets[[preset_id]]
      shiny::updateTextInput(session, "fwd_primer", value = p$fwd %||% "")
      shiny::updateTextInput(session, "rev_primer", value = p$rev %||% "")
    }

    # Remembers the current custom fwd/rev pair (only called when the dropdown
    # is on "Custom") and refreshes the dropdown's "My saved primers" group
    # in place, without disturbing the current selection.
    remember_if_custom <- function() {
      if (!identical(input$preset, "custom")) return(invisible())
      fwd <- trimws(input$fwd_primer %||% "")
      rev <- trimws(input$rev_primer %||% "")
      if (!nzchar(fwd) || !nzchar(rev)) return(invisible())
      saved_sets(save_primer_set(fwd, rev))
      shiny::updateSelectInput(session, "preset",
        choices = build_preset_choices(presets, saved_sets()), selected = "custom")
    }

    initialized <- shiny::reactiveVal(FALSE)
    shiny::observe({
      shiny::req(!initialized())
      st <- sample_table()
      shiny::req(nrow(st) > 0, rv$project_dir)
      len <- detect_read_length(file.path(rv$project_dir, st$fwd[1]))
      pid <- default_preset_id(len, presets)
      shiny::updateSelectInput(session, "preset", selected = pid)
      apply_preset(pid)
      initialized(TRUE)
    })

    shiny::observeEvent(input$preset, {
      apply_preset(input$preset)
    })

    output$read_length_caption <- shiny::renderUI({
      st <- sample_table()
      if (is.null(rv$project_dir) || nrow(st) == 0) {
        return(shiny::div(class = "text-muted small mb-2", "No paired samples detected yet -- finish Setup first."))
      }
      len <- detect_read_length(file.path(rv$project_dir, st$fwd[1]))
      shiny::div(
        class = "text-muted small mb-2",
        if (is.na(len)) "Could not detect read length from the first sample."
        else sprintf("Detected read length: %d bp (from %s)", len, st$fwd[1])
      )
    })

    # A plain <a target="_blank"> gets hijacked by RStudio's Viewer pane into
    # another RStudio window instead of the system browser. actionLink +
    # server-side utils::browseURL() sidesteps that entirely -- browseURL()
    # always launches the OS-registered default browser, Viewer or not.
    MAX_PRESET_DOIS <- 2
    output$citation_caption <- shiny::renderUI({
      if (identical(input$preset, "custom")) return(NULL)
      p <- presets[[input$preset]]
      shiny::req(p)
      n_dois <- min(length(p$dois), MAX_PRESET_DOIS)
      links <- if (n_dois > 0) {
        shiny::tagList(lapply(seq_len(n_dois), function(i) {
          label <- if (n_dois == 1) "[paper]" else sprintf("[paper %d]", i)
          shiny::tagList(shiny::actionLink(ns(paste0("open_doi_", i)), label), " ")
        }))
      } else NULL
      shiny::tagList(
        shiny::div(class = "text-muted small mb-2 fst-italic", paste0(p$region, " · ", p$platform, " — ", p$citation), " ", links)
      )
    })

    lapply(seq_len(MAX_PRESET_DOIS), function(i) {
      shiny::observeEvent(input[[paste0("open_doi_", i)]], {
        p <- presets[[input$preset]]
        shiny::req(p, length(p$dois) >= i)
        utils::browseURL(paste0("https://doi.org/", p$dois[[i]]))
      })
    })

    # --- cutadapt readiness: ask once, validate or install, remember the choice ---

    validate_cutadapt <- function(path) {
      tryCatch({
        out <- system2(path, "--version", stdout = TRUE, stderr = TRUE)
        if (!is.null(attr(out, "status")) && attr(out, "status") != 0) return(NULL)
        trimws(out[1])
      }, error = function(e) NULL)
    }

    cutadapt_invocation <- shiny::reactiveVal(NULL)  # shell-ready command prefix once resolved
    cutadapt_version <- shiny::reactiveVal(NULL)

    shiny::isolate({
      cfg <- read_app_config()
      if (!is.null(cfg$cutadapt_path)) {
        v <- validate_cutadapt(cfg$cutadapt_path)
        if (!is.null(v)) {
          cutadapt_invocation(shQuote(cfg$cutadapt_path))
          cutadapt_version(v)
        }
      }
    })

    install_status <- shiny::reactiveVal("IDLE")  # IDLE, RUNNING, SUCCESS, ERROR
    install_log <- shiny::reactiveVal(character(0))
    install_handle <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$validate_path, {
      shiny::req(nzchar(input$cutadapt_path))
      v <- validate_cutadapt(input$cutadapt_path)
      if (is.null(v)) {
        output$path_validation <- shiny::renderUI(
          shiny::div(class = "text-danger small mt-1", "Could not run that path with --version. Check it points to the cutadapt executable.")
        )
      } else {
        write_app_config(utils::modifyList(read_app_config(), list(cutadapt_path = input$cutadapt_path)))
        cutadapt_invocation(shQuote(input$cutadapt_path))
        cutadapt_version(v)
        output$path_validation <- shiny::renderUI(NULL)
      }
    })

    shiny::observeEvent(input$install_cutadapt, {
      shiny::req(identical(install_status(), "IDLE") || identical(install_status(), "ERROR"))
      h <- processx::process$new("python3", c("-m", "pip", "install", "--user", "cutadapt"), stdout = "|", stderr = "|")
      install_handle(h)
      install_status("RUNNING")
      install_log(character(0))
    })

    shiny::observe({
      h <- install_handle()
      shiny::req(!is.null(h), identical(install_status(), "RUNNING"))
      shiny::invalidateLater(750, session)

      new_out <- tryCatch(h$read_output_lines(), error = function(e) character(0))
      new_err <- tryCatch(h$read_error_lines(), error = function(e) character(0))
      if (length(c(new_out, new_err)) > 0) install_log(c(install_log(), new_out, new_err))

      if (!h$is_alive()) {
        if (identical(h$get_exit_status(), 0L)) {
          resolved_path <- Sys.which("cutadapt")
          inv <- if (nzchar(resolved_path)) shQuote(resolved_path) else "python3 -m cutadapt"
          v <- validate_cutadapt(if (nzchar(resolved_path)) resolved_path else "cutadapt")
          cutadapt_invocation(inv)
          cutadapt_version(v %||% "installed")
          write_app_config(utils::modifyList(read_app_config(), list(
            cutadapt_path = if (nzchar(resolved_path)) resolved_path else NA
          )))
          install_status("SUCCESS")
        } else {
          install_status("ERROR")
        }
      }
    })

    output$cutadapt_setup <- shiny::renderUI({
      if (!is.null(cutadapt_invocation())) {
        return(shiny::div(
          class = "alert alert-success small mb-3 d-flex justify-content-between align-items-center",
          shiny::span(paste("cutadapt ready:", cutadapt_version())),
          shiny::actionLink(ns("reset_cutadapt"), "Change")
        ))
      }

      shiny::tagList(
        shiny::div(
          class = "border rounded p-3 mb-3",
          shiny::h6("Do you have cutadapt installed?"),
          shiny::selectInput(
            ns("has_cutadapt"), NULL,
            choices = c("Choose one" = "", "Yes, it's installed" = "yes", "No, install it for me" = "no"),
            width = "260px"
          ),
          if (identical(input$has_cutadapt, "yes")) {
            shiny::tagList(
              shiny::fluidRow(
                shiny::column(9, shiny::textInput(ns("cutadapt_path"), "Path to cutadapt", value = Sys.which("cutadapt"), width = "100%")),
                shiny::column(3, shiny::actionButton(ns("validate_path"), "Validate & save", class = "btn-outline-primary w-100 mt-4"))
              ),
              shiny::uiOutput(ns("path_validation"))
            )
          } else if (identical(input$has_cutadapt, "no")) {
            shiny::tagList(
              shiny::div(class = "text-muted small mb-2", "Installs via: python3 -m pip install --user cutadapt"),
              shiny::actionButton(ns("install_cutadapt"), "Install cutadapt", class = "btn-primary btn-sm mb-2"),
              shiny::div(class = "mb-2", status_badge(install_status())),
              if (length(install_log()) > 0) {
                shiny::tags$div(
                  class = "amf-log-console",
                  shiny::HTML(paste(htmltools::htmlEscape(install_log()), collapse = "<br/>"))
                )
              }
            )
          }
        )
      )
    })

    shiny::observeEvent(input$reset_cutadapt, {
      cutadapt_invocation(NULL)
      cutadapt_version(NULL)
    })

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    # --- primer detection: "does this pair show up in the data?" ---
    # Self-contained (not part of rv's step FSM) since it's a diagnostic
    # probe, not a pipeline step -- it can run before OR after the real
    # trim, and neither counts toward pipeline progress.

    # err_rate is threaded through so "found in R1/R2" reflects the same
    # match tolerance the real Run uses -- testing at cutadapt's own default
    # (-e 0.1) regardless of what error rate was actually used to trim would
    # make the "after" check report leftover primer that a stricter Run
    # correctly left untouched.
    start_primer_check <- function(status_rv, handle_rv, log_rv, results_rv, files_df, fwd_primer, rev_primer, err_rate) {
      if (nrow(files_df) == 0) return(invisible())
      script <- build_primer_test_script(cutadapt_invocation(), files_df, fwd_primer, rev_primer, err_rate)
      h <- processx::process$new("bash", c("-c", script), stdout = "|", stderr = "|")
      handle_rv(h)
      status_rv("RUNNING")
      log_rv(character(0))
      results_rv(NULL)
    }

    poll_primer_check <- function(status_rv, handle_rv, log_rv, results_rv) {
      h <- handle_rv()
      shiny::req(!is.null(h), identical(status_rv(), "RUNNING"))
      shiny::invalidateLater(750, session)

      new_out <- tryCatch(h$read_output_lines(), error = function(e) character(0))
      new_err <- tryCatch(h$read_error_lines(), error = function(e) character(0))
      if (length(c(new_out, new_err)) > 0) log_rv(c(log_rv(), new_out, new_err))

      if (!h$is_alive()) {
        ok <- identical(h$get_exit_status(), 0L)
        status_rv(if (ok) "SUCCESS" else "ERROR")
        if (ok) results_rv(parse_primer_results(log_rv()))
      }
    }

    # "Before" test, against raw input -- triggered by the Test button.
    before_status <- shiny::reactiveVal("IDLE")
    before_handle <- shiny::reactiveVal(NULL)
    before_log <- shiny::reactiveVal(character(0))
    before_results <- shiny::reactiveVal(NULL)

    # Why Test/Run can't start, or NULL if they can.
    blocked_reason <- function(fwd_list, rev_list) {
      if (is.null(cutadapt_invocation())) return("Set up cutadapt first (top of this page).")
      if (nrow(sample_table()) == 0) return("No paired samples detected -- finish Setup first.")
      invalid_primer_msg(fwd_list, rev_list)
    }

    output$primer_feedback <- shiny::renderUI({
      fwd <- parse_primer_list(input$fwd_primer)
      rev <- parse_primer_list(input$rev_primer)
      if (length(fwd) == 0 && length(rev) == 0) return(NULL)
      msg <- invalid_primer_msg(fwd, rev)
      if (is.null(msg)) return(NULL)
      shiny::div(class = "text-danger small mb-1", msg)
    })

    shiny::observeEvent(input$test_primers, {
      shiny::req(!any_running(rv))
      st <- sample_table()
      fwd_list <- parse_primer_list(input$fwd_primer)
      rev_list <- parse_primer_list(input$rev_primer)
      msg <- blocked_reason(fwd_list, rev_list)
      if (!is.null(msg)) return(shiny::showNotification(msg, type = "error"))
      remember_if_custom()
      start_primer_check(before_status, before_handle, before_log, before_results,
                          raw_files_df(rv$project_dir, st), fwd_list, rev_list, input$err_rate)
    })

    shiny::observe(poll_primer_check(before_status, before_handle, before_log, before_results))

    show_all_before <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(input$toggle_test_table, show_all_before(!show_all_before()))
    shiny::observeEvent(before_results(), show_all_before(FALSE))

    output$test_badge <- shiny::renderUI(status_badge(before_status()))
    output$test_results_table <- DT::renderDataTable(primer_results_table(before_results(), show_all_before()))
    output$test_results_toggle_ui <- shiny::renderUI(
      primer_table_toggle_ui(ns, "toggle_test_table", before_results(), show_all_before())
    )

    output$test_results_download_ui <- shiny::renderUI({
      shiny::req(before_results(), nrow(before_results()) > 0)
      shiny::downloadButton(ns("download_test_results"), "Download as CSV", class = "btn-outline-secondary btn-sm mb-2")
    })

    output$download_test_results <- shiny::downloadHandler(
      filename = function() sprintf("primer_test_results_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        utils::write.csv(before_results(), file, row.names = FALSE)
      }
    )

    # "After" check, against trimmed output -- triggered automatically once
    # the real Run finishes successfully, using the same primers it used.
    after_status <- shiny::reactiveVal("IDLE")
    after_handle <- shiny::reactiveVal(NULL)
    after_log <- shiny::reactiveVal(character(0))
    after_results <- shiny::reactiveVal(NULL)
    run_primers <- shiny::reactiveVal(NULL)  # primers + err_rate used by the last real Run

    step_status <- shiny::reactive(rv$status[[step_id]])
    shiny::observeEvent(step_status(), {
      shiny::req(identical(step_status(), "SUCCESS"), !is.null(run_primers()))
      st <- sample_table()
      tdf <- trimmed_files_df(rv$project_dir, st)
      start_primer_check(after_status, after_handle, after_log, after_results,
                          tdf, run_primers()$fwd, run_primers()$rev, run_primers()$err_rate)
    }, ignoreInit = TRUE)

    shiny::observe(poll_primer_check(after_status, after_handle, after_log, after_results))

    output$after_results_section <- shiny::renderUI({
      shiny::req(!identical(after_status(), "IDLE"))
      df <- after_results()
      verdict <- NULL
      if (identical(after_status(), "SUCCESS") && !is.null(df) && nrow(df) > 0) {
        clean <- all(df$r1_count == 0 & df$r2_count == 0)
        verdict <- shiny::div(
          class = paste("alert small mb-2", if (clean) "alert-success" else "alert-warning"),
          if (clean) "Primers fully removed -- count is zero in every sample."
          else "Primer still detected in trimmed output -- consider a higher error rate or double-checking the sequence."
        )
      }
      result_card(shiny::tagList("Primer check after removal ", status_badge(after_status())),
        shiny::div(class = "amf-param-hint mb-2", "The same test, re-run on the trimmed reads -- counts should now be zero."),
        verdict,
        dt_output(ns("after_results_table")),
        shiny::uiOutput(ns("after_results_toggle_ui"))
      )
    })

    show_all_after <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(input$toggle_after_table, show_all_after(!show_all_after()))
    shiny::observeEvent(after_results(), show_all_after(FALSE))

    output$after_results_table <- DT::renderDataTable(primer_results_table(after_results(), show_all_after()))
    output$after_results_toggle_ui <- shiny::renderUI(
      primer_table_toggle_ui(ns, "toggle_after_table", after_results(), show_all_after())
    )

    # Dual-primer removal: strips the 5' primer AND the other primer's
    # revcomp read-through at the 3' end (R1_flags/R2_flags below), run
    # twice per read (-n 2) so a read-through match doesn't hide behind the
    # first pass. Still enforces --discard-untrimmed / -e
    # on top, so reads without a real primer match are dropped. -m 1 drops
    # pairs trimmed down to nothing: zero-length reads crash DADA2's
    # plotQualityProfile (Filter and Trim's preview). The real minimum read
    # length is enforced downstream by Filter and Trim (minLen).
    #
    # fwd_primers/rev_primers may each hold multiple candidate variants
    # (comma-separated in the UI) -- cutadapt accepts repeated -g/-G flags
    # natively and tries each in order, using whichever matches, so this is
    # a straight extension of the single-primer case rather than needing
    # separate cutadapt invocations per variant.
    build_cutadapt_script <- function(cutadapt_cmd, st, fwd_primers, rev_primers, err_rate) {
      trimmed_dir <- trimmed_dir_path(rv$project_dir)
      lines <- c("set -e", sprintf("mkdir -p %s", shQuote(trimmed_dir)))

      fwd_rc <- vapply(fwd_primers, revcomp, character(1))
      rev_rc <- vapply(rev_primers, revcomp, character(1))
      r1_flags <- paste(c(
        sprintf("-g %s", vapply(fwd_primers, shQuote, character(1))),
        sprintf("-a %s", vapply(rev_rc, shQuote, character(1)))
      ), collapse = " ")
      r2_flags <- paste(c(
        sprintf("-G %s", vapply(rev_primers, shQuote, character(1))),
        sprintf("-A %s", vapply(fwd_rc, shQuote, character(1)))
      ), collapse = " ")

      for (i in seq_len(nrow(st))) {
        s <- st$sample[i]
        in1 <- file.path(rv$project_dir, st$fwd[i])
        in2 <- file.path(rv$project_dir, st$rev[i])
        out_paths <- trimmed_fastq_paths(rv$project_dir, s)
        out1 <- out_paths$fwd
        out2 <- out_paths$rev
        lines <- c(
          lines,
          sprintf("echo '--- %s ---'", s),
          sprintf(
            "%s %s %s -n 2 --discard-untrimmed -m 1 -e %s -o %s -p %s %s %s",
            cutadapt_cmd, r1_flags, r2_flags, err_rate,
            shQuote(out1), shQuote(out2), shQuote(in1), shQuote(in2)
          )
        )
      }
      lines <- c(lines, "echo 'All samples trimmed.'")
      paste(lines, collapse = "\n")
    }

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv))
      st <- sample_table()
      fwd_list <- parse_primer_list(input$fwd_primer)
      rev_list <- parse_primer_list(input$rev_primer)
      msg <- blocked_reason(fwd_list, rev_list)
      if (!is.null(msg)) return(shiny::showNotification(msg, type = "error"))
      remember_if_custom()

      run_primers(list(fwd = fwd_list, rev = rev_list, err_rate = input$err_rate))
      script <- build_cutadapt_script(cutadapt_invocation(), st, fwd_list, rev_list, input$err_rate)
      launch_step(rv, step_id, "bash", c("-c", script))
    })

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    mod_logpanel_server("log", shiny::reactive(rv$log[[step_id]] %||% character(0)))
  })
}
