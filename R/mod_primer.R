# Step 4 -- Primer Removal (cutadapt). Presets come from inst/primer_presets.yaml
# (each with its literature citation); the dropdown auto-selects the pair
# that matches the detected read length.

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

# Background job (callr) counting, per sample, the reads in which each
# primer orientation occurs -- primerHits() from the DADA2 ITS workflow:
# vcountPattern(fixed = FALSE) honors IUPAC codes and matches anywhere in
# the read, exactly (no mismatches). R1 and R2 are counted independently,
# not assuming a fixed F/R pairing or orientation, which is what lets this
# surface "this primer actually shows up reverse-complemented in R2".
# Returns one row per (sample, primer, read file) with the four orientation
# counts as columns. Must be self-contained: it runs in a bare R session.
primer_check_job <- function(files_df, variants, n_jobs) {
  primerHits <- function(primer, reads) sum(Biostrings::vcountPattern(primer, reads, fixed = FALSE) > 0)
  primers <- unique(variants[c("role", "primer")])
  count_sample <- function(i) {
    reads <- list(
      "Forward reads" = ShortRead::sread(ShortRead::readFastq(files_df$fwd[i])),
      "Reverse reads" = ShortRead::sread(ShortRead::readFastq(files_df$rev[i]))
    )
    rows <- list()
    for (j in seq_len(nrow(primers))) {
      v <- variants[variants$role == primers$role[j] & variants$primer == primers$primer[j], ]
      for (rn in names(reads)) {
        hits <- vapply(v$sequence, primerHits, numeric(1), reads = reads[[rn]])
        rows[[length(rows) + 1]] <- data.frame(
          sample = files_df$sample[i], role = primers$role[j], primer = primers$primer[j],
          reads = rn, t(stats::setNames(hits, v$orientation)),
          stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows)
  }
  do.call(rbind, parallel::mclapply(seq_len(nrow(files_df)), count_sample, mc.cores = n_jobs))
}

PRIMER_ORIENTATIONS <- c("Forward", "Complement", "Reverse", "RevComp")

PRIMER_TABLE_PREVIEW_ROWS <- 10

# show_all = FALSE truncates to the first PRIMER_TABLE_PREVIEW_ROWS rows, so
# a long sample list doesn't turn the page into a scroll marathon; callers
# pair this with a "Show full table" button (see primer_table_toggle_ui()).
primer_results_table <- function(df, show_all = FALSE) {
  if (is.null(df) || nrow(df) == 0) {
    empty_cols <- data.frame(Sample = character(0), Primer = character(0), Reads = character(0), check.names = FALSE)
    for (o in PRIMER_ORIENTATIONS) empty_cols[[o]] <- numeric(0)
    return(DT::datatable(empty_cols, rownames = FALSE, options = list(dom = "t")))
  }
  out <- data.frame(
    Sample = df$sample, Primer = sprintf("%s (%s)", df$primer, df$role), Reads = df$reads,
    df[PRIMER_ORIENTATIONS],
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
    title = "4. Primer Removal",
    description = "Trims the PCR primer sequences off the ends of every read with cutadapt, so primer bases are not mistaken for biological variation. Test your primers on the N-filtered reads first if you are unsure which were used.",
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
          shiny::span(class = "amf-param-hint m-0", "Counts reads containing each primer orientation (exact match, IUPAC-aware), without writing any output.")
        ),
        dt_output(ns("test_results_table")),
        shiny::uiOutput(ns("test_results_toggle_ui")),
        shiny::uiOutput(ns("test_results_download_ui"))
      ),
      result_card("Remove primers with cutadapt",
        run_row(ns, "Run Cutadapt"),
        shiny::uiOutput(ns("last_run_caption")),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("log"))
      ),
      shiny::uiOutput(ns("after_results_section"))
    )
  )
}

# --- cutadapt discovery and self-install ---

# An app launched from RStudio or Finder doesn't inherit the login shell's
# PATH, so Sys.which() misses Homebrew/conda installs; look here too.
COMMON_BIN_DIRS <- c(
  "/opt/homebrew/bin", "/usr/local/bin",
  file.path(path.expand("~"), c("miniforge3", "mambaforge", "miniconda3", "anaconda3"), "bin"),
  "/usr/bin"
)

# The app's own install lives in a private venv, so pip never touches (or is
# refused by, per PEP 668) the user's own Python.
cutadapt_venv_dir <- function() file.path(path.expand("~"), ".arbusca", "cutadapt-venv")

find_executables <- function(name) {
  cands <- c(Sys.which(name), file.path(COMMON_BIN_DIRS, name))
  unique(cands[nzchar(cands) & file.exists(cands)])
}

validate_cutadapt <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(NULL)
  tryCatch({
    out <- system2(path, "--version", stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0) return(NULL)
    trimws(out[1])
  }, error = function(e) NULL)
}

# First working cutadapt among: the saved path, the app's own venv, then PATH
# and common install dirs. Returns list(path, version) or NULL.
locate_cutadapt <- function(saved_path = NULL) {
  cands <- c(saved_path, file.path(cutadapt_venv_dir(), "bin", "cutadapt"), find_executables("cutadapt"))
  for (p in unique(cands)) {
    v <- validate_cutadapt(p)
    if (!is.null(v)) return(list(path = p, version = v))
  }
  NULL
}

# A python3 that can create a venv (Debian/Ubuntu ship python3 without
# ensurepip unless python3-venv is installed). NULL if none found.
find_python3 <- function() {
  for (p in find_executables("python3")) {
    ok <- tryCatch({
      out <- system2(p, c("-c", shQuote("import venv, ensurepip")), stdout = TRUE, stderr = TRUE)
      is.null(attr(out, "status")) || attr(out, "status") == 0
    }, error = function(e) FALSE)
    if (ok) return(p)
  }
  NULL
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

    # --- cutadapt readiness: find it silently, else ask once and validate or install ---

    cutadapt_invocation <- shiny::reactiveVal(NULL)  # shell-ready command prefix once resolved
    cutadapt_version <- shiny::reactiveVal(NULL)

    shiny::isolate({
      found <- locate_cutadapt(read_app_config()$cutadapt_path)
      if (!is.null(found)) {
        cutadapt_invocation(shQuote(found$path))
        cutadapt_version(found$version)
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
      py <- find_python3()
      if (is.null(py)) {
        install_log("No usable Python 3 was found, so cutadapt can't be installed automatically. Install Python 3 (https://www.python.org/downloads/) and click Install again, or choose \"Yes\" above and point to an existing cutadapt.")
        install_status("ERROR")
        return()
      }
      venv <- cutadapt_venv_dir()
      script <- paste(
        "set -e",
        sprintf("echo 'Creating a private Python environment in %s'", venv),
        sprintf("%s -m venv --clear %s", shQuote(py), shQuote(venv)),
        sprintf("%s -m pip install --disable-pip-version-check cutadapt", shQuote(file.path(venv, "bin", "python"))),
        sep = "\n"
      )
      h <- processx::process$new("bash", c("-c", script), stdout = "|", stderr = "|")
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
        exe <- file.path(cutadapt_venv_dir(), "bin", "cutadapt")
        v <- if (identical(h$get_exit_status(), 0L)) validate_cutadapt(exe)
        if (!is.null(v)) {
          # No config write needed: locate_cutadapt() checks the venv on every launch.
          cutadapt_invocation(shQuote(exe))
          cutadapt_version(v)
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
              shiny::div(class = "text-muted small mb-2", sprintf("Installs cutadapt into a private Python environment in %s (about a minute).", cutadapt_venv_dir())),
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

    # Samples run in parallel, leaving one core free so the app itself stays
    # responsive.
    start_primer_check <- function(status_rv, handle_rv, results_rv, files_df, fwd_primer, rev_primer) {
      if (nrow(files_df) == 0) return(invisible())
      cores <- parallel::detectCores()
      n_jobs <- if (is.na(cores)) 1 else max(1, cores - 1)
      variants <- rbind(
        expand_primer_orientations(fwd_primer, "Forward"),
        expand_primer_orientations(rev_primer, "Reverse")
      )
      h <- callr::r_bg(primer_check_job, args = list(files_df = files_df, variants = variants, n_jobs = n_jobs),
                       stdout = "|", stderr = "|", supervise = TRUE)
      handle_rv(h)
      status_rv("RUNNING")
      results_rv(NULL)
    }

    poll_primer_check <- function(status_rv, handle_rv, results_rv) {
      h <- handle_rv()
      shiny::req(!is.null(h), identical(status_rv(), "RUNNING"))
      shiny::invalidateLater(750, session)
      if (!h$is_alive()) {
        result <- tryCatch(h$get_result(), error = function(e) e)
        if (inherits(result, "error")) {
          shiny::showNotification(paste("Primer check failed:", conditionMessage(result)), type = "error")
          status_rv("ERROR")
        } else {
          results_rv(result)
          status_rv("SUCCESS")
        }
      }
    }

    # "Before" test, against the N-filtered reads -- triggered by the Test button.
    before_status <- shiny::reactiveVal("IDLE")
    before_handle <- shiny::reactiveVal(NULL)
    before_results <- shiny::reactiveVal(NULL)

    # Step 3's output, re-read once it finishes.
    filtn_files <- shiny::reactive({
      shiny::req(rv$project_dir)
      rv$status$filtn
      filtn_files_df(rv$project_dir, sample_table())
    })

    # Why Test (or, with need_cutadapt, Run) can't start, or NULL if they can.
    blocked_reason <- function(fwd_list, rev_list, need_cutadapt = FALSE) {
      if (need_cutadapt && is.null(cutadapt_invocation())) return("Set up cutadapt first (top of this page).")
      if (nrow(sample_table()) == 0) return("No paired samples detected -- finish Setup first.")
      if (nrow(filtn_files()) < nrow(sample_table())) return("No N-filtered reads yet -- finish Remove Ambiguous Bases first.")
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
      fwd_list <- parse_primer_list(input$fwd_primer)
      rev_list <- parse_primer_list(input$rev_primer)
      msg <- blocked_reason(fwd_list, rev_list)
      if (!is.null(msg)) return(shiny::showNotification(msg, type = "error"))
      remember_if_custom()
      start_primer_check(before_status, before_handle, before_results, filtn_files(), fwd_list, rev_list)
    })

    shiny::observe(poll_primer_check(before_status, before_handle, before_results))

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
    after_results <- shiny::reactiveVal(NULL)
    run_primers <- shiny::reactiveVal(NULL)  # primers + err_rate used by the last real Run
    run_record_saved <- shiny::reactiveVal(0)  # bumped when primer_run.yaml is (re)written

    step_status <- shiny::reactive(rv$status[[step_id]])
    shiny::observeEvent(step_status(), {
      shiny::req(identical(step_status(), "SUCCESS"), !is.null(run_primers()))
      yaml::write_yaml(c(run_primers(), list(run_at = format(Sys.time(), "%Y-%m-%d %H:%M"))),
                       file.path(trimmed_dir_path(rv$project_dir), "primer_run.yaml"))
      run_record_saved(run_record_saved() + 1)
      st <- sample_table()
      tdf <- trimmed_files_df(rv$project_dir, st)
      start_primer_check(after_status, after_handle, after_results, tdf, run_primers()$fwd, run_primers()$rev)
    }, ignoreInit = TRUE)

    shiny::observe(poll_primer_check(after_status, after_handle, after_results))

    # Which primers produced the trimmed files currently on disk -- recorded
    # in primer_run.yaml by the observer above, so it survives reopening the
    # project (when the primer fields just show the auto-detected preset).
    output$last_run_caption <- shiny::renderUI({
      shiny::req(identical(step_status(), "SUCCESS"), rv$project_dir)
      run_record_saved()
      path <- file.path(trimmed_dir_path(rv$project_dir), "primer_run.yaml")
      if (!file.exists(path)) {
        return(shiny::div(class = "alert alert-warning small mt-2 mb-0",
          "No record of which primers produced the trimmed reads on disk (they predate run tracking). Re-run to be sure."))
      }
      run <- yaml::read_yaml(path)
      shiny::div(class = "alert alert-info small mt-2 mb-0",
        shiny::strong(sprintf("Trimmed reads on disk were made on %s with:", run$run_at)),
        shiny::tags$br(), "Forward: ", shiny::tags$code(paste(run$fwd, collapse = ", ")),
        shiny::tags$br(), "Reverse: ", shiny::tags$code(paste(run$rev, collapse = ", ")),
        shiny::tags$br(), sprintf("Error rate: %s", run$err_rate))
    })

    output$after_results_section <- shiny::renderUI({
      shiny::req(!identical(after_status(), "IDLE"))
      df <- after_results()
      verdict <- NULL
      if (identical(after_status(), "SUCCESS") && !is.null(df) && nrow(df) > 0) {
        clean <- all(df[PRIMER_ORIENTATIONS] == 0)
        verdict <- shiny::div(
          class = paste("alert small mb-2", if (clean) "alert-success" else "alert-warning"),
          if (clean) "Primers fully removed -- count is zero in every sample."
          else "Primer still detected in trimmed output -- consider a higher error rate or double-checking the sequence."
        )
      }
      result_card(shiny::tagList("Primer check after removal ", status_badge(after_status())),
        shiny::div(class = "amf-param-hint mb-2", "The same test, re-run on the trimmed reads. Counts should now be zero."),
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
    build_cutadapt_script <- function(cutadapt_cmd, files_df, fwd_primers, rev_primers, err_rate) {
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

      for (i in seq_len(nrow(files_df))) {
        s <- files_df$sample[i]
        in1 <- files_df$fwd[i]
        in2 <- files_df$rev[i]
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
      fwd_list <- parse_primer_list(input$fwd_primer)
      rev_list <- parse_primer_list(input$rev_primer)
      msg <- blocked_reason(fwd_list, rev_list, need_cutadapt = TRUE)
      if (!is.null(msg)) return(shiny::showNotification(msg, type = "error"))
      remember_if_custom()

      run_primers(list(fwd = fwd_list, rev = rev_list, err_rate = input$err_rate))
      unlink(file.path(trimmed_dir_path(rv$project_dir), "primer_run.yaml"))
      script <- build_cutadapt_script(cutadapt_invocation(), filtn_files(), fwd_list, rev_list, input$err_rate)
      launch_step(rv, step_id, "bash", c("-c", script))
    })

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    mod_logpanel_server("log", shiny::reactive(rv$log[[step_id]] %||% character(0)))
  })
}
