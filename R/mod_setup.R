# Dependency readiness: install/load what the app can fix itself, silently.
# Nothing about this is shown to the user unless something is broken that
# they -- not the app -- have to fix (a missing system binary).

R_PKGS <- c("shiny", "bslib", "DT", "processx", "callr",
            "dada2", "Biostrings", "ShortRead", "fs", "yaml", "jsonlite")
BIOC_PKGS <- c("dada2", "Biostrings", "ShortRead")

BIN_SPECS <- list(
  cutadapt    = "install via 'pip install cutadapt' or 'conda install -c bioconda cutadapt'",
  makeblastdb = "install via 'conda install -c bioconda blast'",
  blastn      = "install via 'conda install -c bioconda blast'"
)

# Installs any missing R packages (quietly). Called once at app startup.
# Returns character() of packages that are still missing (should be empty in
# the normal case). Only checks that each package is installed -- loading
# dada2/Biostrings/ShortRead here would add seconds to every app launch, and
# the steps that need them load them on first use anyway.
ensure_r_packages <- function(pkgs = R_PKGS) {
  installed <- rownames(installed.packages())
  missing <- setdiff(pkgs, installed)

  if (length(missing) > 0) {
    cran_missing <- setdiff(missing, BIOC_PKGS)
    bioc_missing <- intersect(missing, BIOC_PKGS)

    if (length(cran_missing) > 0) {
      try(suppressMessages(install.packages(cran_missing, quiet = TRUE)), silent = TRUE)
    }
    if (length(bioc_missing) > 0) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        try(suppressMessages(install.packages("BiocManager", quiet = TRUE)), silent = TRUE)
      }
      try(suppressMessages(BiocManager::install(bioc_missing, update = FALSE, ask = FALSE, quiet = TRUE)), silent = TRUE)
    }
  }

  pkgs[!nzchar(vapply(pkgs, function(p) system.file(package = p), character(1)))]
}

# Binaries the app cannot install itself -- system-level tools. Returns a
# data.frame of the ones that are missing, with an install hint each.
missing_binaries <- function() {
  rows <- lapply(names(BIN_SPECS), function(b) {
    if (nzchar(Sys.which(b))) return(NULL)
    data.frame(name = b, hint = BIN_SPECS[[b]], stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

mod_setup_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "1. Setup & Raw Reads",
    description = "Choose the folder that holds your paired-end FASTQ files. The app detects the samples, and all results are written into subfolders of this project directory.",
    params = shiny::tagList(
      shiny::uiOutput(ns("readiness_alert")),
      param_group("Project directory",
        shiny::uiOutput(ns("last_dir_shortcut")),
        shiny::div(
          class = "amf-inline-row",
          shiny::textInput(ns("project_dir"), NULL, placeholder = "/path/to/project", width = "100%"),
          shiny::actionButton(ns("browse_dir"), "Browse", icon = bsicons::bs_icon("folder2-open"), class = "btn-primary")
        ),
        shiny::uiOutput(ns("project_status"))
      )
    ),
    results = shiny::uiOutput(ns("files_section"))
  )
}

mod_setup_server <- function(id, rv) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    missing_pkgs <- shiny::isolate(ensure_r_packages())
    # Binaries (cutadapt, blastn, makeblastdb) are checked lazily by the
    # specific step that needs them, not here -- no point warning about
    # BLAST+ before the user has even reached the taxonomy step.

    output$readiness_alert <- shiny::renderUI({
      if (length(missing_pkgs) == 0) return(NULL)
      shiny::div(
        class = "alert alert-warning small mb-3",
        sprintf("R package(s) could not be installed automatically: %s.", paste(missing_pkgs, collapse = ", "))
      )
    })

    set_project_dir <- function(dir) {
      shiny::req(nzchar(dir))
      if (!fs::dir_exists(dir)) {
        ok <- tryCatch({ fs::dir_create(dir); TRUE }, error = function(e) FALSE)
        if (!ok) {
          rv$project_dir <- NULL
          output$project_status <- shiny::renderUI(
            shiny::div(class = "text-danger small mt-1", "Could not create that directory.")
          )
          return(invisible())
        }
      }
      rv$project_dir <- fs::path_abs(dir)
      output$project_status <- shiny::renderUI(
        shiny::div(class = "text-success small mt-1", paste("Using:", rv$project_dir))
      )
      write_app_config(utils::modifyList(read_app_config(), list(last_project_dir = rv$project_dir)))
    }

    # --- last-used directory shortcut ---

    last_dir <- shiny::isolate({
      cfg <- read_app_config()
      if (!is.null(cfg$last_project_dir) && fs::dir_exists(cfg$last_project_dir)) cfg$last_project_dir else NULL
    })
    last_dir_dismissed <- shiny::reactiveVal(FALSE)

    output$last_dir_shortcut <- shiny::renderUI({
      shiny::req(!is.null(last_dir), !last_dir_dismissed(), is.null(rv$project_dir))
      shiny::div(
        class = "alert alert-light border d-flex justify-content-between align-items-center py-2 mb-2",
        shiny::span(
          shiny::span(class = "text-muted small", "Last used: "),
          shiny::span(class = "font-monospace small", last_dir)
        ),
        shiny::actionButton(ns("use_last_dir"), "Use this", class = "btn-outline-primary btn-sm")
      )
    })

    # Typed/pasted paths: applied once typing pauses, and only for a folder
    # that already exists (so half-typed paths don't create stray folders --
    # Browse remains the way to pick or create one).
    typed_dir <- shiny::debounce(shiny::reactive(trimws(input$project_dir %||% "")), 800)
    shiny::observeEvent(typed_dir(), {
      dir <- typed_dir()
      shiny::req(nzchar(dir), fs::dir_exists(dir))
      shiny::req(!identical(as.character(fs::path_abs(dir)), as.character(rv$project_dir %||% "")))
      set_project_dir(dir)
    })

    shiny::observeEvent(input$use_last_dir, {
      shiny::updateTextInput(session, "project_dir", value = last_dir)
      set_project_dir(last_dir)
      last_dir_dismissed(TRUE)
    })

    # --- custom directory browser (no dependency on shinyFiles' modal JS) ---

    browse_path <- shiny::reactiveVal(NULL)

    list_subdirs <- function(path) {
      entries <- tryCatch(fs::dir_ls(path, type = "directory"), error = function(e) character(0))
      sort(fs::path_file(entries))
    }

    # Quick-access shortcuts, like a Finder/Explorer sidebar: Home + common
    # folders, plus every cloud-storage mount under ~/Library/CloudStorage
    # (OneDrive, iCloud Drive, Dropbox, ...) which otherwise sit three
    # folders deep and are easy to miss when browsing from Home.
    quick_access_paths <- function() {
      candidates <- c(
        Home = fs::path_home(),
        Desktop = fs::path_home("Desktop"),
        Documents = fs::path_home("Documents")
      )
      cloud_dir <- fs::path_home("Library", "CloudStorage")
      if (fs::dir_exists(cloud_dir)) {
        clouds <- tryCatch(fs::dir_ls(cloud_dir, type = "directory"), error = function(e) character(0))
        if (length(clouds) > 0) {
          names(clouds) <- fs::path_file(clouds)
          candidates <- c(candidates, clouds)
        }
      }
      candidates[fs::dir_exists(candidates)]
    }

    render_browser_ui <- function(path) {
      subdirs <- list_subdirs(path)
      quick <- quick_access_paths()
      shiny::tagList(
        shiny::div(
          class = "d-flex flex-wrap gap-1 mb-2",
          lapply(names(quick), function(nm) {
            shiny::actionLink(
              inputId = ns(paste0("quick_", make.names(nm))),
              label = nm,
              class = "btn btn-outline-secondary btn-sm",
              onclick = sprintf(
                "Shiny.setInputValue('%s', '%s', {priority: 'event'})",
                ns("browse_quick"), htmltools::htmlEscape(quick[[nm]], attribute = TRUE)
              )
            )
          })
        ),
        shiny::div(class = "font-monospace small text-muted mb-2 text-truncate", path),
        shiny::div(
          class = "d-flex gap-2 mb-2",
          shiny::actionButton(ns("browse_up"), "↑ Up", class = "btn-outline-secondary btn-sm"),
          shiny::actionButton(ns("browse_select"), "Select this directory", class = "btn-success btn-sm")
        ),
        if (length(subdirs) == 0) {
          shiny::div(class = "text-muted small fst-italic", "No subdirectories")
        } else {
          shiny::div(
            class = "list-group",
            style = "max-height: 320px; overflow-y: auto;",
            lapply(seq_along(subdirs), function(i) {
              name <- subdirs[[i]]
              shiny::actionLink(
                inputId = ns(paste0("enter_slot_", i)),
                label = shiny::tagList(shiny::icon("folder"), " ", name),
                class = "list-group-item list-group-item-action",
                onclick = sprintf(
                  "Shiny.setInputValue('%s', '%s', {priority: 'event'})",
                  ns("browse_enter"), htmltools::htmlEscape(name, attribute = TRUE)
                )
              )
            })
          )
        }
      )
    }

    shiny::observeEvent(input$browse_dir, {
      start <- if (!is.null(rv$project_dir) && fs::dir_exists(rv$project_dir)) rv$project_dir else fs::path_home()
      browse_path(fs::path_abs(start))
      shiny::showModal(shiny::modalDialog(
        title = "Choose project directory",
        shiny::uiOutput(ns("browser_body")),
        easyClose = TRUE,
        footer = shiny::modalButton("Cancel")
      ))
    })

    output$browser_body <- shiny::renderUI({
      shiny::req(browse_path())
      render_browser_ui(browse_path())
    })

    shiny::observeEvent(input$browse_up, {
      parent <- fs::path_dir(browse_path())
      browse_path(parent)
    })

    shiny::observeEvent(input$browse_enter, {
      new_path <- fs::path(browse_path(), input$browse_enter)
      if (fs::dir_exists(new_path)) browse_path(new_path)
    })

    shiny::observeEvent(input$browse_quick, {
      if (fs::dir_exists(input$browse_quick)) browse_path(fs::path_abs(input$browse_quick))
    })

    shiny::observeEvent(input$browse_select, {
      dir <- browse_path()
      shiny::removeModal()
      shiny::updateTextInput(session, "project_dir", value = dir)
      set_project_dir(dir)
    })

    # --- file listing + delimiter-based sample-name extraction + R1/R2 pairing ---

    files_in_dir <- shiny::reactive({
      shiny::req(rv$project_dir)
      sort(fs::path_file(fs::dir_ls(rv$project_dir, type = c("file", "symlink"))))
    })

    # Sample name = tokens[first_idx:last_idx] of filename split on delim,
    # rejoined with delim. Indices are clamped to the token count so a
    # too-large "last index" (e.g. default 1 on a longer filename) doesn't error.
    extract_sample <- function(filename, delim, first_idx, last_idx) {
      if (!nzchar(delim)) return(filename)
      tokens <- strsplit(filename, delim, fixed = TRUE)[[1]]
      n <- length(tokens)
      first_idx <- max(1, min(first_idx, n))
      last_idx <- max(first_idx, min(last_idx, n))
      paste(tokens[first_idx:last_idx], collapse = delim)
    }

    # The naming settings aren't saved, so when reopening a project that was
    # already run, recover them from the sample names its 02_trimmed/ files
    # were written under -- otherwise the defaults produce different names,
    # no checkpoint matches, and nothing resumes. NULL if there's nothing to
    # match (fresh project) or no candidate reproduces those names exactly.
    guess_sample_naming <- function(project_dir, files) {
      if (!fs::dir_exists(trimmed_dir_path(project_dir))) return(NULL)
      trimmed <- fs::path_file(fs::dir_ls(trimmed_dir_path(project_dir), regexp = "_R1\\.trimmed\\.fastq\\.gz$"))
      if (length(trimmed) == 0) return(NULL)
      target <- sort(sub("_R1\\.trimmed\\.fastq\\.gz$", "", trimmed))
      fwd_files <- files[grepl("(^|[._-])R1([._-]|$)", files)]
      for (delim in c("_", "-", ".")) for (first in 1:3) for (last in first:6) {
        samples <- vapply(fwd_files, extract_sample, character(1), delim = delim, first_idx = first, last_idx = last, USE.NAMES = FALSE)
        if (identical(sort(unique(samples)), target)) return(list(delim = delim, first = first, last = last))
      }
      NULL
    }

    output$files_section <- shiny::renderUI({
      if (is.null(rv$project_dir)) return(NULL)
      naming <- guess_sample_naming(rv$project_dir, files_in_dir()) %||% list(delim = "_", first = 1, last = 1)
      result_card("Samples detected",
        shiny::uiOutput(ns("sample_count")),
        shiny::fluidRow(
          shiny::column(4, shiny::textInput(ns("sample_delim"), "Sample name delimiter", value = naming$delim, width = "100%")),
          shiny::column(4, shiny::numericInput(ns("token_first"), "First index", value = naming$first, min = 1, width = "100%")),
          shiny::column(4, shiny::numericInput(ns("token_last"), "Last index", value = naming$last, min = 1, width = "100%"))
        ),
        shiny::div(
          class = "amf-param-hint mb-2",
          "The sample name is built by cutting each filename at the delimiter into numbered pieces, then keeping pieces First index through Last index. ",
          "E.g. with delimiter _, Sample1_S145_L001_R1_001.fastq.gz splits into Sample1 (1), S145 (2), L001 (3), R1 (4), 001.fastq.gz (5), so indices 1 to 1 give Sample1. ",
          "Forward and reverse files are paired by their R1/R2 tag."
        ),
        dt_output(ns("file_table"))
      )
    })

    output$sample_count <- shiny::renderUI({
      n_files <- length(files_in_dir())
      n <- nrow(sample_table())
      if (n == 0) {
        return(shiny::div(class = "alert alert-warning small mb-2", if (n_files == 0)
          "This folder has no files. Choose the folder that contains your R1/R2 FASTQ files."
          else "No R1/R2 read pairs could be matched. Check the delimiter and index settings below."))
      }
      shiny::div(class = "text-muted small mb-2", sprintf("%d paired sample%s found.", n, if (n == 1) "" else "s"))
    })

    # Sample sheet: one row per sample with its forward/reverse basenames
    # (single match each; multi-/zero-match samples are dropped here and
    # would need tighter delimiter/index settings to resolve). Shared with
    # downstream steps (e.g. primer removal) via the returned reactive.
    sample_table <- shiny::reactive({
      files <- files_in_dir()
      empty <- data.frame(sample = character(0), fwd = character(0), rev = character(0), stringsAsFactors = FALSE)
      if (length(files) == 0) return(empty)

      delim <- input$sample_delim %||% "_"
      first_idx <- input$token_first %||% 1
      last_idx <- input$token_last %||% 1

      samples <- vapply(files, extract_sample, character(1), delim = delim, first_idx = first_idx, last_idx = last_idx, USE.NAMES = FALSE)
      # R1/R2 as a separate token (_R1_001, .R1., -R1.fastq), not any "R1"
      # substring -- a sample named e.g. "PR1-root" must not read as forward.
      is_fwd <- grepl("(^|[._-])R1([._-]|$)", files)
      is_rev <- grepl("(^|[._-])R2([._-]|$)", files)

      rows <- lapply(unique(samples), function(s) {
        idx <- samples == s
        fwd <- files[idx & is_fwd]
        rev <- files[idx & is_rev]
        if (length(fwd) != 1 || length(rev) != 1) return(NULL)
        data.frame(sample = s, fwd = fwd, rev = rev, stringsAsFactors = FALSE)
      })
      rows <- Filter(Negate(is.null), rows)
      if (length(rows) == 0) return(empty)
      do.call(rbind, rows)
    })

    output$file_table <- DT::renderDataTable({
      files <- files_in_dir()
      st <- sample_table()

      if (length(files) == 0) {
        df <- data.frame(`Sample name` = character(0), `Forward file` = character(0), `Reverse file` = character(0), check.names = FALSE)
      } else {
        df <- data.frame(`Sample name` = st$sample, `Forward file` = st$fwd, `Reverse file` = st$rev, check.names = FALSE)
      }
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
    })

    all_ok <- shiny::reactive({
      length(missing_pkgs) == 0
    })

    list(all_ok = all_ok, sample_table = sample_table)
  })
}
