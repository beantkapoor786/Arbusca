# Pipeline step order + FSM gating helpers.
# rv$status[[step_id]] is one of IDLE, QUEUED, RUNNING, SUCCESS, ERROR.

STEP_ORDER <- c("filtn", "primer", "qc", "denoise", "merge", "chimera", "taxonomy", "output")

# TRUE if step_id's upstream dependency (if any) is SUCCESS, so its Run
# button should be enabled. The first step has no upstream.
upstream_ok <- function(rv, step_id) {
  idx <- match(step_id, STEP_ORDER)
  if (is.na(idx) || idx == 1) return(TRUE)
  identical(rv$status[[STEP_ORDER[idx - 1]]], "SUCCESS")
}

# On-disk output(s) owned by a single step -- used by reset_downstream() to
# delete exactly the stale files a re-run invalidates. denoise/merge/chimera
# all write into the shared 04_denoise/ dir, so each lists only its own
# file(s) there rather than the whole directory -- e.g. re-running "merge"
# alone must not delete denoise's still-valid errF.rds/dadaF.rds.
step_output_paths <- function(project_dir, step_id) {
  denoise_dir <- denoise_dir_path(project_dir)
  switch(step_id,
    primer   = trimmed_dir_path(project_dir),
    qc       = filtered_dir_path(project_dir),
    denoise  = file.path(denoise_dir, c("errF.rds", "errR.rds", "dadaF.rds", "dadaR.rds")),
    merge    = file.path(denoise_dir, "seqtab.rds"),
    chimera  = file.path(denoise_dir, c("seqtab_nochim.rds", "track.rds")),
    taxonomy = taxonomy_dir_path(project_dir),
    output   = file.path(project_dir, "06_output"),
    character(0)
  )
}

delete_step_outputs <- function(project_dir, step_id) {
  for (p in step_output_paths(project_dir, step_id)) {
    tryCatch({
      if (fs::is_dir(p)) fs::dir_delete(p)
      else if (fs::file_exists(p)) fs::file_delete(p)
    }, error = function(e) NULL)
  }
}

# Every step after step_id in STEP_ORDER back to IDLE, artifacts/log
# cleared, and its on-disk output deleted -- their outputs were built from
# step_id's *previous* run, so they're stale the moment step_id is launched
# again with (possibly) different parameters. Deleting the files (not just
# resetting rv$status) matters because scan_checkpoints() marks a step
# SUCCESS purely by checking whether its expected output exists on disk --
# leaving stale files behind would let a later project-dir reselect silently
# resurrect a SUCCESS badge for output that no longer matches the new
# upstream run. Called from launch_step()/launch_step_callr() so every real
# pipeline run gets this for free, not just re-runs after a param change (a
# first run has nothing downstream yet to reset or delete anyway).
reset_downstream <- function(rv, step_id) {
  idx <- match(step_id, STEP_ORDER)
  if (is.na(idx)) return(invisible())
  downstream <- STEP_ORDER[-seq_len(idx)]
  for (s in downstream) {
    rv$status[[s]] <- "IDLE"
    rv$artifacts[[s]] <- NULL
    rv$log[[s]] <- character(0)
    if (!is.null(rv$project_dir)) delete_step_outputs(rv$project_dir, s)
  }
  invisible()
}

# Analysis steps (Figures/PERMANOVA/ANCOM-BC2) read 06_output/ directly, so
# any rebuild or re-transform of the phyloseq object makes their on-screen
# results stale. mod_output calls this after either; the analysis modules
# watch rv$phyloseq_rev and clear what they're showing.
invalidate_analyses <- function(rv) {
  rv$phyloseq_rev <- (shiny::isolate(rv$phyloseq_rev) %||% 0) + 1
  for (s in c("permanova", "ancombc")) {
    rv$status[[s]] <- "IDLE"
    rv$artifacts[[s]] <- NULL
    rv$log[[s]] <- character(0)
  }
  invisible()
}

# TRUE while any step has a live process handle -- used to gate the global
# polling loop and to disable Run buttons on other steps.
any_running <- function(rv) {
  !is.null(rv$running_step)
}

# Checkpoint scan (DESIGN.md section 4.4): a step counts as SUCCESS only if
# its *complete* expected output set is already on disk for every sample --
# a partially-written step (some samples missing) is left alone so the user
# resumes from the first genuinely incomplete step, not a half-finished one.
# Only ever upgrades IDLE -> SUCCESS; never touches a RUNNING/ERROR/already-
# SUCCESS status, so this is safe to call any time project_dir changes.
scan_checkpoints <- function(rv, st) {
  if (is.null(rv$project_dir) || nrow(st) == 0) return(invisible())

  filtn <- filtn_files_df(rv$project_dir, st)
  if (identical(rv$status$filtn, "IDLE") && nrow(filtn) == nrow(st)) {
    rv$status$filtn <- "SUCCESS"
    if (is.null(rv$artifacts$filtn)) {
      raw <- raw_files_df(rv$project_dir, st)
      rv$artifacts$filtn <- data.frame(
        sample = st$sample,
        reads_in = vapply(raw$fwd, count_fastq_reads, numeric(1)),
        reads_out = vapply(filtn$fwd, count_fastq_reads, numeric(1)),
        row.names = NULL
      )
    }
  }

  trimmed <- trimmed_files_df(rv$project_dir, st)
  if (identical(rv$status$primer, "IDLE") && nrow(trimmed) == nrow(st) && nrow(filtn) == nrow(st)) {
    rv$status$primer <- "SUCCESS"
  }

  filtered <- filtered_files_df(rv$project_dir, st)
  if (identical(rv$status$qc, "IDLE") && nrow(filtered) == nrow(st) && nrow(trimmed) == nrow(st)) {
    rv$status$qc <- "SUCCESS"
    if (is.null(rv$artifacts$qc)) {
      reads_in <- vapply(trimmed$fwd, count_fastq_reads, numeric(1))
      reads_out <- vapply(filtered$fwd, count_fastq_reads, numeric(1))
      rv$artifacts$qc <- data.frame(sample = st$sample, reads_in = reads_in, reads_out = reads_out, row.names = NULL)
    }
  }

  # errF.rds/errR.rds (6a) are checked directly by mod_denoise's own
  # "resume mid-step" logic, not here -- this only covers 6b's completion,
  # which is what actually gates the pipeline/Proceed button.
  denoise_dir <- denoise_dir_path(rv$project_dir)
  dada_files <- file.path(denoise_dir, c("dadaF.rds", "dadaR.rds"))
  if (identical(rv$status$denoise, "IDLE") && nrow(filtered) == nrow(st) && all(fs::file_exists(dada_files))) {
    rv$status$denoise <- "SUCCESS"
    if (is.null(rv$artifacts$denoise)) {
      getN <- function(x) sum(dada2::getUniques(x))
      dadaFs <- readRDS(dada_files[1])
      dadaRs <- readRDS(dada_files[2])
      rv$artifacts$denoise <- list(summary = data.frame(
        sample = st$sample,
        denoisedF = as.numeric(sapply(dadaFs, getN)[st$sample]),
        denoisedR = as.numeric(sapply(dadaRs, getN)[st$sample]),
        row.names = NULL
      ))
    }
  }

  seqtab_path <- file.path(denoise_dir, "seqtab.rds")
  if (identical(rv$status$merge, "IDLE") && identical(rv$status$denoise, "SUCCESS") && fs::file_exists(seqtab_path)) {
    rv$status$merge <- "SUCCESS"
    if (is.null(rv$artifacts$merge)) {
      seqtab <- readRDS(seqtab_path)
      merged <- rowSums(seqtab)
      rv$artifacts$merge <- list(
        merged = data.frame(sample = st$sample, merged = as.numeric(merged[st$sample]), row.names = NULL),
        asv_lengths = nchar(colnames(seqtab)),
        n_asvs = ncol(seqtab)
      )
    }
  }

  chimera_files <- file.path(denoise_dir, c("seqtab_nochim.rds", "track.rds"))
  if (identical(rv$status$chimera, "IDLE") && identical(rv$status$merge, "SUCCESS") && all(fs::file_exists(chimera_files))) {
    rv$status$chimera <- "SUCCESS"
    if (is.null(rv$artifacts$chimera)) {
      seqtab <- readRDS(seqtab_path)
      seqtab_nochim <- readRDS(chimera_files[1])
      rv$artifacts$chimera <- list(
        track = readRDS(chimera_files[2]),
        asv_lengths = nchar(colnames(seqtab_nochim)),
        n_asvs = ncol(seqtab_nochim),
        n_asvs_before = ncol(seqtab)
      )
    }
  }

  vt_path <- file.path(taxonomy_dir_path(rv$project_dir), "vt_assignments.rds")
  if (identical(rv$status$taxonomy, "IDLE") && identical(rv$status$chimera, "SUCCESS") && fs::file_exists(vt_path)) {
    rv$status$taxonomy <- "SUCCESS"
    if (is.null(rv$artifacts$taxonomy)) {
      rv$artifacts$taxonomy <- readRDS(vt_path)
    }
  }

  phyloseq_path <- file.path(rv$project_dir, "06_output", "phyloseq.rds")
  if (identical(rv$status$output, "IDLE") && identical(rv$status$taxonomy, "SUCCESS") &&
      fs::file_exists(phyloseq_path) && requireNamespace("phyloseq", quietly = TRUE)) {
    rv$status$output <- "SUCCESS"
    if (is.null(rv$artifacts$output)) {
      ps <- readRDS(phyloseq_path)
      rv$artifacts$output <- list(phyloseq = ps, n_samples = phyloseq::nsamples(ps), n_taxa = phyloseq::ntaxa(ps), unmatched_samples = character(0))
    }
  }
}

# First pipeline step (STEP_ORDER) that is not yet SUCCESS -- where a resumed
# session should land. Once the whole pipeline is done, that's Figures.
# (Quality Check and the analysis steps carry no rv$status, so they're
# never a resume point.)
first_incomplete_step <- function(rv) {
  for (step_id in STEP_ORDER) {
    if (!identical(rv$status[[step_id]], "SUCCESS")) return(step_id)
  }
  "figures"
}
