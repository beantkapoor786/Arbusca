# Per-project directory layout, shared across pipeline steps: each stage
# writes to <project>/NN_stage/.

# Sample table (sample, full fwd path, full rev path) for the raw read pairs
# sitting directly in project_dir, as named by mod_setup's sample_table()
# (basenames only) -- before any trimming or filtering has happened.
raw_files_df <- function(project_dir, st) {
  fwd <- file.path(project_dir, st$fwd)
  rev <- file.path(project_dir, st$rev)
  keep <- fs::file_exists(fwd) & fs::file_exists(rev)
  data.frame(sample = st$sample[keep], fwd = fwd[keep], rev = rev[keep], stringsAsFactors = FALSE)
}

filtn_dir_path <- function(project_dir) file.path(project_dir, "01_filtN")

filtn_fastq_paths <- function(project_dir, sample) {
  dir <- filtn_dir_path(project_dir)
  list(
    fwd = file.path(dir, paste0(sample, "_R1.filtN.fastq.gz")),
    rev = file.path(dir, paste0(sample, "_R2.filtN.fastq.gz"))
  )
}

# Sample table restricted to samples whose N-filtered pair exists on disk.
filtn_files_df <- function(project_dir, st) {
  paths <- filtn_fastq_paths(project_dir, st$sample)
  keep <- fs::file_exists(paths$fwd) & fs::file_exists(paths$rev)
  data.frame(sample = st$sample[keep], fwd = paths$fwd[keep], rev = paths$rev[keep], stringsAsFactors = FALSE)
}

trimmed_dir_path <- function(project_dir) file.path(project_dir, "02_trimmed")

trimmed_fastq_paths <- function(project_dir, sample) {
  dir <- trimmed_dir_path(project_dir)
  list(
    fwd = file.path(dir, paste0(sample, "_R1.trimmed.fastq.gz")),
    rev = file.path(dir, paste0(sample, "_R2.trimmed.fastq.gz"))
  )
}

# Sample table (sample, full fwd path, full rev path) restricted to samples
# whose trimmed pair actually exists on disk.
trimmed_files_df <- function(project_dir, st) {
  paths <- trimmed_fastq_paths(project_dir, st$sample)
  keep <- fs::file_exists(paths$fwd) & fs::file_exists(paths$rev)
  data.frame(sample = st$sample[keep], fwd = paths$fwd[keep], rev = paths$rev[keep], stringsAsFactors = FALSE)
}

filtered_dir_path <- function(project_dir) file.path(project_dir, "03_filtered")

filtered_fastq_paths <- function(project_dir, sample) {
  dir <- filtered_dir_path(project_dir)
  list(
    fwd = file.path(dir, paste0(sample, "_R1.filt.fastq.gz")),
    rev = file.path(dir, paste0(sample, "_R2.filt.fastq.gz"))
  )
}

filtered_files_df <- function(project_dir, st) {
  paths <- filtered_fastq_paths(project_dir, st$sample)
  keep <- fs::file_exists(paths$fwd) & fs::file_exists(paths$rev)
  data.frame(sample = st$sample[keep], fwd = paths$fwd[keep], rev = paths$rev[keep], stringsAsFactors = FALSE)
}

denoise_dir_path <- function(project_dir) file.path(project_dir, "04_denoise")

taxonomy_dir_path <- function(project_dir) file.path(project_dir, "05_taxonomy")

# Counts FASTQ records in a (optionally gzipped) file, for reconstructing a
# read-tracking table when resuming a step whose original in-memory result
# wasn't persisted across sessions.
count_fastq_reads <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  total <- 0
  repeat {
    lines <- readLines(con, n = 4000)
    if (length(lines) == 0) break
    total <- total + length(lines)
  }
  total / 4
}
