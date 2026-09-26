# Installs R packages + checks external binaries needed by Arbusca.
# Run once: Rscript install_deps.R

cran_pkgs <- c(
  "shiny", "bslib", "bsicons", "DT", "plotly", "ggplot2",
  "processx", "callr", "fs", "yaml", "jsonlite", "readr", "vegan"
)

bioc_pkgs <- c("dada2", "Biostrings", "ShortRead", "phyloseq", "ANCOMBC")

installed <- rownames(installed.packages())

missing_cran <- setdiff(cran_pkgs, installed)
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

# ANCOMBC 2.12 imports CVXR::solve, which CVXR >= 1.9 no longer exports --
# pin the last 1.0.x release (from the CRAN archive) until ANCOMBC catches up.
if (!"CVXR" %in% installed || packageVersion("CVXR") >= "1.9") {
  install.packages(c("ECOSolveR", "scs", "osqp", "gmp", "Rmpfr", "bit64"))
  install.packages("https://cran.r-project.org/src/contrib/Archive/CVXR/CVXR_1.0-15.tar.gz", repos = NULL, type = "source")
}

missing_bioc <- setdiff(bioc_pkgs, installed)
if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

cat("\nR packages OK.\n")

check_binary <- function(name, hint) {
  path <- Sys.which(name)
  if (nzchar(path)) {
    cat(sprintf("[ok] %s -> %s\n", name, path))
  } else {
    cat(sprintf("[missing] %s -- %s\n", name, hint))
  }
}

cat("\nExternal binaries:\n")
check_binary("cutadapt", "install via 'pip install cutadapt' or 'conda install -c bioconda cutadapt'")
check_binary("makeblastdb", "install via 'conda install -c bioconda blast' or from NCBI BLAST+ downloads")
check_binary("blastn", "install via 'conda install -c bioconda blast' or from NCBI BLAST+ downloads")
