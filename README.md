# Arbusca

A point-and-click Shiny app for processing arbuscular mycorrhizal fungi (AMF) 18S/SSU amplicon data with DADA2: from raw paired-end reads to ASVs, MaarjAM taxonomy, a phyloseq object, figures and statistics.

Each step runs in the background with a live log, and results are saved to the project folder, so you can close the app and resume where you left off.

## Pipeline

| # | Step | What it does |
|---|------|--------------|
| 1 | Setup & raw reads | Pick the folder of R1/R2 FASTQ files; samples are detected from filenames |
| 2 | Quality check | Quality profiles of the raw reads (nothing written to disk) |
| 3 | Remove ambiguous bases | `dada2::filterAndTrim(maxN = 0)` drops read pairs containing an N; shows raw vs. filtered read counts |
| 4 | Primer removal | `cutadapt`, with AMF primer presets (see below) or custom primers; primer counts before and after trimming (`Biostrings::vcountPattern`, all four orientations) |
| 5 | Filter and trim | `dada2::filterAndTrim` with truncation lengths chosen from quality profiles |
| 6 | Learn error rate & denoise | `dada2::learnErrors` + `dada2::dada` to infer ASVs |
| 7 | Merge pairs | Overlap merge, forward-only, or concatenate (for amplicons too long to overlap) |
| 8 | Remove chimeras | `dada2::removeBimeraDenovo` (consensus, pooled or per-sample) |
| 9 | Taxonomy | `blastn` against a MaarjAM database (defaults: ≥97% identity, ≥95% coverage, e-value < 1e-50) |
| 10 | Phyloseq object | Combines ASVs, taxonomy and sample metadata; optional rarefaction, relative abundance or CLR transform |
| 11 | Figures | Rarefaction curves; alpha diversity (Observed, Shannon, Simpson); NMDS and PCoA ordinations (Bray-Curtis, or Aitchison for CLR data) |
| 12 | PERMANOVA | `vegan::adonis2`, plus a `betadisper` dispersion test and Bonferroni-corrected pairwise PERMANOVA |
| 13 | Differential abundance | `ANCOMBC::ancombc2` at family, genus, virtual taxon or ASV level |

### Primer presets

Defined in [`inst/primer_presets.yaml`](inst/primer_presets.yaml), each with its literature citation:

- NS31 / AML2
- nu-SSU-0595 / nu-SSU-0948 (default for 2x250)
- nu-SSU-0450 / nu-SSU-0899 (default for 2x300)
- AMV4.5NF / AMDGR
- WANDA / AML2
- AML1 / AML2 (nested PCR / PacBio)

## Requirements

**R packages** (CRAN and Bioconductor), installed by `install_deps.R`:
shiny, bslib, bsicons, DT, plotly, ggplot2, processx, callr, fs, yaml, jsonlite, readr, vegan, dada2, Biostrings, ShortRead, phyloseq, ANCOMBC.

> `install_deps.R` pins CVXR to 1.0-15, because current ANCOMBC releases don't work with CVXR ≥ 1.9.

**Command-line tools:**

- [cutadapt](https://cutadapt.readthedocs.io/): used automatically if already installed. If not, the Primer removal step can install it for you into a private Python environment (`~/.arbusca/cutadapt-venv`); this only needs Python 3.
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (`blastn`, `makeblastdb`, `blastdbcmd`), which must be on your `PATH`: `conda install -c bioconda blast`

**Reference database:** a MaarjAM SSU FASTA file or a BLAST database built from it. The app does not download it for you. Point it at an existing BLAST database, or give it the FASTA and it will build one with `makeblastdb`.

## Installation

```bash
git clone https://github.com/beantkapoor786/Arbusca.git
cd Arbusca
Rscript install_deps.R
```

`install_deps.R` installs any missing R packages and reports whether `cutadapt` and BLAST+ are found.

## Running the app

From the repository folder, in R:

```r
shiny::runApp("Arbusca.R")
```

Or open `Arbusca.R` in RStudio and click **Run App**.

## Input data

- **Reads:** one folder containing paired-end FASTQ files (R1/R2). In Setup, you choose the delimiter and which filename parts make up the sample name. Forward and reverse files are paired by an `R1`/`R2` tag that must stand on its own in the filename (e.g. `Sample1_S1_L001_R1_001.fastq.gz`).
- **Sample metadata** (needed from step 10): a `.csv` file or a tab-delimited `.txt` file with one row per sample. One column must contain the same sample names used in Setup; you choose which column in the app. For example:

  ```csv
  sample_id,site,treatment
  1-827-1-1-1-Root-AMF,FieldA,Control
  10-405-1-5-1-Root-AMF,FieldA,Treated
  ```

## Output

Everything is written into subfolders of the project (reads) folder:

```
<project>/
├── 01_filtN/      reads with no ambiguous bases (N)
├── 02_trimmed/    primer-trimmed reads
├── 03_filtered/   filtered and trimmed reads
├── 04_denoise/    error models, denoised reads, sequence tables, read tracking
├── 05_taxonomy/   asv.fasta, raw and best-hit BLAST tables, taxonomy assignments
└── 06_output/     phyloseq.rds, phyloseq_transformed.rds, asv_table_transformed.csv
```

Figures, PERMANOVA and ANCOM-BC2 results are not written to the project folder; use the download buttons in each step to save them (PNG at 300 dpi, or CSV).

### Resuming and re-running

Selecting a project folder again marks each step whose output files already exist as done, so you can pick up where you left off. Re-running a step deletes the output of every step after it, because that output was built from the previous run.

App settings that persist across sessions, such as a custom `cutadapt` path, the BLAST database location and saved custom primers, are kept in `~/.arbusca/config.yaml`.

## Citation

If you use Arbusca, please cite its Zenodo record for the version you used.

Please also cite the tools Arbusca runs:

- **DADA2:** Callahan et al. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods* 13:581–583. [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- **cutadapt:** Martin (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal* 17:10–12. [doi:10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)
- **BLAST+:** Camacho et al. (2009) BLAST+: architecture and applications. *BMC Bioinformatics* 10:421. [doi:10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)
- **MaarjAM:** Öpik et al. (2010) The online database MaarjAM reveals global and ecosystemic distribution patterns in arbuscular mycorrhizal fungi (Glomeromycota). *New Phytologist* 188:223–241. [doi:10.1111/j.1469-8137.2010.03334.x](https://doi.org/10.1111/j.1469-8137.2010.03334.x)
- **phyloseq:** McMurdie & Holmes (2013) phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE* 8:e61217. [doi:10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)
- **vegan** (Steps 11–12): cite the version you used; run `citation("vegan")` in R.
- **ANCOM-BC2** (Step 13): Lin & Peddada (2024) Multigroup analysis of compositions of microbiomes with covariate adjustment and repeated measures. *Nature Methods* 21:83–91. [doi:10.1038/s41592-023-02092-7](https://doi.org/10.1038/s41592-023-02092-7)

The primer preset you choose also shows its original citation in the app.

## License

Arbusca is released under the [MIT License](LICENSE).
