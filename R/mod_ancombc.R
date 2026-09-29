# Step 13 -- ANCOM-BC2 differential abundance (ANCOMBC::ancombc2) on Step
# 10's raw phyloseq object: ANCOM-BC2 does its own sample/taxon bias
# correction, so it is never fed the transformed object.
#
# Same design as Step 12 (PERMANOVA): build_ancombc_script() renders the R
# script from the current inputs, the preview shows it, and the callr job
# evaluates exactly that text. Analysis-only; nothing depends on it.

ANCOMBC_TAX_LEVELS <- c("Family" = "family", "Genus" = "genus", "Virtual taxon (VT)" = "vt", "ASV" = "asv")
ANCOMBC_SEED <- 123

# Result tables ancombc2 can return, in display order: element name in the
# ancombc2 output, card title, one-line explanation, and whether the table is
# "wide" (lfc_/q_/diff_... columns per covariate or comparison) or "long"
# (one W/p/q/diff_abn per taxon).
ANCOMBC_TABLES <- list(
  primary  = list(el = "res",       title = "Primary analysis",          wide = TRUE,
                  desc = "Log fold changes (LFC) per covariate in the fixed-effects formula. For a categorical covariate, each column compares one level against the reference level."),
  global   = list(el = "res_global", title = "Global test",              wide = FALSE,
                  desc = "Taxa whose abundance differs across any of the groups."),
  pairwise = list(el = "res_pair",   title = "Pairwise directional test", wide = TRUE,
                  desc = "Every pair of groups, with mdFDR control. A comparison named like 'GroupB' is B vs. the reference level; 'GroupC_GroupB' is C vs. B. Positive LFC = higher in the first-named group."),
  dunnett  = list(el = "res_dunn",   title = "Dunnett's type test",       wide = TRUE,
                  desc = "Each group compared against the reference (control) level, with mdFDR control."),
  trend    = list(el = "res_trend",  title = "Trend test",                wide = FALSE,
                  desc = "Taxa showing a monotone increasing or decreasing pattern across the groups, in the level order you set.")
)

# Monotone contrasts over the K-1 non-reference coefficients (ANCOMBC
# vignette pattern): increasing = 0 <= b2 <= b3 <= ...; decreasing = its
# negation. Returned as the literal matrix(c(...)) text for the script.
ancombc_trend_contrast_text <- function(n_levels) {
  k <- n_levels - 1
  inc <- diag(k)
  inc[cbind(seq_len(k)[-1], seq_len(k - 1))] <- -1
  as_text <- function(m) sprintf("matrix(%s, nrow = %d, byrow = TRUE)", deparse1(c(t(m))), k)
  list(increasing = as_text(inc), decreasing = as_text(-inc), node = k)
}

# p: named list of validated settings (see ancombc_params in the server).
build_ancombc_script <- function(ps_path, p) {
  q <- function(x) deparse1(x)
  vec <- function(x) if (length(x) == 1) q(x) else sprintf("c(%s)", paste(vapply(x, q, ""), collapse = ", "))
  has_group <- !is.null(p$group)

  lines <- c(
    "library(phyloseq)",
    "library(ANCOMBC)",
    sprintf("set.seed(%d)", ANCOMBC_SEED),
    "",
    sprintf("ps <- readRDS(%s)", q(ps_path)),
    "",
    "# Drop samples with a missing value in any model variable",
    'meta <- as(sample_data(ps), "data.frame")',
    sprintf("ps <- prune_samples(complete.cases(meta[, %s, drop = FALSE]), ps)", vec(p$vars))
  )

  if (has_group) {
    lines <- c(lines, "",
      sprintf("# Level order of %s -- the first level is the reference / control", p$group),
      sprintf("sample_data(ps)[[%s]] <- factor(sample_data(ps)[[%s]], levels = %s)", q(p$group), q(p$group), vec(p$levels))
    )
  }

  if (!identical(p$tax_level, "asv")) {
    lines <- c(lines, "",
      sprintf("# Agglomerate to %s level; ASVs without a %s assignment are dropped", p$tax_level, p$tax_level),
      sprintf("ps <- tax_glom(ps, %s, NArm = TRUE)", q(p$tax_level)),
      sprintf("taxa_names(ps) <- make.unique(as.character(tax_table(ps)[, %s]))", q(p$tax_level))
    )
  }

  if (isTRUE(p$trend)) {
    tc <- ancombc_trend_contrast_text(length(p$levels))
    lines <- c(lines, "",
      "# Trend test contrasts (monotone increasing / decreasing across the level order)",
      "trend_contrast <- list(",
      sprintf("  increasing = %s,", tc$increasing),
      sprintf("  decreasing = %s", tc$decreasing),
      ")"
    )
  }

  args <- c(
    "data = ps",
    "tax_level = NULL",
    sprintf("fix_formula = %s", q(p$fix)),
    sprintf("rand_formula = %s", if (is.null(p$rand)) "NULL" else q(p$rand)),
    sprintf("p_adj_method = %s", q(p$p_adj)),
    sprintf("pseudo = %s", q(p$pseudo)),
    sprintf("pseudo_sens = %s", q(p$pseudo_sens)),
    sprintf("prv_cut = %s", q(p$prv_cut)),
    sprintf("lib_cut = %s", q(p$lib_cut)),
    sprintf("s0_perc = %s", q(p$s0_perc)),
    sprintf("group = %s", if (has_group) q(p$group) else "NULL"),
    sprintf("struc_zero = %s", q(p$struc_zero)),
    sprintf("neg_lb = %s", q(p$neg_lb)),
    sprintf("alpha = %s", q(p$alpha)),
    sprintf("n_cl = %s", q(p$n_cl)),
    "verbose = TRUE",
    sprintf("global = %s", q(p$global)),
    sprintf("pairwise = %s", q(p$pairwise)),
    sprintf("dunnet = %s", q(p$dunnet)),
    sprintf("trend = %s", q(p$trend)),
    sprintf("iter_control = list(tol = %s, max_iter = %s, verbose = TRUE)", q(p$iter_tol), q(p$iter_max)),
    sprintf("em_control = list(tol = %s, max_iter = %s)", q(p$em_tol), q(p$em_max)),
    sprintf("mdfdr_control = list(fwer_ctrl_method = %s, B = %s)", q(p$fwer), q(p$mdfdr_B))
  )
  if (isTRUE(p$trend)) {
    args <- c(args, sprintf("trend_control = list(contrast = trend_contrast, node = list(%d, %d), solver = %s, B = %s)",
                            tc$node, tc$node, q(p$solver), q(p$trend_B)))
  }

  lines <- c(lines, "",
    "# ANCOM-BC2",
    "out <- ancombc2(",
    paste0("  ", args, c(rep(",", length(args) - 1), "")),
    ")",
    "print(head(out$res))"
  )

  paste(lines, collapse = "\n")
}

# A diff_*/diff_abn flag counts as significant only if the taxon also passed
# the pseudocount sensitivity analysis (passed_ss); a missing passed_ss column
# (pseudo_sens = FALSE) counts as passed. NA flags count as not significant.
ancombc_is_sig <- function(diff, passed) (diff %in% TRUE) & (is.null(passed) | (passed %in% TRUE))

# Wide table -> one row per significant (taxon, covariate/comparison).
ancombc_significant_long <- function(df) {
  terms <- sub("^lfc_", "", grep("^lfc_", names(df), value = TRUE))
  terms <- setdiff(terms, "(Intercept)")
  rows <- lapply(terms, function(s) {
    sig <- ancombc_is_sig(df[[paste0("diff_", s)]], df[[paste0("passed_ss_", s)]])
    if (!any(sig)) return(NULL)
    lfc <- df[[paste0("lfc_", s)]][sig]
    data.frame(
      Taxon = df$taxon[sig], Comparison = s, LFC = lfc,
      SE = df[[paste0("se_", s)]][sig], q = df[[paste0("q_", s)]][sig],
      Direction = ifelse(lfc > 0, "Increased", "Decreased"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame(Taxon = character(0), Comparison = character(0), LFC = numeric(0),
                      SE = numeric(0), q = numeric(0), Direction = character(0)))
  }
  out[order(out$Comparison, out$q), , drop = FALSE]
}

# Readable name and direction wording for one primary-analysis term.
# levels_map: categorical fixed-effect variable -> its levels, reference
# first (so "StateWashington" reads "Washington vs Texas (State)");
# numeric_vars: the numeric ones. Anything else (e.g. interactions) keeps
# its raw name.
ancombc_term_label <- function(term, levels_map, numeric_vars) {
  if (term %in% numeric_vars) {
    return(list(label = term, up = sprintf("increase with %s", term), down = sprintf("decrease with %s", term)))
  }
  for (v in names(levels_map)) {
    lv <- levels_map[[v]]
    hit <- lv[-1][paste0(v, lv[-1]) == term]
    if (length(hit) == 1) {
      return(list(label = sprintf("%s vs %s (%s)", hit, lv[1], v),
                  up = sprintf("higher in %s", hit), down = sprintf("lower in %s", hit)))
    }
    # Pairwise-test name "<v><B>_<v><A>": B vs A.
    for (a in lv) for (b in lv) {
      if (identical(paste0(v, b, "_", v, a), term)) {
        return(list(label = sprintf("%s vs %s (%s)", b, a, v),
                    up = sprintf("higher in %s", b), down = sprintf("lower in %s", b)))
      }
    }
  }
  list(label = term, up = "positive LFC", down = "negative LFC")
}

# One alert per non-intercept term, like the PERMANOVA summary: how many
# taxa differ and which way, listing the strongest (largest |LFC|) first.
ancombc_primary_summary <- function(df, levels_map, numeric_vars, max_listed = 10) {
  terms <- setdiff(sub("^lfc_", "", grep("^lfc_", names(df), value = TRUE)), "(Intercept)")
  sig <- ancombc_significant_long(df)
  lapply(terms, function(term) {
    lab <- ancombc_term_label(term, levels_map, numeric_vars)
    s <- sig[sig$Comparison == term, , drop = FALSE]
    if (nrow(s) == 0) {
      return(shiny::div(class = "alert alert-secondary small mb-2",
        sprintf("%s: none of the %d taxa tested differ significantly.", lab$label, nrow(df))))
    }
    s <- s[order(-abs(s$LFC)), , drop = FALSE]
    shown <- utils::head(s, max_listed)
    shiny::div(class = "alert alert-success small mb-2",
      sprintf("%s: %d of %d taxa differ significantly (%d %s, %d %s):", lab$label, nrow(s), nrow(df),
              sum(s$LFC > 0), lab$up, sum(s$LFC < 0), lab$down),
      shiny::tags$ul(class = "mb-0", lapply(seq_len(nrow(shown)), function(i) shiny::tags$li(sprintf(
        "%s: %s (LFC = %.2f, q = %s)", shown$Taxon[i], if (shown$LFC[i] > 0) lab$up else lab$down, shown$LFC[i], signif(shown$q[i], 3))))),
      if (nrow(s) > max_listed) shiny::div(sprintf("...and %d more, see the table.", nrow(s) - max_listed))
    )
  })
}

# Pairwise test: often dozens of comparisons, so list only those with
# significant taxa (up to max_listed taxa each) and count the rest.
ancombc_pairwise_summary <- function(df, levels_map, numeric_vars, max_listed = 5) {
  terms <- setdiff(sub("^lfc_", "", grep("^lfc_", names(df), value = TRUE)), "(Intercept)")
  sig <- ancombc_significant_long(df)
  hit_terms <- terms[terms %in% sig$Comparison]
  if (length(hit_terms) == 0) {
    return(shiny::div(class = "alert alert-secondary small mb-2",
      sprintf("None of the %d pairwise comparisons show a significantly different taxon (%d taxa tested).", length(terms), nrow(df))))
  }
  items <- lapply(hit_terms, function(term) {
    lab <- ancombc_term_label(term, levels_map, numeric_vars)
    s <- sig[sig$Comparison == term, , drop = FALSE]
    s <- s[order(-abs(s$LFC)), , drop = FALSE]
    shown <- utils::head(s, max_listed)
    shiny::tags$li(sprintf("%s: %d taxa (%s%s)", lab$label, nrow(s),
      paste(sprintf("%s %s", shown$Taxon, ifelse(shown$LFC > 0, lab$up, lab$down)), collapse = ", "),
      if (nrow(s) > max_listed) sprintf(", and %d more", nrow(s) - max_listed) else ""))
  })
  shiny::div(class = "alert alert-success small mb-2",
    sprintf("%d of %d pairwise comparisons have significantly different taxa:", length(hit_terms), length(terms)),
    shiny::tags$ul(class = "mb-0", items),
    if (length(hit_terms) < length(terms)) shiny::div(sprintf("The other %d comparisons show none.", length(terms) - length(hit_terms)))
  )
}

# Global test: taxa whose abundance differs across the groups at all,
# smallest q first. (Sorted by q, not W: ancombc2's global W doesn't rank
# significance -- significant taxa can have W near 0.)
ancombc_global_summary <- function(df, group, max_listed = 10) {
  s <- df[ancombc_is_sig(df$diff_abn, df$passed_ss), , drop = FALSE]
  across <- if (is.null(group)) "the groups" else sprintf("the %s groups", group)
  if (nrow(s) == 0) {
    return(shiny::div(class = "alert alert-secondary small mb-2",
      sprintf("None of the %d taxa tested differ significantly across %s.", nrow(df), across)))
  }
  s <- s[order(s$q_val), , drop = FALSE]
  shown <- utils::head(s, max_listed)
  shiny::div(class = "alert alert-success small mb-2",
    sprintf("%d of %d taxa differ significantly across %s:", nrow(s), nrow(df), across),
    shiny::tags$ul(class = "mb-0", lapply(seq_len(nrow(shown)), function(i) shiny::tags$li(sprintf(
      "%s (q = %s)", shown$taxon[i], signif(shown$q_val[i], 3))))),
    if (nrow(s) > max_listed) shiny::div(sprintf("...and %d more, see the table.", nrow(s) - max_listed)),
    shiny::div(class = "mt-1", "The global test says only that a taxon differs somewhere among the groups; the pairwise test shows which groups.")
  )
}

ancombc_round <- function(df) {
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(x) signif(x, 4))
  df
}

mod_ancombc_ui <- function(id) {
  ns <- shiny::NS(id)
  num <- function(inputId, label, value, step = NULL, min = NULL) shiny::numericInput(ns(inputId), label, value = value, step = step, min = min, width = "100%")
  chk <- function(inputId, label, value) shiny::checkboxInput(ns(inputId), label, value = value)

  step_card(
    stacked = TRUE,
    title = "13. Differential Abundance (ANCOM-BC2)",
    description = "Identifies taxa whose absolute abundance differs between groups or along covariates, correcting for sample- and taxon-specific biases.",
    params = shiny::tagList(
      shiny::div(class = "alert alert-info small mb-3",
        shiny::strong("Runs on raw counts, not transformed data. "),
        "ANCOM-BC2 log-transforms the counts and corrects for sample-specific sequencing bias itself, so any transformation applied in the Phyloseq Object step is ignored here."),
      shiny::uiOutput(ns("gate_alert")),
      shiny::uiOutput(ns("input_status")),
      shiny::fluidRow(
        shiny::column(4,
          param_group("Model",
            shiny::textInput(ns("fix_formula"),
              help_label("Fixed effects", "Right-hand side of the fixed-effects model, e.g. Treatment + Age. Differential abundance is reported for each term."),
              placeholder = "Treatment + Age", width = "100%"),
            shiny::textInput(ns("rand_formula"),
              help_label("Random effects (optional)", "lme4-style random effects, e.g. (1 | Site) for repeated sampling within sites. Leave empty for none."),
              placeholder = "(1 | Site)", width = "100%"),
            shiny::uiOutput(ns("formula_feedback")),
            shiny::uiOutput(ns("variable_list"))
          )
        ),
        shiny::column(4,
          param_group("Groups and taxa",
            shiny::uiOutput(ns("group_picker")),
            shiny::uiOutput(ns("level_picker")),
            shiny::selectInput(ns("tax_level"), "Taxonomic level", choices = ANCOMBC_TAX_LEVELS, selected = "genus", width = "100%"),
            shiny::uiOutput(ns("group_tests"))
          )
        ),
        shiny::column(4,
          param_group("Filtering and correction",
            shiny::selectInput(ns("p_adj"), "P-value adjustment", choices = stats::p.adjust.methods, selected = "holm", width = "100%"),
            num("prv_cut", help_label("Prevalence filter", "Taxa present in fewer than this fraction of samples are excluded (0.10 = 10%)."), 0.10, step = 0.05, min = 0),
            num("lib_cut", help_label("Library size cutoff", "Samples with fewer total reads than this are excluded. 0 keeps all samples."), 0, step = 500, min = 0)
          )
        )
      ),
      bslib::accordion(
        open = FALSE,
        bslib::accordion_panel(
          "Advanced settings",
          shiny::fluidRow(
            shiny::column(4,
              param_group("Pseudocount",
                num("pseudo", help_label("Pseudocount", "Added to counts before log-transforming. 0 (default) lets ANCOM-BC2 handle zeros itself."), 0, step = 0.5, min = 0),
                chk("pseudo_sens", "Pseudocount sensitivity analysis", TRUE)
              ),
              param_group("Zeros",
                chk("struc_zero", "Detect structural zeros (needs a group variable)", FALSE),
                chk("neg_lb", "Use negative lower bound for structural zeros", FALSE)
              )
            ),
            shiny::column(4,
              param_group("Regularization and significance",
                num("s0_perc", help_label("s0 percentile", "Percentile of standard errors added as a small constant to stabilize the test statistic for taxa with tiny variance."), 0.05, step = 0.01, min = 0),
                num("alpha", "Significance level (alpha)", 0.05, step = 0.01, min = 0)
              ),
              param_group("Iteration (bias estimation)",
                num("iter_tol", "Tolerance", 0.01, step = 0.005, min = 0),
                num("iter_max", "Max iterations", 20, step = 5, min = 1)
              ),
              param_group("EM algorithm",
                num("em_tol", "Tolerance", 1e-5, step = 1e-5, min = 0),
                num("em_max", "Max iterations", 100, step = 10, min = 1)
              )
            ),
            shiny::column(4,
              param_group("mdFDR control (pairwise / Dunnett)",
                shiny::selectInput(ns("fwer"), "FWER control method", choices = stats::p.adjust.methods, selected = "holm", width = "100%"),
                num("mdfdr_B", "Bootstrap samples (B)", 100, step = 50, min = 1)
              ),
              param_group("Additional group tests",
                chk("dunnet", "Dunnett's type test (each group vs. the reference level)", FALSE),
                chk("trend", "Trend test (monotone across the level order)", FALSE),
                shiny::selectInput(ns("solver"), "Trend solver", choices = c("ECOS", "SCS"), width = "100%"),
                num("trend_B", "Trend bootstrap samples (B)", 100, step = 50, min = 1)
              ),
              param_group("Performance",
                num("n_cl", "CPU cores", 1, step = 1, min = 1)
              )
            )
          )
        )
      ),
      param_group("Command preview",
        shiny::div(class = "amf-param-hint mb-2", "These are the exact R commands that run when you click Run ANCOM-BC2."),
        shiny::verbatimTextOutput(ns("script_preview"))
      )
    ),
    results = shiny::tagList(
      result_card("Run ANCOM-BC2",
        run_row(ns, "Run ANCOM-BC2"),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("log"))
      ),
      shiny::uiOutput(ns("error_alert")),
      shiny::uiOutput(ns("results_section")),
      results_placeholder("Run ANCOM-BC2 to see differentially abundant taxa.")
    )
  )
}

# sample_table unused directly (kept for signature consistency with other
# step modules) -- this step reads the phyloseq object straight off disk.
mod_ancombc_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "ancombc"

    ps_path <- shiny::reactive(file.path(rv$project_dir, "06_output", "phyloseq.rds"))
    inputs_ready <- shiny::reactive({
      shiny::req(rv$project_dir)
      rv$status$output
      fs::file_exists(ps_path())
    })

    pkgs_ready <- shiny::reactive(requireNamespace("phyloseq", quietly = TRUE) && requireNamespace("ANCOMBC", quietly = TRUE))
    output$gate_alert <- shiny::renderUI({
      if (isTRUE(pkgs_ready())) return(NULL)
      shiny::div(class = "alert alert-warning small mb-3",
                 "phyloseq and ANCOMBC are both required for this step -- install any missing via BiocManager::install().")
    })

    output$input_status <- shiny::renderUI({
      if (!isTRUE(inputs_ready())) {
        return(shiny::div(class = "alert alert-warning small", "No phyloseq object found yet -- finish the Phyloseq Object step first."))
      }
      NULL
    })

    meta <- shiny::reactive({
      shiny::req(isTRUE(inputs_ready()), requireNamespace("phyloseq", quietly = TRUE))
      as(phyloseq::sample_data(readRDS(ps_path())), "data.frame")
    })

    output$variable_list <- shiny::renderUI(metadata_variable_list(meta()))

    # Pre-fill the fixed formula with the first categorical variable.
    shiny::observeEvent(meta(), {
      if (nzchar(input$fix_formula %||% "")) return()
      m <- meta()
      cat_vars <- names(m)[!vapply(m, is.numeric, logical(1))]
      if (length(cat_vars) > 0) shiny::updateTextInput(session, "fix_formula", value = cat_vars[[1]])
    })

    # Parse one formula box. Returns list(ok, msg) or list(ok, rhs, vars);
    # rhs is rebuilt from the parsed formula so only a single well-formed
    # expression ever reaches the script.
    parse_rhs <- function(txt, m, label) {
      txt <- trimws(sub("^\\s*~", "", txt %||% ""))
      f <- tryCatch(stats::as.formula(paste("~", txt)), error = function(e) NULL)
      if (is.null(f) || length(f) != 2) return(list(ok = FALSE, msg = sprintf("Could not parse the %s formula.", label)))
      vars <- all.vars(f)
      missing <- setdiff(vars, names(m))
      if (length(missing) > 0) {
        return(list(ok = FALSE, msg = sprintf("%s formula -- not in the sample metadata: %s.", label, paste(missing, collapse = ", "))))
      }
      list(ok = TRUE, rhs = deparse1(f[[2]]), vars = vars)
    }

    formulas <- shiny::reactive({
      m <- meta()
      if (!nzchar(trimws(input$fix_formula %||% ""))) return(list(ok = FALSE, msg = "Enter at least one fixed-effect variable."))
      fix <- parse_rhs(input$fix_formula, m, "Fixed-effects")
      if (!fix$ok) return(fix)
      rand <- NULL
      if (nzchar(trimws(input$rand_formula %||% ""))) {
        rand <- parse_rhs(input$rand_formula, m, "Random-effects")
        if (!rand$ok) return(rand)
      }
      # ancombc2 can't estimate a term the others fully determine (e.g. State
      # when every Location lies in one State) and fails only after its
      # checks, with a cryptic "Estimation failed" -- catch it here instead.
      mm_data <- m[stats::complete.cases(m[fix$vars]), fix$vars, drop = FALSE]
      mm <- tryCatch(stats::model.matrix(stats::as.formula(paste("~", fix$rhs)), mm_data), error = function(e) e)
      if (inherits(mm, "error")) {
        # e.g. a variable with a single level: nothing to compare.
        return(list(ok = FALSE, msg = sprintf("Can't build the fixed-effects model: %s. Check that each categorical variable has at least two levels.", conditionMessage(mm))))
      }
      qr_mm <- qr(mm)
      if (qr_mm$rank < ncol(mm)) {
        term_labels <- attr(stats::terms(stats::as.formula(paste("~", fix$rhs))), "term.labels")
        aliased <- unique(term_labels[attr(mm, "assign")[qr_mm$pivot[(qr_mm$rank + 1):ncol(mm)]]])
        return(list(ok = FALSE, msg = sprintf(
          "%s is confounded with %s: one is completely determined by the other (e.g. every site lies in a single region), so ANCOM-BC2 can't separate their effects. Remove one of them from the fixed-effects formula.",
          paste(aliased, collapse = ", "), paste(setdiff(term_labels, aliased), collapse = ", "))))
      }
      cat_fix <- fix$vars[!vapply(m[fix$vars], is.numeric, logical(1))]
      list(ok = TRUE, fix = fix$rhs, rand = rand$rhs, vars = unique(c(fix$vars, rand$vars)), cat_fix = cat_fix)
    })

    output$formula_feedback <- shiny::renderUI({
      fm <- formulas()
      if (!isTRUE(fm$ok)) return(shiny::div(class = "text-danger small mb-2", fm$msg))
      NULL
    })

    # Group variable: a categorical term of the fixed formula (ancombc2
    # requires the group variable to be in fix_formula).
    output$group_picker <- shiny::renderUI({
      fm <- formulas()
      choices <- if (isTRUE(fm$ok)) fm$cat_fix else character(0)
      current <- shiny::isolate(input$group)
      selected <- if (!is.null(current) && current %in% choices) current else if (length(choices) > 0) choices[[1]] else ""
      shiny::selectInput(ns("group"),
        help_label("Group variable", "Categorical fixed-effect variable used for structural zeros and the global, pairwise, Dunnett's and trend tests. Only categorical terms of the fixed formula are offered."),
        choices = c("None" = "", choices), selected = selected, width = "100%")
    })

    group_levels <- shiny::reactive({
      g <- input$group %||% ""
      shiny::req(nzchar(g), g %in% names(meta()))
      sort(unique(as.character(stats::na.omit(meta()[[g]]))))
    })

    output$level_picker <- shiny::renderUI({
      lv <- group_levels()
      shiny::selectizeInput(ns("levels"),
        help_label("Level order (first = reference / control)", "The first level is the reference for the primary analysis and the control for Dunnett's test. The order also defines the trend test's direction. To reorder, remove levels and re-add them in the order you want."),
        choices = lv, selected = lv, multiple = TRUE, width = "100%")
    })

    n_group_levels <- shiny::reactive(if (nzchar(input$group %||% "")) length(group_levels()) else 0L)

    output$group_tests <- shiny::renderUI({
      if (n_group_levels() < 3) {
        return(shiny::div(class = "amf-param-hint", "Global, pairwise, Dunnett's and trend tests need a group variable with 3 or more levels."))
      }
      shiny::tagList(
        shiny::checkboxInput(ns("global"), "Global test", value = shiny::isolate(input$global) %||% TRUE),
        shiny::checkboxInput(ns("pairwise"), "Pairwise directional test", value = shiny::isolate(input$pairwise) %||% TRUE)
      )
    })

    num_in <- function(x, default) {
      v <- suppressWarnings(as.numeric(x))
      if (length(v) != 1 || is.na(v)) default else v
    }

    params <- shiny::reactive({
      fm <- formulas()
      if (!isTRUE(fm$ok)) return(fm)
      group <- input$group %||% ""
      has_group <- nzchar(group)
      levels <- NULL
      if (has_group) {
        levels <- input$levels
        if (!setequal(levels, group_levels()) || length(levels) != length(group_levels())) {
          return(list(ok = FALSE, msg = "Level order must include every level of the group variable exactly once."))
        }
      }
      multi <- has_group && length(levels) >= 3
      list(
        ok = TRUE,
        fix = fm$fix, rand = fm$rand, vars = fm$vars,
        group = if (has_group) group else NULL, levels = levels,
        tax_level = input$tax_level %||% "genus",
        p_adj = input$p_adj %||% "holm",
        prv_cut = num_in(input$prv_cut, 0.10), lib_cut = num_in(input$lib_cut, 0),
        pseudo = num_in(input$pseudo, 0), pseudo_sens = isTRUE(input$pseudo_sens %||% TRUE),
        s0_perc = num_in(input$s0_perc, 0.05),
        struc_zero = has_group && isTRUE(input$struc_zero), neg_lb = has_group && isTRUE(input$neg_lb),
        alpha = num_in(input$alpha, 0.05), n_cl = max(1, round(num_in(input$n_cl, 1))),
        global = multi && isTRUE(input$global %||% TRUE),
        pairwise = multi && isTRUE(input$pairwise %||% TRUE),
        dunnet = multi && isTRUE(input$dunnet),
        trend = multi && isTRUE(input$trend),
        iter_tol = num_in(input$iter_tol, 0.01), iter_max = max(1, round(num_in(input$iter_max, 20))),
        em_tol = num_in(input$em_tol, 1e-5), em_max = max(1, round(num_in(input$em_max, 100))),
        fwer = input$fwer %||% "holm", mdfdr_B = max(1, round(num_in(input$mdfdr_B, 100))),
        solver = input$solver %||% "ECOS", trend_B = max(1, round(num_in(input$trend_B, 100)))
      )
    })

    script <- shiny::reactive({
      p <- params()
      shiny::validate(shiny::need(isTRUE(p$ok), p$msg))
      build_ancombc_script(ps_path(), p)
    })

    output$script_preview <- shiny::renderText(script())

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv), isTRUE(pkgs_ready()), isTRUE(inputs_ready()))
      p <- params()
      if (!isTRUE(p$ok)) {
        shiny::showNotification(p$msg, type = "error")
        return()
      }

      # Warnings (e.g. too few taxa to estimate the sample biases) matter for
      # interpretation, so they're collected and shown above the results
      # rather than left in the log.
      ancombc_job <- function(code, table_els) {
        env <- new.env()
        warns <- character(0)
        withCallingHandlers(
          eval(parse(text = code), envir = env),
          warning = function(w) warns <<- c(warns, conditionMessage(w))
        )
        out <- env$out
        list(tables = Filter(Negate(is.null), out[table_els]), warnings = unique(warns))
      }

      launch_step_callr(rv, step_id, ancombc_job,
                        list(code = script(), table_els = vapply(ANCOMBC_TABLES, `[[`, "", "el")))

      # Levels as the model saw them (group: the chosen order; others: R's
      # default sorted factor levels), for naming comparisons in the summary.
      m <- meta()
      fix_vars <- all.vars(stats::as.formula(paste("~", p$fix)))
      is_num <- vapply(m[fix_vars], is.numeric, logical(1))
      last_run(list(
        group = p$group,
        numeric_vars = fix_vars[is_num],
        levels_map = sapply(fix_vars[!is_num], function(v) {
          if (identical(v, p$group)) p$levels else sort(unique(as.character(stats::na.omit(m[[v]]))))
        }, simplify = FALSE)
      ))
    })
    last_run <- shiny::reactiveVal(list(numeric_vars = character(0), levels_map = list()))

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    # A failed run leaves no results, so surface the error itself (the last
    # "ERROR:" line proc_runner appends to the log) where results would be.
    output$error_alert <- shiny::renderUI({
      shiny::req(identical(rv$status[[step_id]], "ERROR"))
      # Keep only the underlying cause, not callr's "! in callr subprocess." wrapper.
      err <- utils::tail(grep("^ERROR:", rv$log[[step_id]], value = TRUE), 1)
      msg <- if (length(err) > 0) trimws(gsub("(^|\n)!\\s*", "\\1", sub("(?s)^.*Caused by error:\\s*", "", sub("^ERROR:\\s*", "", err), perl = TRUE)))
      shiny::div(class = "alert alert-danger small", style = "white-space: pre-line;",
        shiny::strong("ANCOM-BC2 failed. "),
        if (length(msg) > 0 && nzchar(msg)) msg else "See the log above for details.")
    })

    result_table <- function(key) {
      res <- rv$artifacts[[step_id]]
      shiny::req(res)
      res$tables[[ANCOMBC_TABLES[[key]]$el]]
    }

    output$results_section <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res))
      primary <- res$tables$res
      n_sig <- nrow(ancombc_significant_long(primary))

      cards <- lapply(names(ANCOMBC_TABLES), function(key) {
        spec <- ANCOMBC_TABLES[[key]]
        if (is.null(res$tables[[spec$el]])) return(NULL)
        result_card(spec$title,
          shiny::div(class = "text-muted small mb-2", spec$desc),
          shiny::checkboxInput(ns(paste0("sig_", key)), "Significant taxa only (q < alpha and passed sensitivity analysis)", value = TRUE, width = "100%"),
          dt_output(ns(paste0("table_", key))),
          shiny::div(class = "mt-2", switch(key,
            primary = ancombc_primary_summary(primary, last_run()$levels_map, last_run()$numeric_vars),
            global = ancombc_global_summary(res$tables[[spec$el]], last_run()$group),
            pairwise = ancombc_pairwise_summary(res$tables[[spec$el]], last_run()$levels_map, last_run()$numeric_vars)
          )),
          actions = shiny::downloadButton(ns(paste0("download_", key)), "Download CSV", class = "btn-outline-secondary btn-sm")
        )
      })

      shiny::tagList(
        if (length(res$warnings) > 0) shiny::div(class = "alert alert-warning small",
          shiny::strong("ANCOM-BC2 reported: "),
          shiny::tags$ul(class = "mb-0", lapply(res$warnings, shiny::tags$li))),
        stat_boxes(list(
          "Taxa tested" = format(nrow(primary), big.mark = ","),
          "Significant findings (primary)" = format(n_sig, big.mark = ",")
        )),
        cards
      )
    })

    lapply(names(ANCOMBC_TABLES), function(key) {
      spec <- ANCOMBC_TABLES[[key]]

      output[[paste0("table_", key)]] <- DT::renderDataTable({
        df <- result_table(key)
        sig_only <- isTRUE(input[[paste0("sig_", key)]])
        shown <- if (!sig_only) {
          df
        } else if (spec$wide) {
          ancombc_significant_long(df)
        } else {
          df[ancombc_is_sig(df$diff_abn, df$passed_ss), , drop = FALSE]
        }
        DT::datatable(ancombc_round(shown), rownames = FALSE,
                      options = list(scrollX = TRUE, pageLength = 15, dom = "ftip"))
      })

      output[[paste0("download_", key)]] <- shiny::downloadHandler(
        filename = function() sprintf("ancombc2_%s_%s.csv", key, format(Sys.time(), "%Y%m%d_%H%M%S")),
        content = function(file) utils::write.csv(result_table(key), file, row.names = FALSE)
      )
    })

    mod_logpanel_server("log", shiny::reactive(rv$log[[step_id]] %||% character(0)))
  })
}
