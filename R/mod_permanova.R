# Step 12 -- PERMANOVA (vegan::adonis2) on Step 10's phyloseq object, plus a
# betadisper homogeneity-of-dispersion test for every categorical formula
# term and Bonferroni-corrected pairwise PERMANOVA for categorical terms with
# more than two groups.
#
# The command preview is not a description of what runs -- it IS what runs:
# build_permanova_script() renders the R script from the current inputs, and
# the callr background job evaluates exactly that text. Like Figures, this is
# analysis-only: nothing in STEP_ORDER depends on it.

PERMANOVA_DISTANCES <- c("Bray-Curtis" = "bray", "Jaccard (presence/absence)" = "jaccard", "Euclidean" = "euclidean")
PERMANOVA_BY <- c("Sequential (terms)" = "terms", "Marginal (margin)" = "margin", "Whole model (NULL)" = "none")
PERMANOVA_SEED <- 123

# rhs: formula right-hand side, already normalized by deparse (see
# formula_check in the server). disp_vars: categorical terms that get
# betadisper; pair_vars: the subset with > 2 groups that get pairwise tests.
build_permanova_script <- function(ps_path, rhs, vars, method, permutations, by, disp_vars, pair_vars) {
  q <- function(x) deparse(x)
  vec <- function(x) if (length(x) == 1) q(x) else sprintf("c(%s)", paste(vapply(x, q, ""), collapse = ", "))
  dist_call <- switch(method,
    bray      = 'vegdist(otu, method = "bray")',
    jaccard   = 'vegdist(otu, method = "jaccard", binary = TRUE)',
    euclidean = 'vegdist(otu, method = "euclidean")'
  )
  by_arg <- if (identical(by, "none")) "NULL" else q(by)

  lines <- c(
    "library(phyloseq)",
    "library(vegan)",
    sprintf("set.seed(%d)", PERMANOVA_SEED),
    "",
    sprintf("ps <- readRDS(%s)", q(ps_path)),
    'otu <- as(otu_table(ps), "matrix")',
    "if (taxa_are_rows(ps)) otu <- t(otu)",
    'meta <- as(sample_data(ps), "data.frame")',
    "",
    "# Drop samples with a missing value in any formula variable",
    sprintf("meta <- meta[complete.cases(meta[, %s, drop = FALSE]), , drop = FALSE]", vec(vars)),
    "otu <- otu[rownames(meta), , drop = FALSE]",
    "",
    sprintf("dist_mat <- %s", dist_call),
    "",
    "# PERMANOVA",
    sprintf("permanova <- adonis2(dist_mat ~ %s, data = meta, permutations = %d, by = %s)", rhs, permutations, by_arg),
    "print(permanova)"
  )

  if (length(disp_vars) > 0) {
    lines <- c(lines, "", "# Homogeneity of multivariate dispersion", "betadisper_fit <- list()", "dispersion_test <- list()")
    for (v in disp_vars) {
      lines <- c(lines,
        sprintf("betadisper_fit[[%s]] <- betadisper(dist_mat, meta[[%s]])", q(v), q(v)),
        sprintf("dispersion_test[[%s]] <- permutest(betadisper_fit[[%s]], permutations = %d)", q(v), q(v), permutations),
        sprintf("print(dispersion_test[[%s]])", q(v))
      )
    }
  }

  if (length(pair_vars) > 0) {
    lines <- c(lines, "",
      "# Pairwise PERMANOVA, Bonferroni-corrected",
      "pairwise_permanova <- function(dist_mat, groups, permutations) {",
      "  groups <- factor(groups)",
      "  pairs <- combn(levels(groups), 2, simplify = FALSE)",
      "  res <- do.call(rbind, lapply(pairs, function(p) {",
      "    keep <- groups %in% p",
      "    sub_dist <- as.dist(as.matrix(dist_mat)[keep, keep])",
      "    fit <- adonis2(sub_dist ~ group, data = data.frame(group = droplevels(groups[keep])), permutations = permutations)",
      '    data.frame(group1 = p[1], group2 = p[2], F = fit$F[1], R2 = fit$R2[1], p = fit[["Pr(>F)"]][1])',
      "  }))",
      '  res$p_bonferroni <- p.adjust(res$p, method = "bonferroni")',
      "  res",
      "}",
      "pairwise <- list()"
    )
    for (v in pair_vars) {
      lines <- c(lines,
        sprintf("pairwise[[%s]] <- pairwise_permanova(dist_mat, meta[[%s]], %d)", q(v), q(v), permutations),
        sprintf("print(pairwise[[%s]])", q(v))
      )
    }
  }

  paste(lines, collapse = "\n")
}

# Rounds numeric columns for display; with first_col, row names (e.g. the
# adonis2 term names) become a named first column.
permanova_display_table <- function(df, first_col = NULL) {
  out <- as.data.frame(df, check.names = FALSE)
  if (!is.null(first_col)) {
    out <- cbind(stats::setNames(data.frame(rownames(out), check.names = FALSE), first_col), out, row.names = NULL)
  }
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], function(x) signif(x, 4))
  out
}

# "Available variables" hint under a formula box: each sample_data column
# with its type. Shared with the ANCOM-BC2 step.
metadata_variable_list <- function(m) {
  items <- lapply(names(m), function(v) {
    x <- m[[v]]
    kind <- if (is.numeric(x)) "numeric" else sprintf("%d groups", length(unique(stats::na.omit(x))))
    shiny::tags$li(shiny::tags$code(v), shiny::span(class = "text-muted", sprintf(" (%s)", kind)))
  })
  shiny::div(class = "amf-param-hint", "Available variables:", shiny::tags$ul(class = "mb-0", items))
}

permanova_dt <- function(df) {
  # height = "auto": inside renderUI a widget otherwise reserves 400px.
  DT::datatable(df, rownames = FALSE, height = "auto", options = list(dom = "t", paging = FALSE, ordering = FALSE))
}

mod_permanova_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "12. PERMANOVA",
    description = "Tests whether community composition differs between groups of samples using PERMANOVA (vegan's adonis2), checks the equal-dispersion assumption with betadisper, and runs Bonferroni-corrected pairwise comparisons when a variable has more than two groups.",
    params = shiny::tagList(
      shiny::uiOutput(ns("gate_alert")),
      shiny::uiOutput(ns("input_status")),
      shiny::fluidRow(
        shiny::column(6,
          param_group("Model",
            shiny::textInput(ns("formula"),
              help_label("Formula (metadata variables)", "Right-hand side of the model, e.g. Treatment + Age. Use + for main effects and * or : for interactions. Numeric variables are treated as continuous."),
              placeholder = "Treatment + Age", width = "100%"),
            shiny::uiOutput(ns("formula_feedback")),
            shiny::uiOutput(ns("variable_list"))
          ),
          param_group("Data",
            shiny::uiOutput(ns("data_source_picker"))
          )
        ),
        shiny::column(6,
          param_group("Test settings",
            shiny::selectInput(ns("distance"),
              help_label("Distance method", "Bray-Curtis: abundance-weighted, the usual choice for community data. Jaccard: presence/absence only (binary = TRUE). Euclidean: suited to transformed data such as CLR/Hellinger."),
              choices = PERMANOVA_DISTANCES, width = "100%"),
            shiny::uiOutput(ns("distance_warning")),
            shiny::numericInput(ns("permutations"), "Permutations", value = 9999, min = 99, step = 1000, width = "100%"),
            shiny::selectInput(ns("by"),
              help_label("Test terms (by =)", "Sequential: each term is tested after the ones before it, so formula order matters. Marginal: each term is tested after all the others. Whole model: one test for the full formula."),
              choices = PERMANOVA_BY, width = "100%")
          )
        )
      ),
      param_group("Command preview",
        shiny::div(class = "amf-param-hint mb-2", "These are the exact R commands that run when you click Run PERMANOVA."),
        shiny::verbatimTextOutput(ns("script_preview"))
      )
    ),
    results = shiny::tagList(
      result_card("Run PERMANOVA",
        run_row(ns, "Run PERMANOVA"),
        shiny::h6("Log", class = "mt-3"),
        mod_logpanel_ui(ns("log"))
      ),
      shiny::uiOutput(ns("results_section")),
      results_placeholder("Run PERMANOVA to see the adonis2 table, dispersion tests and pairwise comparisons.")
    )
  )
}

# sample_table unused directly (kept for signature consistency with other
# step modules) -- this step reads phyloseq object(s) straight off disk.
mod_permanova_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    step_id <- "permanova"

    raw_path <- shiny::reactive(file.path(rv$project_dir, "06_output", "phyloseq.rds"))
    transformed_path <- shiny::reactive(file.path(rv$project_dir, "06_output", "phyloseq_transformed.rds"))

    inputs_ready <- shiny::reactive({
      shiny::req(rv$project_dir)
      rv$status$output
      fs::file_exists(raw_path())
    })
    transformed_available <- shiny::reactive({
      inputs_ready()
      fs::file_exists(transformed_path())
    })

    pkgs_ready <- shiny::reactive(requireNamespace("phyloseq", quietly = TRUE) && requireNamespace("vegan", quietly = TRUE))
    output$gate_alert <- shiny::renderUI({
      if (isTRUE(pkgs_ready())) return(NULL)
      shiny::div(class = "alert alert-warning small mb-3",
                 "phyloseq and vegan are both required for this step -- install any missing via install.packages()/BiocManager::install().")
    })

    output$input_status <- shiny::renderUI({
      if (!isTRUE(inputs_ready())) {
        return(shiny::div(class = "alert alert-warning small", "No phyloseq object found yet -- finish the Phyloseq Object step first."))
      }
      NULL
    })

    meta <- shiny::reactive({
      shiny::req(isTRUE(inputs_ready()), isTRUE(pkgs_ready()))
      as(phyloseq::sample_data(readRDS(raw_path())), "data.frame")
    })

    output$variable_list <- shiny::renderUI(metadata_variable_list(meta()))

    # Pre-fill the formula with the first categorical variable, once per
    # project, so the preview isn't empty on first visit.
    shiny::observeEvent(meta(), {
      if (nzchar(input$formula %||% "")) return()
      m <- meta()
      cat_vars <- names(m)[!vapply(m, is.numeric, logical(1))]
      if (length(cat_vars) > 0) shiny::updateTextInput(session, "formula", value = cat_vars[[1]])
    })

    output$data_source_picker <- shiny::renderUI({
      shiny::req(isTRUE(inputs_ready()))
      if (!isTRUE(transformed_available())) {
        return(shiny::div(class = "text-muted small mb-2", "Using raw counts (no transformation has been built in the Phyloseq Object step)."))
      }
      method <- attr(readRDS(transformed_path()), "amf_transform_method")
      label <- if (!is.null(method)) sprintf("Transformed (%s)", TRANSFORM_METHOD_LABELS[[method]] %||% method) else "Transformed"
      shiny::radioButtons(ns("data_source"), "Data source",
                           choices = stats::setNames(c("transformed", "raw"), c(label, "Raw counts")),
                           selected = "transformed", inline = TRUE)
    })

    ps_path <- shiny::reactive({
      if (isTRUE(transformed_available()) && identical(input$data_source %||% "transformed", "transformed")) transformed_path() else raw_path()
    })

    output$distance_warning <- shiny::renderUI({
      shiny::req(identical(ps_path(), transformed_path()), !identical(input$distance, "euclidean"))
      shiny::req(identical(attr(readRDS(transformed_path()), "amf_transform_method"), "clr"))
      shiny::div(class = "text-danger small mb-2", "CLR data contain negative values, so Bray-Curtis and Jaccard are not meaningful here. Use Euclidean (the Aitchison distance) or switch to raw counts.")
    })

    # Parse + validate the formula. rhs is rebuilt from the parsed formula
    # (not the raw text) so only a single well-formed expression ever gets
    # spliced into the script.
    formula_check <- shiny::reactive({
      m <- meta()
      txt <- trimws(sub("^\\s*~", "", input$formula %||% ""))
      if (!nzchar(txt)) return(list(ok = FALSE, msg = "Enter at least one metadata variable."))
      f <- tryCatch(stats::as.formula(paste("~", txt)), error = function(e) NULL)
      if (is.null(f) || length(f) != 2) return(list(ok = FALSE, msg = "Could not parse the formula."))
      vars <- all.vars(f)
      missing <- setdiff(vars, names(m))
      if (length(missing) > 0) {
        return(list(ok = FALSE, msg = sprintf("Not in the sample metadata: %s.", paste(missing, collapse = ", "))))
      }
      n_groups <- vapply(vars, function(v) length(unique(stats::na.omit(m[[v]]))), numeric(1))
      categorical <- vars[!vapply(m[vars], is.numeric, logical(1)) & n_groups >= 2]
      list(
        ok = TRUE, rhs = deparse1(f[[2]]), vars = vars,
        disp_vars = categorical,
        pair_vars = categorical[n_groups[categorical] > 2]
      )
    })

    output$formula_feedback <- shiny::renderUI({
      fc <- formula_check()
      if (!isTRUE(fc$ok)) return(shiny::div(class = "text-danger small mb-2", fc$msg))
      NULL
    })

    permutations <- shiny::reactive({
      n <- suppressWarnings(as.integer(input$permutations))
      if (is.na(n) || n < 99) 99L else n
    })

    script <- shiny::reactive({
      fc <- formula_check()
      shiny::req(isTRUE(fc$ok))
      build_permanova_script(ps_path(), fc$rhs, fc$vars, input$distance %||% "bray", permutations(),
                             input$by %||% "terms", fc$disp_vars, fc$pair_vars)
    })

    output$script_preview <- shiny::renderText(script())

    output$badge <- shiny::renderUI(status_badge(rv$status[[step_id]]))

    # Settings of the last launched run, so the results header describes the
    # run that produced them, not whatever the inputs say now.
    last_settings <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$run, {
      shiny::req(!any_running(rv), isTRUE(pkgs_ready()), isTRUE(inputs_ready()))
      fc <- formula_check()
      if (!isTRUE(fc$ok)) {
        shiny::showNotification(fc$msg, type = "error")
        return()
      }

      permanova_job <- function(code) {
        env <- new.env()
        eval(parse(text = code), envir = env)
        list(
          permanova = env$permanova,
          betadisper = env$betadisper_fit,
          dispersion = env$dispersion_test,
          pairwise = env$pairwise,
          n_samples = nrow(env$meta)
        )
      }

      settings <- list(
        formula = fc$rhs,
        distance = names(PERMANOVA_DISTANCES)[PERMANOVA_DISTANCES == (input$distance %||% "bray")],
        permutations = permutations(),
        by = names(PERMANOVA_BY)[PERMANOVA_BY == (input$by %||% "terms")]
      )
      launch_step_callr(rv, step_id, permanova_job, list(code = script()))
      last_settings(settings)
    })

    shiny::observeEvent(input$cancel, {
      cancel_step(rv, step_id)
    })

    output$results_section <- shiny::renderUI({
      res <- rv$artifacts[[step_id]]
      shiny::req(!is.null(res))
      s <- last_settings()

      perm_tab <- permanova_display_table(res$permanova, "Term")

      # One sentence per tested term (every row but Residual/Total): R2 as the
      # share of variation explained, plus a pointer back to betadisper when
      # that term's dispersion also differed.
      perm <- as.data.frame(res$permanova, check.names = FALSE)
      terms <- setdiff(rownames(perm), c("Residual", "Total"))
      perm_interp <- lapply(terms, function(term) {
        p <- perm[term, "Pr(>F)"]
        r2 <- sprintf("%.1f%%", 100 * perm[term, "R2"])
        label <- if (identical(term, "Model")) "The model as a whole" else term
        disp_p <- if (term %in% names(res$dispersion)) res$dispersion[[term]]$tab[1, "Pr(>F)"] else NA
        if (p < 0.05) {
          shiny::div(class = if (isTRUE(disp_p < 0.05)) "alert alert-warning small mb-2" else "alert alert-success small mb-2",
            sprintf("%s significantly explains community composition (R\u00b2 = %s of the variation, p = %s).", label, r2, signif(p, 3)),
            if (isTRUE(disp_p < 0.05)) " Its groups also differ in dispersion (see betadisper below), so part of this effect may reflect unequal spread rather than a shift in composition."
          )
        } else {
          shiny::div(class = "alert alert-secondary small mb-2",
            sprintf("%s has no significant effect on community composition (R\u00b2 = %s, p = %s).", label, r2, signif(p, 3)))
        }
      })
      if (length(terms) > 1 && identical(s$by, names(PERMANOVA_BY)[PERMANOVA_BY == "terms"])) {
        perm_interp <- c(perm_interp, list(shiny::div(class = "text-muted small",
          "Sequential tests: each term is assessed after the ones listed above it, so reordering the formula can change these results. Use Marginal to test each term after all the others.")))
      }

      disp_ui <- lapply(names(res$dispersion), function(v) {
        tab <- res$dispersion[[v]]$tab
        p <- tab[1, "Pr(>F)"]
        interp <- if (p < 0.05) {
          shiny::div(class = "alert alert-warning small",
            sprintf("Dispersion differs among %s groups (p = %s). PERMANOVA is sensitive to differences in spread, so a significant %s effect may partly reflect unequal within-group variability rather than a shift in community composition. Check the boxplot and ordination before interpreting it.",
                    v, signif(p, 3), v))
        } else {
          shiny::div(class = "alert alert-success small",
            sprintf("No significant difference in dispersion among %s groups (p = %s). The equal-dispersion assumption is not violated, so a significant PERMANOVA result for %s can be read as a difference in community composition (group centroids).",
                    v, signif(p, 3), v))
        }
        shiny::tagList(
          shiny::h6(sprintf("betadisper: %s", v), class = "mt-3"),
          permanova_dt(permanova_display_table(tab, "Source")),
          shiny::div(class = "mt-2", interp)
        )
      })

      # Which pairs differ after Bonferroni correction, strongest (largest R2)
      # first.
      pair_ui <- lapply(names(res$pairwise), function(v) {
        pw <- res$pairwise[[v]]
        sig <- pw[pw$p_bonferroni < 0.05, , drop = FALSE]
        sig <- sig[order(-sig$R2), , drop = FALSE]
        interp <- if (nrow(sig) == 0) {
          shiny::div(class = "alert alert-secondary small",
            sprintf("None of the %d %s pairs differ significantly after Bonferroni correction.", nrow(pw), v))
        } else {
          shiny::div(class = "alert alert-success small",
            sprintf("%d of %d %s pairs differ significantly after Bonferroni correction:", nrow(sig), nrow(pw), v),
            shiny::tags$ul(class = "mb-0", lapply(seq_len(nrow(sig)), function(i) shiny::tags$li(sprintf(
              "%s vs %s (R² = %.1f%%, adjusted p = %s)", sig$group1[i], sig$group2[i], 100 * sig$R2[i], signif(sig$p_bonferroni[i], 3)))))
          )
        }
        shiny::tagList(
          shiny::h6(sprintf("Pairwise PERMANOVA: %s", v), class = "mt-3"),
          permanova_dt(permanova_display_table(pw)),
          shiny::div(class = "mt-2", interp)
        )
      })

      shiny::tagList(
        result_card("PERMANOVA (adonis2)",
          if (!is.null(s)) shiny::div(class = "text-muted small mb-2", sprintf(
            "~ %s | %s | %s permutations | %s | %d samples", s$formula, s$distance,
            format(s$permutations, big.mark = ","), s$by, res$n_samples)),
          permanova_dt(perm_tab),
          shiny::div(class = "mt-2", perm_interp),
          actions = shiny::downloadButton(ns("download_permanova"), "Download CSV", class = "btn-outline-secondary btn-sm")
        ),
        result_card("Homogeneity of dispersion (betadisper)",
          if (length(disp_ui) == 0) shiny::div(class = "text-muted small", "No categorical variable in the formula, so there are no groups to compare dispersion between.")
          else shiny::tagList(disp_ui, shiny::plotOutput(ns("dispersion_plot"), height = "320px")),
          actions = if (length(disp_ui) > 0) shiny::downloadButton(ns("download_dispersion"), "Download CSV", class = "btn-outline-secondary btn-sm")
        ),
        result_card("Pairwise PERMANOVA (Bonferroni)",
          if (length(pair_ui) == 0) shiny::div(class = "text-muted small", "No categorical variable in the formula has more than two groups, so pairwise tests aren't needed.")
          else pair_ui,
          actions = if (length(pair_ui) > 0) shiny::downloadButton(ns("download_pairwise"), "Download CSV", class = "btn-outline-secondary btn-sm")
        )
      )
    })

    output$dispersion_plot <- shiny::renderPlot({
      res <- rv$artifacts[[step_id]]
      shiny::req(res, length(res$betadisper) > 0)
      df <- do.call(rbind, lapply(names(res$betadisper), function(v) {
        bd <- res$betadisper[[v]]
        data.frame(Variable = v, Group = as.character(bd$group), Distance = unname(bd$distances))
      }))
      ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Distance, fill = Group)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
        ggplot2::facet_wrap(~Variable, scales = "free_x") +
        ggplot2::labs(x = NULL, y = "Distance to group centroid") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    })

    stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
    output$download_permanova <- shiny::downloadHandler(
      filename = function() sprintf("permanova_adonis2_%s.csv", stamp()),
      content = function(file) {
        res <- rv$artifacts[[step_id]]
        utils::write.csv(permanova_display_table(res$permanova, "Term"), file, row.names = FALSE)
      }
    )
    output$download_dispersion <- shiny::downloadHandler(
      filename = function() sprintf("permanova_betadisper_%s.csv", stamp()),
      content = function(file) {
        res <- rv$artifacts[[step_id]]
        df <- do.call(rbind, lapply(names(res$dispersion), function(v) {
          cbind(Variable = v, permanova_display_table(res$dispersion[[v]]$tab, "Source"))
        }))
        utils::write.csv(df, file, row.names = FALSE)
      }
    )
    output$download_pairwise <- shiny::downloadHandler(
      filename = function() sprintf("permanova_pairwise_%s.csv", stamp()),
      content = function(file) {
        res <- rv$artifacts[[step_id]]
        df <- do.call(rbind, lapply(names(res$pairwise), function(v) cbind(Variable = v, res$pairwise[[v]])))
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    mod_logpanel_server("log", shiny::reactive(rv$log[[step_id]] %||% character(0)))
  })
}
