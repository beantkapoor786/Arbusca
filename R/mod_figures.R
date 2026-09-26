# Step 10 -- Figures. Basic microbial ecology plots built on top of Step 9's
# phyloseq object(s): rarefaction curves + alpha diversity always use the raw
# (untransformed) object (diversity indices assume integer counts); beta
# diversity ordinations let the user pick raw vs. transformed (Bray-Curtis
# on counts/proportions, Aitchison distance on CLR -- see beta_distance()).
#
# Pure visualization, like Step 2's Quality Check -- not part of the
# rv$status pipeline FSM (it's the last step; nothing downstream to gate).

NO_GROUP <- "__none__"

# samples-x-taxa count/abundance matrix regardless of the phyloseq object's
# own row/column orientation.
otu_matrix <- function(ps) {
  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) mat <- t(mat)
  mat
}

sample_group_values <- function(ps, group_col) {
  if (is.null(group_col) || identical(group_col, NO_GROUP)) {
    return(stats::setNames(rep("All samples", phyloseq::nsamples(ps)), phyloseq::sample_names(ps)))
  }
  meta <- as(phyloseq::sample_data(ps), "data.frame")
  stats::setNames(as.character(meta[[group_col]]), rownames(meta))
}

alpha_diversity_df <- function(ps, group_col) {
  rich <- phyloseq::estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
  groups <- sample_group_values(ps, group_col)
  # estimate_richness() mangles sample names via make.names() (e.g.
  # "1-827-Root" -> "X1.827.Root"), so take them from ps -- same row order.
  rich$Sample <- phyloseq::sample_names(ps)
  rich$Group <- unname(groups[rich$Sample])
  do.call(rbind, lapply(c("Observed", "Shannon", "Simpson"), function(m) {
    data.frame(Sample = rich$Sample, Group = rich$Group, Measure = m, Value = rich[[m]], stringsAsFactors = FALSE)
  }))
}

alpha_diversity_plot <- function(df) {
  df$Measure <- factor(df$Measure, levels = c("Observed", "Shannon", "Simpson"))
  ggplot2::ggplot(df, ggplot2::aes(x = Group, y = Value, fill = Group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
    ggplot2::facet_wrap(~Measure, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Diversity index") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = if (length(unique(df$Group)) > 1) "right" else "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

# Bray-Curtis for counts/proportions; Euclidean (= Aitchison distance) for
# CLR data, whose negative values make Bray-Curtis meaningless.
beta_distance <- function(ps) {
  if (identical(attr(ps, "amf_transform_method"), "clr")) c(euclidean = "Aitchison") else c(bray = "Bray-Curtis")
}

# method: "NMDS" or "PCoA", against beta_distance(ps).
beta_ordination_plot <- function(ps, method, group_col) {
  dist <- beta_distance(ps)
  ord <- phyloseq::ordinate(ps, method = method, distance = names(dist))
  color <- if (!is.null(group_col) && !identical(group_col, NO_GROUP)) group_col else NULL
  p <- phyloseq::plot_ordination(ps, ord, color = color) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(sprintf("%s (%s)", method, dist))
  if (identical(method, "NMDS") && !is.null(ord$stress)) {
    p <- p + ggplot2::labs(subtitle = sprintf("Stress = %.3f", ord$stress))
  }
  list(plot = p, ord = ord)
}

mod_figures_ui <- function(id) {
  ns <- shiny::NS(id)
  gen_btn <- function(id, label) shiny::actionButton(ns(id), label, icon = bsicons::bs_icon("play-fill"), class = "btn-primary")
  hint <- function(text) shiny::div(class = "amf-param-hint mb-2", text)
  # One self-contained card per figure: what it is, its controls, Generate,
  # then the plot; downloads sit in the card header once a plot exists.
  step_card(
    title = "10. Figures",
    description = "Generates rarefaction curves and alpha- and beta-diversity plots from the phyloseq object. Each figure is created on demand and can be downloaded as a 300 dpi PNG.",
    results = shiny::tagList(
      shiny::uiOutput(ns("gate_alert")),
      shiny::uiOutput(ns("input_status")),
      result_card("Rarefaction curves",
        hint("Always computed from raw (untransformed) ASV counts -- rarefaction curves only make sense on actual read counts."),
        shiny::div(class = "amf-run-row mb-2", gen_btn("gen_rarecurve", "Generate rarefaction curves"), shiny::uiOutput(ns("rarecurve_badge"), inline = TRUE)),
        shiny::uiOutput(ns("rarecurve_error_alert")),
        shiny::plotOutput(ns("rarecurve_plot"), height = "380px"),
        actions = shiny::uiOutput(ns("rarecurve_download_ui"))
      ),
      result_card("Alpha diversity",
        hint("Observed richness, Shannon, and Simpson indices, computed from raw ASV counts."),
        shiny::fluidRow(shiny::column(4, shiny::uiOutput(ns("group_picker")))),
        shiny::div(class = "amf-run-row mb-2", gen_btn("gen_alpha", "Generate alpha diversity plots"), shiny::uiOutput(ns("alpha_badge"), inline = TRUE)),
        shiny::uiOutput(ns("alpha_error_alert")),
        shiny::plotOutput(ns("alpha_plot"), height = "420px"),
        actions = shiny::uiOutput(ns("alpha_download_ui"))
      ),
      result_card("Beta diversity",
        hint("NMDS and PCoA ordinations: Bray-Curtis dissimilarity, or Aitchison (Euclidean) distance for CLR-transformed data."),
        shiny::fluidRow(
          shiny::column(4, shiny::uiOutput(ns("beta_group_picker"))),
          shiny::column(8, shiny::uiOutput(ns("data_source_picker")))
        ),
        shiny::div(class = "amf-run-row mb-2", gen_btn("gen_beta", "Generate beta diversity plots"), shiny::uiOutput(ns("beta_badge"), inline = TRUE)),
        shiny::uiOutput(ns("beta_error_alert")),
        shiny::fluidRow(
          shiny::column(6, shiny::plotOutput(ns("nmds_plot"), height = "380px")),
          shiny::column(6, shiny::plotOutput(ns("pcoa_plot"), height = "380px"))
        ),
        actions = shiny::uiOutput(ns("beta_download_ui"))
      )
    )
  )
}

# sample_table unused directly (kept for signature consistency with other
# step modules) -- this step reads phyloseq object(s) straight off disk.
mod_figures_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw_path <- shiny::reactive(file.path(rv$project_dir, "05_output", "phyloseq.rds"))
    transformed_path <- shiny::reactive(file.path(rv$project_dir, "05_output", "phyloseq_transformed.rds"))

    inputs_ready <- shiny::reactive({
      shiny::req(rv$project_dir)
      rv$status$output
      fs::file_exists(raw_path())
    })

    output$input_status <- shiny::renderUI({
      if (!isTRUE(inputs_ready())) {
        return(shiny::div(class = "alert alert-warning small", "No phyloseq object found yet -- finish the Phyloseq Object step first."))
      }
      NULL
    })

    figures_ready <- shiny::reactive(
      requireNamespace("phyloseq", quietly = TRUE) && requireNamespace("vegan", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
    )
    output$gate_alert <- shiny::renderUI({
      if (isTRUE(figures_ready())) return(NULL)
      shiny::div(class = "alert alert-warning small mb-3",
                 "phyloseq, vegan, and ggplot2 are all required for this step -- install any missing via install.packages()/BiocManager::install().")
    })

    ps_raw <- shiny::reactive({
      shiny::req(isTRUE(inputs_ready()))
      readRDS(raw_path())
    })
    ps_transformed <- shiny::reactive({
      shiny::req(fs::file_exists(transformed_path()))
      readRDS(transformed_path())
    })
    transformed_available <- shiny::reactive({
      inputs_ready()  # re-check whenever output's status changes
      fs::file_exists(transformed_path())
    })

    group_choices <- shiny::reactive({
      shiny::req(isTRUE(inputs_ready()))
      cols <- names(as(phyloseq::sample_data(ps_raw()), "data.frame"))
      stats::setNames(c(NO_GROUP, cols), c("None", cols))
    })
    output$group_picker <- shiny::renderUI(
      shiny::selectInput(ns("group_col"), "Group by (optional)", choices = group_choices(), width = "100%")
    )
    output$beta_group_picker <- shiny::renderUI(
      shiny::selectInput(ns("beta_group_col"), "Color by (optional)", choices = group_choices(), width = "100%")
    )

    output$data_source_picker <- shiny::renderUI({
      shiny::req(isTRUE(inputs_ready()))
      if (!isTRUE(transformed_available())) {
        return(shiny::div(class = "text-muted small mb-2", "Using raw counts (no transformation has been built in the Phyloseq Object step)."))
      }
      method <- attr(ps_transformed(), "amf_transform_method")
      label <- if (!is.null(method)) sprintf("Transformed (%s)", TRANSFORM_METHOD_LABELS[[method]] %||% method) else "Transformed"
      shiny::radioButtons(ns("data_source"), "Data source for ordination",
                           choices = stats::setNames(c("transformed", "raw"), c(label, "Raw counts")),
                           selected = "transformed", inline = TRUE)
    })

    beta_ps <- shiny::reactive({
      if (isTRUE(transformed_available()) && identical(input$data_source %||% "transformed", "transformed")) ps_transformed() else ps_raw()
    })

    # --- Rarefaction curves ---
    rarecurve_status <- shiny::reactiveVal("IDLE")
    rarecurve_error <- shiny::reactiveVal(NULL)
    # A counter, not a boolean: reactiveVal only invalidates dependents when
    # the new value differs from the old one, so a plain TRUE flag would stop
    # refreshing the plot on a second click once it was already TRUE (e.g.
    # after rebuilding the phyloseq object and re-generating).
    rarecurve_done <- shiny::reactiveVal(0)

    shiny::observeEvent(input$gen_rarecurve, {
      shiny::req(isTRUE(figures_ready()), isTRUE(inputs_ready()))
      rarecurve_status("RUNNING")
      ok <- tryCatch({ rarecurve_done(rarecurve_done() + 1); TRUE }, error = function(e) { rarecurve_error(conditionMessage(e)); FALSE })
      rarecurve_status(if (ok) "SUCCESS" else "ERROR")
    })

    output$rarecurve_badge <- shiny::renderUI(status_badge(rarecurve_status()))
    output$rarecurve_error_alert <- shiny::renderUI({
      shiny::req(identical(rarecurve_status(), "ERROR"), rarecurve_error())
      shiny::div(class = "alert alert-danger small mb-3", rarecurve_error())
    })

    output$rarecurve_plot <- shiny::renderPlot({
      shiny::req(rarecurve_done() > 0)
      mat <- otu_matrix(ps_raw())
      vegan::rarecurve(mat, step = 100, label = TRUE, xlab = "Sequencing depth", ylab = "ASVs observed")
    })

    output$rarecurve_download_ui <- shiny::renderUI({
      shiny::req(rarecurve_done() > 0)
      shiny::downloadButton(ns("download_rarecurve"), "Download plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm mb-2")
    })
    output$download_rarecurve <- shiny::downloadHandler(
      filename = function() "rarefaction_curves.png",
      content = function(file) {
        grDevices::png(file, width = 8, height = 6, units = "in", res = 300)
        on.exit(grDevices::dev.off())
        mat <- otu_matrix(ps_raw())
        vegan::rarecurve(mat, step = 100, label = TRUE, xlab = "Sequencing depth", ylab = "ASVs observed")
      }
    )

    # --- Alpha diversity ---
    alpha_status <- shiny::reactiveVal("IDLE")
    alpha_error <- shiny::reactiveVal(NULL)
    alpha_result <- shiny::reactiveVal(NULL)  # list(df=, plot=)

    shiny::observeEvent(input$gen_alpha, {
      shiny::req(isTRUE(figures_ready()), isTRUE(inputs_ready()))
      alpha_status("RUNNING")
      out <- tryCatch({
        df <- alpha_diversity_df(ps_raw(), input$group_col)
        list(df = df, plot = alpha_diversity_plot(df))
      }, error = function(e) e)
      if (inherits(out, "error")) {
        alpha_status("ERROR")
        alpha_error(conditionMessage(out))
      } else {
        alpha_error(NULL)
        alpha_result(out)
        alpha_status("SUCCESS")
      }
    })

    output$alpha_badge <- shiny::renderUI(status_badge(alpha_status()))
    output$alpha_error_alert <- shiny::renderUI({
      shiny::req(identical(alpha_status(), "ERROR"), alpha_error())
      shiny::div(class = "alert alert-danger small mb-3", alpha_error())
    })
    output$alpha_plot <- shiny::renderPlot({
      shiny::req(alpha_result())
      alpha_result()$plot
    })
    output$alpha_download_ui <- shiny::renderUI({
      shiny::req(alpha_result())
      shiny::div(
        class = "d-flex gap-2 mt-2 mb-2",
        shiny::downloadButton(ns("download_alpha_plot"), "Download plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm"),
        shiny::downloadButton(ns("download_alpha_table"), "Download values as CSV", class = "btn-outline-secondary btn-sm")
      )
    })
    output$download_alpha_plot <- shiny::downloadHandler(
      filename = function() "alpha_diversity.png",
      content = function(file) ggplot2::ggsave(file, plot = alpha_result()$plot, dpi = 300, width = 9, height = 5, units = "in")
    )
    output$download_alpha_table <- shiny::downloadHandler(
      filename = function() sprintf("alpha_diversity_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) utils::write.csv(alpha_result()$df, file, row.names = FALSE)
    )

    # --- Beta diversity ---
    beta_status <- shiny::reactiveVal("IDLE")
    beta_error <- shiny::reactiveVal(NULL)
    beta_result <- shiny::reactiveVal(NULL)  # list(nmds=, pcoa=) each list(plot=, ord=)

    shiny::observeEvent(input$gen_beta, {
      shiny::req(isTRUE(figures_ready()), isTRUE(inputs_ready()))
      beta_status("RUNNING")
      out <- tryCatch({
        ps <- beta_ps()
        list(
          nmds = beta_ordination_plot(ps, "NMDS", input$beta_group_col),
          pcoa = beta_ordination_plot(ps, "PCoA", input$beta_group_col)
        )
      }, error = function(e) e)
      if (inherits(out, "error")) {
        beta_status("ERROR")
        beta_error(conditionMessage(out))
      } else {
        beta_error(NULL)
        beta_result(out)
        beta_status("SUCCESS")
      }
    })

    output$beta_badge <- shiny::renderUI(status_badge(beta_status()))
    output$beta_error_alert <- shiny::renderUI({
      shiny::req(identical(beta_status(), "ERROR"), beta_error())
      shiny::div(class = "alert alert-danger small mb-3", beta_error())
    })
    output$nmds_plot <- shiny::renderPlot({
      shiny::req(beta_result())
      beta_result()$nmds$plot
    })
    output$pcoa_plot <- shiny::renderPlot({
      shiny::req(beta_result())
      beta_result()$pcoa$plot
    })
    output$beta_download_ui <- shiny::renderUI({
      shiny::req(beta_result())
      shiny::div(
        class = "d-flex gap-2 mt-2 mb-2",
        shiny::downloadButton(ns("download_nmds_plot"), "Download NMDS plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm"),
        shiny::downloadButton(ns("download_pcoa_plot"), "Download PCoA plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm")
      )
    })
    output$download_nmds_plot <- shiny::downloadHandler(
      filename = function() "beta_nmds.png",
      content = function(file) ggplot2::ggsave(file, plot = beta_result()$nmds$plot, dpi = 300, width = 7, height = 6, units = "in")
    )
    output$download_pcoa_plot <- shiny::downloadHandler(
      filename = function() "beta_pcoa.png",
      content = function(file) ggplot2::ggsave(file, plot = beta_result()$pcoa$plot, dpi = 300, width = 7, height = 6, units = "in")
    )

    # A rebuilt/re-transformed phyloseq object makes every figure stale.
    shiny::observeEvent(rv$phyloseq_rev, {
      rarecurve_status("IDLE"); rarecurve_done(0)
      alpha_status("IDLE"); alpha_result(NULL)
      beta_status("IDLE"); beta_result(NULL)
    }, ignoreInit = TRUE)

    # Exposed so the stepper badge can show SUCCESS once at least one figure
    # has been generated -- same pattern as Step 2 (Quality Check).
    done <- shiny::reactive(
      identical(rarecurve_status(), "SUCCESS") || identical(alpha_status(), "SUCCESS") || identical(beta_status(), "SUCCESS")
    )

    list(done = done)
  })
}
