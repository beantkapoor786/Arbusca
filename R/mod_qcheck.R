# Step 2 -- Quality Check. Pure visualization of the raw (untrimmed) reads'
# quality profiles, so the user can eyeball run quality before anything is
# processed -- distinct from Step 5's "Quality Control" (mod_qc.R), which
# filters/trims and happens to preview the *trimmed* reads. Not a pipeline
# stage with artifacts, so it carries no rv$status entry and is always
# treated as complete (see Arbusca.R's setup/qcheck special-casing).

# Interactive version of a dada2::plotQualityProfile() ggplot (also used by
# Step 5's preview), rebuilt from the data inside it rather than converted
# with ggplotly(): ggplotly splits the heatmap into hundreds of partial
# heatmaps that plotly stretches across the gaps between tiles, drawing solid
# grey blocks where the real plot has sparse cells. Same marks and colors as
# DADA2's plot, one panel per file (or one when aggregated, where the panel
# column is `label` and the read count is a fixed label, not per-file data);
# hover shows the values.
quality_plotly <- function(p) {
  panel_col <- names(p$facet$params$facets)[1]
  tiles <- p$data
  lines <- Filter(function(l) inherits(l$geom, "GeomLine"), p$layers)
  stats <- lines[[1]]$data
  has_cum <- any(vapply(lines, function(l) identical(rlang::as_label(l$mapping$y), "Cum"), logical(1)))
  text_layer <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)[[1]]
  in_panel <- function(df, f) if (is.null(df[[panel_col]])) df else df[df[[panel_col]] == f, , drop = FALSE]
  reads_label <- function(f) text_layer$aes_params$label %||% in_panel(text_layer$data, f)$rclabel
  files <- if (is.factor(tiles[[panel_col]])) levels(tiles[[panel_col]]) else sort(unique(tiles[[panel_col]]))
  n <- length(files)
  n_col <- ceiling(sqrt(n))  # facet_wrap's default grid
  x_range <- range(tiles$Cycle) + c(-1, 1)
  y_range <- c(-1, max(tiles$Score) + 1)
  z_range <- range(tiles$Count)

  panels <- lapply(files, function(f) {
    t <- in_panel(tiles, f)
    s <- in_panel(stats, f)
    cyc <- sort(unique(t$Cycle))
    sc <- seq(min(t$Score), max(t$Score))
    z <- matrix(NA_real_, length(sc), length(cyc))  # full grid: empty cells stay blank
    z[cbind(match(t$Score, sc), match(t$Cycle, cyc))] <- t$Count
    add_stat <- function(pl, y, color, width, dash, label) {
      plotly::add_lines(pl, x = s$Cycle, y = y, line = list(color = color, width = width, dash = dash),
                        hovertemplate = paste0("Cycle %{x}: ", label, " Q%{y:.1f}<extra></extra>"))
    }
    pl <- plotly::plot_ly() |>
      plotly::add_heatmap(x = cyc, y = sc, z = z, zmin = z_range[1], zmax = z_range[2],
                          colorscale = list(c(0, "#F5F5F5"), c(1, "black")), showscale = FALSE,
                          hovertemplate = "Cycle %{x}, Q%{y}: %{z} reads<extra></extra>") |>
      add_stat(s$Mean, "#66C2A5", 1.5, "solid", "mean") |>
      add_stat(s$Q25, "#FC8D62", 0.8, "dash", "25th percentile") |>
      add_stat(s$Q50, "#FC8D62", 0.8, "solid", "median") |>
      add_stat(s$Q75, "#FC8D62", 0.8, "dash", "75th percentile")
    if (has_cum) {
      # Cum is the fraction of reads at least this long, scaled x10 to share the Q axis.
      pl <- plotly::add_lines(pl, x = s$Cycle, y = s$Cum, customdata = 10 * s$Cum, line = list(color = "red", width = 0.8),
                              hovertemplate = "Cycle %{x}: %{customdata:.1f}% of reads this long<extra></extra>")
    }
    pl
  })

  out <- plotly::subplot(panels, nrows = ceiling(n / n_col), margin = c(0.02, 0.02, 0.05, 0.05))
  axes <- list()
  notes <- list()
  for (i in seq_len(n)) {
    k <- if (i == 1) "" else i
    axis <- list(showline = TRUE, mirror = TRUE, linecolor = "#999999", zeroline = FALSE, showgrid = FALSE)
    axes[[paste0("xaxis", k)]] <- c(axis, list(range = x_range, title = list(text = if (i + n_col > n) "Cycle" else "")))
    axes[[paste0("yaxis", k)]] <- c(axis, list(range = y_range, title = list(text = if ((i - 1) %% n_col == 0) "Quality Score" else "")))
    notes[[length(notes) + 1]] <- list(text = files[i], xref = paste0("x", k, " domain"), yref = paste0("y", k, " domain"),
                                        x = 0.5, y = 1, yanchor = "bottom", showarrow = FALSE, font = list(size = 11))
    notes[[length(notes) + 1]] <- list(text = reads_label(files[i]), xref = paste0("x", k), yref = paste0("y", k),
                                        x = x_range[1], y = 0, xanchor = "left", yanchor = "bottom", showarrow = FALSE,
                                        font = list(color = "red"))
  }
  out <- do.call(plotly::layout, c(list(out), axes, list(
    annotations = notes, showlegend = FALSE, hovermode = "closest",
    paper_bgcolor = "white", plot_bgcolor = "white", margin = list(t = 30)
  )))
  plotly::config(out, displaylogo = FALSE, modeBarButtonsToRemove = c("lasso2d", "select2d"))
}

mod_qcheck_ui <- function(id) {
  ns <- shiny::NS(id)
  step_card(
    stacked = TRUE,
    title = "2. Quality Check",
    description = "Plots the base-quality profile of your raw reads so you can judge run quality before anything is trimmed. Nothing is written to disk in this step.",
    params = shiny::tagList(
      shiny::uiOutput(ns("input_status")),
      param_group("Profile options",
        shiny::fluidRow(
          shiny::column(4, shiny::numericInput(ns("n_samples"), "Number of samples to profile", value = 2, min = 1, width = "100%")),
          shiny::column(8,
            shiny::radioButtons(ns("view_mode"), "View", choices = c("Per sample" = "per_sample", "Aggregated" = "aggregate"), selected = "per_sample", inline = TRUE),
            shiny::div(class = "amf-param-hint", "Per sample facets each sample separately; Aggregated pools them into one combined profile.")
          )
        )
      ),
      shiny::div(
        class = "amf-run-row",
        shiny::actionButton(ns("generate"), "Generate quality profiles", icon = bsicons::bs_icon("play-fill"), class = "btn-primary"),
        shiny::uiOutput(ns("badge"), inline = TRUE)
      )
    ),
    results = shiny::tagList(
      result_card("Raw read quality profiles",
        shiny::uiOutput(ns("profile_error_alert")),
        shiny::uiOutput(ns("plots")),
        actions = shiny::uiOutput(ns("downloads"))
      )
    )
  )
}

# sample_table: reactive() -> data.frame(sample, fwd, rev) of basenames from
# mod_setup. This step reads the raw files directly in project_dir, before
# any primer trimming or filtering has happened.
mod_qcheck_server <- function(id, rv, sample_table) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw_files <- shiny::reactive({
      st <- sample_table()
      shiny::req(rv$project_dir)
      raw_files_df(rv$project_dir, st)
    })

    output$input_status <- shiny::renderUI({
      files <- raw_files()
      if (nrow(files) == 0) {
        return(shiny::div(class = "alert alert-warning small", "No raw read pairs found yet -- finish Setup & Raw Reads first."))
      }
      shiny::div(class = "text-muted small mb-2", sprintf("%d raw sample pair(s) available.", nrow(files)))
    })

    # Rendered as a background callr job (like Step 6's "Learn error rate")
    # so the RUNNING badge actually shows while plotQualityProfile() -- which
    # can take a while over many raw reads -- is still computing, instead of
    # freezing the UI with no feedback until it's done.
    profile_status <- shiny::reactiveVal("IDLE")
    profile_handle <- shiny::reactiveVal(NULL)
    profile_result <- shiny::reactiveVal(NULL)  # list(fwd=, rev=) ggplot objects
    profile_error <- shiny::reactiveVal(NULL)
    profile_n <- shiny::reactiveVal(1)  # facets per plot (1 when aggregated)
    output$profile_error_alert <- shiny::renderUI({
      shiny::req(identical(profile_status(), "ERROR"), profile_error())
      shiny::div(class = "alert alert-danger small mb-2", profile_error())
    })

    shiny::observeEvent(input$generate, {
      shiny::req(!identical(profile_status(), "RUNNING"))
      files <- raw_files()
      shiny::req(nrow(files) > 0)
      n <- max(1, min(input$n_samples, nrow(files)))
      aggregate <- identical(input$view_mode, "aggregate")

      profile_job <- function(fwd_files, rev_files, aggregate) {
        library(dada2)
        list(
          fwd = plotQualityProfile(fwd_files, aggregate = aggregate),
          rev = plotQualityProfile(rev_files, aggregate = aggregate)
        )
      }

      h <- callr::r_bg(
        func = profile_job,
        args = list(fwd_files = files$fwd[seq_len(n)], rev_files = files$rev[seq_len(n)], aggregate = aggregate),
        stdout = "|", stderr = "|", supervise = TRUE
      )
      profile_handle(h)
      profile_n(if (aggregate) 1 else n)
      profile_status("RUNNING")
      profile_result(NULL)
    })

    shiny::observe({
      h <- profile_handle()
      shiny::req(!is.null(h), identical(profile_status(), "RUNNING"))
      shiny::invalidateLater(500, session)
      if (!h$is_alive()) {
        result <- tryCatch(h$get_result(), error = function(e) e)
        if (inherits(result, "error")) {
          profile_error(conditionMessage(result))
          profile_status("ERROR")
        } else {
          profile_error(NULL)
          profile_result(result)
          profile_status("SUCCESS")
        }
      }
    })

    output$badge <- shiny::renderUI(status_badge(profile_status()))

    # plotQualityProfile facets with facet_wrap's default grid
    # (ceiling(sqrt(n)) columns), so size the plot to its facet rows. Beyond
    # 2 samples the facets get cramped in half-width columns, so Forward and
    # Reverse are stacked full-width instead of side by side.
    plot_height <- shiny::reactive({
      n <- profile_n()
      max(360, ceiling(n / ceiling(sqrt(n))) * 300)
    })
    output$plots <- shiny::renderUI({
      width <- if (profile_n() > 2) 12 else 6
      shiny::fluidRow(
        shiny::column(width, shiny::h6("Forward"), plotly::plotlyOutput(ns("quality_plot_fwd"), height = paste0(plot_height(), "px"))),
        shiny::column(width, shiny::h6("Reverse"), plotly::plotlyOutput(ns("quality_plot_rev"), height = paste0(plot_height(), "px")))
      )
    })
    output$quality_plot_fwd <- plotly::renderPlotly({
      shiny::req(profile_result())
      quality_plotly(profile_result()$fwd)
    })
    output$quality_plot_rev <- plotly::renderPlotly({
      shiny::req(profile_result())
      quality_plotly(profile_result()$rev)
    })

    # Only offer downloads once profiles actually exist for this session.
    output$downloads <- shiny::renderUI({
      shiny::req(profile_result())
      shiny::div(
        class = "d-flex gap-2 mt-2",
        shiny::downloadButton(ns("dl_fwd"), "Download forward plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm"),
        shiny::downloadButton(ns("dl_rev"), "Download reverse plot (300 dpi PNG)", class = "btn-outline-secondary btn-sm")
      )
    })

    output$dl_fwd <- shiny::downloadHandler(
      filename = function() "quality_profile_forward.png",
      content = function(file) ggplot2::ggsave(file, plot = profile_result()$fwd, dpi = 300, width = 8, height = 6, units = "in")
    )
    output$dl_rev <- shiny::downloadHandler(
      filename = function() "quality_profile_reverse.png",
      content = function(file) ggplot2::ggsave(file, plot = profile_result()$rev, dpi = 300, width = 8, height = 6, units = "in")
    )

    # Exposed so the stepper badge can show SUCCESS only once profiles have
    # actually been generated -- not on launch, before anything has run.
    done <- shiny::reactive(identical(profile_status(), "SUCCESS"))

    list(done = done)
  })
}
