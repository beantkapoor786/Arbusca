# Shared live-log console: monospace, auto-scrolling, reads from a reactive
# character-vector source (rv$log[[step_id]]).

mod_logpanel_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tags$div(
    class = "amf-log-console",
    id = ns("console"),
    shiny::uiOutput(ns("log_out"))
  )
}

# log_lines: a reactive() returning character() of log lines for one step.
mod_logpanel_server <- function(id, log_lines) {
  shiny::moduleServer(id, function(input, output, session) {
    output$log_out <- shiny::renderUI({
      lines <- log_lines()
      if (length(lines) == 0) {
        return(shiny::tags$span(class = "text-muted", "(no output yet)"))
      }
      shiny::HTML(paste(htmltools::htmlEscape(lines), collapse = "<br/>"))
    })

    shiny::observe({
      log_lines()
      session$sendCustomMessage("amf-scroll-log", session$ns("console"))
    })
  })
}
