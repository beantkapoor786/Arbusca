# Launch/poll/cancel a non-blocking job per step, streaming its
# stdout/stderr into rv$log[[step_id]]. Two runner backends share this same
# interface:
#   - launch_step()      processx, for external binaries (cutadapt, blastn)
#   - launch_step_callr() callr::r_bg, for heavy in-R work (DADA2)
# callr's r_bg handle is itself a processx handle under the hood, so
# poll_step()/cancel_step() work on either without knowing which backs a
# given step -- they only special-case r_process to pull back its R return
# value into rv$artifacts[[step_id]].

launch_step <- function(rv, step_id, cmd, args = character(), workdir = NULL) {
  reset_downstream(rv, step_id)
  handle <- processx::process$new(cmd, args, stdout = "|", stderr = "|", wd = workdir)
  rv$handles[[step_id]] <- handle
  rv$status[[step_id]] <- "RUNNING"
  rv$log[[step_id]] <- character(0)
  rv$running_step <- step_id
  invisible(handle)
}

# func: a function run in a fresh background R session; its return value
# ends up in rv$artifacts[[step_id]] on success. func must load any packages
# it needs itself (the background session starts bare).
launch_step_callr <- function(rv, step_id, func, args = list()) {
  reset_downstream(rv, step_id)
  handle <- callr::r_bg(func = func, args = args, stdout = "|", stderr = "|", supervise = TRUE)
  rv$handles[[step_id]] <- handle
  rv$status[[step_id]] <- "RUNNING"
  rv$log[[step_id]] <- character(0)
  rv$running_step <- step_id
  invisible(handle)
}

poll_step <- function(rv, step_id) {
  handle <- rv$handles[[step_id]]
  if (is.null(handle)) return(invisible())

  new_out <- tryCatch(handle$read_output_lines(), error = function(e) character(0))
  new_err <- tryCatch(handle$read_error_lines(), error = function(e) character(0))
  new_lines <- c(new_out, new_err)
  if (length(new_lines) > 0) {
    rv$log[[step_id]] <- c(rv$log[[step_id]], new_lines)
  }

  if (!handle$is_alive()) {
    if (inherits(handle, "r_process")) {
      # A callr background job can exit 0 yet have raised an R error inside
      # func(); get_result() re-raises that error here rather than returning.
      result <- tryCatch(handle$get_result(), error = function(e) e)
      if (inherits(result, "error")) {
        rv$status[[step_id]] <- "ERROR"
        rv$log[[step_id]] <- c(rv$log[[step_id]], paste("ERROR:", conditionMessage(result)))
      } else {
        rv$status[[step_id]] <- "SUCCESS"
        rv$artifacts[[step_id]] <- result
      }
    } else {
      exit_status <- handle$get_exit_status()
      rv$status[[step_id]] <- if (identical(exit_status, 0L)) "SUCCESS" else "ERROR"
    }
    rv$handles[[step_id]] <- NULL
    if (identical(rv$running_step, step_id)) rv$running_step <- NULL
  }
}

cancel_step <- function(rv, step_id) {
  handle <- rv$handles[[step_id]]
  if (!is.null(handle) && handle$is_alive()) handle$kill()
  rv$status[[step_id]] <- "IDLE"
  rv$log[[step_id]] <- character(0)
  rv$handles[[step_id]] <- NULL
  if (identical(rv$running_step, step_id)) rv$running_step <- NULL
}
