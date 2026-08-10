## Resolve editable defaults without losing caller-supplied runtime overrides.
##
## Analysis scripts source analysis_setup.R more than once.  Remember whether a
## value came from this file so an edited default takes effect when setup is
## sourced again, while a value explicitly assigned by the caller remains an
## override throughout the workflow.
resolve_runtime_setting <- function(name, default, envir = parent.frame()) {
  state_name <- paste0(".analysis_setup_state_", name)
  previous_state <- get0(state_name, envir = envir, inherits = FALSE)
  has_value <- exists(name, envir = envir, inherits = TRUE)

  use_default <- !has_value || (
    is.list(previous_state) &&
      isTRUE(previous_state$from_default) &&
      identical(get(name, envir = envir, inherits = TRUE), previous_state$value)
  )
  value <- if (use_default) default else get(name, envir = envir, inherits = TRUE)

  assign(name, value, envir = envir)
  assign(
    state_name,
    list(value = value, from_default = use_default),
    envir = envir
  )
  invisible(value)
}
