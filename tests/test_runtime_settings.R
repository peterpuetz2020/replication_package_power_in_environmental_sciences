source(file.path("scripts", "runtime_settings.R"))

settings <- new.env(parent = emptyenv())

resolve_runtime_setting("n_iterations", 1000L, settings)
stopifnot(identical(settings$n_iterations, 1000L))

## Simulate editing the default in analysis_setup.R and sourcing it again in
## the same R session.
resolve_runtime_setting("n_iterations", 10L, settings)
stopifnot(identical(settings$n_iterations, 10L))

## An explicit caller override must survive every subsequent setup source.
settings$n_iterations <- 25L
resolve_runtime_setting("n_iterations", 1000L, settings)
resolve_runtime_setting("n_iterations", 1000L, settings)
stopifnot(identical(settings$n_iterations, 25L))

cat("Runtime setting precedence checks passed.\n")
