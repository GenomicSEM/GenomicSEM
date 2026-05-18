#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lavaan)
  library(lavaanrust)
})

parse_args <- function(args) {
  values <- list(
    n_vars = c(4L, 8L, 12L, 16L, 20L, 24L, 28L, 32L),
    fit_repeats = 3L,
    surface_repeats = 5L,
    surface_iterations = 1000L,
    seed = 20260507L
  )

  for (arg in args) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) {
      stop("Arguments must have the form key=value.", call. = FALSE)
    }

    key <- parts[[1L]]
    value <- parts[[2L]]
    if (key == "n_vars") {
      values$n_vars <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
    } else if (key == "fit_repeats") {
      values$fit_repeats <- as.integer(value)
    } else if (key == "surface_repeats") {
      values$surface_repeats <- as.integer(value)
    } else if (key == "surface_iterations") {
      values$surface_iterations <- as.integer(value)
    } else if (key == "seed") {
      values$seed <- as.integer(value)
    } else {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
  }

  values
}

make_full_covariance_problem <- function(n_vars, seed) {
  variable_names <- paste0("V", seq_len(n_vars))
  pair_grid <- t(utils::combn(variable_names, 2L))
  model_lines <- c(
    sprintf("%s ~~ %s", variable_names, variable_names),
    sprintf("%s ~~ %s", pair_grid[, 1L], pair_grid[, 2L])
  )

  set.seed(seed + n_vars)
  raw <- matrix(stats::rnorm(n_vars * n_vars), nrow = n_vars)
  sample_cov <- crossprod(raw) / n_vars
  diag(sample_cov) <- diag(sample_cov) + seq(1, 2, length.out = n_vars)
  dimnames(sample_cov) <- list(variable_names, variable_names)

  list(
    n_vars = n_vars,
    n_free = n_vars * (n_vars + 1L) / 2L,
    model = paste(model_lines, collapse = "\n"),
    sample_cov = sample_cov,
    wls_v = diag(n_vars * (n_vars + 1L) / 2L)
  )
}

fit_lavaan <- function(problem) {
  suppressWarnings(lavaan::sem(
    problem$model,
    sample.cov = problem$sample_cov,
    sample.nobs = 1000L,
    estimator = "DWLS",
    WLS.V = problem$wls_v,
    fixed.x = FALSE,
    se = "none",
    test = "none"
  ))
}

fit_rust <- function(problem) {
  lavaanrust::sem_rust(
    problem$model,
    sample.cov = problem$sample_cov,
    sample.nobs = 1000L,
    estimator = "DWLS",
    WLS.V = problem$wls_v,
    fixed.x = FALSE,
    se = "none",
    test = "none"
  )
}

timed_fit <- function(problem, backend) {
  gc()
  system.time({
    if (backend == "lavaan") {
      fit_lavaan(problem)
    } else {
      fit_rust(problem)
    }
  })[["elapsed"]]
}

compile_rust_problem <- function(problem) {
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(
    problem$model,
    problem$sample_cov,
    std.lv = FALSE
  )
  lavaanrust:::.lavaan_fast_compile_par_table(parsed, colnames(problem$sample_cov))
}

timed_rust_surface <- function(compiled, iterations) {
  free_values <- compiled$default_free_values
  gc()
  elapsed <- system.time({
    for (idx in seq_len(iterations)) {
      lavaanrust:::.lavaan_fast_implied_surfaces_rust_flat(compiled, free_values)
    }
  })[["elapsed"]]

  elapsed / iterations
}

summarize_scaling <- function(frame, group_columns) {
  medians <- stats::aggregate(
    elapsed_sec ~ .,
    data = frame[c(group_columns, "elapsed_sec")],
    FUN = stats::median
  )
  fits <- do.call(rbind, lapply(split(medians, medians$backend), function(group) {
    fit <- stats::lm(log(elapsed_sec) ~ log(n_free), data = group)
    data.frame(
      backend = group$backend[[1L]],
      intercept = unname(stats::coef(fit)[[1L]]),
      exponent = unname(stats::coef(fit)[[2L]]),
      r_squared = summary(fit)$r.squared,
      stringsAsFactors = FALSE
    )
  }))

  list(medians = medians, fits = fits)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
problems <- lapply(args$n_vars, make_full_covariance_problem, seed = args$seed)

equivalence <- do.call(rbind, lapply(problems, function(problem) {
  old <- fit_lavaan(problem)
  rust <- fit_rust(problem)
  data.frame(
    n_vars = problem$n_vars,
    n_free = problem$n_free,
    lavaan_free = max(lavaan::parTable(old)$free),
    rust_free = max(lavaanrust::parTable_rust(rust)$free),
    max_abs_cov_diff = max(abs(lavaan::fitted(old)$cov - lavaanrust::fitted_rust(rust)$cov)),
    stringsAsFactors = FALSE
  )
}))

fit_timings <- do.call(rbind, lapply(seq_len(args$fit_repeats), function(repeat_id) {
  do.call(rbind, lapply(problems, function(problem) {
    do.call(rbind, lapply(c("lavaan", "lavaanrust"), function(backend) {
      data.frame(
        repeat_id = repeat_id,
        backend = backend,
        n_vars = problem$n_vars,
        n_free = problem$n_free,
        elapsed_sec = timed_fit(problem, backend),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))

compiled_problems <- lapply(problems, compile_rust_problem)
surface_timings <- do.call(rbind, lapply(seq_len(args$surface_repeats), function(repeat_id) {
  do.call(rbind, lapply(seq_along(problems), function(idx) {
    problem <- problems[[idx]]
    data.frame(
      repeat_id = repeat_id,
      backend = "lavaanrust_surface",
      n_vars = problem$n_vars,
      n_free = problem$n_free,
      iterations = args$surface_iterations,
      elapsed_sec = timed_rust_surface(compiled_problems[[idx]], args$surface_iterations),
      stringsAsFactors = FALSE
    )
  }))
}))

fit_summary <- summarize_scaling(fit_timings, c("backend", "n_vars", "n_free"))
surface_summary <- summarize_scaling(surface_timings, c("backend", "n_vars", "n_free", "iterations"))

cat("equivalence_check\n")
print(equivalence, row.names = FALSE)
cat("fit_timings\n")
print(fit_timings, row.names = FALSE)
cat("fit_median_timings\n")
print(fit_summary$medians, row.names = FALSE)
cat("fit_log_log_scaling\n")
print(fit_summary$fits, row.names = FALSE)
cat("surface_timings\n")
print(surface_timings, row.names = FALSE)
cat("surface_median_timings\n")
print(surface_summary$medians, row.names = FALSE)
cat("surface_log_log_scaling\n")
print(surface_summary$fits, row.names = FALSE)
