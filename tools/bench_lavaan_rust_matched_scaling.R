#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lavaan)
  library(lavaanrust)
})

parse_args <- function(args) {
  values <- list(
    variable_n_vars = c(6L, 8L, 12L, 16L, 20L, 24L, 28L, 32L, 40L, 48L),
    parameter_n_vars = 24L,
    parameter_extra_covars = c(0L, 24L, 48L, 72L, 96L, 120L, 144L, 176L),
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
    if (key == "variable_n_vars") {
      values$variable_n_vars <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
    } else if (key == "parameter_n_vars") {
      values$parameter_n_vars <- as.integer(value)
    } else if (key == "parameter_extra_covars") {
      values$parameter_extra_covars <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
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

lower_triangle_size <- function(n_vars) {
  n_vars * (n_vars + 1L) / 2L
}

make_pair_order <- function(n_vars, seed) {
  pair_grid <- utils::combn(seq_len(n_vars), 2L)
  set.seed(seed + 1000L + n_vars)
  pair_grid[, sample.int(ncol(pair_grid)), drop = FALSE]
}

make_two_factor_sample_cov <- function(n_vars, seed) {
  if (n_vars %% 2L != 0L || n_vars < 6L) {
    stop("n_vars must be even and at least 6.", call. = FALSE)
  }

  half <- n_vars / 2L
  loadings <- matrix(0, nrow = n_vars, ncol = 2L)
  loadings[seq_len(half), 1L] <- seq(0.62, 0.82, length.out = half)
  loadings[half + seq_len(half), 2L] <- seq(0.58, 0.78, length.out = half)

  latent_cov <- matrix(c(1, 0.30, 0.30, 1), nrow = 2L)
  residual_cov <- diag(seq(0.42, 0.74, length.out = n_vars))
  base_cov <- loadings %*% latent_cov %*% t(loadings) + residual_cov

  set.seed(seed + n_vars)
  noise <- matrix(stats::rnorm(n_vars * n_vars), nrow = n_vars)
  sample_cov <- base_cov + 0.02 * crossprod(noise) / n_vars

  variable_names <- paste0("V", seq_len(n_vars))
  dimnames(sample_cov) <- list(variable_names, variable_names)
  sample_cov
}

make_two_factor_model <- function(n_vars, extra_covars, pair_order) {
  variable_names <- paste0("V", seq_len(n_vars))
  half <- n_vars / 2L
  loadings_f1 <- paste(variable_names[seq_len(half)], collapse = " + ")
  loadings_f2 <- paste(variable_names[half + seq_len(half)], collapse = " + ")

  residual_covars <- character()
  if (extra_covars > 0L) {
    selected_pairs <- pair_order[, seq_len(extra_covars), drop = FALSE]
    residual_covars <- apply(
      selected_pairs,
      2L,
      function(pair) sprintf("%s ~~ %s", variable_names[pair[[1L]]], variable_names[pair[[2L]]])
    )
  }

  paste(
    c(
      sprintf("F1 =~ %s", loadings_f1),
      sprintf("F2 =~ %s", loadings_f2),
      "F1 ~~ F2",
      residual_covars
    ),
    collapse = "\n"
  )
}

make_problem <- function(family, n_vars, extra_covars, seed) {
  pair_order <- make_pair_order(n_vars, seed)
  max_covars <- ncol(pair_order)
  if (extra_covars > max_covars) {
    stop("extra_covars exceeds the number of observed-variable pairs.", call. = FALSE)
  }

  list(
    family = family,
    n_vars = n_vars,
    extra_covars = extra_covars,
    model = make_two_factor_model(n_vars, extra_covars, pair_order),
    sample_cov = make_two_factor_sample_cov(n_vars, seed),
    wls_v = diag(lower_triangle_size(n_vars))
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
    std.lv = TRUE,
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
    std.lv = TRUE,
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
    std.lv = TRUE
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

problem_frame <- function(problem, lavaan_fit, rust_fit) {
  data.frame(
    family = problem$family,
    n_vars = problem$n_vars,
    n_stats = lower_triangle_size(problem$n_vars),
    extra_covars = problem$extra_covars,
    n_free = max(lavaan::parTable(lavaan_fit)$free),
    lavaan_free = max(lavaan::parTable(lavaan_fit)$free),
    rust_free = max(lavaanrust::parTable_rust(rust_fit)$free),
    rust_model_kind = rust_fit@model_kind,
    lavaan_converged = lavaan::inspect(lavaan_fit, "converged"),
    rust_converged = rust_fit@converged,
    max_abs_cov_diff = max(abs(lavaan::fitted(lavaan_fit)$cov - lavaanrust::fitted_rust(rust_fit)$cov)),
    stringsAsFactors = FALSE
  )
}

summarize_timings <- function(frame, group_columns) {
  stats::aggregate(
    elapsed_sec ~ .,
    data = frame[c(group_columns, "elapsed_sec")],
    FUN = stats::median
  )
}

fit_scaling <- function(frame, axis, family) {
  subset <- frame[frame$family == family, , drop = FALSE]
  do.call(rbind, lapply(split(subset, subset$backend), function(group) {
    group <- group[is.finite(group$elapsed_sec) & group$elapsed_sec > 0, , drop = FALSE]
    if (nrow(group) < 2L || length(unique(group[[axis]])) < 2L) {
      return(data.frame(
        family = family,
        axis = axis,
        backend = if (nrow(group)) group$backend[[1L]] else NA_character_,
        intercept = NA_real_,
        exponent = NA_real_,
        r_squared = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    fit_data <- data.frame(
      elapsed_sec = group$elapsed_sec,
      axis_value = group[[axis]]
    )
    fit <- stats::lm(log(elapsed_sec) ~ log(axis_value), data = fit_data)
    data.frame(
      family = family,
      axis = axis,
      backend = group$backend[[1L]],
      intercept = unname(stats::coef(fit)[[1L]]),
      exponent = unname(stats::coef(fit)[[2L]]),
      r_squared = summary(fit)$r.squared,
      stringsAsFactors = FALSE
    )
  }))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
variable_problems <- lapply(
  args$variable_n_vars,
  function(n_vars) make_problem("variable_scaling", n_vars, 0L, args$seed)
)
parameter_problems <- lapply(
  args$parameter_extra_covars,
  function(extra_covars) make_problem("parameter_scaling", args$parameter_n_vars, extra_covars, args$seed)
)
problems <- c(variable_problems, parameter_problems)

equivalence <- do.call(rbind, lapply(problems, function(problem) {
  problem_frame(problem, fit_lavaan(problem), fit_rust(problem))
}))

fit_timings <- do.call(rbind, lapply(seq_len(args$fit_repeats), function(repeat_id) {
  do.call(rbind, lapply(problems, function(problem) {
    do.call(rbind, lapply(c("lavaan", "lavaanrust"), function(backend) {
      data.frame(
        repeat_id = repeat_id,
        family = problem$family,
        backend = backend,
        n_vars = problem$n_vars,
        n_stats = lower_triangle_size(problem$n_vars),
        extra_covars = problem$extra_covars,
        n_free = if (problem$family == "variable_scaling") {
          2L * problem$n_vars + 1L
        } else {
          2L * problem$n_vars + 1L + problem$extra_covars
        },
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
      family = problem$family,
      backend = "lavaanrust_surface",
      n_vars = problem$n_vars,
      n_stats = lower_triangle_size(problem$n_vars),
      extra_covars = problem$extra_covars,
      n_free = if (problem$family == "variable_scaling") {
        2L * problem$n_vars + 1L
      } else {
        2L * problem$n_vars + 1L + problem$extra_covars
      },
      iterations = args$surface_iterations,
      elapsed_sec = timed_rust_surface(compiled_problems[[idx]], args$surface_iterations),
      stringsAsFactors = FALSE
    )
  }))
}))

fit_medians <- summarize_timings(
  fit_timings,
  c("family", "backend", "n_vars", "n_stats", "extra_covars", "n_free")
)
surface_medians <- summarize_timings(
  surface_timings,
  c("family", "backend", "n_vars", "n_stats", "extra_covars", "n_free", "iterations")
)
fit_scaling_summary <- rbind(
  fit_scaling(fit_medians, "n_vars", "variable_scaling"),
  fit_scaling(fit_medians, "n_free", "variable_scaling"),
  fit_scaling(fit_medians, "n_free", "parameter_scaling")
)
surface_scaling_summary <- rbind(
  fit_scaling(surface_medians, "n_vars", "variable_scaling"),
  fit_scaling(surface_medians, "n_free", "variable_scaling"),
  fit_scaling(surface_medians, "n_free", "parameter_scaling")
)

cat("equivalence_check\n")
print(equivalence, row.names = FALSE)
cat("fit_timings\n")
print(fit_timings, row.names = FALSE)
cat("fit_median_timings\n")
print(fit_medians, row.names = FALSE)
cat("fit_scaling\n")
print(fit_scaling_summary, row.names = FALSE)
cat("surface_timings\n")
print(surface_timings, row.names = FALSE)
cat("surface_median_timings\n")
print(surface_medians, row.names = FALSE)
cat("surface_scaling\n")
print(surface_scaling_summary, row.names = FALSE)
