#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lavaan)
  library(lavaanrust)
})

source("lavaanrust/tests/testthat/helper-fixtures-user-gwas.R")

parse_counts <- function(value) {
  counts <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (!length(counts) || any(!is.finite(counts)) || any(counts <= 0L)) {
    stop("counts must be a comma-separated list of positive integers.", call. = FALSE)
  }

  counts
}

parse_args <- function(args) {
  values <- list(
    counts = c(1L, 10L, 100L, 1000L, 10000L),
    repeats = 3L
  )

  for (arg in args) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) {
      stop("Arguments must have the form key=value.", call. = FALSE)
    }

    key <- parts[[1L]]
    value <- parts[[2L]]
    if (key == "counts") {
      values$counts <- parse_counts(value)
    } else if (key == "repeats") {
      values$repeats <- as.integer(value)
    } else {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
  }

  if (!is.finite(values$repeats) || values$repeats <= 0L) {
    stop("repeats must be a positive integer.", call. = FALSE)
  }

  values
}

elapsed <- function(expr) {
  unname(system.time(force(expr))[["elapsed"]])
}

make_sample_covs <- function(base, n_fit) {
  stopifnot(identical(colnames(base), c("SNP", "A", "B", "C")))
  lapply(seq_len(n_fit), function(idx) {
    sample_cov <- base
    sample_cov["SNP", "SNP"] <- 0.42 + 0.004 * cos(idx / 17)
    snp_cov <- c(
      0.05 + 0.004 * sin(idx / 13),
      0.04 + 0.003 * cos(idx / 11),
      0.03 + 0.002 * sin(idx / 7)
    )
    sample_cov[c("A", "B", "C"), "SNP"] <- snp_cov
    sample_cov["SNP", c("A", "B", "C")] <- snp_cov

    sample_cov
  })
}

lavaan_reuse_fit <- function(base_fit, sample_cov, wls_v) {
  suppressWarnings(lavaan::lavaan(
    sample.cov = sample_cov,
    WLS.V = wls_v,
    ordered = NULL,
    sampling.weights = NULL,
    se = "standard",
    sample.mean = NULL,
    sample.th = NULL,
    sample.nobs = 2,
    group = NULL,
    cluster = NULL,
    constraints = "",
    NACOV = NULL,
    slotOptions = base_fit@Options,
    slotParTable = base_fit@ParTable,
    slotSampleStats = NULL,
    slotData = base_fit@Data,
    slotModel = base_fit@Model,
    slotCache = NULL,
    sloth1 = NULL
  ))
}

run_lavaan <- function(base_fit, sample_covs, wls_v) {
  out <- NULL
  for (sample_cov in sample_covs) {
    out <- lavaan_reuse_fit(base_fit, sample_cov, wls_v)
  }

  out
}

run_lavaanrust_recompiled <- function(base_fit, sample_covs, wls_v) {
  out <- NULL
  for (sample_cov in sample_covs) {
    out <- lavaanrust:::.fit_lavaan_fast_ram_model(
      base_fit@Model@par_table,
      sample_cov,
      wls_v
    )
  }

  out
}

run_lavaanrust_plan_reuse <- function(base_fit, sample_covs, wls_v) {
  out <- NULL
  for (sample_cov in sample_covs) {
    out <- lavaanrust::lavaan_rust(
      sample.cov = sample_cov,
      WLS.V = wls_v,
      slotOptions = base_fit@Options,
      slotParTable = base_fit@ParTable,
      slotData = base_fit@Data,
      slotModel = base_fit@Model
    )
  }

  out
}

cov_diff <- function(old_fit, rust_fit) {
  max(abs(lavaan::fitted(old_fit)$cov - lavaanrust::fitted_rust(rust_fit)$cov))
}

linear_summary <- function(summary) {
  counts <- summary$n_fit >= 10L
  if (sum(counts) < 2L) {
    return(data.frame(
      backend = c("lavaan", "lavaanrust_recompiled", "lavaanrust_plan_reuse"),
      intercept_ms = NA_real_,
      slope_ms_per_fit = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    backend = c("lavaan", "lavaanrust_recompiled", "lavaanrust_plan_reuse"),
    intercept_ms = c(
      coef(stats::lm(lavaan_s ~ n_fit, data = summary[counts, ]))[[1L]] * 1000,
      coef(stats::lm(lavaanrust_recompiled_s ~ n_fit, data = summary[counts, ]))[[1L]] * 1000,
      coef(stats::lm(lavaanrust_plan_reuse_s ~ n_fit, data = summary[counts, ]))[[1L]] * 1000
    ),
    slope_ms_per_fit = c(
      coef(stats::lm(lavaan_s ~ n_fit, data = summary[counts, ]))[[2L]] * 1000,
      coef(stats::lm(lavaanrust_recompiled_s ~ n_fit, data = summary[counts, ]))[[2L]] * 1000,
      coef(stats::lm(lavaanrust_plan_reuse_s ~ n_fit, data = summary[counts, ]))[[2L]] * 1000
    ),
    stringsAsFactors = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
fixture <- user_gwas_fixture()
model <- paste(
  "F1 =~ NA*A + l2*B + l3*C",
  "F1 ~~ 1*F1",
  "F1 ~ gamma*SNP",
  "A ~ direct*SNP",
  "combo := gamma * direct",
  sep = "\n"
)
sample_covs <- make_sample_covs(fixture$sample_cov, max(args$counts))
wls_v <- fixture$wls_v
base_lavaan <- suppressWarnings(lavaan::sem(
  model,
  sample.cov = sample_covs[[1L]],
  estimator = "DWLS",
  se = "standard",
  WLS.V = wls_v,
  sample.nobs = 2,
  optim.dx.tol = 0.01
))
base_lavaanrust <- lavaanrust::sem_rust(
  model,
  sample.cov = sample_covs[[1L]],
  estimator = "DWLS",
  WLS.V = wls_v
)

# Warm each backend once before timing.
invisible(run_lavaan(base_lavaan, sample_covs[1L], wls_v))
invisible(run_lavaanrust_recompiled(base_lavaanrust, sample_covs[1L], wls_v))
invisible(run_lavaanrust_plan_reuse(base_lavaanrust, sample_covs[1L], wls_v))

timings <- do.call(rbind, lapply(args$counts, function(n_fit) {
  inputs <- sample_covs[seq_len(n_fit)]
  do.call(rbind, lapply(seq_len(args$repeats), function(repeat_id) {
    old_fit <- NULL
    rust_recompiled_fit <- NULL
    rust_plan_fit <- NULL
    old_elapsed <- elapsed({
      old_fit <- run_lavaan(base_lavaan, inputs, wls_v)
    })
    recompiled_elapsed <- elapsed({
      rust_recompiled_fit <- run_lavaanrust_recompiled(base_lavaanrust, inputs, wls_v)
    })
    plan_elapsed <- elapsed({
      rust_plan_fit <- run_lavaanrust_plan_reuse(base_lavaanrust, inputs, wls_v)
    })

    data.frame(
      n_fit = n_fit,
      repeat_id = repeat_id,
      lavaan_s = old_elapsed,
      lavaanrust_recompiled_s = recompiled_elapsed,
      lavaanrust_plan_reuse_s = plan_elapsed,
      max_abs_cov_diff_lavaan_vs_recompiled = cov_diff(old_fit, rust_recompiled_fit),
      max_abs_cov_diff_lavaan_vs_plan = cov_diff(old_fit, rust_plan_fit),
      max_abs_est_diff_recompiled_vs_plan = max(abs(
        lavaanrust::parTable_rust(rust_recompiled_fit)$est -
          lavaanrust::parTable_rust(rust_plan_fit)$est
      )),
      stringsAsFactors = FALSE
    )
  }))
}))

summary <- aggregate(
  cbind(
    lavaan_s,
    lavaanrust_recompiled_s,
    lavaanrust_plan_reuse_s,
    max_abs_cov_diff_lavaan_vs_recompiled,
    max_abs_cov_diff_lavaan_vs_plan,
    max_abs_est_diff_recompiled_vs_plan
  ) ~ n_fit,
  data = timings,
  FUN = stats::median
)
summary$recompiled_speedup_vs_lavaan <- summary$lavaan_s / summary$lavaanrust_recompiled_s
summary$plan_speedup_vs_lavaan <- summary$lavaan_s / summary$lavaanrust_plan_reuse_s
summary$plan_speedup_vs_recompiled <- summary$lavaanrust_recompiled_s / summary$lavaanrust_plan_reuse_s

cat("timings\n")
print(timings, row.names = FALSE)
cat("\nsummary\n")
print(summary, row.names = FALSE)
cat("\nlinear_summary\n")
print(linear_summary(summary), row.names = FALSE)
