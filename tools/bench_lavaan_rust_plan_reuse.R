#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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
    counts = c(1L, 10L, 100L, 1000L),
    repeats = 5L
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

  values
}

elapsed <- function(expr) {
  unname(system.time(force(expr))[["elapsed"]])
}

run_recompiled <- function(model, sample_cov, wls_v, n_fit) {
  out <- NULL
  for (idx in seq_len(n_fit)) {
    out <- lavaanrust:::.fit_lavaan_fast_ram_model(model, sample_cov, wls_v)
  }

  out
}

run_native_plan <- function(base_fit, sample_cov, wls_v, n_fit) {
  out <- NULL
  for (idx in seq_len(n_fit)) {
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

args <- parse_args(commandArgs(trailingOnly = TRUE))
fixture <- user_gwas_fixture()
model <- paste(
  "F1 =~ NA*A + l2*B + l3*C",
  "F1 ~~ 1*F1",
  "F1 ~ gamma*SNP",
  "A ~ direct*SNP",
  sep = "\n"
)
base_fit <- lavaanrust::sem_rust(
  model,
  sample.cov = fixture$sample_cov,
  estimator = "DWLS",
  WLS.V = fixture$wls_v
)

timings <- do.call(rbind, lapply(args$counts, function(n_fit) {
  do.call(rbind, lapply(seq_len(args$repeats), function(repeat_id) {
    recompiled <- NULL
    native_plan <- NULL
    recompiled_elapsed <- elapsed({
      recompiled <- run_recompiled(
        base_fit@Model@par_table,
        fixture$sample_cov,
        fixture$wls_v,
        n_fit
      )
    })
    native_plan_elapsed <- elapsed({
      native_plan <- run_native_plan(
        base_fit,
        fixture$sample_cov,
        fixture$wls_v,
        n_fit
      )
    })

    data.frame(
      n_fit = n_fit,
      repeat_id = repeat_id,
      recompiled_s = recompiled_elapsed,
      native_plan_s = native_plan_elapsed,
      max_abs_est_diff = max(abs(
        lavaanrust::parTable_rust(recompiled)$est -
          lavaanrust::parTable_rust(native_plan)$est
      )),
      stringsAsFactors = FALSE
    )
  }))
}))

summary <- aggregate(
  cbind(recompiled_s, native_plan_s, max_abs_est_diff) ~ n_fit,
  data = timings,
  FUN = stats::median
)
summary$speedup <- summary$recompiled_s / summary$native_plan_s

cat("timings\n")
print(timings, row.names = FALSE)
cat("\nsummary\n")
print(summary, row.names = FALSE)
