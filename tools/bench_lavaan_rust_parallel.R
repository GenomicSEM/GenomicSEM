#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
})

parse_args <- function(args) {
  values <- list(
    n_snp = 1000L,
    k = 5L,
    cores = c(1L, 2L, 4L, 8L, 16L),
    repeats = 1L,
    seed = 20260504L,
    fix_measurement = c(TRUE, FALSE),
    q_snp = c(TRUE, FALSE)
  )

  for (arg in args) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    if (length(parts) != 2L) {
      stop("Arguments must have the form key=value.", call. = FALSE)
    }

    key <- parts[[1]]
    value <- parts[[2]]

    if (key == "n_snp") {
      values$n_snp <- as.integer(value)
    } else if (key == "k") {
      values$k <- as.integer(value)
    } else if (key == "cores") {
      values$cores <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
    } else if (key == "repeats") {
      values$repeats <- as.integer(value)
    } else if (key == "seed") {
      values$seed <- as.integer(value)
    } else if (key == "fix_measurement") {
      values$fix_measurement <- strsplit(value, ",", fixed = TRUE)[[1]] %in% c("TRUE", "true", "1", "yes")
    } else if (key == "q_snp") {
      values$q_snp <- strsplit(value, ",", fixed = TRUE)[[1]] %in% c("TRUE", "true", "1", "yes")
    } else {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
  }

  values$cores <- sort(unique(values$cores[values$cores > 0L]))
  values
}

make_inputs <- function(n_snp, k, seed) {
  set.seed(seed)
  traits <- paste0("T", seq_len(k))

  s_ld <- matrix(0.35, k, k)
  diag(s_ld) <- 1
  dimnames(s_ld) <- list(traits, traits)

  v_ld <- diag(1e-4, k * (k + 1L) / 2L)

  i_ld <- matrix(0.02, k, k)
  diag(i_ld) <- 1.05
  dimnames(i_ld) <- list(traits, traits)

  snps <- data.frame(
    SNP = paste0("rs", seq_len(n_snp)),
    CHR = rep(1L, n_snp),
    BP = seq_len(n_snp),
    MAF = runif(n_snp, 0.05, 0.5),
    A1 = rep("A", n_snp),
    A2 = rep("G", n_snp),
    check.names = FALSE
  )

  beta <- as.data.frame(
    matrix(rnorm(n_snp * k, sd = 0.02), n_snp, k),
    check.names = FALSE
  )
  names(beta) <- paste0("beta.", traits)

  se <- as.data.frame(
    matrix(runif(n_snp * k, min = 0.03, max = 0.08), n_snp, k),
    check.names = FALSE
  )
  names(se) <- paste0("se.", traits)

  model <- paste(
    paste0("F1 =~ ", paste(traits, collapse = " + ")),
    "F1 ~ SNP",
    sep = "\n"
  )

  list(
    covstruc = list(V = v_ld, S = s_ld, I = i_ld),
    SNPs = cbind(snps, beta, se),
    model = model
  )
}

bind_usergwas <- function(result) {
  do.call(rbind, result)
}

numeric_max_abs_diff <- function(left, right, columns) {
  max(vapply(
    columns,
    function(column) {
      max(abs(as.numeric(left[[column]]) - as.numeric(right[[column]])), na.rm = TRUE)
    },
    numeric(1)
  ))
}

usergwas_numeric_columns <- function(q_snp) {
  columns <- c("est", "SE", "chisq")
  if (q_snp) {
    columns <- c(columns, "Q_SNP", "Q_SNP_df", "Q_SNP_pval")
  }

  columns
}

run_commonfactor <- function(inputs, backend, parallel, cores) {
  fun <- if (backend == "lavaan") commonfactorGWAS else commonfactorGWAS_rust

  suppressWarnings(fun(
    covstruc = inputs$covstruc,
    SNPs = inputs$SNPs,
    estimation = "DWLS",
    parallel = parallel,
    cores = if (parallel) cores else NULL,
    GC = "standard"
  ))
}

run_usergwas <- function(inputs, backend, parallel, cores, fix_measurement, q_snp) {
  fun <- if (backend == "lavaan") userGWAS else userGWAS_rust

  suppressWarnings(fun(
    covstruc = inputs$covstruc,
    SNPs = inputs$SNPs,
    model = inputs$model,
    estimation = "DWLS",
    parallel = parallel,
    cores = if (parallel) cores else NULL,
    GC = "standard",
    fix_measurement = fix_measurement,
    Q_SNP = q_snp,
    printwarn = TRUE
  ))
}

benchmark_once <- function(workflow, backend, parallel, cores, inputs, fix_measurement = NA, q_snp = NA) {
  gc()
  elapsed <- system.time({
    result <- if (workflow == "commonfactorGWAS") {
      run_commonfactor(inputs, backend, parallel, cores)
    } else {
      run_usergwas(inputs, backend, parallel, cores, fix_measurement, q_snp)
    }
  })[["elapsed"]]

  data.frame(
    workflow = workflow,
    backend = backend,
    mode = if (parallel) "parallel" else "sequential",
    cores = if (parallel) cores else 1L,
    fix_measurement = fix_measurement,
    Q_SNP = q_snp,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

check_equivalence <- function(inputs, user_matrix) {
  commonfactor_seq <- run_commonfactor(inputs, "lavaanrust", FALSE, 1L)
  commonfactor_par <- run_commonfactor(inputs, "lavaanrust", TRUE, 2L)

  commonfactor_row <- data.frame(
    workflow = "commonfactorGWAS",
    fix_measurement = NA,
    Q_SNP = NA,
    comparison = "rust_sequential_vs_parallel",
    max_abs_diff = numeric_max_abs_diff(
      commonfactor_seq,
      commonfactor_par,
      c("est", "se_c", "Q", "Z_Estimate", "Pval_Estimate")
    ),
    stringsAsFactors = FALSE
  )

  user_rows <- do.call(rbind, lapply(seq_len(nrow(user_matrix)), function(row_idx) {
    fix_measurement <- user_matrix$fix_measurement[[row_idx]]
    q_snp <- user_matrix$Q_SNP[[row_idx]]
    columns <- usergwas_numeric_columns(q_snp)

    old_seq <- bind_usergwas(run_usergwas(inputs, "lavaan", FALSE, 1L, fix_measurement, q_snp))
    rust_seq <- bind_usergwas(run_usergwas(inputs, "lavaanrust", FALSE, 1L, fix_measurement, q_snp))
    rust_par <- bind_usergwas(run_usergwas(inputs, "lavaanrust", TRUE, 2L, fix_measurement, q_snp))

    data.frame(
      workflow = "userGWAS",
      fix_measurement = fix_measurement,
      Q_SNP = q_snp,
      comparison = c("lavaan_vs_rust_sequential", "rust_sequential_vs_parallel"),
      max_abs_diff = c(
        numeric_max_abs_diff(old_seq, rust_seq, columns),
        numeric_max_abs_diff(rust_seq, rust_par, columns)
      ),
      stringsAsFactors = FALSE
    )
  }))

  rbind(commonfactor_row, user_rows)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
inputs <- make_inputs(args$n_snp, args$k, args$seed)
user_matrix <- expand.grid(
  fix_measurement = args$fix_measurement,
  Q_SNP = args$q_snp,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

cat("equivalence_check\n")
print(check_equivalence(inputs, user_matrix), row.names = FALSE)

timings <- do.call(rbind, lapply(seq_len(args$repeats), function(repeat_id) {
  commonfactor_rows <- do.call(rbind, lapply(c("lavaan", "lavaanrust"), function(backend) {
    rows <- list(
      benchmark_once("commonfactorGWAS", backend, FALSE, 1L, inputs)
    )
    rows <- c(rows, lapply(args$cores, function(core_count) {
      benchmark_once("commonfactorGWAS", backend, TRUE, core_count, inputs)
    }))

    out <- do.call(rbind, rows)
    out$repeat_id <- repeat_id
    out
  }))

  user_rows <- do.call(rbind, lapply(seq_len(nrow(user_matrix)), function(row_idx) {
    fix_measurement <- user_matrix$fix_measurement[[row_idx]]
    q_snp <- user_matrix$Q_SNP[[row_idx]]

    do.call(rbind, lapply(c("lavaan", "lavaanrust"), function(backend) {
      rows <- list(
        benchmark_once("userGWAS", backend, FALSE, 1L, inputs, fix_measurement, q_snp)
      )
      rows <- c(rows, lapply(args$cores, function(core_count) {
        benchmark_once("userGWAS", backend, TRUE, core_count, inputs, fix_measurement, q_snp)
      }))

      out <- do.call(rbind, rows)
      out$repeat_id <- repeat_id
      out
    }))
  }))

  rbind(commonfactor_rows, user_rows)
}))

timings$n_snp <- args$n_snp
timings$k <- args$k
timings$seed <- args$seed

cat("timings\n")
print(timings, row.names = FALSE)
