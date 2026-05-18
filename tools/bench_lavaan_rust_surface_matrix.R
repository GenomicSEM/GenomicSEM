#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
})

parse_args <- function(args) {
  values <- list(
    n_snp = 100L,
    repeats = 1L,
    seed = 20260506L
  )

  for (arg in args) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) {
      stop("Arguments must have the form key=value.", call. = FALSE)
    }

    key <- parts[[1L]]
    value <- parts[[2L]]
    if (key == "n_snp") {
      values$n_snp <- as.integer(value)
    } else if (key == "repeats") {
      values$repeats <- as.integer(value)
    } else if (key == "seed") {
      values$seed <- as.integer(value)
    } else {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
  }

  values
}

make_snps <- function(traits, n_snp, seed) {
  set.seed(seed)
  snps <- data.frame(
    SNP = paste0("rs", seq_len(n_snp)),
    CHR = rep(1L, n_snp),
    BP = seq_len(n_snp),
    MAF = runif(n_snp, 0.05, 0.5),
    A1 = rep("A", n_snp),
    A2 = rep("G", n_snp),
    check.names = FALSE
  )
  beta <- as.data.frame(outer(
    seq_len(n_snp),
    seq_along(traits),
    function(snp_idx, trait_idx) {
      c(0.02, -0.01, 0.015, -0.005)[((snp_idx + trait_idx - 2L) %% 4L) + 1L]
    }
  ))
  names(beta) <- paste0("beta.", traits)
  se <- as.data.frame(outer(
    seq_len(n_snp),
    seq_along(traits),
    function(snp_idx, trait_idx) {
      0.045 + 0.003 * ((snp_idx - 1L) %% 4L) + 0.002 * (trait_idx - 1L)
    }
  ))
  names(se) <- paste0("se.", traits)
  cbind(snps, beta, se)
}

make_one_factor_inputs <- function(n_snp, seed) {
  traits <- c("A", "B", "C")
  s_ld <- matrix(
    c(
      1.00, 0.45, 0.35,
      0.45, 1.10, 0.40,
      0.35, 0.40, 0.95
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(traits, traits)
  )
  i_ld <- matrix(0.02, nrow = 3L, ncol = 3L, dimnames = list(traits, traits))
  diag(i_ld) <- 1.05

  list(
    covstruc = list(V = diag(6L), S = s_ld, I = i_ld),
    SNPs = make_snps(traits, n_snp, seed),
    simple_model = paste("F1 =~ A + B + C", "F1 ~ SNP", sep = "\n"),
    flexible_model = paste(
      "F1 =~ NA*A + l2*B + l3*C",
      "F1 ~~ 1*F1",
      "F1 ~ gamma*SNP",
      "A ~ direct*SNP",
      "combo := gamma * direct",
      sep = "\n"
    ),
    usermodel = "F1 =~ A + B + C"
  )
}

make_two_factor_inputs <- function(n_snp, seed) {
  traits <- c("A", "B", "C", "D", "E", "F")
  s_ld <- matrix(0.18, nrow = 6L, ncol = 6L, dimnames = list(traits, traits))
  s_ld[1:3, 1:3] <- 0.45
  s_ld[4:6, 4:6] <- 0.40
  diag(s_ld) <- 1
  i_ld <- matrix(0.02, nrow = 6L, ncol = 6L, dimnames = list(traits, traits))
  diag(i_ld) <- 1.05

  list(
    covstruc = list(V = diag(21L), S = s_ld, I = i_ld),
    SNPs = make_snps(traits, n_snp, seed + 1L),
    gwas_model = paste(
      "F1 =~ A + B + C",
      "F2 =~ D + E + F",
      "F1 ~~ F2",
      "F1 ~ SNP",
      "F2 ~ SNP",
      sep = "\n"
    ),
    usermodel = paste(
      "F1 =~ A + B + C",
      "F2 =~ D + E + F",
      "F1 ~~ F2",
      sep = "\n"
    )
  )
}

bind_usergwas <- function(result) {
  do.call(rbind, result)
}

numeric_max_abs_diff <- function(left, right) {
  columns <- intersect(
    names(left)[vapply(left, is.numeric, logical(1L))],
    names(right)[vapply(right, is.numeric, logical(1L))]
  )
  max(vapply(
    columns,
    function(column) max(abs(as.numeric(left[[column]]) - as.numeric(right[[column]])), na.rm = TRUE),
    numeric(1L)
  ))
}

bind_rows_fill <- function(...) {
  frames <- list(...)
  all_names <- unique(unlist(lapply(frames, names), use.names = FALSE))
  frames <- lapply(frames, function(frame) {
    missing <- setdiff(all_names, names(frame))
    for (column in missing) {
      frame[[column]] <- NA
    }
    frame[all_names]
  })
  do.call(rbind, frames)
}

run_usermodel <- function(inputs, model, backend, std_lv) {
  fun <- if (backend == "lavaan") usermodel else usermodel_rust
  result <- NULL
  utils::capture.output({
    result <- suppressWarnings(fun(
      covstruc = inputs$covstruc,
      model = model,
      estimation = "DWLS",
      std.lv = std_lv,
      imp_cov = FALSE,
      fix_resid = TRUE
    ))
  })
  result$results
}

run_usergwas <- function(inputs, model, backend, std_lv, fix_measurement, q_snp) {
  fun <- if (backend == "lavaan") userGWAS else userGWAS_rust
  result <- NULL
  utils::capture.output({
    result <- suppressWarnings(fun(
      covstruc = inputs$covstruc,
      SNPs = inputs$SNPs,
      model = model,
      estimation = "DWLS",
      parallel = FALSE,
      GC = "standard",
      std.lv = std_lv,
      fix_measurement = fix_measurement,
      Q_SNP = q_snp,
      printwarn = TRUE
    ))
  })
  bind_usergwas(result)
}

benchmark_case <- function(case, backend, runner) {
  gc()
  elapsed <- system.time(runner(backend))[["elapsed"]]
  transform(case, backend = backend, elapsed_sec = elapsed)
}

select_gwas_model <- function(inputs, case_name) {
  if (case_name == "flexible_one_factor") {
    return(inputs$flexible_model)
  }
  if (!is.null(inputs$gwas_model)) {
    return(inputs$gwas_model)
  }

  inputs$simple_model
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
one_factor <- make_one_factor_inputs(args$n_snp, args$seed)
two_factor <- make_two_factor_inputs(args$n_snp, args$seed)

usermodel_cases <- rbind(
  data.frame(case = "one_factor", model_family = "one_factor", std_lv = c(FALSE, TRUE)),
  data.frame(case = "two_factor", model_family = "two_factor", std_lv = c(FALSE, TRUE))
)
usergwas_cases <- rbind(
  transform(expand.grid(
    std_lv = c(FALSE, TRUE),
    fix_measurement = c(FALSE, TRUE),
    q_snp = c(FALSE, TRUE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), case = "simple_one_factor", model_family = "one_factor"),
  transform(expand.grid(
    std_lv = TRUE,
    fix_measurement = c(FALSE, TRUE),
    q_snp = FALSE,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), case = "flexible_one_factor", model_family = "one_factor"),
  transform(expand.grid(
    std_lv = c(FALSE, TRUE),
    fix_measurement = c(FALSE, TRUE),
    q_snp = c(FALSE, TRUE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), case = "two_factor", model_family = "two_factor")
)

cat("equivalence_check\n")
usermodel_equivalence <- do.call(rbind, lapply(seq_len(nrow(usermodel_cases)), function(row_idx) {
  case <- usermodel_cases[row_idx, , drop = FALSE]
  inputs <- if (case$model_family == "one_factor") one_factor else two_factor
  old <- run_usermodel(inputs, inputs$usermodel, "lavaan", case$std_lv)
  rust <- run_usermodel(inputs, inputs$usermodel, "lavaanrust", case$std_lv)
  transform(case, workflow = "usermodel", max_abs_diff = numeric_max_abs_diff(old, rust))
}))
usergwas_equivalence <- do.call(rbind, lapply(seq_len(nrow(usergwas_cases)), function(row_idx) {
  case <- usergwas_cases[row_idx, , drop = FALSE]
  inputs <- if (case$model_family == "one_factor") one_factor else two_factor
  model <- select_gwas_model(inputs, case$case)
  old <- run_usergwas(inputs, model, "lavaan", case$std_lv, case$fix_measurement, case$q_snp)
  rust <- run_usergwas(inputs, model, "lavaanrust", case$std_lv, case$fix_measurement, case$q_snp)
  transform(case, workflow = "userGWAS", max_abs_diff = numeric_max_abs_diff(old, rust))
}))
print(bind_rows_fill(usermodel_equivalence, usergwas_equivalence), row.names = FALSE)

timings <- do.call(rbind, lapply(seq_len(args$repeats), function(repeat_id) {
  usermodel_rows <- do.call(rbind, lapply(seq_len(nrow(usermodel_cases)), function(row_idx) {
    case <- usermodel_cases[row_idx, , drop = FALSE]
    inputs <- if (case$model_family == "one_factor") one_factor else two_factor
    rows <- lapply(c("lavaan", "lavaanrust"), function(backend) {
      benchmark_case(
        transform(case, workflow = "usermodel"),
        backend,
        function(chosen_backend) run_usermodel(inputs, inputs$usermodel, chosen_backend, case$std_lv)
      )
    })
    do.call(rbind, rows)
  }))
  usergwas_rows <- do.call(rbind, lapply(seq_len(nrow(usergwas_cases)), function(row_idx) {
    case <- usergwas_cases[row_idx, , drop = FALSE]
    inputs <- if (case$model_family == "one_factor") one_factor else two_factor
    model <- select_gwas_model(inputs, case$case)
    rows <- lapply(c("lavaan", "lavaanrust"), function(backend) {
      benchmark_case(
        transform(case, workflow = "userGWAS"),
        backend,
        function(chosen_backend) run_usergwas(
          inputs,
          model,
          chosen_backend,
          case$std_lv,
          case$fix_measurement,
          case$q_snp
        )
      )
    })
    do.call(rbind, rows)
  }))
  transform(bind_rows_fill(usermodel_rows, usergwas_rows), repeat_id = repeat_id)
}))

timings$n_snp <- args$n_snp
timings$seed <- args$seed
cat("timings\n")
print(timings, row.names = FALSE)
