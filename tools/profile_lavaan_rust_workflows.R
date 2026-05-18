#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicSEM)
})

source("tests/testthat/helper-user-gwas-rust.R")

parse_args <- function(args) {
  values <- list(
    n_snp = 200L,
    interval = 0.001,
    output_dir = tempfile("lavaan-rust-profile-"),
    case = "usergwas_flexible_one_factor"
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
    } else if (key == "interval") {
      values$interval <- as.numeric(value)
    } else if (key == "output_dir") {
      values$output_dir <- value
    } else if (key == "case") {
      values$case <- value
    } else {
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    }
  }

  values
}

make_snps <- function(snps, n_snp) {
  snps[rep(seq_len(nrow(snps)), length.out = n_snp), , drop = FALSE]
}

profile_case <- function(name, output_dir, interval, runner) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  profile_path <- file.path(output_dir, paste0(name, ".out"))
  invisible(runner())
  Rprof(profile_path, interval = interval)
  on.exit(Rprof(NULL), add = TRUE)
  invisible(runner())
  Rprof(NULL)

  summary <- summaryRprof(profile_path)
  self <- utils::head(summary$by.self, 12L)
  total <- utils::head(summary$by.total, 12L)
  list(
    case = name,
    profile_path = profile_path,
    self = self,
    total = total
  )
}

print_profile <- function(profile) {
  cat("\ncase:", profile$case, "\n")
  cat("by_self\n")
  print(profile$self)
  cat("by_total\n")
  print(profile$total)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
one_factor <- user_gwas_wrapper_fixture()
one_factor$SNPs <- make_snps(one_factor$SNPs, args$n_snp)
two_factor <- two_factor_wrapper_fixture()
two_factor$SNPs <- make_snps(two_factor$SNPs, args$n_snp)

cases <- list(
  usermodel_bounded = function() run_user_model_wrapper(
    usermodel_rust,
    one_factor,
    one_factor$bounded_usermodel,
    std_lv = FALSE
  ),
  usergwas_flexible_one_factor = function() run_user_gwas_wrapper(
    userGWAS_rust,
    one_factor,
    one_factor$flexible_model,
    std_lv = TRUE,
    fix_measurement = FALSE,
    q_snp = FALSE
  ),
  usergwas_two_factor = function() run_user_gwas_wrapper(
    userGWAS_rust,
    two_factor,
    two_factor$gwas_model,
    std_lv = FALSE,
    fix_measurement = FALSE,
    q_snp = FALSE
  )
)
if (!args$case %in% names(cases)) {
  stop(
    sprintf("Unknown case '%s'. Available cases: %s", args$case, paste(names(cases), collapse = ", ")),
    call. = FALSE
  )
}
profile <- profile_case(args$case, args$output_dir, args$interval, cases[[args$case]])

cat("profile_dir\n")
cat(args$output_dir, "\n")
print_profile(profile)
