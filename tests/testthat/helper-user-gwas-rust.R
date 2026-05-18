user_gwas_wrapper_fixture <- function() {
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
  i_ld <- matrix(
    c(
      1.05, 0.02, 0.02,
      0.02, 1.05, 0.02,
      0.02, 0.02, 1.05
    ),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(traits, traits)
  )
  snps <- data.frame(
    SNP = c("rs1", "rs2"),
    CHR = c(1L, 1L),
    BP = c(1L, 2L),
    MAF = c(0.30, 0.20),
    A1 = c("A", "A"),
    A2 = c("G", "G"),
    beta.A = c(0.02, -0.01),
    se.A = c(0.05, 0.06),
    beta.B = c(0.03, 0.01),
    se.B = c(0.05, 0.06),
    beta.C = c(0.01, -0.02),
    se.C = c(0.05, 0.06),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  list(
    covstruc = list(
      V = diag(6L),
      S = s_ld,
      I = i_ld
    ),
    SNPs = snps,
    simple_model = paste(
      "F1 =~ A + B + C",
      "F1 ~ SNP",
      sep = "\n"
    ),
    bounded_model = paste(
      "F1 =~ A + l2*B + l3*C",
      "F1 ~ gamma*SNP",
      "l2 > .1",
      sep = "\n"
    ),
    flexible_model = paste(
      "F1 =~ NA*A + l2*B + l3*C",
      "F1 ~~ 1*F1",
      "F1 ~ gamma*SNP",
      "A ~ direct*SNP",
      "combo := gamma * direct",
      sep = "\n"
    ),
    bounded_usermodel = paste(
      "F1 =~ A + B + C",
      "A ~~ rvA*A",
      "rvA > .001",
      sep = "\n"
    )
  )
}

two_factor_wrapper_fixture <- function() {
  traits <- c("A", "B", "C", "D", "E", "F")
  s_ld <- matrix(0.18, nrow = 6L, ncol = 6L, dimnames = list(traits, traits))
  s_ld[1:3, 1:3] <- 0.45
  s_ld[4:6, 4:6] <- 0.40
  diag(s_ld) <- 1
  i_ld <- matrix(0.02, nrow = 6L, ncol = 6L, dimnames = list(traits, traits))
  diag(i_ld) <- 1.05
  snps <- data.frame(
    SNP = c("rs1", "rs2"),
    CHR = c(1L, 1L),
    BP = c(1L, 2L),
    MAF = c(0.30, 0.20),
    A1 = c("A", "A"),
    A2 = c("G", "G"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  for (trait in traits) {
    snps[[paste0("beta.", trait)]] <- c(0.02, -0.01)
    snps[[paste0("se.", trait)]] <- c(0.05, 0.06)
  }

  list(
    covstruc = list(
      V = diag(21L),
      S = s_ld,
      I = i_ld
    ),
    SNPs = snps,
    model = paste(
      "F1 =~ A + B + C",
      "F2 =~ D + E + F",
      "F1 ~~ F2",
      sep = "\n"
    ),
    gwas_model = paste(
      "F1 =~ A + B + C",
      "F2 =~ D + E + F",
      "F1 ~~ F2",
      "F1 ~ SNP",
      "F2 ~ SNP",
      sep = "\n"
    )
  )
}

run_user_gwas_wrapper <- function(fun, fixture, model, std_lv, fix_measurement, q_snp) {
  result <- NULL
  utils::capture.output({
    result <- suppressWarnings(fun(
      covstruc = fixture$covstruc,
      SNPs = fixture$SNPs,
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

  do.call(rbind, result)
}

run_user_model_wrapper <- function(fun, fixture, model, std_lv) {
  result <- NULL
  utils::capture.output({
    result <- suppressWarnings(fun(
      covstruc = fixture$covstruc,
      model = model,
      estimation = "DWLS",
      std.lv = std_lv,
      imp_cov = FALSE,
      fix_resid = TRUE
    ))
  })

  result$results
}

user_gwas_row_keys <- function(result) {
  paste(result$lhs, result$op, result$rhs, sep = "|")
}

user_gwas_numeric_columns <- function(result) {
  names(result)[vapply(result, is.numeric, logical(1L))]
}

user_gwas_numeric_frame <- function(result, columns = user_gwas_numeric_columns(result)) {
  out <- result[columns]
  rownames(out) <- NULL
  out
}

user_gwas_structural_rows <- function(result) {
  result[!result$op %in% c("==", "<", ">", "<=", ">="), , drop = FALSE]
}
