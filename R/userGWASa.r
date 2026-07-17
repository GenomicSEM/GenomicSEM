#' userGWASa: Ultra-fast multivariate GWAS with flexible analytic estimation
#'
#' Runs a multivariate GWAS across a set of
#' GWAS summary statistics and a user-specified factor model. Factor-specific
#' SNP effects (betas, SEs, Z-statistics, p-values) and an omnibus
#' heterogeneity statistic (Q_omnibus) are computed.
#'
#' @param sumstats A \code{data.frame} of merged GWAS summary statistics,
#'   as produced by the \code{sumstats()} function in GenomicSEM. Must contain
#'   columns \code{SNP}, \code{A1}, \code{A2}, \code{MAF}, \code{N}, and
#'   trait-specific \code{beta.*} and \code{se.*} columns.
#' @param LDSCoutput A list object returned by the \code{ldsc()} function.
#' @param model A character string specifying the factor model in
#'   \code{lavaan}-style syntax. Ignored if
#'   \code{usermod} is provided.
#' @param usermod Optional. A pre-fitted no-SNP model results data frame
#'   (the \code{$results} element from a \code{usermodel()} call). When
#'   supplied, the function skips fitting the no-SNP model internally and uses
#'   these parameter estimates directly to extract lambda coefficients.
#'   Default is \code{NULL}.
#' @param batch_size Integer. Number of SNPs to process per batch. Larger
#'   values increase memory use but reduce overhead. Default is \code{100000}.
#'
#' @return A \code{data.frame} with one row per SNP and the following columns:
#'   \itemize{
#'     \item The first 6 columns from \code{sumstats} (SNP identifiers and
#'           allele information).
#'     \item \code{beta_<factor>}: GLS-estimated SNP effect on each factor.
#'     \item \code{SE_<factor>}: Sandwich-corrected standard error of the
#'           factor beta.
#'     \item \code{Z_beta_<factor>}: Z-statistic for the factor beta.
#'     \item \code{p_val_<factor>}: Two-sided p-value for the factor beta.
#'     \item \code{Q_omnibus}, \code{Q_omnibus_df}, \code{Q_omnibus_pval}:
#'           Omnibus Q_SNP statistic across all traits, its degrees of freedom,
#'           and p-value.
#'   }
#'
#' @details
#' The function implements a two-stage approach. First, a no-SNP factor model
#' is fit using \code{\link[GenomicSEM]{usermodel}} with DWLS estimation to
#' obtain factor loading estimates (lambdas). Second, for each batch of SNPs,
#' SNP-to-factor betas are estimated via GLS using the diagonal of the
#' SNP-specific sampling covariance matrix as weights, with a sandwich
#' variance estimator for the standard errors.
#'
#' The Q_omnibus statistic tests whether the observed SNP-trait association
#' vector is consistent with the implied factor model.
#'
#' @seealso \code{\link[GenomicSEM]{usermodel}}, \code{\link[GenomicSEM]{userGWAS}}
#' @noRd
#'
#' @examples
#' \dontrun{
#' load("LDSC_PSYCH.RData")
#' sumstats <- data.table::fread("Psych_sumstats_4GLS.txt", data.table = FALSE)
#'
#' model <- '
#'   Psych =~ a*SCZ + a*BIP
#'   Neuro =~ ADHD + MDD + ASD
#'   Psych ~~ Neuro
#'   Psych ~ SNP
#'   Neuro ~ SNP
#' '
#'
#' results <- userGWASa(
#'   sumstats   = sumstats,
#'   LDSCoutput = LDSC_P,
#'   model      = model,
#'   batch_size = 50000
#' )
#' }
#'
#' @keywords internal
.userGWASa <- function(sumstats, LDSCoutput, model, usermod = NULL, batch_size = 100000) {

  # Helper function: V' M V
  VMV <- function(V1, M, V2) { V1 %*% M %*% V2 }

  start_time <- Sys.time()
  cat("userGWASa started at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

  # Coerce to plain data.frame to ensure consistent column subsetting
  # regardless of whether input is a data.table, tibble, or other tabular class
  sumstats <- as.data.frame(sumstats)

  # ── No-SNP model ──────────────────────────────────────────────────────────────
  if (is.character(model) & is.null(usermod)) {
    model_lines  <- strsplit(model, "\n")[[1]]
    snp_lines    <- grep("~.*\\bSNP\\b", model_lines, value = TRUE)
    nosnp_model  <- paste(grep("~.*\\bSNP\\b", model_lines, value = TRUE, invert = TRUE),
                          collapse = "\n")

    captured_output <- capture.output(
    suppressWarnings(suppressMessages({
        nosnpmod <- usermodel(
        LDSCoutput, estimation = "DWLS", model = nosnp_model,
        CFIcalc = FALSE, std.lv = FALSE, imp_cov = FALSE
        )
    }))
    )
    nosnpmod <- nosnpmod$results

    # Re-emit smoothing warning if it occurred
    if (any(grepl("smoothed", captured_output))) {
    warning("The S matrix was smoothed prior to model estimation. ")
    }

  } else {
    nosnpmod <- usermod
  }

  # ── Extract lambda coefficients ───────────────────────────────────────────────
  factors    <- unique(nosnpmod$lhs[nosnpmod$op == "=~"])
  traits     <- colnames(LDSCoutput$S)
  num_traits  <- ncol(LDSCoutput$S)
  num_factors <- length(factors)

  combinations <- expand.grid(traits = traits, factors = factors)
  column_names <- paste0("lambda.", combinations$traits, "_", combinations$factors)

  extract_lambdas <- function(df, factors, traits, num_traits, num_factors) {
    lambdas <- rep(0, num_traits * num_factors)
    for (factor_idx in seq_along(factors)) {
      factor <- factors[factor_idx]
      for (trait_idx in seq_along(traits)) {
        trait <- traits[trait_idx]
        row   <- df[df$lhs == factor & df$rhs == trait & df$op == "=~", ]
        lambda_value <- if (nrow(row) > 0) row$Unstand_Est else 0
        lambdas[(factor_idx - 1) * num_traits + trait_idx] <- lambda_value
      }
    }
    lambdas_df <- data.frame(matrix(lambdas, nrow = 1, byrow = TRUE))
    colnames(lambdas_df) <- column_names
    return(lambdas_df)
  }

  lambdas <- extract_lambdas(nosnpmod, factors, traits, num_traits, num_factors)

  # ── Initialise output data frame ──────────────────────────────────────────────
  GLS_mGWAS_results <- sumstats[, 1:6]

  for (j in factors) {
    GLS_mGWAS_results[[paste0("beta_",   j)]] <- NA_real_
    GLS_mGWAS_results[[paste0("SE_",     j)]] <- NA_real_
    GLS_mGWAS_results[[paste0("Z_beta_", j)]] <- NA_real_
    GLS_mGWAS_results[[paste0("p_val_",  j)]] <- NA_real_
  }

  GLS_mGWAS_results <- GLS_mGWAS_results %>%
    mutate(Q_omnibus = NA_real_, Q_omnibus_df = NA_real_, Q_omnibus_pval = NA_real_)

  # ── Batch loop ────────────────────────────────────────────────────────────────
  total_batches <- ceiling(nrow(sumstats) / batch_size)
  pb <- txtProgressBar(min = 0, max = total_batches, style = 3)

  for (batch_num in seq_len(total_batches)) {

    i             <- (batch_num - 1) * batch_size + 1
    batch_end     <- min(i + batch_size - 1, nrow(sumstats))
    snp_batch     <- sumstats[i:batch_end, ]
    batch_indices <- i:batch_end

    betas <- snp_batch %>% select(contains("beta."))
    SEs   <- snp_batch %>% select(contains("se."))

    # Lambda matrix
    lambdas_snp <- as.numeric(lambdas)
    R_SNP       <- LDSCoutput$I
    diag(R_SNP)[diag(R_SNP) < 1] <- 1
    X      <- matrix(lambdas_snp, nrow = num_traits, ncol = num_factors)
    colnames(X) <- factors

    # SNP-trait beta and SE lists
    beta_l   <- lapply(transpose(betas), function(x) as.numeric(unlist(x)))
    se_snp   <- lapply(transpose(SEs),   function(x) as.numeric(unlist(x)))

    # V_SNP list and its diagonalised inverse
    V_SNP_list <- apply(SEs, 1, function(se) {
    lavaan::lav_cor2cov(R = as.matrix(R_SNP), sds = as.numeric(se))
    }, simplify = FALSE)

    V_d_list_inv <- lapply(V_SNP_list, function(V_SNP) diag(1 / diag(V_SNP)))

    # GLS factor betas (sandwich SE)
    Beta_list <- unname(Map(function(V_d_inv, beta) {
      solve(t(X) %*% V_d_inv %*% X) %*% t(X) %*% V_d_inv %*% beta
    }, V_d_list_inv, beta_l))

    SE_parallel_list <- unname(Map(function(V_d_inv, V_SNP) {
      bread    <- solve(t(X) %*% V_d_inv %*% X)
      meat     <- t(X) %*% V_d_inv %*% V_SNP %*% V_d_inv %*% X
      sandwich <- bread %*% meat %*% bread
      sqrt(diag(sandwich))
    }, V_d_list_inv, V_SNP_list))

    # Convert to data frames
    SE_parallel_df <- as.data.frame(do.call(rbind, SE_parallel_list))
    colnames(SE_parallel_df) <- paste0("SE_", factors)

    Beta_parallel_df <- as.data.frame(do.call(rbind, lapply(Beta_list, as.numeric)))
    colnames(Beta_parallel_df) <- paste0("Beta_", factors)

    Z_df_parallel <- Beta_parallel_df / SE_parallel_df
    colnames(Z_df_parallel) <- paste0("Z_", factors)

    # ── Q_omnibus ───────────────────────────────────────────────────────────────
    solveI <- solve(LDSCoutput$I)

    beta_hats_list <- lapply(Beta_list, function(beta) {
      as.matrix(beta)[, 1] %*% t(X)
    })

    Resid_parallel_list <- mapply(function(beta, beta_hat) {
      unname(as.numeric(as.vector(beta) - beta_hat))
    }, beta_l, beta_hats_list, SIMPLIFY = FALSE)

    inside_list <- lapply(se_snp, function(SEs_snp) {
      SEs_snp <- as.numeric(SEs_snp)
      solveI / (SEs_snp %*% t(SEs_snp))
    })

    Q_Omnibus_parallel <- mapply(Resid_parallel_list, inside_list, Resid_parallel_list,
                                 FUN = VMV)

    # ── Write batch results ──────────────────────────────────────────────────────
    for (j in factors) {
      GLS_mGWAS_results[batch_indices, paste0("beta_",   j)] <- Beta_parallel_df[, paste0("Beta_", j)]
      GLS_mGWAS_results[batch_indices, paste0("SE_",     j)] <- SE_parallel_df[,   paste0("SE_",   j)]
      GLS_mGWAS_results[batch_indices, paste0("Z_beta_", j)] <- Z_df_parallel[,    paste0("Z_",    j)]
      GLS_mGWAS_results[batch_indices, paste0("p_val_",  j)] <- 2 * pnorm(-abs(
        GLS_mGWAS_results[batch_indices, paste0("Z_beta_", j)]))
    }

    GLS_mGWAS_results[batch_indices, "Q_omnibus"]    <- Q_Omnibus_parallel
    GLS_mGWAS_results[batch_indices, "Q_omnibus_df"] <- length(colnames(betas)) - ncol(Beta_parallel_df)
    GLS_mGWAS_results[batch_indices, "Q_omnibus_pval"] <- pchisq(
      GLS_mGWAS_results[batch_indices, "Q_omnibus"],
      df = GLS_mGWAS_results[batch_indices, "Q_omnibus_df"],
      lower.tail = FALSE
    )

    GLS_mGWAS_results <- GLS_mGWAS_results %>% select(where(~ !all(is.na(.))))

    setTxtProgressBar(pb, batch_num)
    rm(snp_batch)
    if (batch_num %% 10 == 0) gc()
  }

  close(pb)

  end_time     <- Sys.time()
  elapsed_time <- end_time - start_time
  cat("Finished at:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total time elapsed:", round(elapsed_time, 2), attr(elapsed_time, "units"), "\n")

  return(GLS_mGWAS_results)
}