test_that("commonfactor GWAS DWLS slice matches frozen lavaan estimates", {
  fixture <- commonfactor_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  par_table <- lavaanrust::parTable_rust(fit)

  expect_s4_class(fit, "lavaan_rust_fit")
  expect_equal(par_table$est, fixture$est, tolerance = 2e-6)
  expect_equal(dim(lavaanrust::lavInspect_rust(fit, "delta")), c(10L, 8L))
  expect_equal(fit@Model$model_kind, "commonfactor_gwas_dwls")
})

test_that("commonfactor GWAS base-model reuse stays native", {
  fixture <- commonfactor_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  refit <- lavaanrust::lavaan_rust(
    sample.cov = fixture$sample_cov,
    WLS.V = fixture$wls_v,
    slotOptions = fit@Options,
    slotParTable = fit@ParTable,
    slotData = fit@Data,
    slotModel = fit@Model
  )

  expect_equal(
    lavaanrust::parTable_rust(refit)$est,
    lavaanrust::parTable_rust(fit)$est,
    tolerance = 2e-6
  )
})

test_that("commonfactor GWAS Q-model refit matches frozen lavaan estimates", {
  fixture <- commonfactor_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  q_table <- lavaanrust::parTable_rust(fit)
  k <- 3L
  q_table$free <- c(rep(0L, k + 1L), seq_len(k * 2L), 0L, 0L)
  q_table$ustart <- q_table$est
  snp_resid <- lavaanrust::resid_rust(fit)$cov[k + 1L, seq_len(k)]

  for (row_idx in seq_len(nrow(q_table))) {
    if (q_table$free[[row_idx]] > 0L && q_table$free[[row_idx]] <= k) {
      q_table$ustart[[row_idx]] <- snp_resid[[q_table$free[[row_idx]]]]
    }
  }

  q_fit <- lavaanrust::sem_rust(
    q_table,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )

  expect_equal(lavaanrust::parTable_rust(q_fit)$est, fixture$q_est, tolerance = 2e-6)
  expect_equal(dim(lavaanrust::lavInspect_rust(q_fit, "delta")), c(10L, 6L))
})

test_that("unsupported model reuse still errors instead of falling back", {
  fixture <- commonfactor_gwas_fixture()
  expect_error(
    lavaanrust::lavaan_rust(
      sample.cov = fixture$sample_cov,
      WLS.V = fixture$wls_v,
      slotModel = list(model_kind = "unknown")
    ),
    "does not fall back to lavaan"
  )
})
