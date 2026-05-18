test_that("userGWAS first-stage DWLS slice matches frozen lavaan estimates", {
  fixture <- user_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )

  expect_s4_class(fit, "lavaan_rust_fit")
  expect_equal(lavaanrust::parTable_rust(fit)$est, fixture$first_stage_est, tolerance = 2e-6)
  expect_equal(dim(lavaanrust::lavInspect_rust(fit, "delta")), c(10L, 8L))
  expect_equal(fit@Model$model_kind, "user_gwas_dwls")
})

test_that("userGWAS fixed-measurement refit matches frozen lavaan estimates", {
  fixture <- user_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$fixed_model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  par_table <- lavaanrust::parTable_rust(fit)

  expect_equal(par_table$est, fixture$fixed_est, tolerance = 2e-6)
  expect_equal(par_table$free, c(0L, 0L, 0L, seq_len(6L)))
  expect_equal(dim(lavaanrust::lavInspect_rust(fit, "delta")), c(10L, 6L))
  expect_equal(
    colnames(lavaanrust::lavInspect_rust(fit, "delta")),
    c("A~~A", "B~~B", "C~~C", "F1~~F1", "F1~SNP", "SNP~~SNP")
  )
})

test_that("userGWAS unrestricted base-model reuse stays native", {
  fixture <- user_gwas_fixture()
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

test_that("userGWAS fixed-measurement base-model reuse stays native", {
  fixture <- user_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$fixed_model,
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
