test_that("one-factor DWLS slice matches frozen lavaan fixtures", {
  fixture <- one_factor_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  par_table <- lavaanrust::parTable_rust(fit)

  expect_s4_class(fit, "lavaan_rust_fit")
  expect_equal(par_table$est[1:3], fixture$loadings, tolerance = 2e-7)
  expect_equal(par_table$est[5:7], fixture$residuals, tolerance = 2e-7)
  expect_equal(lavaanrust::fitted_rust(fit)$cov, fixture$implied, tolerance = 2e-7)
  expect_equal(lavaanrust::lavInspect_rust(fit, "delta"), fixture$delta, tolerance = 2e-7)
  expect_true(lavaanrust::lavInspect_rust(fit, "converged"))
})

test_that("unsupported syntax errors instead of falling back to lavaan", {
  fixture <- one_factor_fixture()
  model <- paste(
    "F1 =~ A + B + C",
    "F1 == F1",
    sep = "\n"
  )
  expect_error(
    lavaanrust::sem_rust(
      model,
      sample.cov = fixture$sample_cov,
      estimator = "DWLS",
      WLS.V = fixture$wls_v,
      se = "standard",
      sample.nobs = 2
    ),
    "does not fall back to lavaan"
  )
})

test_that("commonfactor null-model slice fits diagonal covariance structure", {
  fixture <- one_factor_fixture()
  model <- paste(
    "A ~~ A",
    "B ~~ B",
    "C ~~ C",
    "VF1 =~ 1*A",
    "VF2 =~ 1*B",
    "VF3 =~ 1*C",
    "VF1 ~~ 0*VF2 + 0*VF3",
    "VF2 ~~ 0*VF3",
    "A ~~ 0*B + 0*C",
    "B ~~ 0*C",
    "VF1 ~~ 0*VF1",
    "VF2 ~~ 0*VF2",
    "VF3 ~~ 0*VF3",
    sep = "\n"
  )

  fit <- lavaanrust::sem_rust(
    model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  expected <- diag(diag(fixture$sample_cov))
  dimnames(expected) <- dimnames(fixture$sample_cov)

  expect_s4_class(fit, "lavaan_rust_fit")
  expect_equal(
    lavaanrust::fitted_rust(fit)$cov,
    expected,
    tolerance = 1e-12
  )
  expect_equal(dim(lavaanrust::lavInspect_rust(fit, "delta")), c(6L, 3L))
})

test_that("parameter-table refit recovers saturated observed covariance model", {
  fixture <- one_factor_fixture()
  model <- paste(
    "A ~~ A",
    "B ~~ B",
    "C ~~ C",
    "VF1 =~ 1*A",
    "VF2 =~ 1*B",
    "VF3 =~ 1*C",
    "VF1 ~~ 0*VF2 + 0*VF3",
    "VF2 ~~ 0*VF3",
    "A ~~ 0*B + 0*C",
    "B ~~ 0*C",
    "VF1 ~~ 0*VF1",
    "VF2 ~~ 0*VF2",
    "VF3 ~~ 0*VF3",
    sep = "\n"
  )
  null_fit <- lavaanrust::sem_rust(
    model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  q_table <- lavaanrust::parTable_rust(null_fit)
  z <- 6L
  p2 <- length(q_table$free) - z
  q_table$free <- c(rep(0L, p2), seq_len(z))
  q_table$ustart <- q_table$est

  refit <- lavaanrust::sem_rust(
    q_table,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )

  expect_equal(lavaanrust::fitted_rust(refit)$cov, fixture$sample_cov, tolerance = 1e-12)
  expect_equal(dim(lavaanrust::lavInspect_rust(refit, "delta")), c(6L, 6L))
})
