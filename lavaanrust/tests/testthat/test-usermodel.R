test_that("marker-scaled one-factor usermodel slice matches frozen lavaan estimates", {
  fixture <- usermodel_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v,
    std.lv = FALSE
  )

  expect_s4_class(fit, "lavaan_rust_fit")
  expect_equal(lavaanrust::parTable_rust(fit)$est, fixture$marker_est, tolerance = 2e-6)
  expect_equal(dim(lavaanrust::lavInspect_rust(fit, "delta")), c(6L, 6L))
  expect_equal(unname(lavaanrust::lavInspect_rust(fit, "cor.lv")), matrix(fixture$marker_est[[7L]]), tolerance = 2e-6)
})

test_that("simple one-factor usermodel slice also supports std.lv scaling", {
  fixture <- usermodel_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v,
    std.lv = TRUE
  )

  expect_equal(lavaanrust::parTable_rust(fit)$est, fixture$std_lv_est, tolerance = 2e-6)
  expect_equal(dim(lavaanrust::lavInspect_rust(fit, "delta")), c(6L, 6L))
})

test_that("standardizedSolution_rust exposes the one-factor standardized rows usermodel needs", {
  fixture <- usermodel_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v,
    std.lv = FALSE
  )
  standardized <- lavaanrust::standardizedSolution_rust(fit)

  expect_equal(standardized$est.std, fixture$standardized_est, tolerance = 2e-6)
  expect_equal(is.na(standardized$pvalue), c(rep(FALSE, 6L), TRUE))
})

test_that("generic RAM standardized output re-evaluates defined parameters", {
  fixture <- usermodel_fixture()
  model <- paste(
    "F1 =~ NA*A + l2*B + l3*C",
    "F1 ~~ 1*F1",
    "omega := l2 * l3",
    sep = "\n"
  )
  fit <- lavaanrust::sem_rust(
    model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v,
    std.lv = FALSE
  )
  standardized <- lavaanrust::standardizedSolution_rust(fit)

  expect_equal(
    standardized$est.std,
    unname(c(
      fixture$std_lv_est[1:3] / sqrt(diag(fixture$sample_cov)),
      1,
      fixture$standardized_est[4:6],
      standardized$est.std[[2L]] * standardized$est.std[[3L]]
    )),
    tolerance = 2e-6
  )
  expect_equal(is.na(standardized$pvalue), c(rep(FALSE, 3L), TRUE, rep(FALSE, 4L)))
})
