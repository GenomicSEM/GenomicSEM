test_that("userGWAS_rust matches lavaan across the supported std.lv matrix", {
  fixture <- user_gwas_wrapper_fixture()
  matrix <- expand.grid(
    std_lv = c(FALSE, TRUE),
    fix_measurement = c(FALSE, TRUE),
    q_snp = c(FALSE, TRUE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  for (row_idx in seq_len(nrow(matrix))) {
    args <- matrix[row_idx, , drop = FALSE]
    old <- run_user_gwas_wrapper(
      userGWAS,
      fixture,
      fixture$simple_model,
      args$std_lv,
      args$fix_measurement,
      args$q_snp
    )
    rust <- run_user_gwas_wrapper(
      userGWAS_rust,
      fixture,
      fixture$simple_model,
      args$std_lv,
      args$fix_measurement,
      args$q_snp
    )

    expect_equal(
      user_gwas_row_keys(rust),
      user_gwas_row_keys(old),
      info = sprintf(
        "std.lv=%s fix_measurement=%s Q_SNP=%s",
        args$std_lv,
        args$fix_measurement,
        args$q_snp
      )
    )
    expect_equal(
      user_gwas_numeric_frame(rust, user_gwas_numeric_columns(old)),
      user_gwas_numeric_frame(old, user_gwas_numeric_columns(old)),
      tolerance = 2e-6,
      info = sprintf(
        "std.lv=%s fix_measurement=%s Q_SNP=%s",
        args$std_lv,
        args$fix_measurement,
        args$q_snp
      )
    )
  }
})

test_that("userGWAS_rust handles labeled direct effects and defined parameters", {
  fixture <- user_gwas_wrapper_fixture()

  for (fix_measurement in c(FALSE, TRUE)) {
    old <- run_user_gwas_wrapper(
      userGWAS,
      fixture,
      fixture$flexible_model,
      std_lv = TRUE,
      fix_measurement = fix_measurement,
      q_snp = FALSE
    )
    rust <- run_user_gwas_wrapper(
      userGWAS_rust,
      fixture,
      fixture$flexible_model,
      std_lv = TRUE,
      fix_measurement = fix_measurement,
      q_snp = FALSE
    )

    expect_true(any(rust$lhs == "combo" & rust$op == ":="))
    expect_equal(user_gwas_row_keys(rust), user_gwas_row_keys(old))
    expect_equal(
      user_gwas_numeric_frame(rust, user_gwas_numeric_columns(old)),
      user_gwas_numeric_frame(old, user_gwas_numeric_columns(old)),
      tolerance = 1e-5
    )
  }
})

test_that("userGWAS_rust matches lavaan on bounded generic models", {
  fixture <- user_gwas_wrapper_fixture()
  old <- run_user_gwas_wrapper(
    userGWAS,
    fixture,
    fixture$bounded_model,
    std_lv = FALSE,
    fix_measurement = FALSE,
    q_snp = FALSE
  )
  rust <- run_user_gwas_wrapper(
    userGWAS_rust,
    fixture,
    fixture$bounded_model,
    std_lv = FALSE,
    fix_measurement = FALSE,
    q_snp = FALSE
  )
  old <- user_gwas_structural_rows(old)
  rust <- user_gwas_structural_rows(rust)

  expect_equal(user_gwas_row_keys(rust), user_gwas_row_keys(old))
  expect_equal(
    user_gwas_numeric_frame(rust, user_gwas_numeric_columns(old)),
    user_gwas_numeric_frame(old, user_gwas_numeric_columns(old)),
    tolerance = 1e-5
  )
})

test_that("userGWAS_rust matches lavaan across two-factor option surfaces", {
  fixture <- two_factor_wrapper_fixture()
  matrix <- expand.grid(
    std_lv = c(FALSE, TRUE),
    fix_measurement = c(FALSE, TRUE),
    q_snp = c(FALSE, TRUE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  for (row_idx in seq_len(nrow(matrix))) {
    args <- matrix[row_idx, , drop = FALSE]
    old <- run_user_gwas_wrapper(
      userGWAS,
      fixture,
      fixture$gwas_model,
      args$std_lv,
      args$fix_measurement,
      args$q_snp
    )
    rust <- run_user_gwas_wrapper(
      userGWAS_rust,
      fixture,
      fixture$gwas_model,
      args$std_lv,
      args$fix_measurement,
      args$q_snp
    )

    expect_equal(
      user_gwas_row_keys(rust),
      user_gwas_row_keys(old),
      info = sprintf(
        "std.lv=%s fix_measurement=%s Q_SNP=%s",
        args$std_lv,
        args$fix_measurement,
        args$q_snp
      )
    )
    expect_equal(
      user_gwas_numeric_frame(rust, user_gwas_numeric_columns(old)),
      user_gwas_numeric_frame(old, user_gwas_numeric_columns(old)),
      tolerance = 2e-6,
      info = sprintf(
        "std.lv=%s fix_measurement=%s Q_SNP=%s",
        args$std_lv,
        args$fix_measurement,
        args$q_snp
      )
    )
  }
})

test_that("usermodel_rust matches lavaan on two-factor models", {
  fixture <- two_factor_wrapper_fixture()

  for (std_lv in c(FALSE, TRUE)) {
    old <- run_user_model_wrapper(usermodel, fixture, fixture$model, std_lv)
    rust <- run_user_model_wrapper(usermodel_rust, fixture, fixture$model, std_lv)

    expect_equal(user_gwas_row_keys(rust), user_gwas_row_keys(old))
    expect_equal(
      user_gwas_numeric_frame(rust, user_gwas_numeric_columns(old)),
      user_gwas_numeric_frame(old, user_gwas_numeric_columns(old)),
      tolerance = 2e-6
    )
  }
})

test_that("usermodel_rust matches lavaan on bounded generic models", {
  fixture <- user_gwas_wrapper_fixture()
  old <- run_user_model_wrapper(
    usermodel,
    fixture,
    fixture$bounded_usermodel,
    std_lv = FALSE
  )
  rust <- run_user_model_wrapper(
    usermodel_rust,
    fixture,
    fixture$bounded_usermodel,
    std_lv = FALSE
  )

  expect_equal(user_gwas_row_keys(rust), user_gwas_row_keys(old))
  expect_equal(
    user_gwas_numeric_frame(rust, user_gwas_numeric_columns(old)),
    user_gwas_numeric_frame(old, user_gwas_numeric_columns(old)),
    tolerance = 2e-6
  )
})
