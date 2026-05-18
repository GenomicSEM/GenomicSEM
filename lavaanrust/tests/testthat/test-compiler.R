test_that("lavaan_fast compiler reproduces unrestricted userGWAS implied covariance", {
  fixture <- user_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    lavaanrust::parTable_rust(fit),
    colnames(fixture$sample_cov)
  )

  expect_s3_class(compiled, "lavaan_fast_compiled")
  expect_equal(compiled$variable_names, c("SNP", "A", "B", "C", "F1"))
  expect_equal(compiled$free_labels, lavaanrust:::.lavaan_fast_free_labels(compiled))
  expect_equal(compiled$stat_names, lavaanrust:::.stat_names(compiled$observed_names))
  free_row_counts <- lengths(lapply(seq_along(compiled$free_ids), function(free_position) {
    which(compiled$free_index == free_position)
  }))
  expect_equal(compiled$free_row_offsets, c(0L, cumsum(free_row_counts)))
  expect_equal(
    compiled$free_row_indices,
    as.integer(unlist(lapply(seq_along(compiled$free_ids), function(free_position) {
      which(compiled$free_index == free_position)
    }), use.names = FALSE))
  )
  expect_equal(
    lavaanrust:::.lavaan_fast_implied_covariance(compiled),
    lavaanrust::fitted_rust(fit)$cov,
    tolerance = 1e-10
  )
  expect_equal(
    lavaanrust:::.lavaan_fast_implied_jacobian(compiled),
    lavaanrust::lavInspect_rust(fit, "delta"),
    tolerance = 1e-10
  )
  rust_surfaces <- lavaanrust:::.lavaan_fast_implied_surfaces_rust(compiled)
  rust_surfaces_flat <- lavaanrust:::.lavaan_fast_implied_surfaces_rust_flat(compiled)
  expect_equal(rust_surfaces$implied, lavaanrust::fitted_rust(fit)$cov, tolerance = 1e-10)
  expect_equal(rust_surfaces$delta, lavaanrust::lavInspect_rust(fit, "delta"), tolerance = 1e-10)
  expect_equal(
    matrix(rust_surfaces_flat$implied, nrow = compiled$n_observed),
    unname(rust_surfaces$implied),
    tolerance = 1e-10
  )
  expect_equal(
    matrix(rust_surfaces_flat$delta, nrow = compiled$n_stats),
    unname(rust_surfaces$delta),
    tolerance = 1e-10
  )
})

test_that("lavaan_fast compiler preserves fixed-measurement userGWAS structure", {
  fixture <- user_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$fixed_model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    lavaanrust::parTable_rust(fit),
    colnames(fixture$sample_cov)
  )

  expect_equal(
    lavaanrust:::.lavaan_fast_implied_covariance(compiled),
    lavaanrust::fitted_rust(fit)$cov,
    tolerance = 1e-10
  )
  expect_equal(
    lavaanrust:::.lavaan_fast_implied_jacobian(compiled),
    lavaanrust::lavInspect_rust(fit, "delta"),
    tolerance = 1e-10
  )
  rust_surfaces <- lavaanrust:::.lavaan_fast_implied_surfaces_rust(compiled)
  expect_equal(rust_surfaces$implied, lavaanrust::fitted_rust(fit)$cov, tolerance = 1e-10)
  expect_equal(rust_surfaces$delta, lavaanrust::lavInspect_rust(fit, "delta"), tolerance = 1e-10)
  expect_equal(compiled$free_ids, seq_len(6L))
})

test_that("lavaan_fast string parser handles modifiers and repeated terms", {
  fixture <- usermodel_fixture()
  model <- paste(
    "F1 =~ NA*A + start(.1)*A + start(1.1)*l2*B + l3*C",
    "F1 ~~ 1*F1",
    "A ~~ rvA*A",
    "B ~~ rvB*B",
    "C ~~ rvC*C",
    "rvA > .001",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, fixture$sample_cov)

  expect_equal(
    paste0(parsed$lhs, parsed$op, parsed$rhs),
    c("F1=~A", "F1=~B", "F1=~C", "F1~~F1", "A~~A", "B~~B", "C~~C")
  )
  expect_equal(parsed$free, c(1L, 2L, 3L, 0L, 4L, 5L, 6L))
  expect_equal(parsed$label, c("", "l2", "l3", "", "rvA", "rvB", "rvC"))
  expect_equal(parsed$ustart, c(0.1, 1.1, NA, 1, NA, NA, NA))
  expect_equal(parsed$lower, c(NA, NA, NA, NA, 0.001, NA, NA))
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(parsed, colnames(fixture$sample_cov))
  expect_equal(compiled$free_lower_bounds, c(-Inf, -Inf, -Inf, 0.001, 1e-10, 1e-10))
})

test_that("lavaan_fast parser supports box bounds and explicit label equalities", {
  fixture <- usermodel_fixture()
  model <- paste(
    "F1 =~ NA*A + l2*B + l3*C",
    "F1 ~~ 1*F1",
    "l2 > .1",
    "l3 < 2",
    "l2 == l3",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, fixture$sample_cov)
  loading_rows <- which(parsed$op == "=~")
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    parsed,
    colnames(fixture$sample_cov)
  )

  expect_equal(parsed$free[loading_rows], c(1L, 2L, 2L))
  expect_equal(parsed$lower[loading_rows], c(NA, 0.1, 0.1))
  expect_equal(parsed$upper[loading_rows], c(NA, 2, 2))
  expect_equal(parsed$lhs[parsed$op == "=="], "l2")
  expect_equal(parsed$rhs[parsed$op == "=="], "l3")
  expect_equal(compiled$free_lower_bounds[[2L]], 0.1)
  expect_equal(compiled$free_upper_bounds[[2L]], 2)
})

test_that("lavaan_fast generic string path fits labeled direct-effect RAM models", {
  fixture <- user_gwas_fixture()
  model <- paste(
    "F1 =~ 1*A + start(1.14285707631494)*l2*B + start(0.888888778086765)*l3*C",
    "F1 ~ start(0.094108558629525)*gamma*SNP",
    "A ~ start(0.02)*direct*SNP",
    "B ~ direct*SNP",
    "A ~~ start(0.60624995424363)*rvA*A",
    "B ~~ start(0.585714285882577)*rvB*B",
    "C ~~ start(0.638888930371178)*rvC*C",
    "F1 ~~ start(0.390030348982832)*psi*F1",
    "SNP ~~ 0.42*SNP",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, fixture$sample_cov)
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    parsed,
    colnames(fixture$sample_cov)
  )
  sample_cov <- lavaanrust:::.lavaan_fast_implied_covariance(compiled)
  fit <- lavaanrust::sem_rust(
    model,
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )

  expect_equal(fit@Model$model_kind, "ram_dwls_generic")
  expect_equal(lavaanrust::fitted_rust(fit)$cov, sample_cov, tolerance = 1e-8)
  expect_equal(
    lavaanrust::parTable_rust(fit)$free[match(c("A~SNP", "B~SNP"), paste0(
      lavaanrust::parTable_rust(fit)$lhs,
      lavaanrust::parTable_rust(fit)$op,
      lavaanrust::parTable_rust(fit)$rhs
    ))],
    rep(lavaanrust::parTable_rust(fit)$free[which(lavaanrust::parTable_rust(fit)$label == "direct")[[1L]]], 2L)
  )
})

test_that("lavaan_fast generic string path enforces simple lower bounds", {
  sample_cov <- matrix(1, nrow = 1L, dimnames = list("A", "A"))
  fit <- lavaanrust::sem_rust(
    paste("A ~~ rvA*A", "rvA > 1.5", sep = "\n"),
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = matrix(1, nrow = 1L)
  )

  expect_equal(lavaanrust::parTable_rust(fit)$est, 1.5, tolerance = 1e-10)
})

test_that("lavaan_fast generic string path enforces upper bounds", {
  sample_cov <- matrix(4, nrow = 1L, dimnames = list("A", "A"))
  fit <- lavaanrust::sem_rust(
    paste("A ~~ rvA*A", "rvA < 1.5", sep = "\n"),
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = matrix(1, nrow = 1L)
  )

  expect_equal(lavaanrust::parTable_rust(fit)$est[[1L]], 1.5, tolerance = 1e-10)
})

test_that("lavaan_fast parser rejects incompatible box bounds", {
  sample_cov <- matrix(1, nrow = 1L, dimnames = list("A", "A"))
  model <- paste("A ~~ rvA*A", "rvA > 2", "rvA < 1", sep = "\n")

  expect_null(lavaanrust:::.lavaan_fast_parse_model_string(model, sample_cov))
})

test_that("lavaan_fast generic parser auto-expands marker-scaled shorthand", {
  sample_cov <- diag(c(1, 1.2, 0.9, 1.1, 0.8, 1.3))
  dimnames(sample_cov) <- list(c("A", "B", "C", "D", "E", "F"), c("A", "B", "C", "D", "E", "F"))
  model <- paste(
    "F1 =~ A + B + C",
    "F2 =~ D + E + F",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, sample_cov, std.lv = FALSE)

  expect_equal(
    paste0(parsed$lhs, parsed$op, parsed$rhs),
    c(
      "F1=~A", "F1=~B", "F1=~C", "F2=~D", "F2=~E", "F2=~F",
      "A~~A", "B~~B", "C~~C", "D~~D", "E~~E", "F~~F",
      "F1~~F1", "F2~~F2", "F1~~F2"
    )
  )
  expect_equal(parsed$user, c(rep(1L, 6L), rep(0L, 9L)))
  expect_equal(parsed$free[1:6], c(0L, 1L, 2L, 0L, 3L, 4L))
  expect_equal(parsed$ustart[c(1L, 4L, 13L, 14L)], c(1, 1, 0.05, 0.05))
  expect_equal(parsed$start[7:12], unname(diag(sample_cov)) / 2)
  expect_equal(parsed$start[13:15], c(0.05, 0.05, 0))
})

test_that("lavaan_fast generic parser supports std.lv auto-identification", {
  sample_cov <- diag(c(1, 1.2, 0.9, 1.1, 0.8, 1.3))
  dimnames(sample_cov) <- list(c("A", "B", "C", "D", "E", "F"), c("A", "B", "C", "D", "E", "F"))
  model <- paste(
    "F1 =~ A + B + C",
    "F2 =~ D + E + F",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, sample_cov, std.lv = TRUE)
  latent_var_rows <- which(parsed$op == "~~" & parsed$lhs %in% c("F1", "F2") & parsed$lhs == parsed$rhs)

  expect_equal(parsed$free[1:6], seq_len(6L))
  expect_equal(parsed$free[latent_var_rows], c(0L, 0L))
  expect_equal(parsed$ustart[latent_var_rows], c(1, 1))
})

test_that("lavaan_fast generic parser auto-expands exogenous observed covariance rows", {
  sample_cov <- diag(c(1, 1.2, 0.9, 0.4, 0.6))
  dimnames(sample_cov) <- list(c("A", "B", "C", "X1", "X2"), c("A", "B", "C", "X1", "X2"))
  sample_cov["X1", "X2"] <- 0.12
  sample_cov["X2", "X1"] <- 0.12
  model <- paste(
    "F1 =~ A + B + C",
    "F1 ~ X1 + X2",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, sample_cov, std.lv = FALSE)
  row_names <- paste0(parsed$lhs, parsed$op, parsed$rhs)

  expect_true(all(c("X1~~X1", "X2~~X2", "X1~~X2") %in% row_names))
  expect_equal(parsed$start[match(c("X1~~X1", "X2~~X2", "X1~~X2"), row_names)], c(0.4, 0.6, 0.12))
})

test_that("lavaan_fast generic parser keeps lavaan auto-row order under std.lv", {
  sample_cov <- diag(c(1, 1.2, 0.9, 0.4, 0.6))
  dimnames(sample_cov) <- list(c("A", "B", "C", "X1", "X2"), c("A", "B", "C", "X1", "X2"))
  sample_cov["X1", "X2"] <- 0.12
  sample_cov["X2", "X1"] <- 0.12
  model <- paste(
    "F1 =~ A + B + C",
    "F1 ~ X1 + X2",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, sample_cov, std.lv = TRUE)

  expect_equal(
    paste0(parsed$lhs, parsed$op, parsed$rhs),
    c(
      "F1=~A", "F1=~B", "F1=~C", "F1~X1", "F1~X2",
      "A~~A", "B~~B", "C~~C", "F1~~F1",
      "X1~~X1", "X1~~X2", "X2~~X2"
    )
  )
})

test_that("lavaan_fast normalizes stale unlabeled free ids like lavaan", {
  fixture <- user_gwas_fixture()
  normalized <- lavaanrust:::.lavaan_fast_normalize_free_ids(fixture$fixed_model)
  labeled <- fixture$fixed_model
  labeled$label[c(4L, 8L)] <- "shared"
  labeled_normalized <- lavaanrust:::.lavaan_fast_normalize_free_ids(labeled)

  expect_equal(normalized$free, c(0L, 0L, 0L, seq_len(6L)))
  expect_equal(labeled_normalized$free[c(4L, 8L)], c(1L, 1L))
})

test_that("parTable_rust adds bound columns to specialized fit tables", {
  fixture <- user_gwas_fixture()
  fit <- lavaanrust::sem_rust(
    fixture$model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  par_table <- lavaanrust::parTable_rust(fit)

  expect_true("lower" %in% names(par_table))
  expect_true("upper" %in% names(par_table))
  expect_equal(names(par_table)[match("label", names(par_table)) + 1L], "lower")
  expect_equal(names(par_table)[match("lower", names(par_table)) + 1L], "upper")
  expect_true(all(is.na(par_table$lower)))
  expect_true(all(is.na(par_table$upper)))
})

test_that("lavaan_fast generic parser tolerates row-only covariance names", {
  sample_cov <- diag(c(1, 1.2, 0.9))
  rownames(sample_cov) <- c("A", "B", "C")
  model <- "F1 =~ A + B + C"
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, sample_cov, std.lv = FALSE)

  expect_equal(parsed$start[4:6], c(0.5, 0.6, 0.45))
})

test_that("lavaan_fast generic fitter can skip naive SE inversion", {
  fixture <- usermodel_fixture()
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(
    fixture$model,
    fixture$sample_cov,
    std.lv = TRUE
  )
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(parsed, colnames(fixture$sample_cov))

  with_se <- lavaanrust:::.lavaan_fast_fit_dwls_rust(
    compiled,
    fixture$sample_cov,
    fixture$wls_v
  )
  without_se <- lavaanrust:::.lavaan_fast_fit_dwls_rust(
    compiled,
    fixture$sample_cov,
    fixture$wls_v,
    compute_se = FALSE
  )

  expect_equal(without_se$estimates, with_se$estimates, tolerance = 1e-10)
  expect_equal(without_se$implied, with_se$implied, tolerance = 1e-10)
  expect_equal(without_se$naive_se, rep(0, length(without_se$estimates)))
})

test_that("lavaan_fast generic shorthand strings fit through the native RAM path", {
  sample_cov <- diag(c(1, 1.2, 0.9, 1.1, 0.8, 1.3))
  dimnames(sample_cov) <- list(c("A", "B", "C", "D", "E", "F"), c("A", "B", "C", "D", "E", "F"))
  model <- paste(
    "F1 =~ A + B + C",
    "F2 =~ D + E + F",
    sep = "\n"
  )
  parsed <- lavaanrust:::.lavaan_fast_parse_model_string(model, sample_cov, std.lv = FALSE)
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(parsed, colnames(sample_cov))
  target_cov <- lavaanrust:::.lavaan_fast_implied_covariance(compiled)
  fit <- lavaanrust::sem_rust(
    model,
    sample.cov = target_cov,
    estimator = "DWLS",
    WLS.V = diag(21)
  )

  expect_equal(fit@Model$model_kind, "ram_dwls_generic")
  expect_equal(lavaanrust::fitted_rust(fit)$cov, target_cov, tolerance = 1e-8)
  expect_equal(dim(lavaanrust::lavInspect_rust(fit, "cor.lv")), c(2L, 2L))
})

test_that("lavaan_fast generic string path evaluates defined parameters", {
  sample_cov <- matrix(
    c(0.999, 0.3996, 0.3996, 0.99884),
    nrow = 2L,
    dimnames = list(c("X", "Y"), c("X", "Y"))
  )
  model <- paste(
    "Y ~ start(.4)*b*X",
    "X ~~ start(.999)*vx*X",
    "Y ~~ start(.839)*vy*Y",
    "ind := b * vx",
    "ratio := ind / vy",
    sep = "\n"
  )
  fit <- lavaanrust::sem_rust(
    model,
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = diag(3)
  )
  refit <- lavaanrust::lavaan_rust(
    sample.cov = sample_cov,
    WLS.V = diag(3),
    slotOptions = fit@Options,
    slotParTable = fit@ParTable,
    slotData = fit@Data,
    slotModel = fit@Model
  )
  par_table <- lavaanrust::parTable_rust(fit)
  defined_rows <- which(par_table$op == ":=")
  free_values <- lavaanrust::lav_model_get_parameters_rust(fit@Model, type = "free")
  jacobian <- lavaanrust::lav_func_jacobian_complex_rust(fit@Model@def.function, free_values)

  expect_s4_class(fit@Model, "lavaan_rust_model")
  expect_equal(par_table$lhs[defined_rows], c("ind", "ratio"))
  expect_equal(par_table$rhs[defined_rows], c("b*vx", "ind/vy"))
  expect_equal(
    par_table$est[defined_rows],
    c(
      free_values[[1L]] * free_values[[2L]],
      free_values[[1L]] * free_values[[2L]] / free_values[[3L]]
    ),
    tolerance = 1e-10
  )
  expect_equal(
    unname(jacobian),
    rbind(
      c(free_values[[2L]], free_values[[1L]], 0),
      c(
        free_values[[2L]] / free_values[[3L]],
        free_values[[1L]] / free_values[[3L]],
        -(free_values[[1L]] * free_values[[2L]]) / free_values[[3L]]^2
      )
    ),
    tolerance = 1e-10
  )
  expect_equal(
    lavaanrust::parTable_rust(refit)$est[defined_rows],
    par_table$est[defined_rows],
    tolerance = 1e-10
  )
})

test_that("lavaan_fast generic parser rejects unsupported defined expressions", {
  sample_cov <- diag(c(1, 1))
  dimnames(sample_cov) <- list(c("X", "Y"), c("X", "Y"))
  model <- paste(
    "Y ~ b*X",
    "X ~~ vx*X",
    "Y ~~ vy*Y",
    "bad := system('echo nope')",
    sep = "\n"
  )

  expect_null(lavaanrust:::.lavaan_fast_parse_model_string(model, sample_cov))
})

test_that("lavaan_fast compiler rejects unsupported operators", {
  fixture <- user_gwas_fixture()
  par_table <- fixture$fixed_model
  par_table <- rbind(
    par_table,
    transform(
      par_table[1L, , drop = FALSE],
      lhs = "ghost",
      op = "<~",
      rhs = "F1~SNP",
      free = 0L
    )
  )

  expect_error(
    lavaanrust:::.lavaan_fast_compile_par_table(par_table, colnames(fixture$sample_cov)),
    "does not yet support operators"
  )
})

test_that("lavaan_fast generic DWLS optimizer reproduces fixed-measurement userGWAS fit", {
  fixture <- user_gwas_fixture()
  specialized <- lavaanrust::sem_rust(
    fixture$fixed_model,
    sample.cov = fixture$sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  par_table <- fixture$fixed_model
  free_rows <- which(par_table$free > 0L)
  par_table$free[] <- 0L
  par_table$free[free_rows] <- seq_along(free_rows)
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    par_table,
    colnames(fixture$sample_cov)
  )
  generic <- lavaanrust:::.lavaan_fast_fit_dwls_rust(
    compiled,
    fixture$sample_cov,
    fixture$wls_v
  )

  expect_equal(generic$estimates, fixture$fixed_est[free_rows], tolerance = 2e-6)
  expect_equal(generic$implied, lavaanrust::fitted_rust(specialized)$cov, tolerance = 1e-8)
  expect_equal(generic$delta, lavaanrust::lavInspect_rust(specialized, "delta"), tolerance = 1e-8)
})

test_that("lavaan_fast generic parameter-table path supports direct SNP effect", {
  fixture <- user_gwas_fixture()
  par_table <- fixture$fixed_model
  free_rows <- which(par_table$free > 0L)
  par_table$free[] <- 0L
  par_table$free[free_rows] <- seq_along(free_rows)
  extra <- transform(
    par_table[1L, , drop = FALSE],
    id = max(par_table$id) + 1L,
    lhs = "A",
    op = "~",
    rhs = "SNP",
    user = 1L,
    free = max(par_table$free) + 1L,
    ustart = NA_real_,
    plabel = paste0(".p", max(par_table$id) + 1L, "."),
    start = 0.02,
    est = 0.02,
    se = 0
  )
  par_table <- rbind(par_table, extra)
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    par_table,
    colnames(fixture$sample_cov)
  )
  sample_cov <- lavaanrust:::.lavaan_fast_implied_covariance(compiled)
  fit <- lavaanrust::sem_rust(
    par_table,
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = fixture$wls_v
  )
  refit <- lavaanrust::lavaan_rust(
    sample.cov = sample_cov,
    WLS.V = fixture$wls_v,
    slotOptions = fit@Options,
    slotParTable = fit@ParTable,
    slotData = fit@Data,
    slotModel = fit@Model
  )

  expect_equal(fit@Model$model_kind, "ram_dwls_generic")
  expect_equal(lavaanrust::fitted_rust(fit)$cov, sample_cov, tolerance = 1e-8)
  expect_equal(lavaanrust::parTable_rust(fit)$est, par_table$est, tolerance = 1e-8)
  expect_equal(lavaanrust::fitted_rust(refit)$cov, sample_cov, tolerance = 1e-8)
})

test_that("lavaan_fast generic model reuse keeps a native plan", {
  sample_cov <- matrix(
    c(1, 0.25, 0.25, 1.1),
    nrow = 2L,
    dimnames = list(c("X", "Y"), c("X", "Y"))
  )
  model <- paste(
    "Y ~ start(.1)*b*X",
    "X ~~ start(.8)*vx*X",
    "Y ~~ start(.8)*vy*Y",
    sep = "\n"
  )
  fit <- lavaanrust::sem_rust(
    model,
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = diag(3)
  )
  refit <- lavaanrust::lavaan_rust(
    sample.cov = sample_cov,
    WLS.V = diag(3),
    slotOptions = fit@Options,
    slotParTable = fit@ParTable,
    slotData = fit@Data,
    slotModel = fit@Model
  )
  direct_refit <- lavaanrust:::.fit_lavaan_fast_ram_model(
    fit@Model@par_table,
    sample_cov,
    diag(3)
  )

  expect_true(inherits(fit@Model@compiled, "lavaan_fast_compiled"))
  expect_identical(typeof(fit@Model@plan), "externalptr")
  expect_true(lavaanrust:::ram_dwls_plan_is_valid(fit@Model@plan))
  expect_identical(refit@Model@plan, fit@Model@plan)
  expect_equal(lavaanrust::parTable_rust(refit)$est, lavaanrust::parTable_rust(direct_refit)$est)
  expect_equal(lavaanrust::fitted_rust(refit)$cov, lavaanrust::fitted_rust(direct_refit)$cov)
})

test_that("lavaan_fast generic model reuse rebuilds serialized native plans", {
  sample_cov <- matrix(
    c(1, 0.25, 0.25, 1.1),
    nrow = 2L,
    dimnames = list(c("X", "Y"), c("X", "Y"))
  )
  fit <- lavaanrust::sem_rust(
    paste("Y ~ X", "X ~~ X", "Y ~~ Y", sep = "\n"),
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = diag(3)
  )
  restored_model <- unserialize(serialize(fit@Model, NULL))
  refit <- lavaanrust::lavaan_rust(
    sample.cov = sample_cov,
    WLS.V = diag(3),
    slotOptions = fit@Options,
    slotParTable = fit@ParTable,
    slotData = fit@Data,
    slotModel = restored_model
  )

  expect_false(lavaanrust:::ram_dwls_plan_is_valid(restored_model@plan))
  expect_true(lavaanrust:::ram_dwls_plan_is_valid(refit@Model@plan))
  expect_equal(lavaanrust::fitted_rust(refit)$cov, sample_cov, tolerance = 1e-8)
})

test_that("lavaan_fast generic model reuse rejects changed structural rows", {
  sample_cov <- matrix(
    c(1, 0.25, 0.25, 1.1),
    nrow = 2L,
    dimnames = list(c("X", "Y"), c("X", "Y"))
  )
  fit <- lavaanrust::sem_rust(
    paste("Y ~ X", "X ~~ X", "Y ~~ Y", sep = "\n"),
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = diag(3)
  )
  mutated_model <- fit@Model
  mutated_model@par_table$rhs[[1L]] <- "Z"

  expect_error(
    lavaanrust::lavaan_rust(
      sample.cov = sample_cov,
      WLS.V = diag(3),
      slotOptions = fit@Options,
      slotParTable = fit@ParTable,
      slotData = fit@Data,
      slotModel = mutated_model
    ),
    "cached model structure changed after compilation"
  )
})

test_that("lavaan_fast generic model reuse rejects changed fixed structural values", {
  sample_cov <- matrix(
    c(1, 0.25, 0.25, 1.1),
    nrow = 2L,
    dimnames = list(c("X", "Y"), c("X", "Y"))
  )
  fit <- lavaanrust::sem_rust(
    paste("Y ~ X", "X ~~ 1*X", "Y ~~ Y", sep = "\n"),
    sample.cov = sample_cov,
    estimator = "DWLS",
    WLS.V = diag(3)
  )
  mutated_model <- fit@Model
  fixed_row <- which(
    mutated_model@par_table$lhs == "X" &
      mutated_model@par_table$op == "~~" &
      mutated_model@par_table$rhs == "X"
  )
  mutated_model@par_table$est[[fixed_row]] <- 2

  expect_error(
    lavaanrust::lavaan_rust(
      sample.cov = sample_cov,
      WLS.V = diag(3),
      slotOptions = fit@Options,
      slotParTable = fit@ParTable,
      slotData = fit@Data,
      slotModel = mutated_model
    ),
    "cached model structure changed after compilation"
  )
})

test_that("lavaan_fast native Jacobian sums shared directed-edge parameters", {
  fixture <- user_gwas_fixture()
  par_table <- fixture$fixed_model
  free_rows <- which(par_table$free > 0L)
  par_table$free[] <- 0L
  par_table$free[free_rows] <- seq_along(free_rows)
  shared_free <- max(par_table$free) + 1L
  direct_rows <- rbind(
    transform(
      par_table[1L, , drop = FALSE],
      id = max(par_table$id) + 1L,
      lhs = "A",
      op = "~",
      rhs = "SNP",
      user = 1L,
      free = shared_free,
      ustart = NA_real_,
      plabel = paste0(".p", max(par_table$id) + 1L, "."),
      start = 0.02,
      est = 0.02,
      se = 0
    ),
    transform(
      par_table[1L, , drop = FALSE],
      id = max(par_table$id) + 2L,
      lhs = "B",
      op = "~",
      rhs = "SNP",
      user = 1L,
      free = shared_free,
      ustart = NA_real_,
      plabel = paste0(".p", max(par_table$id) + 2L, "."),
      start = 0.02,
      est = 0.02,
      se = 0
    )
  )
  compiled <- lavaanrust:::.lavaan_fast_compile_par_table(
    rbind(par_table, direct_rows),
    colnames(fixture$sample_cov)
  )
  rust_surfaces <- lavaanrust:::.lavaan_fast_implied_surfaces_rust(compiled)

  expect_equal(
    rust_surfaces$delta,
    lavaanrust:::.lavaan_fast_implied_jacobian(compiled),
    tolerance = 1e-10
  )
})
