.with_lavaan_rust_backend <- function(fun, helper_names = character()) {
  if (!requireNamespace("lavaanrust", quietly = TRUE)) {
    stop(
      "The experimental lavaanrust package must be installed before using *_rust() wrappers.",
      call. = FALSE
    )
  }

  rust_env <- new.env(parent = environment(fun))
  rust_env$sem <- lavaanrust::sem_rust
  rust_env$lavaan <- lavaanrust::lavaan_rust
  rust_env$lavInspect <- lavaanrust::lavInspect_rust
  rust_env$inspect <- lavaanrust::inspect_rust
  rust_env$parTable <- lavaanrust::parTable_rust
  rust_env$fitted <- lavaanrust::fitted_rust
  rust_env$resid <- lavaanrust::resid_rust
  rust_env$standardizedSolution <- lavaanrust::standardizedSolution_rust
  rust_env$lav_model_get_parameters <- lavaanrust::lav_model_get_parameters_rust
  rust_env$lav_func_jacobian_complex <- lavaanrust::lav_func_jacobian_complex_rust
  rust_env$.parallel_worker_packages <- "lavaanrust"
  rust_env$.parallel_worker_packages_windows <- c("lavaanrust", "gdata")
  rust_env$class <- function(x) {
    if (methods::is(x, "lavaan_rust_fit")) {
      return("lavaan")
    }

    base::class(x)
  }

  for (helper_name in helper_names) {
    helper <- get(helper_name, envir = environment(fun), inherits = TRUE)
    environment(helper) <- rust_env
    assign(helper_name, helper, envir = rust_env)
  }

  environment(fun) <- rust_env
  fun
}

# This reuses the exact original function body while rebinding only the lavaan
# surface underneath it.
commonfactor_rust <- .with_lavaan_rust_backend(commonfactor)
usermodel_rust <- function(...) {
  rust_fun <- .with_lavaan_rust_backend(
    usermodel,
    helper_names = ".rearrange"
  )

  rust_fun(...)
}
commonfactorGWAS_rust <- function(...) {
  # Supported native slice:
  # - the current one-factor commonfactorGWAS model family
  # - either `parallel = TRUE/FALSE`
  # - the DWLS SNP path backed by the strict rust sem/lavaan surface
  rust_fun <- .with_lavaan_rust_backend(
    commonfactorGWAS,
    helper_names = c(".commonfactorGWAS_main", ".rearrange")
  )

  rust_fun(...)
}

userGWAS_rust <- function(...) {
  dots <- list(...)
  parallel <- if ("parallel" %in% names(dots)) dots$parallel else formals(userGWAS)$parallel
  estimation <- if ("estimation" %in% names(dots)) dots$estimation else formals(userGWAS)$estimation
  twas <- if ("TWAS" %in% names(dots)) dots$TWAS else formals(userGWAS)$TWAS

  # Supported native slice:
  # - generic DWLS RAM models accepted by `lavaanrust`, including the current
  #   one-factor and two-factor workflow fixtures
  # - either `parallel = TRUE/FALSE`
  # - `estimation = "DWLS"`
  # - `TWAS = FALSE`
  # - either `std.lv = TRUE/FALSE`
  # - both `fix_measurement = TRUE/FALSE`
  # - both `Q_SNP = TRUE/FALSE`
  #
  # The wrapper stays strict on purpose: unsupported combinations error instead
  # of silently mixing lavaan and rust execution.
  if (!identical(estimation, "DWLS")) {
    stop(
      "userGWAS_rust() currently supports estimation = \"DWLS\" only; unsupported paths do not fall back to lavaan.",
      call. = FALSE
    )
  }

  if (isTRUE(twas)) {
    stop(
      "userGWAS_rust() currently supports TWAS = FALSE only; unsupported paths do not fall back to lavaan.",
      call. = FALSE
    )
  }

  rust_fun <- .with_lavaan_rust_backend(
    userGWAS,
    helper_names = c(".userGWAS_main", ".rearrange")
  )

  rust_fun(...)
}
