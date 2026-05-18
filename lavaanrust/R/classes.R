#' Minimal compatibility fit object for the supported lavaanrust surface.
#'
#' @exportClass lavaan_rust_fit
methods::setClass(
  "lavaan_rust_fit",
  slots = c(
    ParTable = "data.frame",
    observed = "matrix",
    implied = "list",
    delta = "matrix",
    WLS.V = "matrix",
    fit = "numeric",
    cor.lv = "matrix",
    se = "list",
    converged = "logical",
    observed_names = "character",
    latent_name = "character",
    model_kind = "character",
    Options = "list",
    Data = "list",
    Model = "ANY"
  )
)

#' Minimal compatibility model object for generic RAM fits.
#'
#' @exportClass lavaan_rust_model
methods::setClass(
  "lavaan_rust_model",
  slots = c(
    model_kind = "character",
    par_table = "data.frame",
    compiled = "ANY",
    plan = "ANY",
    free_values = "numeric",
    def.function = "function"
  )
)

methods::setMethod(
  "$",
  "lavaan_rust_model",
  function(x, name) {
    if (!name %in% methods::slotNames(x)) {
      stop("Unknown lavaan_rust_model field: ", name, call. = FALSE)
    }

    methods::slot(x, name)
  }
)
