.lav_matrix_vech_rust <- function(matrix) {
  matrix[lower.tri(matrix, diag = TRUE)]
}

.is_one_factor_dwls_model <- function(model) {
  if (!is.character(model) || length(model) != 1L) {
    return(FALSE)
  }

  lines <- trimws(strsplit(model, "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]

  if (length(lines) != 2L) {
    return(FALSE)
  }

  loading_line <- lines[grepl("=~", lines, fixed = TRUE)]
  variance_line <- lines[grepl("~~", lines, fixed = TRUE)]

  if (length(loading_line) != 1L || length(variance_line) != 1L) {
    return(FALSE)
  }

  loading_parts <- strsplit(loading_line, "=~", fixed = TRUE)[[1L]]
  if (length(loading_parts) != 2L) {
    return(FALSE)
  }

  latent <- trimws(loading_parts[[1L]])
  variance_no_space <- gsub("[[:space:]]+", "", variance_line)

  identical(variance_no_space, paste0(latent, "~~1*", latent))
}

.parse_simple_one_factor_model <- function(model) {
  if (!is.character(model) || length(model) != 1L) {
    return(NULL)
  }

  lines <- trimws(strsplit(model, "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]

  if (length(lines) != 1L || !grepl("=~", lines, fixed = TRUE)) {
    return(NULL)
  }

  parts <- trimws(strsplit(lines, "=~", fixed = TRUE)[[1L]])
  if (length(parts) != 2L) {
    return(NULL)
  }

  latent_name <- parts[[1L]]
  observed_names <- trimws(strsplit(parts[[2L]], "+", fixed = TRUE)[[1L]])
  observed_names <- observed_names[nzchar(observed_names)]

  if (
    !nzchar(latent_name) ||
      length(observed_names) < 2L ||
      any(grepl("[*:=~><]", observed_names))
  ) {
    return(NULL)
  }

  list(
    latent_name = latent_name,
    observed_names = observed_names
  )
}

.lavaan_fast_numeric_value <- function(value) {
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) == 1L && is.finite(parsed)) {
    parsed
  } else {
    NA_real_
  }
}

.lavaan_fast_parse_term <- function(term) {
  pieces <- trimws(strsplit(trimws(term), "*", fixed = TRUE)[[1L]])
  pieces <- pieces[nzchar(pieces)]

  if (!length(pieces)) {
    return(NULL)
  }

  rhs <- pieces[[length(pieces)]]
  modifiers <- if (length(pieces) > 1L) pieces[-length(pieces)] else character()
  fixed_value <- NA_real_
  start_value <- NA_real_
  label <- ""
  free_marker <- FALSE

  for (modifier in modifiers) {
    if (identical(toupper(modifier), "NA")) {
      free_marker <- TRUE
      next
    }

    start_match <- regexec("^start\\(([^()]*)\\)$", modifier)
    start_parts <- regmatches(modifier, start_match)[[1L]]
    if (length(start_parts)) {
      parsed_start <- .lavaan_fast_numeric_value(start_parts[[2L]])
      if (!is.finite(parsed_start)) {
        return(NULL)
      }

      if (is.finite(start_value) && !isTRUE(all.equal(start_value, parsed_start))) {
        return(NULL)
      }

      start_value <- parsed_start
      next
    }

    numeric_value <- .lavaan_fast_numeric_value(modifier)
    if (is.finite(numeric_value)) {
      if (is.finite(fixed_value) && !isTRUE(all.equal(fixed_value, numeric_value))) {
        return(NULL)
      }

      fixed_value <- numeric_value
      next
    }

    if (!grepl("^[A-Za-z.][A-Za-z0-9_.]*$", modifier)) {
      return(NULL)
    }

    if (nzchar(label) && !identical(label, modifier)) {
      return(NULL)
    }

    label <- modifier
  }

  if (
    !grepl("^[A-Za-z.][A-Za-z0-9_.]*$", rhs) ||
      (is.finite(fixed_value) && free_marker)
  ) {
    return(NULL)
  }

  list(
    rhs = rhs,
    fixed_value = fixed_value,
    start_value = start_value,
    label = label,
    is_free = !is.finite(fixed_value),
    free_marker = free_marker
  )
}

.lavaan_fast_allowed_definition_calls <- c("+", "-", "*", "/", "^", "(", "sqrt", "exp", "log")

.lavaan_fast_validate_definition_expression <- function(expr, allowed_symbols) {
  if (is.numeric(expr) || is.complex(expr)) {
    return(length(expr) == 1L)
  }

  if (is.name(expr)) {
    return(as.character(expr) %in% allowed_symbols)
  }

  if (!is.call(expr)) {
    return(FALSE)
  }

  call_name <- as.character(expr[[1L]])
  if (!call_name %in% .lavaan_fast_allowed_definition_calls) {
    return(FALSE)
  }

  args <- as.list(expr)[-1L]
  if (call_name %in% c("sqrt", "exp", "log", "(") && length(args) != 1L) {
    return(FALSE)
  }
  if (call_name %in% c("*", "/", "^") && length(args) != 2L) {
    return(FALSE)
  }
  if (call_name %in% c("+", "-") && !length(args) %in% c(1L, 2L)) {
    return(FALSE)
  }

  all(vapply(
    args,
    .lavaan_fast_validate_definition_expression,
    logical(1L),
    allowed_symbols = allowed_symbols
  ))
}

.lavaan_fast_parse_definition <- function(line) {
  definition_match <- regexec(
    "^([A-Za-z.][A-Za-z0-9_.]*)\\s*:=\\s*(.+)$",
    line
  )
  definition_parts <- regmatches(line, definition_match)[[1L]]
  if (!length(definition_parts)) {
    return(NULL)
  }

  expr_text <- gsub("[[:space:]]+", "", definition_parts[[3L]])
  parsed <- tryCatch(parse(text = expr_text, keep.source = FALSE), error = function(e) NULL)
  if (is.null(parsed) || length(parsed) != 1L) {
    return(NULL)
  }

  list(
    lhs = definition_parts[[2L]],
    rhs = expr_text,
    expr = parsed[[1L]]
  )
}

.lavaan_fast_append_defined_rows <- function(par_table, definitions) {
  if (!length(definitions)) {
    return(par_table)
  }

  available_symbols <- unique(par_table$label[nzchar(par_table$label)])
  seen_definitions <- character()
  rows <- vector("list", length(definitions))
  for (idx in seq_along(definitions)) {
    definition <- definitions[[idx]]
    if (
      definition$lhs %in% available_symbols ||
        definition$lhs %in% seen_definitions ||
        !.lavaan_fast_validate_definition_expression(definition$expr, available_symbols)
    ) {
      return(NULL)
    }

    row_idx <- nrow(par_table) + idx
    rows[[idx]] <- data.frame(
      id = row_idx,
      lhs = definition$lhs,
      op = ":=",
      rhs = definition$rhs,
      user = 1L,
      block = 1L,
      group = 1L,
      free = 0L,
      ustart = NA_real_,
      exo = 0L,
      label = definition$lhs,
      lower = NA_real_,
      upper = NA_real_,
      plabel = "",
      start = 0,
      est = 0,
      se = 0,
      stringsAsFactors = FALSE
    )
    available_symbols <- c(available_symbols, definition$lhs)
    seen_definitions <- c(seen_definitions, definition$lhs)
  }

  rbind(par_table, do.call(rbind, rows))
}

.lavaan_fast_definition_plan <- function(par_table, compiled) {
  defined_rows <- which(par_table$op == ":=")
  if (!length(defined_rows)) {
    return(list(
      rows = integer(),
      def.function = function(x) numeric()
    ))
  }

  structural_rows <- which(par_table$op %in% c("=~", "~", "~~"))
  label_rows <- structural_rows[nzchar(par_table$label[structural_rows])]
  label_info <- split(label_rows, par_table$label[label_rows])
  free_positions <- setNames(integer(length(label_info)), names(label_info))
  fixed_values <- setNames(numeric(length(label_info)), names(label_info))
  for (label in names(label_info)) {
    rows <- label_info[[label]]
    free_ids <- unique(par_table$free[rows][par_table$free[rows] > 0L])
    if (length(free_ids) > 1L) {
      stop("Defined-parameter label maps to multiple free parameters: ", label, call. = FALSE)
    }

    if (length(free_ids) == 1L) {
      free_position <- match(free_ids[[1L]], compiled$free_ids)
      if (is.na(free_position)) {
        stop("Defined-parameter label is not present in the compiled free vector: ", label, call. = FALSE)
      }
      free_positions[[label]] <- free_position
      fixed_values[[label]] <- NA_real_
    } else {
      free_positions[[label]] <- 0L
      fixed_values[[label]] <- .lavaan_fast_row_value(par_table, rows[[1L]])
    }
  }

  available_symbols <- names(label_info)
  expressions <- vector("list", length(defined_rows))
  names(expressions) <- par_table$lhs[defined_rows]
  if (anyDuplicated(names(expressions)) || any(names(expressions) %in% available_symbols)) {
    stop("Defined-parameter names must be unique and distinct from structural labels.", call. = FALSE)
  }

  for (idx in seq_along(defined_rows)) {
    row_idx <- defined_rows[[idx]]
    parsed <- tryCatch(parse(text = par_table$rhs[[row_idx]], keep.source = FALSE), error = function(e) NULL)
    if (
      is.null(parsed) ||
        length(parsed) != 1L ||
        !.lavaan_fast_validate_definition_expression(parsed[[1L]], available_symbols)
    ) {
      stop("Unsupported defined-parameter expression: ", par_table$rhs[[row_idx]], call. = FALSE)
    }

    expressions[[idx]] <- parsed[[1L]]
    available_symbols <- c(available_symbols, par_table$lhs[[row_idx]])
  }

  def_function <- local({
    label_free_positions <- free_positions
    label_fixed_values <- fixed_values
    definition_expressions <- expressions
    expected_free <- length(compiled$free_ids)

    function(x) {
      if (length(x) != expected_free) {
        stop("Defined-parameter evaluator received the wrong free-vector length.", call. = FALSE)
      }

      env <- new.env(parent = baseenv())
      for (label in names(label_free_positions)) {
        free_position <- label_free_positions[[label]]
        value <- if (free_position > 0L) x[[free_position]] else label_fixed_values[[label]]
        assign(label, value, envir = env)
      }

      values <- vector("list", length(definition_expressions))
      names(values) <- names(definition_expressions)
      for (idx in seq_along(definition_expressions)) {
        value <- eval(definition_expressions[[idx]], envir = env)
        if (
          length(value) != 1L ||
            !(is.numeric(value) || is.complex(value))
        ) {
          stop("Defined-parameter expression did not evaluate to a scalar number.", call. = FALSE)
        }

        values[[idx]] <- value
        assign(names(definition_expressions)[[idx]], value, envir = env)
      }

      unlist(values, use.names = TRUE)
    }
  })

  list(
    rows = defined_rows,
    def.function = def_function
  )
}

.lavaan_fast_evaluate_defined_from_labels <- function(par_table, label_values) {
  defined_rows <- which(par_table$op == ":=")
  if (!length(defined_rows)) {
    return(numeric())
  }

  available_symbols <- names(label_values)
  env <- list2env(as.list(label_values), parent = baseenv())
  values <- setNames(vector("list", length(defined_rows)), par_table$lhs[defined_rows])
  for (idx in seq_along(defined_rows)) {
    row_idx <- defined_rows[[idx]]
    parsed <- tryCatch(parse(text = par_table$rhs[[row_idx]], keep.source = FALSE), error = function(e) NULL)
    if (
      is.null(parsed) ||
        length(parsed) != 1L ||
        !.lavaan_fast_validate_definition_expression(parsed[[1L]], available_symbols)
    ) {
      stop("Unsupported defined-parameter expression: ", par_table$rhs[[row_idx]], call. = FALSE)
    }

    value <- eval(parsed[[1L]], envir = env)
    if (
      length(value) != 1L ||
        !(is.numeric(value) || is.complex(value))
    ) {
      stop("Defined-parameter expression did not evaluate to a scalar number.", call. = FALSE)
    }

    values[[idx]] <- value
    assign(par_table$lhs[[row_idx]], value, envir = env)
    available_symbols <- c(available_symbols, par_table$lhs[[row_idx]])
  }

  unlist(values, use.names = TRUE)
}

.lavaan_fast_label_groups <- function(labels, equalities = list()) {
  labels <- unique(labels[nzchar(labels)])
  if (!length(labels)) {
    return(character())
  }

  parent <- stats::setNames(labels, labels)
  find_root <- function(label) {
    while (!identical(parent[[label]], label)) {
      parent[[label]] <<- parent[[parent[[label]]]]
      label <- parent[[label]]
    }

    label
  }

  for (equality in equalities) {
    lhs <- equality$lhs
    rhs <- equality$rhs
    if (!lhs %in% labels || !rhs %in% labels) {
      next
    }

    lhs_root <- find_root(lhs)
    rhs_root <- find_root(rhs)
    if (!identical(lhs_root, rhs_root)) {
      parent[[rhs_root]] <- lhs_root
    }
  }

  vapply(labels, find_root, character(1L))
}

.lavaan_fast_normalize_free_ids <- function(par_table) {
  par_table <- as.data.frame(par_table, stringsAsFactors = FALSE)
  if (!nrow(par_table)) {
    return(par_table)
  }

  labels <- par_table$label[nzchar(par_table$label)]
  equalities <- lapply(which(par_table$op == "=="), function(row_idx) {
    list(lhs = par_table$lhs[[row_idx]], rhs = par_table$rhs[[row_idx]])
  })
  label_groups <- .lavaan_fast_label_groups(labels, equalities)

  next_free <- 1L
  label_free <- integer()
  normalized <- integer(nrow(par_table))
  for (row_idx in seq_len(nrow(par_table))) {
    if (par_table$free[[row_idx]] <= 0L) {
      next
    }

    label <- par_table$label[[row_idx]]
    free_key <- if (nzchar(label) && label %in% names(label_groups)) label_groups[[label]] else label
    if (nzchar(free_key) && free_key %in% names(label_free)) {
      normalized[[row_idx]] <- label_free[[free_key]]
      next
    }

    normalized[[row_idx]] <- next_free
    if (nzchar(free_key)) {
      label_free[[free_key]] <- next_free
    }
    next_free <- next_free + 1L
  }

  par_table$free <- normalized
  par_table
}

.lavaan_fast_default_start <- function(op, lhs, rhs, observed_names, sample_cov) {
  if (identical(op, "=~")) {
    return(1)
  }

  if (identical(op, "~~") && identical(lhs, rhs)) {
    if (lhs %in% observed_names) {
      return(sample_cov[lhs, lhs])
    }

    return(1)
  }

  0
}

.lavaan_fast_merge_term_rows <- function(rows) {
  merged <- list()
  positions <- list()

  for (row in rows) {
    key <- paste(row$lhs, row$op, row$rhs, sep = "\r")
    position <- positions[[key]]

    if (is.null(position)) {
      merged[[length(merged) + 1L]] <- row
      positions[[key]] <- length(merged)
      next
    }

    previous <- merged[[position]]
    if (
      isTRUE(previous$is_free != row$is_free) ||
        (
          !previous$is_free &&
            !isTRUE(all.equal(previous$fixed_value, row$fixed_value))
        ) ||
        (
          nzchar(previous$label) &&
            nzchar(row$label) &&
            !identical(previous$label, row$label)
        ) ||
        (
          is.finite(previous$start_value) &&
            is.finite(row$start_value) &&
            !isTRUE(all.equal(previous$start_value, row$start_value))
        )
    ) {
      return(NULL)
    }

    if (!nzchar(previous$label)) {
      previous$label <- row$label
    }

    if (!is.finite(previous$start_value)) {
      previous$start_value <- row$start_value
    }

    merged[[position]] <- previous
  }

  merged
}

.lavaan_fast_edge_exists <- function(rows, lhs, op, rhs) {
  any(vapply(rows, function(row) {
    identical(row$lhs, lhs) &&
      identical(row$op, op) &&
      identical(row$rhs, rhs)
  }, logical(1L)))
}

.lavaan_fast_new_auto_row <- function(lhs, op, rhs, is_free, fixed_value = NA_real_, start_value = NA_real_) {
  list(
    lhs = lhs,
    op = op,
    rhs = rhs,
    fixed_value = fixed_value,
    start_value = start_value,
    label = "",
    is_free = is_free,
    free_marker = FALSE,
    user = 0L
  )
}

.lavaan_fast_apply_identification <- function(rows, std.lv) {
  latent_names <- unique(vapply(rows, function(row) {
    if (identical(row$op, "=~")) row$lhs else ""
  }, character(1L)))
  latent_names <- latent_names[nzchar(latent_names)]

  if (isTRUE(std.lv)) {
    return(rows)
  }

  for (latent_name in latent_names) {
    loading_rows <- which(vapply(rows, function(row) {
      identical(row$op, "=~") && identical(row$lhs, latent_name)
    }, logical(1L)))
    if (!length(loading_rows)) {
      next
    }

    row_idx <- loading_rows[[1L]]
    row <- rows[[row_idx]]
    if (row$is_free && !isTRUE(row$free_marker)) {
      row$is_free <- FALSE
      row$fixed_value <- if (is.finite(row$start_value)) row$start_value else 1
      rows[[row_idx]] <- row
    }
  }

  rows
}

.lavaan_fast_expand_auto_rows <- function(rows, observed_names, sample_cov, std.lv) {
  rows <- .lavaan_fast_apply_identification(rows, std.lv)
  latent_names <- unique(vapply(rows, function(row) {
    if (identical(row$op, "=~")) row$lhs else ""
  }, character(1L)))
  latent_names <- latent_names[nzchar(latent_names)]
  indicator_names <- unique(vapply(rows, function(row) {
    if (identical(row$op, "=~")) row$rhs else ""
  }, character(1L)))
  indicator_names <- indicator_names[nzchar(indicator_names)]
  observed_regression_lhs <- unique(vapply(rows, function(row) {
    if (identical(row$op, "~") && row$lhs %in% observed_names) row$lhs else ""
  }, character(1L)))
  observed_regression_lhs <- observed_regression_lhs[nzchar(observed_regression_lhs)]
  observed_endogenous <- unique(c(indicator_names, observed_regression_lhs))
  observed_exogenous <- setdiff(observed_names, observed_endogenous)

  for (observed_name in observed_endogenous) {
    if (.lavaan_fast_edge_exists(rows, observed_name, "~~", observed_name)) {
      next
    }

    rows[[length(rows) + 1L]] <- .lavaan_fast_new_auto_row(
      observed_name,
      "~~",
      observed_name,
      is_free = TRUE,
      start_value = sample_cov[observed_name, observed_name] / 2
    )
  }

  for (latent_name in latent_names) {
    if (.lavaan_fast_edge_exists(rows, latent_name, "~~", latent_name)) {
      next
    }

    rows[[length(rows) + 1L]] <- if (isTRUE(std.lv)) {
      .lavaan_fast_new_auto_row(
        latent_name,
        "~~",
        latent_name,
        is_free = FALSE,
        fixed_value = 1,
        start_value = 1
      )
    } else {
      .lavaan_fast_new_auto_row(
        latent_name,
        "~~",
        latent_name,
        is_free = TRUE,
        start_value = 0.05
      )
    }
  }

  if (length(observed_exogenous)) {
    for (col in seq_along(observed_exogenous)) {
      lhs <- observed_exogenous[[col]]
      if (!.lavaan_fast_edge_exists(rows, lhs, "~~", lhs)) {
        rows[[length(rows) + 1L]] <- .lavaan_fast_new_auto_row(
          lhs,
          "~~",
          lhs,
          is_free = TRUE,
          start_value = sample_cov[lhs, lhs]
        )
      }

      if (col == length(observed_exogenous)) {
        next
      }

      for (row in (col + 1L):length(observed_exogenous)) {
        rhs <- observed_exogenous[[row]]
        if (.lavaan_fast_edge_exists(rows, lhs, "~~", rhs) || .lavaan_fast_edge_exists(rows, rhs, "~~", lhs)) {
          next
        }

        rows[[length(rows) + 1L]] <- .lavaan_fast_new_auto_row(
          lhs,
          "~~",
          rhs,
          is_free = TRUE,
          start_value = sample_cov[lhs, rhs]
        )
      }
    }
  }

  latent_regression_lhs <- unique(vapply(rows, function(row) {
    if (identical(row$op, "~") && row$lhs %in% latent_names) row$lhs else ""
  }, character(1L)))
  latent_regression_lhs <- latent_regression_lhs[nzchar(latent_regression_lhs)]
  exogenous_latents <- setdiff(latent_names, latent_regression_lhs)

  if (length(exogenous_latents) > 1L) {
    for (col in seq_len(length(exogenous_latents) - 1L)) {
      for (row in (col + 1L):length(exogenous_latents)) {
        lhs <- exogenous_latents[[col]]
        rhs <- exogenous_latents[[row]]
        if (.lavaan_fast_edge_exists(rows, lhs, "~~", rhs) || .lavaan_fast_edge_exists(rows, rhs, "~~", lhs)) {
          next
        }

        rows[[length(rows) + 1L]] <- .lavaan_fast_new_auto_row(
          lhs,
          "~~",
          rhs,
          is_free = TRUE,
          start_value = 0
        )
      }
    }
  }

  rows
}

.lavaan_fast_parse_model_string <- function(model, sample_cov, std.lv = FALSE) {
  if (
    !is.character(model) ||
      length(model) != 1L ||
      !is.matrix(sample_cov)
  ) {
    return(NULL)
  }

  observed_names <- colnames(sample_cov)
  if (is.null(observed_names)) {
    observed_names <- rownames(sample_cov)
  }
  if (is.null(observed_names)) {
    return(NULL)
  }
  if (is.null(colnames(sample_cov))) {
    colnames(sample_cov) <- observed_names
  }
  if (is.null(rownames(sample_cov))) {
    rownames(sample_cov) <- observed_names
  }

  lines <- trimws(strsplit(model, "\n", fixed = TRUE)[[1L]])
  lines <- trimws(sub("#.*$", "", lines))
  lines <- lines[nzchar(lines)]

  if (!length(lines)) {
    return(NULL)
  }

  rows <- list()
  definitions <- list()
  lower_bounds <- numeric()
  upper_bounds <- numeric()
  equalities <- list()
  for (line in lines) {
    definition <- .lavaan_fast_parse_definition(line)
    if (!is.null(definition)) {
      definitions[[length(definitions) + 1L]] <- definition
      next
    }

    if (grepl(":=", line, fixed = TRUE)) {
      return(NULL)
    }

    equality_match <- regexec(
      "^([A-Za-z.][A-Za-z0-9_.]*)\\s*==\\s*([A-Za-z.][A-Za-z0-9_.]*)$",
      line
    )
    equality_parts <- regmatches(line, equality_match)[[1L]]
    if (length(equality_parts)) {
      equalities[[length(equalities) + 1L]] <- list(
        lhs = equality_parts[[2L]],
        rhs = equality_parts[[3L]]
      )
      next
    }

    bound_match <- regexec(
      "^([A-Za-z.][A-Za-z0-9_.]*)\\s*(>|<)\\s*([+-]?(?:[0-9]*\\.?[0-9]+)(?:[eE][+-]?[0-9]+)?)$",
      line
    )
    bound_parts <- regmatches(line, bound_match)[[1L]]
    if (length(bound_parts)) {
      label <- bound_parts[[2L]]
      op <- bound_parts[[3L]]
      value <- .lavaan_fast_numeric_value(bound_parts[[4L]])
      if (!is.finite(value)) {
        return(NULL)
      }

      if (identical(op, ">")) {
        if (label %in% names(lower_bounds) && !isTRUE(all.equal(lower_bounds[[label]], value))) {
          return(NULL)
        }
        lower_bounds[[label]] <- value
      } else {
        if (label %in% names(upper_bounds) && !isTRUE(all.equal(upper_bounds[[label]], value))) {
          return(NULL)
        }
        upper_bounds[[label]] <- value
      }
      next
    }

    if (grepl("(==|>|<)", line)) {
      return(NULL)
    }

    op <- if (grepl("=~", line, fixed = TRUE)) {
      "=~"
    } else if (grepl("~~", line, fixed = TRUE)) {
      "~~"
    } else if (grepl("~", line, fixed = TRUE)) {
      "~"
    } else {
      ""
    }
    if (!nzchar(op)) {
      return(NULL)
    }

    parts <- trimws(strsplit(line, op, fixed = TRUE)[[1L]])
    if (
      length(parts) != 2L ||
        !grepl("^[A-Za-z.][A-Za-z0-9_.]*$", parts[[1L]])
    ) {
      return(NULL)
    }

    lhs <- parts[[1L]]
    rhs_terms <- trimws(strsplit(parts[[2L]], "+", fixed = TRUE)[[1L]])
    rhs_terms <- rhs_terms[nzchar(rhs_terms)]
    if (!length(rhs_terms)) {
      return(NULL)
    }

    for (rhs_term in rhs_terms) {
      parsed_term <- .lavaan_fast_parse_term(rhs_term)
      if (is.null(parsed_term)) {
        return(NULL)
      }

      rows[[length(rows) + 1L]] <- c(
        list(lhs = lhs, op = op, user = 1L),
        parsed_term
      )
    }
  }

  rows <- .lavaan_fast_merge_term_rows(rows)
  if (is.null(rows)) {
    return(NULL)
  }
  rows <- .lavaan_fast_expand_auto_rows(rows, observed_names, sample_cov, std.lv)
  row_labels <- vapply(rows, `[[`, character(1L), "label")
  constrained_labels <- unique(c(
    names(lower_bounds),
    names(upper_bounds),
    unlist(equalities, use.names = FALSE)
  ))
  if (length(constrained_labels) && !all(constrained_labels %in% row_labels[nzchar(row_labels)])) {
    return(NULL)
  }
  label_groups <- .lavaan_fast_label_groups(row_labels, equalities)
  group_lower_bounds <- numeric()
  group_upper_bounds <- numeric()
  for (label in names(label_groups)) {
    group <- label_groups[[label]]
    if (label %in% names(lower_bounds)) {
      previous_lower <- if (group %in% names(group_lower_bounds)) group_lower_bounds[[group]] else -Inf
      group_lower_bounds[[group]] <- max(previous_lower, lower_bounds[[label]])
    }
    if (label %in% names(upper_bounds)) {
      previous_upper <- if (group %in% names(group_upper_bounds)) group_upper_bounds[[group]] else Inf
      group_upper_bounds[[group]] <- min(previous_upper, upper_bounds[[label]])
    }
  }
  overlapping_bounds <- intersect(names(group_lower_bounds), names(group_upper_bounds))
  if (length(overlapping_bounds) && any(group_lower_bounds[overlapping_bounds] > group_upper_bounds[overlapping_bounds])) {
    return(NULL)
  }

  free_labels <- integer()
  next_free <- 1L
  n_rows <- length(rows)
  lhs <- character(n_rows)
  op <- character(n_rows)
  rhs <- character(n_rows)
  user <- integer(n_rows)
  free <- integer(n_rows)
  ustart <- rep(NA_real_, n_rows)
  label <- character(n_rows)
  lower <- rep(NA_real_, n_rows)
  upper <- rep(NA_real_, n_rows)
  start <- numeric(n_rows)
  for (row_idx in seq_along(rows)) {
    row <- rows[[row_idx]]
    free_value <- 0L

    if (row$is_free) {
      free_key <- if (nzchar(row$label)) label_groups[[row$label]] else ""
      if (nzchar(free_key) && free_key %in% names(free_labels)) {
        free_value <- free_labels[[free_key]]
      } else {
        free_value <- next_free
        if (nzchar(free_key)) {
          free_labels[[free_key]] <- free_value
        }
        next_free <- next_free + 1L
      }
    }

    start_value <- if (is.finite(row$fixed_value)) {
      row$fixed_value
    } else if (is.finite(row$start_value)) {
      row$start_value
    } else {
      .lavaan_fast_default_start(row$op, row$lhs, row$rhs, observed_names, sample_cov)
    }

    lhs[[row_idx]] <- row$lhs
    op[[row_idx]] <- row$op
    rhs[[row_idx]] <- row$rhs
    user[[row_idx]] <- row$user
    free[[row_idx]] <- free_value
    ustart[[row_idx]] <- if (row$is_free) row$start_value else row$fixed_value
    label[[row_idx]] <- row$label
    lower[[row_idx]] <- if (nzchar(row$label) && label_groups[[row$label]] %in% names(group_lower_bounds)) {
        group_lower_bounds[[label_groups[[row$label]]]]
      } else {
        NA_real_
      }
    upper[[row_idx]] <- if (nzchar(row$label) && label_groups[[row$label]] %in% names(group_upper_bounds)) {
        group_upper_bounds[[label_groups[[row$label]]]]
      } else {
        NA_real_
      }
    start[[row_idx]] <- start_value
  }

  par_table <- data.frame(
    id = seq_len(n_rows),
    lhs = lhs,
    op = op,
    rhs = rhs,
    user = user,
    block = rep.int(1L, n_rows),
    group = rep.int(1L, n_rows),
    free = free,
    ustart = ustart,
    exo = integer(n_rows),
    label = label,
    lower = lower,
    upper = upper,
    plabel = paste0(".p", seq_len(n_rows), "."),
    start = start,
    est = start,
    se = numeric(n_rows),
    stringsAsFactors = FALSE
  )
  if (length(equalities)) {
    equality_rows <- lapply(seq_along(equalities), function(idx) {
      equality <- equalities[[idx]]
      row_idx <- nrow(par_table) + idx
      data.frame(
        id = row_idx,
        lhs = equality$lhs,
        op = "==",
        rhs = equality$rhs,
        user = 1L,
        block = 1L,
        group = 1L,
        free = 0L,
        ustart = NA_real_,
        exo = 0L,
        label = "",
        lower = NA_real_,
        upper = NA_real_,
        plabel = "",
        start = 0,
        est = 0,
        se = 0,
        stringsAsFactors = FALSE
      )
    })
    par_table <- rbind(par_table, do.call(rbind, equality_rows))
  }

  .lavaan_fast_append_defined_rows(par_table, definitions)
}

.one_factor_latent_name <- function(model) {
  loading_line <- trimws(strsplit(model, "\n", fixed = TRUE)[[1L]])
  loading_line <- loading_line[grepl("=~", loading_line, fixed = TRUE)][[1L]]
  trimws(strsplit(loading_line, "=~", fixed = TRUE)[[1L]][[1L]])
}

.stat_names <- function(observed_names) {
  k <- length(observed_names)
  n_stats <- k * (k + 1L) / 2L
  stat_names <- character(n_stats)
  idx <- 1L

  for (col in seq_len(k)) {
    for (row in col:k) {
      stat_names[[idx]] <- paste0(observed_names[[col]], "~~", observed_names[[row]])
      idx <- idx + 1L
    }
  }

  stat_names
}

.vech_rust <- function(matrix) {
  matrix[lower.tri(matrix, diag = TRUE)]
}

.one_factor_par_table <- function(observed_names, latent_name, fit) {
  k <- length(observed_names)
  ids <- seq_len(2L * k + 1L)
  free <- c(seq_len(k), 0L, k + seq_len(k))
  est <- c(fit$loadings, 1, fit$residuals)
  se <- c(fit$naive_se[seq_len(k)], 0, fit$naive_se[k + seq_len(k)])
  starts <- c(fit$loadings, 1, fit$residuals)

  data.frame(
    id = ids,
    lhs = c(rep(latent_name, k), latent_name, observed_names),
    op = c(rep("=~", k), "~~", rep("~~", k)),
    rhs = c(observed_names, latent_name, observed_names),
    user = c(rep(1L, k + 1L), rep(0L, k)),
    block = rep(1L, 2L * k + 1L),
    group = rep(1L, 2L * k + 1L),
    free = free,
    ustart = c(rep(NA_real_, k), 1, rep(NA_real_, k)),
    exo = rep(0L, 2L * k + 1L),
    label = rep("", 2L * k + 1L),
    plabel = paste0(".p", ids, "."),
    start = starts,
    est = est,
    se = se,
    stringsAsFactors = FALSE
  )
}

.new_one_factor_fit <- function(model, sample.cov, WLS.V, fit) {
  observed_names <- colnames(sample.cov)
  if (is.null(observed_names)) {
    observed_names <- rownames(sample.cov)
  }
  if (is.null(observed_names)) {
    observed_names <- paste0("V", seq_len(ncol(sample.cov)))
  }

  dimnames(sample.cov) <- list(observed_names, observed_names)
  latent_name <- .one_factor_latent_name(model)
  k <- ncol(sample.cov)
  n_stats <- k * (k + 1L) / 2L
  delta <- matrix(fit$delta, nrow = n_stats, ncol = 2L * k)
  implied <- matrix(fit$implied, nrow = k, ncol = k)
  dimnames(implied) <- dimnames(sample.cov)
  dimnames(delta) <- list(
    .stat_names(observed_names),
    c(
      paste0(latent_name, "=~", observed_names),
      paste0(observed_names, "~~", observed_names)
    )
  )
  cor.lv <- matrix(1, nrow = 1L, ncol = 1L, dimnames = list(latent_name, latent_name))
  par_table <- .one_factor_par_table(observed_names, latent_name, fit)
  npar <- 2L * k
  df <- n_stats - npar
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = df,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = sample.cov,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = cor.lv,
    se = list(theta = implied),
    converged = isTRUE(fit$converged),
    observed_names = observed_names,
    latent_name = latent_name,
    model_kind = "one_factor_dwls",
    Options = list(estimator = "DWLS"),
    Data = list(observed_names = observed_names),
    Model = list(model_kind = "one_factor_dwls", model = model)
  )
}

.simple_std_lv_one_factor_par_table <- function(observed_names, latent_name, fit) {
  k <- length(observed_names)
  ids <- seq_len(2L * k + 1L)
  free <- c(seq_len(k), k + seq_len(k), 0L)
  est <- c(fit$loadings, fit$residuals, 1)
  se <- c(fit$naive_se[seq_len(k)], fit$naive_se[k + seq_len(k)], 0)

  data.frame(
    id = ids,
    lhs = c(rep(latent_name, k), observed_names, latent_name),
    op = c(rep("=~", k), rep("~~", k), "~~"),
    rhs = c(observed_names, observed_names, latent_name),
    user = c(rep(1L, k), rep(0L, k + 1L)),
    block = rep(1L, 2L * k + 1L),
    group = rep(1L, 2L * k + 1L),
    free = free,
    ustart = c(rep(NA_real_, 2L * k), 1),
    exo = rep(0L, 2L * k + 1L),
    label = rep("", 2L * k + 1L),
    plabel = paste0(".p", ids, "."),
    start = est,
    est = est,
    se = se,
    stringsAsFactors = FALSE
  )
}

.new_simple_std_lv_one_factor_fit <- function(model, sample.cov, WLS.V, fit) {
  spec <- .parse_simple_one_factor_model(model)
  observed_names <- spec$observed_names
  dimnames(sample.cov) <- list(observed_names, observed_names)
  k <- length(observed_names)
  n_stats <- k * (k + 1L) / 2L
  delta <- matrix(fit$delta, nrow = n_stats, ncol = 2L * k)
  implied <- matrix(fit$implied, nrow = k, ncol = k)
  dimnames(implied) <- dimnames(sample.cov)
  dimnames(delta) <- list(
    .stat_names(observed_names),
    c(
      paste0(spec$latent_name, "=~", observed_names),
      paste0(observed_names, "~~", observed_names)
    )
  )
  par_table <- .simple_std_lv_one_factor_par_table(observed_names, spec$latent_name, fit)
  fit_stats <- c(
    npar = 2L * k,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = n_stats - 2L * k,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = sample.cov,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = matrix(1, nrow = 1L, ncol = 1L, dimnames = list(spec$latent_name, spec$latent_name)),
    se = list(theta = implied),
    converged = isTRUE(fit$converged),
    observed_names = observed_names,
    latent_name = spec$latent_name,
    model_kind = "simple_std_lv_one_factor_dwls",
    Options = list(estimator = "DWLS"),
    Data = list(observed_names = observed_names),
    Model = list(model_kind = "simple_std_lv_one_factor_dwls", model = model)
  )
}

.marker_one_factor_par_table <- function(observed_names, latent_name, fit) {
  k <- length(observed_names)
  ids <- seq_len(2L * k + 1L)
  free <- c(0L, if (k > 1L) seq_len(k - 1L) else integer(), k - 1L + seq_len(k), 2L * k)
  est <- c(fit$loadings, fit$residuals, fit$phi)
  se <- c(0, if (k > 1L) fit$naive_se[seq_len(k - 1L)] else numeric(), fit$naive_se[k - 1L + seq_len(k)], fit$naive_se[[2L * k]])

  data.frame(
    id = ids,
    lhs = c(rep(latent_name, k), observed_names, latent_name),
    op = c(rep("=~", k), rep("~~", k), "~~"),
    rhs = c(observed_names, observed_names, latent_name),
    user = c(rep(1L, k), rep(0L, k + 1L)),
    block = rep(1L, 2L * k + 1L),
    group = rep(1L, 2L * k + 1L),
    free = free,
    ustart = c(1, rep(NA_real_, k - 1L), rep(NA_real_, k + 1L)),
    exo = rep(0L, 2L * k + 1L),
    label = rep("", 2L * k + 1L),
    plabel = paste0(".p", ids, "."),
    start = est,
    est = est,
    se = se,
    stringsAsFactors = FALSE
  )
}

.marker_one_factor_delta <- function(loadings, phi) {
  k <- length(loadings)
  delta <- matrix(0, nrow = k * (k + 1L) / 2L, ncol = 2L * k)
  row_idx <- 1L

  for (col_idx in seq_len(k)) {
    for (row_obs in col_idx:k) {
      if (row_obs > 1L) {
        delta[row_idx, row_obs - 1L] <- delta[row_idx, row_obs - 1L] + phi * loadings[[col_idx]]
      }

      if (col_idx > 1L) {
        delta[row_idx, col_idx - 1L] <- delta[row_idx, col_idx - 1L] + phi * loadings[[row_obs]]
      }

      if (identical(row_obs, col_idx)) {
        delta[row_idx, k - 1L + row_obs] <- 1
      }

      delta[row_idx, 2L * k] <- loadings[[row_obs]] * loadings[[col_idx]]
      row_idx <- row_idx + 1L
    }
  }

  delta
}

.new_marker_one_factor_fit <- function(model, sample.cov, WLS.V, fit) {
  spec <- .parse_simple_one_factor_model(model)
  observed_names <- spec$observed_names
  dimnames(sample.cov) <- list(observed_names, observed_names)
  marker_loading <- fit$loadings[[1L]]

  if (abs(marker_loading) < sqrt(.Machine$double.eps)) {
    stop("Marker-scaled one-factor rust slice requires a non-zero first loading.", call. = FALSE)
  }

  loadings <- fit$loadings / marker_loading
  phi <- marker_loading^2
  implied <- matrix(fit$implied, nrow = length(observed_names), ncol = length(observed_names))
  dimnames(implied) <- dimnames(sample.cov)
  delta <- .marker_one_factor_delta(loadings, phi)
  dimnames(delta) <- list(
    .stat_names(observed_names),
    c(
      paste0(spec$latent_name, "=~", observed_names[-1L]),
      paste0(observed_names, "~~", observed_names),
      paste0(spec$latent_name, "~~", spec$latent_name)
    )
  )
  bread <- solve(t(delta) %*% WLS.V %*% delta)
  marker_fit <- list(
    loadings = loadings,
    residuals = fit$residuals,
    phi = phi,
    naive_se = sqrt(pmax(diag(bread), 0)),
    objective = fit$objective,
    srmr = fit$srmr,
    converged = fit$converged
  )
  par_table <- .marker_one_factor_par_table(observed_names, spec$latent_name, marker_fit)
  npar <- 2L * length(observed_names)
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = nrow(delta) - npar,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = sample.cov,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = matrix(phi, nrow = 1L, ncol = 1L, dimnames = list(spec$latent_name, spec$latent_name)),
    se = list(theta = implied),
    converged = isTRUE(fit$converged),
    observed_names = observed_names,
    latent_name = spec$latent_name,
    model_kind = "marker_one_factor_dwls",
    Options = list(estimator = "DWLS"),
    Data = list(observed_names = observed_names),
    Model = list(model_kind = "marker_one_factor_dwls", model = model)
  )
}

.is_commonfactor_null_model <- function(model) {
  is.character(model) &&
    length(model) == 1L &&
    grepl("VF1 =~ 1*", model, fixed = TRUE) &&
    grepl("~~ 0*", model, fixed = TRUE) &&
    !grepl(":=", model, fixed = TRUE)
}

.commonfactor_null_par_table <- function(observed_names) {
  k <- length(observed_names)
  latent_names <- paste0("VF", seq_len(k))
  rows <- vector("list", 0L)

  add_row <- function(lhs, op, rhs, free, ustart, est) {
    rows[[length(rows) + 1L]] <<- list(
      lhs = lhs,
      op = op,
      rhs = rhs,
      free = free,
      ustart = ustart,
      est = est
    )
  }

  for (idx in seq_len(k)) {
    add_row(observed_names[[idx]], "~~", observed_names[[idx]], idx, NA_real_, NA_real_)
  }

  for (idx in seq_len(k)) {
    add_row(latent_names[[idx]], "=~", observed_names[[idx]], 0L, 1, 1)
  }

  if (k >= 2L) {
    for (lhs in seq_len(k - 1L)) {
      for (rhs in seq.int(lhs + 1L, k)) {
        add_row(latent_names[[lhs]], "~~", latent_names[[rhs]], 0L, 0, 0)
      }
    }

    for (lhs in seq_len(k - 1L)) {
      for (rhs in seq.int(lhs + 1L, k)) {
        add_row(observed_names[[lhs]], "~~", observed_names[[rhs]], 0L, 0, 0)
      }
    }
  }

  for (idx in seq_len(k)) {
    add_row(latent_names[[idx]], "~~", latent_names[[idx]], 0L, 0, 0)
  }

  raw <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  n <- nrow(raw)
  data.frame(
    id = seq_len(n),
    lhs = raw$lhs,
    op = raw$op,
    rhs = raw$rhs,
    user = rep(1L, n),
    block = rep(1L, n),
    group = rep(1L, n),
    free = as.integer(raw$free),
    ustart = as.numeric(raw$ustart),
    exo = rep(0L, n),
    label = rep("", n),
    plabel = paste0(".p", seq_len(n), "."),
    start = as.numeric(raw$est),
    est = as.numeric(raw$est),
    se = rep(0, n),
    stringsAsFactors = FALSE
  )
}

.observed_covariance_spec <- function(par_table, sample.cov) {
  observed_names <- colnames(sample.cov)
  if (is.null(observed_names)) {
    observed_names <- rownames(sample.cov)
  }
  if (is.null(observed_names)) {
    observed_names <- paste0("V", seq_len(ncol(sample.cov)))
  }

  dimnames(sample.cov) <- list(observed_names, observed_names)
  stat_names <- .stat_names(observed_names)
  free_mask <- integer(length(stat_names))
  fixed_values <- numeric(length(stat_names))
  row_lookup <- integer(length(stat_names))
  loading_rows <- par_table[par_table$op == "=~", , drop = FALSE]
  latent_for_observed <- stats::setNames(loading_rows$lhs, loading_rows$rhs)

  par_value <- function(row_idx) {
    if (!length(row_idx)) {
      return(0)
    }

    if (!is.na(par_table$ustart[[row_idx]])) {
      return(par_table$ustart[[row_idx]])
    }

    if (!is.na(par_table$est[[row_idx]])) {
      return(par_table$est[[row_idx]])
    }

    0
  }

  for (idx in seq_along(stat_names)) {
    parts <- strsplit(stat_names[[idx]], "~~", fixed = TRUE)[[1L]]
    direct_row <- which(
      par_table$op == "~~" &
        par_table$lhs == parts[[1L]] &
        par_table$rhs == parts[[2L]]
    )

    if (!length(direct_row)) {
      direct_row <- which(
        par_table$op == "~~" &
          par_table$lhs == parts[[2L]] &
          par_table$rhs == parts[[1L]]
      )
    }

    if (!length(direct_row)) {
      stop("Observed covariance row missing from parameter table: ", stat_names[[idx]], call. = FALSE)
    }

    direct_row <- direct_row[[1L]]

    if (identical(parts[[1L]], parts[[2L]]) && parts[[1L]] %in% names(latent_for_observed)) {
      latent <- latent_for_observed[[parts[[1L]]]]
      latent_var_row <- which(
        par_table$op == "~~" &
          par_table$lhs == latent &
          par_table$rhs == latent
      )
      latent_var_row <- if (length(latent_var_row)) latent_var_row[[1L]] else integer()

      direct_free <- par_table$free[[direct_row]] > 0L
      latent_free <- length(latent_var_row) && par_table$free[[latent_var_row]] > 0L

      if (direct_free && latent_free) {
        stop("Unsupported model: both observed and latent variance rows are free", call. = FALSE)
      }

      if (direct_free) {
        free_mask[[idx]] <- 1L
        row_lookup[[idx]] <- direct_row
        fixed_values[[idx]] <- par_value(latent_var_row)
      } else if (latent_free) {
        free_mask[[idx]] <- 1L
        row_lookup[[idx]] <- latent_var_row
        fixed_values[[idx]] <- par_value(direct_row)
      } else {
        fixed_values[[idx]] <- par_value(direct_row) + par_value(latent_var_row)
      }
    } else {
      row_lookup[[idx]] <- direct_row
      free_mask[[idx]] <- as.integer(par_table$free[[direct_row]] > 0L)

      if (!free_mask[[idx]]) {
        fixed_values[[idx]] <- par_value(direct_row)
      }
    }
  }

  list(
    observed_names = observed_names,
    sample_cov = sample.cov,
    free_mask = free_mask,
    fixed_values = fixed_values,
    row_lookup = row_lookup,
    stat_names = stat_names
  )
}

.new_observed_covariance_fit <- function(par_table, sample.cov, WLS.V, fit, model_kind) {
  spec <- .observed_covariance_spec(par_table, sample.cov)
  implied <- matrix(fit$implied, nrow = nrow(sample.cov), ncol = ncol(sample.cov))
  dimnames(implied) <- list(spec$observed_names, spec$observed_names)
  delta <- matrix(fit$delta, nrow = length(spec$stat_names))
  dimnames(delta) <- list(spec$stat_names, spec$stat_names[spec$free_mask > 0L])
  par_table <- as.data.frame(par_table, stringsAsFactors = FALSE)
  free_rows <- spec$row_lookup[spec$free_mask > 0L]
  par_table$est[free_rows] <- fit$estimates
  par_table$se[] <- 0
  par_table$se[free_rows] <- fit$naive_se
  npar <- sum(spec$free_mask)
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = length(spec$stat_names) - npar,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = spec$sample_cov,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = matrix(numeric(), nrow = 0L, ncol = 0L),
    se = list(theta = matrix(0, nrow = nrow(sample.cov), ncol = ncol(sample.cov))),
    converged = isTRUE(fit$converged),
    observed_names = spec$observed_names,
    latent_name = character(),
    model_kind = model_kind,
    Options = list(estimator = "DWLS"),
    Data = list(observed_names = spec$observed_names),
    Model = list(model_kind = model_kind, par_table = par_table)
  )
}

.fit_observed_covariance_model <- function(par_table, sample.cov, WLS.V, model_kind) {
  spec <- .observed_covariance_spec(par_table, sample.cov)
  fit <- fit_observed_covariance_dwls(
    as.double(spec$sample_cov),
    as.double(WLS.V),
    as.integer(spec$free_mask),
    as.double(spec$fixed_values),
    as.integer(nrow(spec$sample_cov))
  )

  .new_observed_covariance_fit(par_table, spec$sample_cov, WLS.V, fit, model_kind)
}

.parse_commonfactor_gwas_model <- function(model) {
  if (!is.character(model) || length(model) != 1L) {
    return(NULL)
  }

  lines <- trimws(strsplit(model, "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]
  loading_line <- lines[grepl("=~", lines, fixed = TRUE)]
  regression_line <- lines[grepl("~", lines, fixed = TRUE) & !grepl("=~", lines, fixed = TRUE)]

  if (length(loading_line) != 1L || length(regression_line) < 2L) {
    return(NULL)
  }

  loading_parts <- trimws(strsplit(loading_line, "=~", fixed = TRUE)[[1L]])
  if (length(loading_parts) != 2L) {
    return(NULL)
  }

  latent_name <- loading_parts[[1L]]
  traits <- trimws(strsplit(loading_parts[[2L]], "+", fixed = TRUE)[[1L]])
  traits <- traits[nzchar(traits)]
  factor_regression <- regression_line[
    grepl(paste0("^", latent_name, "[[:space:]]*~"), regression_line)
  ]

  if (length(factor_regression) != 1L || !length(traits)) {
    return(NULL)
  }

  predictor <- trimws(strsplit(factor_regression, "~", fixed = TRUE)[[1L]][[2L]])
  expected_direct <- paste0(traits, " ~ 0*", predictor)
  direct_regression <- regression_line[regression_line != factor_regression]

  if (!identical(sort(gsub("[[:space:]]+", "", direct_regression)), sort(gsub("[[:space:]]+", "", expected_direct)))) {
    return(NULL)
  }

  list(
    latent_name = latent_name,
    predictor = predictor,
    traits = traits
  )
}

.commonfactor_gwas_par_table <- function(spec, fit) {
  k <- length(spec$traits)
  n_rows <- 3L * k + 3L
  free <- c(
    0L,
    if (k > 1L) seq_len(k - 1L) else integer(),
    k,
    rep(0L, k),
    k + seq_len(k),
    2L * k + 1L,
    2L * k + 2L
  )
  est <- c(
    fit$loadings,
    fit$gamma,
    rep(0, k),
    fit$residuals,
    fit$psi,
    fit$phi
  )
  se <- c(
    0,
    if (k > 1L) fit$naive_se[seq_len(k - 1L)] else numeric(),
    fit$naive_se[[k]],
    rep(0, k),
    fit$naive_se[k + seq_len(k)],
    fit$naive_se[[2L * k + 1L]],
    fit$naive_se[[2L * k + 2L]]
  )

  data.frame(
    id = seq_len(n_rows),
    lhs = c(
      rep(spec$latent_name, k),
      spec$latent_name,
      spec$traits,
      spec$traits,
      spec$latent_name,
      spec$predictor
    ),
    op = c(
      rep("=~", k),
      "~",
      rep("~", k),
      rep("~~", k),
      "~~",
      "~~"
    ),
    rhs = c(
      spec$traits,
      spec$predictor,
      rep(spec$predictor, k),
      spec$traits,
      spec$latent_name,
      spec$predictor
    ),
    user = c(rep(1L, 2L * k + 1L), rep(0L, k + 2L)),
    block = rep(1L, n_rows),
    group = rep(1L, n_rows),
    free = as.integer(free),
    ustart = c(1, rep(NA_real_, k - 1L), NA_real_, rep(0, k), rep(NA_real_, k + 2L)),
    exo = rep(0L, n_rows),
    label = rep("", n_rows),
    plabel = paste0(".p", seq_len(n_rows), "."),
    start = est,
    est = est,
    se = se,
    stringsAsFactors = FALSE
  )
}

.new_commonfactor_gwas_fit <- function(model, sample.cov, WLS.V, fit) {
  spec <- .parse_commonfactor_gwas_model(model)
  observed_names <- c(spec$predictor, spec$traits)
  dimnames(sample.cov) <- list(observed_names, observed_names)
  delta <- matrix(fit$delta, nrow = length(.stat_names(observed_names)))
  dimnames(delta) <- list(
    .stat_names(observed_names),
    c(
      paste0(spec$latent_name, "=~", spec$traits[-1L]),
      paste0(spec$latent_name, "~", spec$predictor),
      paste0(spec$traits, "~~", spec$traits),
      paste0(spec$latent_name, "~~", spec$latent_name),
      paste0(spec$predictor, "~~", spec$predictor)
    )
  )
  implied <- matrix(fit$implied, nrow = nrow(sample.cov), ncol = ncol(sample.cov))
  dimnames(implied) <- dimnames(sample.cov)
  par_table <- .commonfactor_gwas_par_table(spec, fit)
  npar <- 2L * length(spec$traits) + 2L
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = length(.stat_names(observed_names)) - npar,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = sample.cov,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = matrix(1, nrow = 1L, ncol = 1L, dimnames = list(spec$latent_name, spec$latent_name)),
    se = list(theta = matrix(0, nrow = nrow(sample.cov), ncol = ncol(sample.cov))),
    converged = isTRUE(fit$converged),
    observed_names = observed_names,
    latent_name = spec$latent_name,
    model_kind = "commonfactor_gwas_dwls",
    Options = list(estimator = "DWLS"),
    Data = list(observed_names = observed_names),
    Model = list(model_kind = "commonfactor_gwas_dwls", model = model)
  )
}

.fit_commonfactor_gwas_model <- function(model, sample.cov, WLS.V) {
  spec <- .parse_commonfactor_gwas_model(model)
  observed_names <- c(spec$predictor, spec$traits)

  if (!identical(colnames(sample.cov), observed_names) && !identical(rownames(sample.cov), observed_names)) {
    stop("commonfactor GWAS rust slice requires sample.cov ordered as predictor followed by indicators.", call. = FALSE)
  }

  fit <- fit_commonfactor_gwas_dwls(
    as.double(sample.cov),
    as.double(WLS.V),
    as.integer(length(spec$traits)),
    400L,
    1e-12
  )

  .new_commonfactor_gwas_fit(model, sample.cov, WLS.V, fit)
}

.parse_user_gwas_model <- function(model) {
  if (!is.character(model) || length(model) != 1L) {
    return(NULL)
  }

  lines <- trimws(strsplit(model, "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]
  loading_line <- lines[grepl("=~", lines, fixed = TRUE)]
  regression_line <- lines[grepl("~", lines, fixed = TRUE) & !grepl("=~", lines, fixed = TRUE)]

  if (length(lines) != 2L || length(loading_line) != 1L || length(regression_line) != 1L) {
    return(NULL)
  }

  loading_parts <- trimws(strsplit(loading_line, "=~", fixed = TRUE)[[1L]])
  if (length(loading_parts) != 2L) {
    return(NULL)
  }

  latent_name <- loading_parts[[1L]]
  traits <- trimws(strsplit(loading_parts[[2L]], "+", fixed = TRUE)[[1L]])
  traits <- traits[nzchar(traits)]
  regression_parts <- trimws(strsplit(regression_line, "~", fixed = TRUE)[[1L]])

  if (
    length(regression_parts) != 2L ||
      !identical(regression_parts[[1L]], latent_name) ||
      !nzchar(regression_parts[[2L]]) ||
      length(traits) < 2L ||
      any(grepl("[*:=~><]", traits))
  ) {
    return(NULL)
  }

  list(
    latent_name = latent_name,
    predictor = regression_parts[[2L]],
    traits = traits
  )
}

.user_gwas_par_table <- function(spec, fit) {
  k <- length(spec$traits)
  n_rows <- 2L * k + 3L
  free <- c(
    0L,
    if (k > 1L) seq_len(k - 1L) else integer(),
    k,
    k + seq_len(k),
    2L * k + 1L,
    2L * k + 2L
  )
  est <- c(
    fit$loadings,
    fit$gamma,
    fit$residuals,
    fit$psi,
    fit$phi
  )
  se <- c(
    0,
    if (k > 1L) fit$naive_se[seq_len(k - 1L)] else numeric(),
    fit$naive_se[[k]],
    fit$naive_se[k + seq_len(k)],
    fit$naive_se[[2L * k + 1L]],
    fit$naive_se[[2L * k + 2L]]
  )

  data.frame(
    id = seq_len(n_rows),
    lhs = c(rep(spec$latent_name, k), spec$latent_name, spec$traits, spec$latent_name, spec$predictor),
    op = c(rep("=~", k), "~", rep("~~", k), "~~", "~~"),
    rhs = c(spec$traits, spec$predictor, spec$traits, spec$latent_name, spec$predictor),
    user = c(rep(1L, k + 1L), rep(0L, k + 2L)),
    block = rep(1L, n_rows),
    group = rep(1L, n_rows),
    free = as.integer(free),
    ustart = c(1, rep(NA_real_, k - 1L), rep(NA_real_, k + 3L)),
    exo = rep(0L, n_rows),
    label = rep("", n_rows),
    plabel = paste0(".p", seq_len(n_rows), "."),
    start = est,
    est = est,
    se = se,
    stringsAsFactors = FALSE
  )
}

.new_user_gwas_fit <- function(model, sample.cov, WLS.V, fit) {
  spec <- .parse_user_gwas_model(model)
  observed_names <- c(spec$predictor, spec$traits)
  dimnames(sample.cov) <- list(observed_names, observed_names)
  delta <- matrix(fit$delta, nrow = length(.stat_names(observed_names)))
  dimnames(delta) <- list(
    .stat_names(observed_names),
    c(
      paste0(spec$latent_name, "=~", spec$traits[-1L]),
      paste0(spec$latent_name, "~", spec$predictor),
      paste0(spec$traits, "~~", spec$traits),
      paste0(spec$latent_name, "~~", spec$latent_name),
      paste0(spec$predictor, "~~", spec$predictor)
    )
  )
  implied <- matrix(fit$implied, nrow = nrow(sample.cov), ncol = ncol(sample.cov))
  dimnames(implied) <- dimnames(sample.cov)
  par_table <- .user_gwas_par_table(spec, fit)
  npar <- 2L * length(spec$traits) + 2L
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = length(.stat_names(observed_names)) - npar,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = sample.cov,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = matrix(1, nrow = 1L, ncol = 1L, dimnames = list(spec$latent_name, spec$latent_name)),
    se = list(theta = matrix(0, nrow = nrow(sample.cov), ncol = ncol(sample.cov))),
    converged = isTRUE(fit$converged),
    observed_names = observed_names,
    latent_name = spec$latent_name,
    model_kind = "user_gwas_dwls",
    Options = list(estimator = "DWLS"),
    Data = list(observed_names = observed_names),
    Model = list(model_kind = "user_gwas_dwls", model = model)
  )
}

.fit_user_gwas_model <- function(model, sample.cov, WLS.V) {
  spec <- .parse_user_gwas_model(model)
  observed_names <- c(spec$predictor, spec$traits)

  if (!identical(colnames(sample.cov), observed_names) && !identical(rownames(sample.cov), observed_names)) {
    stop("userGWAS rust slice requires sample.cov ordered as predictor followed by indicators.", call. = FALSE)
  }

  fit <- fit_commonfactor_gwas_dwls(
    as.double(sample.cov),
    as.double(WLS.V),
    as.integer(length(spec$traits)),
    400L,
    1e-12
  )

  .new_user_gwas_fit(model, sample.cov, WLS.V, fit)
}

.par_value_rust <- function(par_table, row_idx) {
  if (!length(row_idx)) {
    return(0)
  }

  row_idx <- row_idx[[1L]]
  if (!is.na(par_table$ustart[[row_idx]])) {
    return(par_table$ustart[[row_idx]])
  }

  if (!is.na(par_table$est[[row_idx]])) {
    return(par_table$est[[row_idx]])
  }

  0
}

.commonfactor_gwas_q_spec <- function(par_table) {
  loading_rows <- which(par_table$op == "=~")
  latent_names <- unique(par_table$lhs[loading_rows])
  if (length(latent_names) != 1L) {
    return(NULL)
  }

  latent_name <- latent_names[[1L]]
  traits <- par_table$rhs[loading_rows]
  factor_regression <- which(
    par_table$op == "~" &
      par_table$lhs == latent_name
  )

  if (length(factor_regression) != 1L) {
    return(NULL)
  }

  predictor <- par_table$rhs[[factor_regression]]
  direct_rows <- vapply(
    traits,
    function(trait) {
      rows <- which(
        par_table$op == "~" &
          par_table$lhs == trait &
          par_table$rhs == predictor
      )
      if (length(rows) != 1L) {
        return(NA_integer_)
      }

      rows[[1L]]
    },
    integer(1L)
  )
  residual_rows <- vapply(
    traits,
    function(trait) {
      rows <- which(
        par_table$op == "~~" &
          par_table$lhs == trait &
          par_table$rhs == trait
      )
      if (length(rows) != 1L) {
        return(NA_integer_)
      }

      rows[[1L]]
    },
    integer(1L)
  )
  psi_row <- which(
    par_table$op == "~~" &
      par_table$lhs == latent_name &
      par_table$rhs == latent_name
  )
  phi_row <- which(
    par_table$op == "~~" &
      par_table$lhs == predictor &
      par_table$rhs == predictor
  )

  if (
    anyNA(direct_rows) ||
      anyNA(residual_rows) ||
      length(psi_row) != 1L ||
      length(phi_row) != 1L ||
      !identical(as.integer(par_table$free[direct_rows]), seq_along(traits)) ||
      !identical(as.integer(par_table$free[residual_rows]), length(traits) + seq_along(traits))
  ) {
    return(NULL)
  }

  list(
    latent_name = latent_name,
    predictor = predictor,
    traits = traits,
    loading_rows = loading_rows,
    factor_regression = factor_regression,
    direct_rows = direct_rows,
    residual_rows = residual_rows,
    psi_row = psi_row[[1L]],
    phi_row = phi_row[[1L]]
  )
}

.new_commonfactor_gwas_q_fit <- function(par_table, sample.cov, WLS.V, spec, fit) {
  observed_names <- c(spec$predictor, spec$traits)
  dimnames(sample.cov) <- list(observed_names, observed_names)
  implied <- matrix(fit$implied, nrow = nrow(sample.cov), ncol = ncol(sample.cov))
  dimnames(implied) <- dimnames(sample.cov)
  delta <- matrix(fit$delta, nrow = length(.stat_names(observed_names)))
  dimnames(delta) <- list(
    .stat_names(observed_names),
    c(
      paste0(spec$traits, "~", spec$predictor),
      paste0(spec$traits, "~~", spec$traits)
    )
  )
  par_table <- as.data.frame(par_table, stringsAsFactors = FALSE)
  par_table$est[spec$direct_rows] <- fit$direct
  par_table$est[spec$residual_rows] <- fit$residuals
  par_table$se[] <- 0
  par_table$se[c(spec$direct_rows, spec$residual_rows)] <- fit$naive_se
  npar <- 2L * length(spec$traits)
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = length(.stat_names(observed_names)) - npar,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = sample.cov,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = matrix(1, nrow = 1L, ncol = 1L, dimnames = list(spec$latent_name, spec$latent_name)),
    se = list(theta = matrix(0, nrow = nrow(sample.cov), ncol = ncol(sample.cov))),
    converged = isTRUE(fit$converged),
    observed_names = observed_names,
    latent_name = spec$latent_name,
    model_kind = "commonfactor_gwas_q_dwls",
    Options = list(estimator = "DWLS"),
    Data = list(observed_names = observed_names),
    Model = list(model_kind = "commonfactor_gwas_q_dwls", par_table = par_table)
  )
}

.fit_commonfactor_gwas_q_model <- function(par_table, sample.cov, WLS.V, spec) {
  loadings <- vapply(par_table$est[spec$loading_rows], as.numeric, numeric(1L))
  gamma <- .par_value_rust(par_table, spec$factor_regression)
  direct <- vapply(spec$direct_rows, function(row_idx) .par_value_rust(par_table, row_idx), numeric(1L))
  residuals <- vapply(spec$residual_rows, function(row_idx) .par_value_rust(par_table, row_idx), numeric(1L))
  psi <- .par_value_rust(par_table, spec$psi_row)
  phi <- .par_value_rust(par_table, spec$phi_row)
  fit <- fit_commonfactor_gwas_q_dwls(
    as.double(sample.cov),
    as.double(WLS.V),
    as.double(loadings),
    as.double(gamma),
    as.double(direct),
    as.double(residuals),
    as.double(psi),
    as.double(phi),
    as.integer(length(spec$traits)),
    400L,
    1e-12
  )

  .new_commonfactor_gwas_q_fit(par_table, sample.cov, WLS.V, spec, fit)
}

.user_gwas_fixed_measurement_spec <- function(par_table) {
  if (!is.data.frame(par_table)) {
    return(NULL)
  }

  loading_rows <- which(par_table$op == "=~")
  latent_names <- unique(par_table$lhs[loading_rows])
  if (length(latent_names) != 1L) {
    return(NULL)
  }

  latent_name <- latent_names[[1L]]
  traits <- par_table$rhs[loading_rows]
  factor_regression <- which(
    par_table$op == "~" &
      par_table$lhs == latent_name
  )
  residual_rows <- vapply(
    traits,
    function(trait) {
      rows <- which(
        par_table$op == "~~" &
          par_table$lhs == trait &
          par_table$rhs == trait
      )
      if (length(rows) != 1L) {
        return(NA_integer_)
      }

      rows[[1L]]
    },
    integer(1L)
  )
  psi_row <- which(
    par_table$op == "~~" &
      par_table$lhs == latent_name &
      par_table$rhs == latent_name
  )

  if (length(factor_regression) != 1L || anyNA(residual_rows) || length(psi_row) != 1L) {
    return(NULL)
  }

  predictor <- par_table$rhs[[factor_regression]]
  phi_row <- which(
    par_table$op == "~~" &
      par_table$lhs == predictor &
      par_table$rhs == predictor
  )
  direct_rows <- which(
    par_table$op == "~" &
      par_table$lhs %in% traits &
      par_table$rhs == predictor
  )
  expected_rows <- unname(sort(c(loading_rows, residual_rows, psi_row, factor_regression, phi_row)))

  if (
    length(phi_row) != 1L ||
      length(direct_rows) > 0L ||
      !identical(expected_rows, seq_len(nrow(par_table))) ||
      any(par_table$free[loading_rows] != 0L) ||
      any(par_table$free[c(residual_rows, psi_row, factor_regression, phi_row)] <= 0L)
  ) {
    return(NULL)
  }

  list(
    latent_name = latent_name,
    predictor = predictor,
    traits = traits,
    loading_rows = loading_rows,
    residual_rows = residual_rows,
    psi_row = psi_row[[1L]],
    factor_regression = factor_regression[[1L]],
    phi_row = phi_row[[1L]]
  )
}

.new_user_gwas_fixed_measurement_fit <- function(par_table, sample.cov, WLS.V, spec, fit) {
  observed_names <- c(spec$predictor, spec$traits)
  dimnames(sample.cov) <- list(observed_names, observed_names)
  implied <- matrix(fit$implied, nrow = nrow(sample.cov), ncol = ncol(sample.cov))
  dimnames(implied) <- dimnames(sample.cov)
  delta <- matrix(fit$delta, nrow = length(.stat_names(observed_names)))
  dimnames(delta) <- list(
    .stat_names(observed_names),
    c(
      paste0(spec$traits, "~~", spec$traits),
      paste0(spec$latent_name, "~~", spec$latent_name),
      paste0(spec$latent_name, "~", spec$predictor),
      paste0(spec$predictor, "~~", spec$predictor)
    )
  )
  par_table <- as.data.frame(par_table, stringsAsFactors = FALSE)
  rownames(par_table) <- NULL
  free_rows <- c(spec$residual_rows, spec$psi_row, spec$factor_regression, spec$phi_row)
  par_table$free[] <- 0L
  par_table$free[free_rows] <- seq_along(free_rows)
  par_table$est[spec$residual_rows] <- fit$residuals
  par_table$est[spec$psi_row] <- fit$psi
  par_table$est[spec$factor_regression] <- fit$gamma
  par_table$est[spec$phi_row] <- fit$phi
  par_table$se[] <- 0
  par_table$se[free_rows] <- fit$naive_se
  npar <- length(free_rows)
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = length(.stat_names(observed_names)) - npar,
    srmr = fit$srmr
  )

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = sample.cov,
    implied = list(cov = implied),
    delta = delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = matrix(1, nrow = 1L, ncol = 1L, dimnames = list(spec$latent_name, spec$latent_name)),
    se = list(theta = matrix(0, nrow = nrow(sample.cov), ncol = ncol(sample.cov))),
    converged = isTRUE(fit$converged),
    observed_names = observed_names,
    latent_name = spec$latent_name,
    model_kind = "user_gwas_fixed_measurement_dwls",
    Options = list(estimator = "DWLS"),
    Data = list(observed_names = observed_names),
    Model = list(model_kind = "user_gwas_fixed_measurement_dwls", par_table = par_table)
  )
}

.fit_user_gwas_fixed_measurement_model <- function(par_table, sample.cov, WLS.V, spec) {
  observed_names <- c(spec$predictor, spec$traits)

  if (!identical(colnames(sample.cov), observed_names) && !identical(rownames(sample.cov), observed_names)) {
    stop("userGWAS fixed-measurement rust slice requires sample.cov ordered as predictor followed by indicators.", call. = FALSE)
  }

  fit <- fit_user_gwas_fixed_measurement_dwls(
    as.double(sample.cov),
    as.double(WLS.V),
    as.double(par_table$est[spec$loading_rows]),
    as.double(par_table$est[spec$residual_rows]),
    as.double(.par_value_rust(par_table, spec$psi_row)),
    as.double(.par_value_rust(par_table, spec$factor_regression)),
    as.double(.par_value_rust(par_table, spec$phi_row)),
    as.integer(length(spec$traits)),
    400L,
    1e-12
  )

  .new_user_gwas_fixed_measurement_fit(par_table, sample.cov, WLS.V, spec, fit)
}

.lavaan_fast_latent_correlation <- function(compiled, free_values) {
  if (!length(compiled$latent_names)) {
    return(matrix(numeric(), nrow = 0L, ncol = 0L))
  }

  matrices <- .lavaan_fast_ram_matrices(compiled, free_values)
  inverse <- solve(diag(nrow(matrices$A)) - matrices$A)
  full_implied <- inverse %*% matrices$S %*% t(inverse)
  latent_idx <- match(compiled$latent_names, compiled$variable_names)
  latent_cov <- full_implied[latent_idx, latent_idx, drop = FALSE]
  latent_sd <- sqrt(diag(latent_cov))
  latent_cor <- latent_cov / outer(latent_sd, latent_sd)
  dimnames(latent_cor) <- list(compiled$latent_names, compiled$latent_names)
  latent_cor
}

.new_lavaan_fast_ram_fit <- function(par_table, sample.cov, WLS.V, compiled, plan, fit) {
  par_table <- as.data.frame(par_table, stringsAsFactors = FALSE)
  rownames(par_table) <- NULL
  free_rows <- which(!is.na(compiled$free_index))
  full_free_rows <- compiled$structural_row_indices[free_rows]
  par_table$est[full_free_rows] <- fit$estimates[compiled$free_index[free_rows]]
  par_table$se[] <- 0
  par_table$se[full_free_rows] <- fit$naive_se[compiled$free_index[free_rows]]
  definition_plan <- .lavaan_fast_definition_plan(par_table, compiled)
  if (length(definition_plan$rows)) {
    par_table$est[definition_plan$rows] <- Re(definition_plan$def.function(fit$estimates))
  }
  npar <- length(compiled$free_ids)
  fit_stats <- c(
    npar = npar,
    fmin = fit$objective,
    chisq = 2 * fit$objective,
    df = length(.stat_names(compiled$observed_names)) - npar,
    srmr = fit$srmr
  )
  cor_lv <- .lavaan_fast_latent_correlation(compiled, fit$estimates)

  methods::new(
    "lavaan_rust_fit",
    ParTable = par_table,
    observed = sample.cov,
    implied = list(cov = fit$implied),
    delta = fit$delta,
    WLS.V = WLS.V,
    fit = fit_stats,
    cor.lv = cor_lv,
    se = list(theta = matrix(0, nrow = nrow(sample.cov), ncol = ncol(sample.cov))),
    converged = isTRUE(fit$converged),
    observed_names = compiled$observed_names,
    latent_name = if (length(compiled$latent_names)) compiled$latent_names[[1L]] else "",
    model_kind = "ram_dwls_generic",
    Options = list(estimator = "DWLS"),
    Data = list(observed_names = compiled$observed_names),
    Model = methods::new(
      "lavaan_rust_model",
      model_kind = "ram_dwls_generic",
      par_table = par_table,
      compiled = compiled,
      plan = plan,
      free_values = as.double(fit$estimates),
      def.function = definition_plan$def.function
    )
  )
}

.lavaan_fast_assert_compatible_data <- function(compiled, sample.cov, WLS.V) {
  observed_names <- colnames(sample.cov)
  if (is.null(observed_names)) {
    observed_names <- rownames(sample.cov)
  }

  if (
    !is.matrix(sample.cov) ||
      !identical(dim(sample.cov), c(compiled$n_observed, compiled$n_observed)) ||
      !identical(observed_names, compiled$observed_names)
  ) {
    stop("lavaan_fast native plan received incompatible observed covariance data.", call. = FALSE)
  }

  if (
    !is.matrix(WLS.V) ||
      !identical(dim(WLS.V), c(compiled$n_stats, compiled$n_stats))
  ) {
    stop("lavaan_fast native plan received an incompatible DWLS weight matrix.", call. = FALSE)
  }
}

.fit_lavaan_fast_compiled_ram_model <- function(
  par_table,
  compiled,
  plan,
  sample.cov,
  WLS.V,
  free_values = compiled$default_free_values,
  compute_se = TRUE
) {
  if (!.lavaan_fast_compiled_model_is_current(par_table, compiled)) {
    stop(
      paste(
        "lavaan_fast cached model structure changed after compilation;",
        "refit with sem_rust() instead of mutating cached model internals."
      ),
      call. = FALSE
    )
  }

  .lavaan_fast_assert_compatible_data(compiled, sample.cov, WLS.V)
  plan <- .lavaan_fast_native_plan(compiled, plan)
  fit <- .lavaan_fast_fit_dwls_plan_rust(
    compiled,
    plan,
    sample.cov,
    WLS.V,
    free_values = free_values,
    compute_se = compute_se
  )
  .new_lavaan_fast_ram_fit(par_table, sample.cov, WLS.V, compiled, plan, fit)
}

.fit_lavaan_fast_ram_model <- function(par_table, sample.cov, WLS.V, compute_se = TRUE) {
  par_table <- .lavaan_fast_normalize_free_ids(par_table)
  observed_names <- colnames(sample.cov)
  if (is.null(observed_names)) {
    observed_names <- rownames(sample.cov)
  }
  if (is.null(observed_names)) {
    stop("lavaan_fast RAM models require named observed variables.", call. = FALSE)
  }

  compiled <- .lavaan_fast_compile_par_table(par_table, observed_names)
  .fit_lavaan_fast_compiled_ram_model(
    par_table = par_table,
    compiled = compiled,
    plan = NULL,
    sample.cov = sample.cov,
    WLS.V = WLS.V,
    compute_se = compute_se
  )
}

#' Rust-backed experimental replacement for the supported `lavaan::sem()` slice.
#'
#' @param model A supported lavaan-style model string or parameter table.
#' @param sample.cov Observed covariance matrix.
#' @param estimator Estimator name. Only `"DWLS"` is currently supported.
#' @param WLS.V DWLS weight matrix.
#' @param ... Additional lavaan-style arguments. Present for signature
#'   compatibility; unsupported paths still error.
#' @export
sem_rust <- function(model, sample.cov, estimator = "ML", WLS.V = NULL, ...) {
  dots <- list(...)
  std.lv <- isTRUE(dots$std.lv)
  compute_se <- !identical(dots$se, "none")

  if (
    identical(estimator, "DWLS") &&
      .is_one_factor_dwls_model(model) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    fit <- fit_one_factor_dwls(
      as.double(sample.cov),
      as.double(WLS.V),
      as.integer(nrow(sample.cov)),
      200L,
      1e-12
    )
    return(.new_one_factor_fit(model, sample.cov, WLS.V, fit))
  }

  if (
    identical(estimator, "DWLS") &&
      !is.null(.parse_simple_one_factor_model(model)) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    fit <- fit_one_factor_dwls(
      as.double(sample.cov),
      as.double(WLS.V),
      as.integer(nrow(sample.cov)),
      200L,
      1e-12
    )

    if (std.lv) {
      return(.new_simple_std_lv_one_factor_fit(model, sample.cov, WLS.V, fit))
    }

    return(.new_marker_one_factor_fit(model, sample.cov, WLS.V, fit))
  }

  if (
    identical(estimator, "DWLS") &&
      !is.null(.parse_commonfactor_gwas_model(model)) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    return(.fit_commonfactor_gwas_model(model, sample.cov, WLS.V))
  }

  if (
    identical(estimator, "DWLS") &&
      !std.lv &&
      !is.null(.parse_user_gwas_model(model)) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    return(.fit_user_gwas_model(model, sample.cov, WLS.V))
  }

  if (
    identical(estimator, "DWLS") &&
      .is_commonfactor_null_model(model) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    observed_names <- colnames(sample.cov)
    if (is.null(observed_names)) {
      observed_names <- rownames(sample.cov)
    }
    if (is.null(observed_names)) {
      observed_names <- paste0("V", seq_len(ncol(sample.cov)))
    }

    return(.fit_observed_covariance_model(
      .commonfactor_null_par_table(observed_names),
      sample.cov,
      WLS.V,
      "commonfactor_null_dwls"
    ))
  }

  if (
    identical(estimator, "DWLS") &&
      is.character(model) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    parsed_model <- .lavaan_fast_parse_model_string(model, sample.cov, std.lv = std.lv)
    if (!is.null(parsed_model)) {
      return(.fit_lavaan_fast_ram_model(parsed_model, sample.cov, WLS.V, compute_se = compute_se))
    }
  }

  if (
    identical(estimator, "DWLS") &&
      is.data.frame(model) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    user_gwas_spec <- .user_gwas_fixed_measurement_spec(model)
    if (!is.null(user_gwas_spec)) {
      return(.fit_user_gwas_fixed_measurement_model(model, sample.cov, WLS.V, user_gwas_spec))
    }

    q_spec <- .commonfactor_gwas_q_spec(model)
    if (!is.null(q_spec)) {
      return(.fit_commonfactor_gwas_q_model(model, sample.cov, WLS.V, q_spec))
    }

    if (any(model$op %in% c("=~", "~") & model$free > 0L)) {
      return(.fit_lavaan_fast_ram_model(model, sample.cov, WLS.V, compute_se = compute_se))
    }

    return(.fit_observed_covariance_model(
      model,
      sample.cov,
      WLS.V,
      "observed_covariance_par_table_dwls"
    ))
  }

  stop(
    "Unsupported sem_rust() model path. The rust wrapper does not fall back to lavaan.",
    call. = FALSE
  )
}

#' Rust-backed experimental replacement for `lavaan::lavaan()`.
#'
#' @param sample.cov Observed covariance matrix.
#' @param WLS.V DWLS weight matrix.
#' @param slotOptions Compatibility payload from a prior rust fit.
#' @param slotParTable Compatibility parameter table from a prior rust fit.
#' @param slotData Compatibility payload from a prior rust fit.
#' @param slotModel Compatibility model payload from a prior rust fit.
#' @param ... Additional lavaan-style arguments. Unsupported paths still error.
#' @export
lavaan_rust <- function(sample.cov, WLS.V = NULL, slotOptions = NULL, slotParTable = NULL, slotData = NULL, slotModel = NULL, ...) {
  if (
    is.list(slotModel) &&
      identical(slotModel$model_kind, "commonfactor_gwas_dwls") &&
      is.character(slotModel$model) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    return(.fit_commonfactor_gwas_model(slotModel$model, sample.cov, WLS.V))
  }

  if (
    is.list(slotModel) &&
      identical(slotModel$model_kind, "user_gwas_dwls") &&
      is.character(slotModel$model) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    return(.fit_user_gwas_model(slotModel$model, sample.cov, WLS.V))
  }

  if (
    is.list(slotModel) &&
      identical(slotModel$model_kind, "user_gwas_fixed_measurement_dwls") &&
      is.data.frame(slotModel$par_table) &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    spec <- .user_gwas_fixed_measurement_spec(slotModel$par_table)
    if (!is.null(spec)) {
      return(.fit_user_gwas_fixed_measurement_model(slotModel$par_table, sample.cov, WLS.V, spec))
    }
  }

  if (
    methods::is(slotModel, "lavaan_rust_model") &&
      identical(slotModel@model_kind, "ram_dwls_generic") &&
      is.data.frame(slotModel@par_table) &&
      inherits(slotModel@compiled, "lavaan_fast_compiled") &&
      is.matrix(sample.cov) &&
      is.matrix(WLS.V)
  ) {
    return(.fit_lavaan_fast_compiled_ram_model(
      par_table = slotModel@par_table,
      compiled = slotModel@compiled,
      plan = slotModel@plan,
      sample.cov = sample.cov,
      WLS.V = WLS.V,
      free_values = slotModel@free_values
    ))
  }

  stop(
    "Unsupported lavaan_rust() model-reuse path. The rust wrapper does not fall back to lavaan.",
    call. = FALSE
  )
}

#' Rust-backed experimental replacement for `lavaan::lavInspect()`.
#'
#' @param object A `lavaan_rust_fit` object.
#' @param what Inspection key.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
lavInspect_rust <- function(object, what, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("lavInspect_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
  }

  switch(
    what,
    "delta" = object@delta,
    "WLS.V" = object@WLS.V,
    "converged" = object@converged,
    "cor.lv" = object@cor.lv,
    "fit" = object@fit,
    stop("Unsupported lavaan_rust inspection key: ", what, call. = FALSE)
  )
}

#' Rust-backed experimental replacement for `lavaan::inspect()`.
#'
#' @param object A `lavaan_rust_fit` object.
#' @param what Inspection key.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
inspect_rust <- function(object, what, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("inspect_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
  }

  if (missing(what)) {
    return(list(object@implied$cov))
  }

  switch(
    what,
    "list" = object@ParTable,
    "se" = object@se,
    stop("Unsupported lavaan_rust inspection key: ", what, call. = FALSE)
  )
}

#' Rust-backed experimental replacement for `lavaan::parTable()`.
#'
#' @param object A `lavaan_rust_fit` object.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
parTable_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("parTable_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
  }

  par_table <- object@ParTable
  if (!"lower" %in% names(par_table)) {
    insert_after <- match("label", names(par_table))
    par_table <- cbind(
      par_table[seq_len(insert_after)],
      lower = rep(NA_real_, nrow(par_table)),
      par_table[-seq_len(insert_after)]
    )
  }
  if (!"upper" %in% names(par_table)) {
    insert_after <- match("lower", names(par_table))
    par_table <- cbind(
      par_table[seq_len(insert_after)],
      upper = rep(NA_real_, nrow(par_table)),
      par_table[-seq_len(insert_after)]
    )
  }

  par_table
}

#' Rust-backed experimental replacement for `stats::fitted()`.
#'
#' @param object A `lavaan_rust_fit` object.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
fitted_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("fitted_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
  }

  object@implied
}

#' Rust-backed experimental replacement for `stats::resid()`.
#'
#' @param object A `lavaan_rust_fit` object.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
resid_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("resid_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
  }

  list(cov = object@observed - object@implied$cov)
}

#' Rust-backed experimental replacement for `lavaan::standardizedSolution()`.
#'
#' @param object A `lavaan_rust_fit` object.
#' @param ... Reserved for lavaan-compatible arguments.
#' @export
standardizedSolution_rust <- function(object, ...) {
  if (!methods::is(object, "lavaan_rust_fit")) {
    stop("standardizedSolution_rust() only supports lavaan_rust_fit objects.", call. = FALSE)
  }

  if (identical(object@model_kind, "ram_dwls_generic")) {
    par_table <- object@ParTable
    compiled <- .lavaan_fast_compile_par_table(par_table, object@observed_names)
    free_values <- lav_model_get_parameters_rust(object@Model, type = "free")
    matrices <- .lavaan_fast_ram_matrices(compiled, free_values)
    inverse <- solve(diag(nrow(matrices$A)) - matrices$A)
    full_implied <- inverse %*% matrices$S %*% t(inverse)
    variable_sd <- sqrt(diag(full_implied))
    names(variable_sd) <- compiled$variable_names

    est.std <- par_table$est
    structural_rows <- which(par_table$op %in% c("=~", "~", "~~"))
    for (row_idx in structural_rows) {
      lhs <- par_table$lhs[[row_idx]]
      rhs <- par_table$rhs[[row_idx]]
      op <- par_table$op[[row_idx]]

      est.std[[row_idx]] <- if (identical(op, "=~")) {
        par_table$est[[row_idx]] * variable_sd[[lhs]] / variable_sd[[rhs]]
      } else if (identical(op, "~")) {
        par_table$est[[row_idx]] * variable_sd[[rhs]] / variable_sd[[lhs]]
      } else {
        par_table$est[[row_idx]] / (variable_sd[[lhs]] * variable_sd[[rhs]])
      }
    }

    label_rows <- structural_rows[nzchar(par_table$label[structural_rows])]
    label_values <- setNames(est.std[label_rows], par_table$label[label_rows])
    label_values <- label_values[!duplicated(names(label_values))]
    defined_rows <- which(par_table$op == ":=")
    if (length(defined_rows)) {
      est.std[defined_rows] <- Re(.lavaan_fast_evaluate_defined_from_labels(par_table, label_values))
    }

    latent_variance_rows <- which(
      par_table$op == "~~" &
        par_table$lhs %in% compiled$latent_names &
        par_table$lhs == par_table$rhs
    )
    pvalue <- rep(0, nrow(par_table))
    pvalue[latent_variance_rows] <- NA_real_

    return(data.frame(
      lhs = par_table$lhs,
      op = par_table$op,
      rhs = par_table$rhs,
      est.std = est.std,
      se = rep(0, nrow(par_table)),
      z = rep(NA_real_, nrow(par_table)),
      pvalue = pvalue,
      ci.lower = est.std,
      ci.upper = est.std,
      stringsAsFactors = FALSE
    ))
  }

  if (!object@model_kind %in% c("marker_one_factor_dwls", "one_factor_dwls", "simple_std_lv_one_factor_dwls")) {
    stop("Unsupported standardizedSolution_rust() model path.", call. = FALSE)
  }

  par_table <- object@ParTable
  observed_sd <- sqrt(diag(object@implied$cov))
  names(observed_sd) <- object@observed_names
  latent_var <- if (object@model_kind == "marker_one_factor_dwls") {
    par_table$est[par_table$lhs == object@latent_name & par_table$op == "~~" & par_table$rhs == object@latent_name][[1L]]
  } else {
    1
  }
  est.std <- par_table$est

  loading_rows <- which(par_table$op == "=~")
  est.std[loading_rows] <- par_table$est[loading_rows] * sqrt(latent_var) / observed_sd[par_table$rhs[loading_rows]]

  residual_rows <- which(par_table$op == "~~" & par_table$lhs %in% object@observed_names & par_table$lhs == par_table$rhs)
  est.std[residual_rows] <- par_table$est[residual_rows] / diag(object@implied$cov)[match(par_table$lhs[residual_rows], object@observed_names)]

  latent_rows <- which(par_table$op == "~~" & par_table$lhs == object@latent_name & par_table$rhs == object@latent_name)
  est.std[latent_rows] <- 1

  pvalue <- rep(0, nrow(par_table))
  pvalue[latent_rows] <- NA_real_

  data.frame(
    lhs = par_table$lhs,
    op = par_table$op,
    rhs = par_table$rhs,
    est.std = est.std,
    se = rep(0, nrow(par_table)),
    z = rep(NA_real_, nrow(par_table)),
    pvalue = pvalue,
    ci.lower = est.std,
    ci.upper = est.std,
    stringsAsFactors = FALSE
  )
}

#' Rust-backed experimental replacement for lavaan's parameter extractor.
#'
#' @param object A `lavaan_rust_model` object.
#' @param type Parameter type. Only `"free"` is currently supported.
#' @export
lav_model_get_parameters_rust <- function(object, type = "free", ...) {
  if (!methods::is(object, "lavaan_rust_model")) {
    stop("lav_model_get_parameters_rust() only supports lavaan_rust_model objects.", call. = FALSE)
  }

  if (!identical(type, "free")) {
    stop("lav_model_get_parameters_rust() currently supports type = \"free\" only.", call. = FALSE)
  }

  object@free_values
}

#' Rust-backed experimental replacement for lavaan's complex-step Jacobian.
#'
#' @param func Defined-parameter function.
#' @param x Free-parameter vector.
#' @export
lav_func_jacobian_complex_rust <- function(func, x, ...) {
  if (!is.function(func) || !is.numeric(x)) {
    stop("lav_func_jacobian_complex_rust() expects a function and numeric vector.", call. = FALSE)
  }

  baseline <- func(x)
  jacobian <- matrix(
    0,
    nrow = length(baseline),
    ncol = length(x),
    dimnames = list(names(baseline), names(x))
  )
  step <- 1e-20
  for (idx in seq_along(x)) {
    perturbed <- as.complex(x)
    perturbed[[idx]] <- perturbed[[idx]] + 1i * step
    jacobian[, idx] <- Im(func(perturbed)) / step
  }

  jacobian
}
