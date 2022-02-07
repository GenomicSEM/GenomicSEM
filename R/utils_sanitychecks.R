
# Sanity check functions
.check_equal_length <- function(left, right) {
  name_left <- deparse(substitute(left))
  name_right <- deparse(substitute(right))
  if (!(is.null(left)) & !(is.null(right))) {
    if (length(left) != length(right)) {
      stop(paste("Length of", name_left,"and",name_right,"should be equal"), call.=FALSE)
    }
  }
}

.check_range <- function(val, min=0, max=1, allowNA=FALSE, inclusive=FALSE, name=NULL) {
  if (is.null(name)) name <- deparse(substitute(val))
  if (length(val) > 1) {
    for (x in val) {.check_range(x, min=min, max=max, allowNA=allowNA, inclusive=inclusive, name=name)}
  } else {
    if (is.na(val)) {
      if (!(allowNA)) stop(paste("Value(s) of", name,"should not be NA"), call.=FALSE)
    } else {
      if (!(is.numeric(val))) stop(paste("Value(s) of", name,"should be numeric"), call.=FALSE)
      if (inclusive) {
        if (val <= min) {
          stop(paste("Value(s) of", name,"should be above",min), call.=FALSE)
        } else if (val >= max) {
          stop(paste("Value(s) of", name,"should be below",max), call.=FALSE)
        }
      } else {
        if (val < min) {
          stop(paste("Value(s) of", name,"should be above",min), call.=FALSE)
        } else if (val > max) {
          stop(paste("Value(s) of", name,"should be below",max), call.=FALSE)
        }
      }
    }
  }
}

.check_file_exists <- function(path, name=NULL) {
  if (is.null(name)) name <- deparse(substitute(path))
  if (length(path) > 1) {
    for (x in path) {.check_file_exists(x, name=name)}
  } else {
    if (!(file.exists(path))) {
      stop(paste("File", path, "passed to", name,"does not exist"), call.=FALSE)
    }
  }
}

.check_boolean <- function(val, name=NULL) {
  if (is.null(name)) name <- deparse(substitute(val))
  if (length(val) > 1) {
    for (x in val) {.check_boolean(x, name=name)}
  } else {
    if (!(val %in% c(TRUE, FALSE))) {
      stop(paste(name, "should be TRUE or FALSE"))
    }
  }
}

.check_one_of <- function(val, one_of) {
  name <- deparse(substitute(val))
  if (sum(one_of == val) != 1) {
    stop(paste(name, "should be one of the following values:\n",paste(one_of, collapse=", ")))
  }
}

