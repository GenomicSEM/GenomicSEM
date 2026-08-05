#' Create an Index Matrix for an LDSC Sampling Covariance Matrix
#'
#' Creates the matrix of index values used to map entries in an LDSC genetic
#' covariance or correlation matrix to the corresponding rows and columns of
#' its sampling covariance matrix.
#'
#' @param LDSC_OBJECT An optional LDSC output object containing an \code{S}
#'   matrix. Supply either \code{LDSC_OBJECT} or \code{MATRIX}.
#' @param MATRIX An optional square numeric matrix. This is typically an LDSC
#'   genetic covariance matrix. Supply either \code{MATRIX} or
#'   \code{LDSC_OBJECT}.
#' @param R Logical. If \code{FALSE} (the default), indices are created for a
#'   covariance matrix and include the diagonal. If \code{TRUE}, indices are
#'   created for a directly estimated genetic correlation matrix and exclude
#'   the diagonal.
#'
#' @details
#' The returned matrix is symmetric. Selecting a row, column, or diagonal block
#' identifies the positions needed to extract matching parameter estimates and
#' sampling covariances with \code{\link{subSV}}.
#'
#' @return A symmetric numeric matrix of index values with the same dimensions
#'   as the supplied matrix.
#'
#' @seealso \code{\link{subSV}}
#' @export


indexS <- function(LDSC_OBJECT = NULL, MATRIX = NULL, R = F){
  if (is.list(LDSC_OBJECT)) {
    S <- LDSC_OBJECT$S
  } else {
    S <- MATRIX
  }

  if (R) {
    k <- nrow(S)
    Snum <- matrix(0, k, k)
    Snum[lower.tri(Snum, diag = F)] <- 1:(k * (k - 1) / 2)
    Snum[upper.tri(Snum, diag = FALSE)] =
      (t(Snum))[upper.tri(Snum, diag = FALSE)]
    Snum
  } else {
    k <- nrow(S)
    Snum <- matrix(0, k, k)
    Snum[lower.tri(Snum, diag = TRUE)] <- 1:(k * (k + 1) / 2)
    Snum[upper.tri(Snum, diag = FALSE)] =
      (t(Snum))[upper.tri(Snum, diag = FALSE)]
    Snum
  }
}
