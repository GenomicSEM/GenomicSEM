indexS <- function(LDSC_OBJECT = NULL, MATRIX = NULL, R = F){
  
  ##put warning here that pops if they have both ldsc object and smatrix
  if (is.list(LDSC_OBJECT)) {
    S <- LDSC_OBJECT$S
  } else {
    S <- MATRIX
  }  
  
  if (R) {
    k <- nrow(S) # number of variables
    Snum <- matrix(0, k, k) 
    Snum[lower.tri(Snum, diag = F)] <- 1:(k * (k - 1) / 2)
    Snum[upper.tri(Snum,diag=FALSE)]=(t(Snum))[upper.tri(Snum,diag=FALSE)]
    Snum
  } else {
    k <- nrow(S)
    Snum<- matrix(0, k, k)
    Snum[lower.tri(Snum,diag=TRUE)] <-1:(k*(k+1)/2)
    Snum[upper.tri(Snum,diag=FALSE)]=(t(Snum))[upper.tri(Snum,diag=FALSE)]
    Snum
  }  