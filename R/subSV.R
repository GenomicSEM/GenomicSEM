subSV <- function(LDSC_OBJECT = NULL, SMATRIX = NULL, VMATRIX = NULL, INDEXVALS, TYPE = "S"){
  
  #WARNINGS:
  #Checks for either an LDSCobject or and S and V matrix
  if (!is.null(LDSC_OBJECT) & (!is.null(SMATRIX) | !is.null(VMATRIX))) {
    stop("You must include either an LDSC object OR an S and V matrix, not both.")
  }
  #If LDSCOBJECT is not provided, check for the presence of SMATRIX and VMATRIX
  if (is.null(LDSC_OBJECT)) {
    #If VMATRIX is provided but SMATRIX is not
    if (!is.null(VMATRIX) & is.null(RMATRIX)) {
      warning("You forgot to enter the SMATRIX!")
    }
    
    # If SMATRIX is provided but VMATRIX is not
    if (!is.null(RMATRIX) & is.null(VMATRIX)) {
      warning("You forgot to enter the VMATRIX!")
    }
  }
  
  if (!is.null(LDSC_OBJECT)) {
    if (TYPE == "S") {
      S <- LDSC_OBJECT$S
      V <- LDSC_OBJECT$V
      
    } else if (TYPE == "S_Stand") {
      S <- LDSC_OBJECT$S_Stand
      V <- LDSC_OBJECT$V_Stand
      
    } else if (TYPE == "R") {
      S <- LDSC_OBJECT$R
      V <- LDSC_OBJECT$V_R
      
    } else {
      # Code to execute if TYPE does not match any of the above
      print("Unknown TYPE")
    }
  } else if (!is.null(SMATRIX) | !is.null(VMATRIX)) {
    S <- SMATRIX
    V <- VMATRIX
  } else {
    return("Must submit either an ldsc object or S/V matrix")
  }
  
  
  ###Creating index values related to rGs of interest
  if (TYPE == "R") {
    k <- nrow(S) # number of variables
    Snum <- matrix(0, k, k) 
    Snum[lower.tri(Snum, diag = F)] <- 1:(k * (k - 1) / 2)
    Snum[upper.tri(Snum,diag=FALSE)]=(t(Snum))[upper.tri(Snum,diag=FALSE)]
    Snum
    ###vectorize lower triangle both Matrices (needs to be lower tri so that the locations and values aren't             repeated)
    S_vec <-  S[lower.tri(S,diag=TRUE)]
    Snum_vec <- Snum[lower.tri(Snum,diag=TRUE)] 
    subS <- S_vec[Snum_vec %in% INDEXVALS]
    subV <- V[INDEXVALS,INDEXVALS]
  } else {
    k <- nrow(S)
    Snum<-matrix(0, k, k)
    Snum[lower.tri(Snum,diag=TRUE)] <-1:(k*(k+1)/2)
    Snum[upper.tri(Snum,diag=FALSE)]=(t(Snum))[upper.tri(Snum,diag=FALSE)]
    Snum
    
    ###vectorize lower triangle both Matrices (needs to be lower tri so that the locations and values aren't repeated)
    S_vec <-  S[lower.tri(S,diag=TRUE)]
    Snum_vec <- Snum[lower.tri(Snum,diag=TRUE)]                      
    subS <- S_vec[Snum_vec %in% INDEXVALS]
    subV <- V[INDEXVALS,INDEXVALS]
  }
  
  ###Turning the subsetted V and S matrices into a list and creating a global variable
  SUBSET_LDSC <- list(subV = subV, subS = subS)
  
  ###Printing the subsetted S and V matrices list
  return(SUBSET_LDSC)
}