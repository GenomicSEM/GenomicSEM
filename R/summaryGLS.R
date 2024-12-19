summaryGLS <- function(OBJECT = NULL, Y = NULL, V_Y = NULL, PREDICTORS, INTERCEPT = T){
  
  #if (length(PREDICTORS[,1]) != length(LDSC_OBJECT$subS)) {
  #   warning("The length of predictors must be the same length as the parameters to be modeled.")
  #}
  
  ###Getting the length of SE estimation predictors
  if (is.vector(PREDICTORS)) {
    # Convert the vector to a matrix with 1 column and the number of rows equal to the length of the vector
    PREDICTORS <- matrix(PREDICTORS, ncol = 1)}
  
  num_predictors <- nrow(PREDICTORS)
  #Default to adding quadratic column to set of predictors
  for (i in 1:ncol(as.matrix(PREDICTORS))) {
    colnames(PREDICTORS)[i] <- paste("b", i, sep = "")
  }
  
  ###Creating object for outcome variables for GLS analysis
  ###Creating matrix of predictors, defaulting to including intercept
  
  if (INTERCEPT) {
    X <- cbind(rep(1, num_predictors), PREDICTORS)
    colnames(X)[1] <- "b0"
  } else {
    X <- as.matrix(PREDICTORS)
  }
  
  ###Creating a matrix of subsetted V matrix
  if (is.null(Y)){
    Ohm <- matrix(unlist(OBJECT[1]), nrow = num_predictors, ncol = num_predictors)
  }  
  
  if (is.null(OBJECT)){
    Ohm <- V_Y
  }  
  
  ###Creating a vector of genetic correlations
  if (is.null(Y)){
    y   <- unlist(OBJECT[2])
  }  
  
  if (is.null(OBJECT)){
    y <- Y
  }  
  
  ###Estimating betas with GLS equation
  BETAS <- solve(t(X) %*% solve(Ohm) %*% as.matrix(X)) %*% t(X) %*% solve(Ohm) %*% y #beta      
  
  ###Grabbing standard errors of 
  SE <- sqrt(diag(solve(t(X) %*% solve(Ohm) %*% as.matrix(X)))) #SE of beta
  Z <- BETAS/SE
  Pvals <- 2*pnorm(abs(BETAS)/SE,lower.tail=FALSE) #p value 
  
  outputGLS <- cbind(BETAS[,1], Pvals, SE, (BETAS/SE))
  colnames(outputGLS) <- c("betas", "pvals", "SE", "Z")
  print(outputGLS)     
}
