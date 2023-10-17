.userGWAS_main <- function(i, cores, k, n, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order, SNPs2, beta_SNP, SE_SNP,
                           varSNP, GC, coords, smooth_check, TWAS, printwarn, toler, estimation, sub, Model1,
                           df, npar, utilfuncs=NULL, basemodel=NULL, returnlavmodel=FALSE) {
  # utilfuncs contains utility functions to enable this code to work on PSOC clusters (for Windows)
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
      assign(j, utilfuncs[[j]], envir=environment())
    }
  }
  V_SNP <- .get_V_SNP(SE_SNP, I_LD, varSNP, GC, coords, k, i)
  
  if (smooth_check) {
    Z_pre <- .get_Z_pre(i, beta_SNP, SE_SNP, I_LD, GC)
  }
  
  V_full <- .get_V_full(k, V_LD, varSNPSE2, V_SNP)
  
  if(eigen(V_full)$values[nrow(V_full)] <= 0){
    V_full <- as.matrix((nearPD(V_full, corr = FALSE))$mat)
    V_smooth <- 1
  }
  
  #reorder sampling covariance matrix based on what lavaan expects given the specified model
  V_full_Reorder <- V_full[order, order]
  u <- nrow(V_full_Reorder)
  W <- diag(u)
  diag(W) <- diag(V_full_Reorder)
  
  ##invert the reordered sampling covariance matrix to create a weight matrix
  W <- solve(W, tol=toler)
  
  #create empty vector for S_SNP
  S_SNP <- vector(mode="numeric",length=k+1)
  
  #enter SNP variance from reference panel as first observation
  S_SNP[1] <- varSNP[i]
  
  #enter SNP covariances (standardized beta * SNP variance from refference panel)
  for (p in 1:k) {
    S_SNP[p+1] <- varSNP[i]*beta_SNP[i, p]
  }
  
  #create shell of the full S (observed covariance) matrix
  S_Fullrun <- diag(k+1)
  
  ##add the LD portion of the S matrix
  S_Fullrun[(2:(k+1)),(2:(k+1))] <- S_LD
  
  ##add in observed SNP variances as first row/column
  S_Fullrun[1:(k+1),1] <- S_SNP
  S_Fullrun[1,1:(k+1)] <- t(S_SNP)
  
  ##smooth to near positive definite if either V or S are non-positive definite
  ks <- nrow(S_Fullrun)
  if(eigen(S_Fullrun)$values[ks] <= 0){
    S_Fullrun <- as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat)
    S_smooth <- 1
  }
  
  if (smooth_check) {
    if(exists("S_smooth") | exists("V_smooth")){
      SE_smooth <- matrix(0, ks, ks)
      SE_smooth[lower.tri(SE_smooth,diag=TRUE)]  <- sqrt(diag(V_full))
      Z_smooth <- (S_Fullrun/SE_smooth)[2:ks,1]
      Z_smooth <- max(abs(Z_smooth-Z_pre))
    }else{
      Z_smooth <- 0
    }
  }
  
  #name the columns
  if(!(TWAS)){
    colnames(S_Fullrun) <- c("SNP", colnames(S_LD))
  } else {
    colnames(S_Fullrun) <- c("Gene", colnames(S_LD))
  }
  
  ##name rows like columns
  rownames(S_Fullrun) <- colnames(S_Fullrun)
  
  ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error.
  if(estimation == "DWLS"){
    if (!is.null(basemodel)) {
      test <- .tryCatch.W.E(Model1_Results <- lavaan(sample.cov = S_Fullrun, WLS.V=W, ordered=NULL, sampling.weights = NULL,
                                                     sample.mean=NULL, sample.th=NULL, sample.nobs=2, group=NULL, cluster= NULL, constraints='', NACOV=NULL,
                                                     slotOptions=basemodel@Options, slotParTable=basemodel@ParTable, slotSampleStats=NULL,
                                                     slotData=basemodel@Data, slotModel=basemodel@Model, slotCache=NULL, sloth1=NULL))
    } else {
      test <- .tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf,std.lv=std.lv))
    }
  } else if(estimation == "ML"){
    if (!is.null(basemodel)) {
      test <- .tryCatch.W.E(Model1_Results <- lavaan(sample.cov = S_Fullrun, WLS.V=NULL, ordered=NULL, sampling.weights = NULL,
                                                     sample.mean=NULL, sample.th=NULL, sample.nobs=200, group=NULL, cluster= NULL, constraints='', NACOV=NULL,
                                                     slotOptions=basemodel@Options, slotParTable=basemodel@ParTable, slotSampleStats=NULL,
                                                     slotData=basemodel@Data, slotModel=basemodel@Model, slotCache=NULL, sloth1=NULL))
    } else {
      test <- .tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML",sample.nobs = 200, optim.dx.tol = +Inf,sample.cov.rescale=FALSE,std.lv=std.lv))
    }
    
  }
  if (returnlavmodel) {
    return(Model1_Results)
  }
  
  test$warning$message[1] <- ifelse(is.null(test$warning$message), test$warning$message[1] <- 0, test$warning$message[1])
  
  if(class(test$value)[1] == "lavaan" & grepl("solution has NOT",  as.character(test$warning)) != TRUE){
    Model_Output <- parTable(Model1_Results)

    #pull the delta matrix (this doesn't depend on N)
    S2.delt <- lavInspect(Model1_Results, "delta")
    
    ##weight matrix from stage 2
    S2.W <- lavInspect(Model1_Results, "WLS.V")
    
    #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
    bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler)
    
    #create the "lettuce" part of the sandwich
    lettuce <- S2.W%*%S2.delt
    
    #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
    Ohtt <- bread %*% t(lettuce)%*%V_full_Reorder%*%lettuce%*%bread
    
    #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
    SE <- as.matrix(sqrt(diag(Ohtt)))
    
    #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
    if(estimation == "DWLS"){
      if(":=" %in% Model_Output$op & !(NA %in% Model_Output$se)){
        #variance-covariance matrix of parameter estimates, q-by-q (this is the naive one)
        vcov <- lavInspect(Model1_Results, "vcov")
        
        #internal lavaan representation of the model
        lavmodel <- Model1_Results@Model
        
        #lavaan representation of the indirect effect
        func <- lavmodel@def.function
        
        #vector of parameter estimates
        x <- lav_model_get_parameters(lavmodel, type = "free")
        
        #vector of indirect effect derivatives evaluated @ parameter estimates
        Jac <- lav_func_jacobian_complex(func = func, x = x)
        
        #replace vcov here with our corrected one. this gives parameter variance
        var.ind <- Jac %*% vcov %*% t(Jac)
        
        #square root of parameter variance = parameter SE.
        se.ghost <- sqrt(diag(var.ind))
        
        #pull the ghost parameter point estiamte
        ghost <- subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
        
        ##combine with delta method SE
        ghost2 <- cbind(ghost,se.ghost)
        colnames(ghost2)[7] <- "SE"
      }else{
        if(":=" %in% Model_Output$op & (NA %in% Model_Output$se)){
          se.ghost <- rep(NA, sum(":=" %in% Model_Output$op))
          warning("SE for ghost parameter not available for ML")
          ghost <- subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
          ghost2 <- cbind(ghost,se.ghost)
          colnames(ghost2)[7] <- "SE"
        }
      }
    }
    
    #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
    if(estimation == "ML"){
      if(":=" %in% Model_Output$op){
        #pull the ghost parameter point estiamte
        ghost <- subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
        se.ghost <- rep(NA, sum(":=" %in% Model_Output$op))
        warning("SE for ghost parameter not available for ML")
        ##combine with delta method SE
        ghost2 <- cbind(ghost,se.ghost)
        colnames(ghost2)[7] <- "SE"
      }else{
        if(":=" %in% Model_Output$op & (NA %in% Model_Output$se)){
          se.ghost <- rep(NA, sum(":=" %in% Model_Output$op))
          warning("SE for ghost parameter not available for ML")
          ghost <- subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
          ghost2 <- cbind(ghost,se.ghost)
          colnames(ghost2)[7] <- "SE"
        }
      }
    }
    
    #calculate model chi-square
    Eig <- as.matrix(eigen(V_full)$values)
    Eig2 <- diag((ncol(S_Fullrun)*(ncol(S_Fullrun)+1))/2)
    diag(Eig2) <- Eig
    
    #Pull P1 (the eigen vectors of V_eta)
    P1 <- eigen(V_full)$vectors
    
    implied <- as.matrix(fitted(Model1_Results))[1]
    implied_order <- colnames(S_Fullrun)
    implied[[1]] <- implied[[1]][implied_order,implied_order]
    implied2 <- S_Fullrun-implied[[1]]
    eta <- as.vector(lowerTriangle(implied2,diag=TRUE))
    Q <- t(eta)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta
    
    ##remove parameter constraints, ghost parameters, and fixed effects from output to merge with SEs
    unstand <- subset(Model_Output, Model_Output$plabel != "" & Model_Output$free > 0)[,c(2:4,8,11,14)]
    
    ##combine ghost parameters with rest of output
    if(exists("ghost2") == "TRUE"){
      unstand2 <- rbind(cbind(unstand,SE),ghost2)
    }else{
      unstand2 <- cbind(unstand,SE)
    }
    
    ##add in fixed effects and parameter constraints to output
    other <- subset(Model_Output, (Model_Output$plabel == "" & Model_Output$op != ":=") | (Model_Output$free == 0 & Model_Output$plabel != ""))[,c(2:4,8,11,14)]
    other$SE <- rep(NA, nrow(other))
    
    ##combine fixed effects and parameter constraints with output if there are any
    if(nrow(other) > 0){
      final <- rbind(unstand2,other)
    }else{
      final <- unstand2
    }
    
    #reorder based on row numbers so it is in order the user provided
    final$index <- as.numeric(row.names(final))
    final <- final[order(final$index), ]
    final$index <- NULL
    
    ##add in p-values
    if(class(final$SE) != "factor"){
      final$Z_Estimate <- final$est/final$SE
      final$Pval_Estimate <- 2*pnorm(abs(final$Z_Estimate),lower.tail=FALSE)
    }else{
      final$SE <- as.character(final$SE)
      final$Z_Estimate <- NA
      final$Pval_Estimate <- NA
    }

    #remove redundant internal representations of model from lavaan
    final<-subset(final, final$op != "da")
  
    ##add in model fit components to each row
    if(!(is.na(Q))){
      final$chisq <- rep(Q,nrow(final))
      final$chisq_df <- df
      final$chisq_pval <- pchisq(final$chisq,final$chisq_df,lower.tail=FALSE)
      final$AIC <- rep(Q + 2*npar,nrow(final))
    }else{final$chisq <- rep(NA, nrow(final))
    final$chisq_df <- rep(NA,nrow(final))
    final$chisq_pval <- rep(NA,nrow(final))
    final$AIC <- rep(NA, nrow(final))}
    
    ##add in error and warning messages
    if(printwarn){
      final$error <- ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
      final$warning <- ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]
    }
    ##combine with rs-id, BP, CHR, etc.
    final2 <- cbind(n + (i-1) * cores,SNPs2[i,],final,row.names=NULL)
    
    if(smooth_check){
      final2 <- cbind(final2,Z_smooth)
    }
    
    if(!(sub[[1]])==FALSE){
      final2 <- subset(final2, paste0(final2$lhs, final2$op, final2$rhs, sep = "") %in% sub)
    }else{##pull results
      final2$est <- ifelse(final2$op == "<" | final2$op == ">" | final2$op == ">=" | final2$op == "<=", final2$est == NA, final2$est)
    }
  }else{
    
     final <- data.frame(t(rep(NA, 13)))
 
    if(printwarn){
      final$error <- ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
      final$warning <- ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]}
    ##combine results with SNP, CHR, BP, A1, A2 for particular model
    final2 <- cbind(n + (i-1) * cores, SNPs2[i,], final, row.names=NULL)
    
    if(smooth_check){
      final2 <- cbind(final2,Z_smooth)
    }
  }
  if(TWAS){
    new_names <- c("i", "Gene","Panel","HSQ", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC","error","warning")
  } else {
    new_names <- c("i", "SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC","error","warning")
  }
  if(smooth_check)
    new_names <- c(new_names, "Z_smooth")
  colnames(final2) <- new_names

  return(final2)
}
