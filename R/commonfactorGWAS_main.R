.commonfactorGWAS_main <- function(i, cores, n, S_LD, V_LD, I_LD, beta_SNP, SE_SNP, varSNP, varSNPSE2, GC, coords, k,
                                   smooth_check, Model1, toler, estimation, order, utilfuncs=NULL, basemodel=NULL, returnlavmodel=FALSE) {
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
      assign(j, utilfuncs[[j]], envir=environment())
    }
  }
  #create empty shell of V_SNP matrix
  V_SNP <- .get_V_SNP(SE_SNP, I_LD, varSNP, GC, coords, k, i)

  if(smooth_check){
    if(GC == "conserv"){
      Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*diag(I_LD))
    } else if(GC == "standard"){
      Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*sqrt(diag(I_LD)))
    } else if(GC == "none"){
      Z_pre<-beta_SNP[i,]/SE_SNP[i,]
    }
  }

  ##create shell of full sampling covariance matrix
  V_Full <- .get_V_full(k, V_LD, varSNPSE2, V_SNP)

  kv<-nrow(V_Full)
  if(eigen(V_Full)$values[kv] <= 0){
    V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat)
    V_smooth<-1
  }

  #reorder sampling covariance matrix based on what lavaan expects given the specified model
  V_Full_Reorder <- V_Full[order, order]
  u<-nrow(V_Full_Reorder)
  W<-diag(u)
  diag(W)<-diag(V_Full_Reorder)

  ##invert the reordered sampling covariance matrix to create a weight matrix
  W <- solve(W,tol=toler)

  #create empty vector for S_SNP
  S_SNP<-vector(mode="numeric",length=k+1)

  #enter SNP variance from reference panel as first observation
  S_SNP[1]<-varSNP[i]

  #enter SNP covariances (standardized beta * SNP variance from refference panel)
  for (p in 1:k) {
    S_SNP[p+1]<-varSNP[i]*beta_SNP[i,p]
  }

  #create shell of the full S (observed covariance) matrix
  S_Fullrun<-diag(k+1)

  ##add the LD portion of the S matrix
  S_Fullrun[(2:(k+1)),(2:(k+1))]<-S_LD

  ##add in observed SNP variances as first row/column
  S_Fullrun[1:(k+1), 1] <- S_SNP
  S_Fullrun[1, 1:(k+1)] <- t(S_SNP)

  colnames(S_Fullrun)<-c("SNP", colnames(S_LD))

  ##name rows like columns
  rownames(S_Fullrun)<-colnames(S_Fullrun)

  ##smooth to near positive definite if either V or S are non-positive definite
  ks<-nrow(S_Fullrun)

  if(eigen(S_Fullrun)$values[ks] <= 0){
    S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat)
    S_smooth<-1
  }

  if(smooth_check){
    if(exists("S_smooth") | exists("V_smooth")){
      SE_smooth<-matrix(0, ks, ks)
      SE_smooth[lower.tri(SE_smooth,diag=TRUE)] <-sqrt(diag(V_Full))
      Z_smooth<-(S_Fullrun/SE_smooth)[2:ks,1]
      Z_smooth<-max(abs(Z_smooth-Z_pre))
    }else{Z_smooth<-0}
  }

  ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error.
  if(estimation == "DWLS"){
    if (!is.null(basemodel)){
      test <- .tryCatch.W.E(Model1_Results <- lavaan(sample.cov = S_Fullrun, WLS.V=W, ordered=NULL, sampling.weights = NULL,
                                                           sample.mean=NULL, sample.th=NULL, sample.nobs=2, group=NULL, cluster= NULL, constraints='', NACOV=NULL,
                                                           slotOptions=basemodel@Options, slotParTable=basemodel@ParTable, slotSampleStats=NULL,
                                                           slotData=basemodel@Data, slotModel=basemodel@Model, slotCache=NULL, sloth1=NULL))
    }
    else {
      test <- .tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf))
    }
  }

  if(estimation == "ML"){
    if (!is.null(basemodel)){
      test <- .tryCatch.W.E(Model1_Results <- lavaan(sample.cov = S_Fullrun, WLS.V=NULL, ordered=NULL, sampling.weights = NULL,
                                                           sample.mean=NULL, sample.th=NULL, sample.nobs=200, group=NULL, cluster= NULL, constraints='', NACOV=NULL,
                                                           slotOptions=basemodel@Options, slotParTable=basemodel@ParTable, slotSampleStats=NULL,
                                                           slotData=basemodel@Data, slotModel=basemodel@Model, slotCache=NULL, sloth1=NULL))
    } else {
      test <- .tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML",sample.nobs = 200, optim.dx.tol = +Inf, sample.cov.rescale=FALSE))
    }

  }
  if (returnlavmodel){
    return(Model1_Results)
  }
  test$warning$message[1]<-ifelse(is.null(test$warning$message), test$warning$message[1]<-0, test$warning$message[1])

  if(class(test$value)[1] == "lavaan" & grepl("solution has NOT",  as.character(test$warning)) != TRUE){
    #pull the delta matrix (this doesn't depend on N)
    S2.delt <- lavInspect(Model1_Results, "delta")

    ##weight matrix from stage 2
    S2.W <- lavInspect(Model1_Results, "WLS.V")

    #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
    bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler)

    #create the "lettuce" part of the sandwich
    lettuce <- S2.W%*%S2.delt

    #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
    Ohtt <- bread %*% t(lettuce)%*%V_Full_Reorder%*%lettuce%*%bread

    #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
    SE <- as.matrix(sqrt(diag(Ohtt)))

    ##pull the corrected SE for SNP effect on P-factor
    se_c<-SE[k,1]

    ##code to estimate Q_SNP##
    #First pull the estimates from Step 1
    ModelQ <- parTable(Model1_Results)
    
    #2023 add: remove additional rows for internal representation of model by lavaan
     ModelQ<-ModelQ[1:((k*3)+3),]

    #fix the indicator loadings from Step 1, free the direct effects of the SNP on the indicators, and fix the factor residual variance
    ModelQ$free <- c(rep(0, k+1), 1:(k*2), 0, 0)

    ##added##
    ModelQ$ustart <- ModelQ$est
    SNPresid<-resid(Model1_Results)$cov[k+1,1:k]

    for(t in 1:nrow(ModelQ)) {
      if(ModelQ$free[t] > 0 & ModelQ$free[t] <= k){
        ModelQ$ustart[t]<-SNPresid[ModelQ$free[t]]} else{}}

    #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
    if(estimation == "DWLS"){
      testQ<-.tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2,  optim.dx.tol = +Inf))
    } else if(estimation == "ML"){
      testQ<-.tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "ML", sample.nobs = 200, optim.dx.tol = +Inf, sample.cov.rescale=FALSE))
    }

    testQ$warning$message[1]<-ifelse(is.null(testQ$warning$message), testQ$warning$message[1]<-"Safe", testQ$warning$message[1])
    testQ$warning$message[1]<-ifelse(exists("ModelQ_Results") == FALSE, testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!", ifelse(is.na(inspect(ModelQ_Results, "se")$theta[1,2]) == TRUE ,testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!" , testQ$warning$message[1]))

    if(as.character(testQ$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
      #pull the delta matrix for Q (this doesn't depend on N)
      S2.delt_Q <- lavInspect(ModelQ_Results, "delta")

      ##weight matrix from stage 2 for Q
      S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V")

      #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
      Q_catch<-.tryCatch.W.E(bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q,tol=toler))

      if(class(Q_catch$value)[1] == "matrix"){

        #create the "lettuce" part of the sandwich
        lettuce_Q <- S2.W_Q%*%S2.delt_Q

        #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
        Ohtt_Q <- bread_Q %*% t(lettuce_Q)%*%V_Full_Reorder%*%lettuce_Q%*%bread_Q

        ##compute diagonal matrix (Ron calls this lambda, we call it Eig) of eigenvalues of the sampling covariance matrix of the model residuals (V_eta)
        V_eta<- Ohtt_Q[1:k,1:k]
        Eig2<-as.matrix(eigen(V_eta)$values)
        Eig<-diag(k)
        diag(Eig)<-Eig2

        #Pull P1 (the eigen vectors of V_eta)
        P1<-eigen(V_eta)$vectors

        ##Pull eta = vector of direct effects of the SNP (Model Residuals)
        eta<-cbind(inspect(ModelQ_Results,"list")[(k+2):(2*k+1),14])

        #Combining all the pieces from above:
        Q<-t(eta)%*%P1%*%solve(Eig,tol=toler)%*%t(P1)%*%eta}} else{Q<-"Not Computed"}

    ##pull all the results into a single row
     if(smooth_check){
      results<-data.frame(n + (i-1) * cores,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value)[1] == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning)[1] == 'NULL', 0, as.character(test$warning$message[1])),Z_smooth,stringsAsFactors = FALSE)
      colnames(results) <- c("i", "lhs", "op", "rhs", "est", "se", "se_c", "Q", "fail", "warning", "Z_smooth")
       } else {
      results<-data.frame(n + (i-1) * cores,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value)[1] == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning)[1] == 'NULL', 0, as.character(test$warning$message[1])),stringsAsFactors = FALSE)
      colnames(results) <- c("i", "lhs", "op", "rhs", "est", "se", "se_c", "Q", "fail", "warning")
        }
    
    
  }else{
    # Reassign i to maintain original order in parallel operation
    if(smooth_check){
      results<-data.frame(n + (i-1) * cores,inspect(Model1_Results,"list")[k+1,-c(1,5:13,15)],t(rep(NA,3)),ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning)[1] == 'NULL', 0, as.character(test$warning$message[1])),Z_smooth,stringsAsFactors = FALSE)
      colnames(results) <- c("i", "lhs", "op", "rhs", "est", "se", "se_c", "Q", "fail", "warning", "Z_smooth")
      } else {
      results<-data.frame(n + (i-1) * cores,inspect(Model1_Results,"list")[k+1,-c(1,5:13,15)],t(rep(NA,3)),ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning)[1] == 'NULL', 0, as.character(test$warning$message[1])),stringsAsFactors = FALSE)
      colnames(results) <- c("i", "lhs", "op", "rhs", "est", "se", "se_c", "Q", "fail", "warning")
       }

  }
  return(results)
}
