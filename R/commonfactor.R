commonfactor <-function(covstruc,estimation="DWLS"){ 
  time<-proc.time()
  
  #function to rearrange the sampling covariance matrix from original order to lavaan's order: 
  # 'k' is the number of variables in the model
  # 'fit' is the fit function of the regression model
  # 'names' is a vector of variable names in the order you used
  rearrange <- function (k, fit, names) {
    order1 <- names
    order2 <- rownames(inspect(fit)[[1]]) #order of variables
    kst <- k*(k+1)/2
    covA <- matrix(NA, k, k)
    covA[lower.tri(covA, diag = TRUE)] <- 1:kst
    covA <- t(covA)
    covA[lower.tri(covA, diag = TRUE)] <- 1:kst 
    colnames(covA) <- rownames(covA) <- order1 #give A actual variable order from lavaan output
    #reorder A by order2
    covA <- covA[order2, order2] #rearrange rows/columns
    vec2 <- lav_matrix_vech(covA) #grab new vectorized order
    return(vec2)
  }
  
  #function to create lavaan syntax for a 1 factor model given k phenotypes
  write.Model1 <- function(k, label = "V", label2 = "VF") {  
    Model1 <- ""
    for (i in 1) {
      linestart <- paste("F1"," =~ NA*",label, i, sep = "")  
      if (k-i > 0) {
        linemid <- ""
        for (j in (i+1):k) {
          linemid <- paste(linemid, " + ", label, j, sep = "")
        }
      } else {linemid <- ""}
      
    }
    Model1 <- paste(Model1, linestart, linemid, " \n ", "F1 ~~ 1*F1", " \n ", sep = "")
    
    ModelsatF<-""
    for (i in 1:(k-1)) {
      linestartc <- paste(label2, i, "~~ 0*", label2, i+1,  sep = "")
      if (k-i >= 2) {
        linemidc <- ""
        for (j in (i+2):k) {
          linemidc <- paste(linemidc, " + 0*", label2, j, sep = "")
        }
      } else {linemidc <- ""}
      ModelsatF <- paste(ModelsatF, linestartc, linemidc, " \n ", sep = "")
    } 
    
    Model1b <- ""
    for (i in 1) {
      linestartb <- paste("F1"," =~ 0*",label2, i, sep = "")  
      if ((k-1)-i > 0) {
        linemidb <- ""
        for (j in (i+1):k) {
          linemidb <- paste(linemidb, " + 0*", label2, j, sep = "")
        }
      } else {linemidb <- ""}
      
    }
    Model1b <- paste(Model1b, linestartb, linemidb, " \n ", sep = "")
    
    Model2<-""
    for (p in 1:k) {
      linestart2 <- paste(label2, p, " =~ 1*", label, p, sep = "")
      Model2<-paste(Model2, linestart2, " \n ", sep = "")}
    
    Model3<-""
    for (p in 1:k) {
      linestart3 <- paste(label, p, " ~~ ", label, p, sep = "")
      Model3<-paste(Model3, linestart3, " \n ", sep = "")}
    
    Model4<-""
    for (p in 1:k) {
      linestart4 <- paste(label2, p, " ~~ 0*", label2, p, sep = "")
      Model4<-paste(Model4, linestart4, " \n ", sep = "")}
    
    Modelsat<-""
    for (i in 1:(k-1)) {
      linestartc <- paste(label, i, "~~ 0*", label, i+1,  sep = "")
      if (k-i >= 2) {
        linemidc <- ""
        for (j in (i+2):k) {
          linemidc <- paste(linemidc, " + 0*", label, j, sep = "")
        }
      } else {linemidc <- ""}
      Modelsat <- paste(Modelsat, linestartc, linemidc, " \n ", sep = "")
    } 
    
    Model5<-paste(Model1, ModelsatF, Model1b, Model2, Model3, Model4, Modelsat, sep = "")
    
    return(Model5)
  } 

  ##code to write null model for calculation of CFI
  write.null<-function(k, label = "V", label2 = "VF") {
    Model3<-""
    for (p in 1:k) {
      linestart3 <- paste(label, p, " ~~ ", label, p, sep = "")
      Model3<-paste(Model3, linestart3, " \n ", sep = "")}
    
    Model2<-""
    for (p in 1:k) {
      linestart2 <- paste(label2, p, " =~ 1*", label, p, sep = "")
      Model2<-paste(Model2, linestart2, " \n ", sep = "")}
    
    Modelsat<-""
    for (i in 1:(k-1)) {
      linestartc <- paste(label, i, "~~ 0*", label, i+1,  sep = "")
      if (k-i >= 2) {
        linemidc <- ""
        for (j in (i+2):k) {
          linemidc <- paste(linemidc, " + 0*", label, j, sep = "")
        }
      } else {linemidc <- ""}
      Modelsat <- paste(Modelsat, linestartc, linemidc, " \n ", sep = "")
    }
    
    ModelsatF<-""
    for (i in 1:(k-1)) {
      linestartc <- paste(" ", label2, i, "~~ 0*", label2, i+1,  sep = "")
      if (k-i >= 2) {
        linemidc <- ""
        for (j in (i+2):k) {
          linemidc <- paste(linemidc, " + 0*", label2, j, sep = "")
        }
      } else {linemidc <- ""}
      ModelsatF <- paste(ModelsatF, linestartc, linemidc, " \n ", sep = "")
    } 
    
    Model4<-""
    for (p in 1:k) {
      linestart4 <- paste(label2, p, " ~~ 0*", label2, p, sep = "")
      Model4<-paste(Model4, linestart4, " \n ", sep = "")}
    
    modelCFI<-paste(Model3, Model2, ModelsatF, Modelsat, Model4)
    return(modelCFI)
  }
  
  #function to creat row/column names for S_LD matrix
  write.names <- function(k, label = "V") {  
    varnames<-vector(mode="character",length=k)
    
    for (i in 1:k){
      varnames[i]<-paste(label,i,sep="")}
    
    return(varnames)
  }
  
  ##modification of trycatch that allows the results of a failed run to still be saved
  tryCatch.W.E <- function(expr)
  {
    .warning <- NULL
    w.handler <- function(w) { # warning handler
      .warning <<- c( .warning, w$message )
      invokeRestart("muffleWarning")
    }
    list(
      value=withCallingHandlers(
        tryCatch( expr, error=function(e) e ),
        warning=w.handler
      ),
      warning=paste( .warning, collapse='\n' )
    )
  }
  
  ##read in the LD portion of the V (sampling covariance) matrix
  V_LD<-as.matrix(covstruc[[1]])
  
  ##read in the LD portion of the S (covariance) matrix
  S_LD<-as.matrix(covstruc[[2]])
  
  ##set the number of observations (typically LD-score bins) behind the covariance matrix
  S_LD_nobs<-covstruc$n  # this is what it should be, but see -> __NOBS_OVERRIDE__
  
  ##k = number of phenotypes in dataset (i.e., number of columns in LD portion of S matrix)
  k<-ncol(S_LD)
  
  ##create the 1 factor model with k # of indicators
  Model1 <- write.Model1(k)
  
  ##create inependence model for calculation of CFI
  modelCFI<-write.null(k)
  
  ##pull the column names specified in the munge function
  traits<-colnames(S_LD)
  
  ##create the names
  S_names<-write.names(k=k)
  
  if(is.null(traits)){
    traits<-S_names} 
  
  ##name the columns and rows of the S matrix
  rownames(S_LD) <- S_names
  colnames(S_LD) <- S_names
  
  ##smooth to near positive definite if either V or S are non-positive definite
  ks<-nrow(S_LD)
  S_LDb<-S_LD
  smooth1<-ifelse(eigen(S_LD)$values[ks] <= 0, S_LD<-as.matrix((nearPD(S_LD, corr = FALSE))$mat), S_LD<-S_LD)
  diff<-(S_LD-S_LDb)
  LD_sdiff<-max(diff)
  
  kv<-nrow(V_LD)
  V_LDb<-V_LD
  smooth2<-ifelse(eigen(V_LD)$values[kv] <= 0, V_LD<-as.matrix((nearPD(V_LD, corr = FALSE))$mat), V_LD<-V_LD)
  diff2<-(V_LD-V_LDb)
  LD_sdiff2<-max(diff2)
  
  ##run model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  #transform V_LD matrix into a diagonal weight matrix: 
  z<-(k*(k+1))/2
 
  ##save the ordering
  order <-(1:nrow(V_LD))

  ##reorder the weight (inverted V_LD) matrix
  V_Reorder<-V_LD[order,order]

  ##run CFI model so it knows the reordering
  #empty2<-tryCatch.W.E(fitCFI <- sem(modelCFI, sample.cov = S_LD, estimator = "DWLS", WLS.V = solve(V_LD),sample.nobs=2, optim.dx.tol = +Inf))
  #orderCFI <- rearrange(k = k, fit =  fitCFI, names =  rownames(S_LD))
  
  ##save the ordering
  orderCFI<-(1:nrow(V_LD))

  ##reorder matrix for independence (i.e., null) model for CFI calculation
  V_Reorder2<-V_LD[orderCFI,orderCFI]
  
  ##set estimation-specific parameters
  S_LD_nobs = 200  # see __NOBS_OVERRIDE__
  W_Reorder = NULL
  W_CFI = NULL
  if ( estimation == "DWLS" ) {
    S_LD_nobs = 2  # see __NOBS_OVERRIDE__
    V_Reorderb<-diag(z)
    diag(V_Reorderb)<-diag(V_Reorder)
    W_Reorder<-solve(V_Reorderb)
    V_Reorder2b<-diag(z)
    diag(V_Reorder2b)<-diag(V_Reorder2)
    W_CFI<-solve(V_Reorder2b)
  }

  ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
  empty <- tryCatch.W.E(
    Model1_Results <- sem(
      Model1,
      sample.cov = S_LD,
      estimator = estimation,
      WLS.V = W_Reorder,
      sample.nobs = S_LD_nobs,
      optim.dx.tol = +Inf
    )
  )
  
  if ( empty$warning == '' ) empty$warning = 0
  
  if ( class(empty$value)[1] == "simpleError" || any( grepl("solution has NOT",  as.character(empty$warning)) ) ) {
    print("The common factor initially failed to converge. A lower bound of 0 on residual variances has been added to try and troubleshoot this.")
    
    write.Model1 <- function(k, label = "V", label2 = "VF") {  
      Model1 <- ""
      for (i in 1) {
        linestart <- paste("F1"," =~ NA*",label, i, sep = "")  
        if (k-i > 0) {
          linemid <- ""
          for (j in (i+1):k) {
            linemid <- paste(linemid, " + ", label, j, sep = "")
          }
        } else {linemid <- ""}
        
      }
      Model1 <- paste(Model1, linestart, linemid, " \n ", "F1 ~~ 1*F1", " \n ", sep = "")
      
      ModelsatF<-""
      for (i in 1:(k-1)) {
        linestartc <- paste(label2, i, "~~ 0*", label2, i+1,  sep = "")
        if (k-i >= 2) {
          linemidc <- ""
          for (j in (i+2):k) {
            linemidc <- paste(linemidc, " + 0*", label2, j, sep = "")
          }
        } else {linemidc <- ""}
        ModelsatF <- paste(ModelsatF, linestartc, linemidc, " \n ", sep = "")
      } 
      
      Model1b <- ""
      for (i in 1) {
        linestartb <- paste("F1"," =~ 0*",label2, i, sep = "")  
        if ((k-1)-i > 0) {
          linemidb <- ""
          for (j in (i+1):k) {
            linemidb <- paste(linemidb, " + 0*", label2, j, sep = "")
          }
        } else {linemidb <- ""}
        
      }
      Model1b <- paste(Model1b, linestartb, linemidb, " \n ", sep = "")
      
      Model2<-""
      for (p in 1:k) {
        linestart2 <- paste(label2, p, " =~ 1*", label, p, sep = "")
        Model2<-paste(Model2, linestart2, " \n ", sep = "")}
      
      Model3<-""
      for (p in 1:k) {
        linestart3a <- paste(label, p, " ~~ ", letters[p], "*", label, p, sep = "")
        linestart3b <- paste(letters[p], " > .001", sep = "")
        Model3<-paste(Model3, linestart3a, " \n ", linestart3b, " \n ", sep = "")}
      
      Model4<-""
      for (p in 1:k) {
        linestart4 <- paste(label2, p, " ~~ 0*", label2, p, sep = "")
        Model4<-paste(Model4, linestart4, " \n ", sep = "")}
      
      Modelsat<-""
      for (i in 1:(k-1)) {
        linestartc <- paste(label, i, "~~ 0*", label, i+1,  sep = "")
        if (k-i >= 2) {
          linemidc <- ""
          for (j in (i+2):k) {
            linemidc <- paste(linemidc, " + 0*", label, j, sep = "")
          }
        } else {linemidc <- ""}
        Modelsat <- paste(Modelsat, linestartc, linemidc, " \n ", sep = "")
      } 
      
      Model5<-paste(Model1, ModelsatF, Model1b, Model2, Model3, Model4, Modelsat, sep = "")
      
      return(Model5)
    }
    Model1 <- write.Model1(k)
    empty <- tryCatch.W.E(
      Model1_Results <- sem(
        Model1,
        sample.cov = S_LD,
        estimator = estimation,
        WLS.V = W_Reorder,
        sample.nobs = S_LD_nobs,
        optim.dx.tol = +Inf
      )
    )
  }
  
  if ( class(empty$value)[1] == "simpleError" ) {
    print("The common factor model failed to converge on a solution. Please try specifying an alternative model using the usermodel function.")
  }
  
  if ( any( grepl("solution has NOT",  as.character(empty$warning)) ) ) {
    print("The common factor model failed to converge on a solution. Please try specifying an alternative model using the usermodel function.")
  }
   
  #pull the delta matrix (this doesn't depend on N)
  S2.delt <- lavInspect(Model1_Results, "delta")
  
  ##weight matrix from stage 2
  S2.W <- lavInspect(Model1_Results, "WLS.V") 
  
  #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
  bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
  
  #create the "lettuce" part of the sandwich
  lettuce <- S2.W%*%S2.delt
  
  #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
  Ohtt <- bread %*% t(lettuce)%*%V_Reorder%*%lettuce%*%bread  
  
  #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
  SE <- as.matrix(sqrt(diag(Ohtt)))
  
  ModelQ <- parTable(Model1_Results)
  
  ModelQ<-subset(ModelQ, ModelQ$plabel != "")
  
  ##fix indicator loadings, residual factor loadings, and residual factor variances from stage 2 estimation
  z<-((k*(k+1))/2)
  f<-length(ModelQ$free)-z
  
  ModelQ$free <- c(rep(0, f), 1:z)
  
  ModelQ$ustart <- ModelQ$est
  ModelQ$ustart<-ifelse(ModelQ$free > 0, .05, ModelQ$ustart)
  
  testQ <- tryCatch.W.E(
    ModelQ_Results <- sem(
      model = ModelQ,
      sample.cov = S_LD,
      estimator = estimation,
      WLS.V = W_Reorder,
      sample.nobs = S_LD_nobs,
      start = ModelQ$ustart,
      optim.dx.tol = +Inf
    )
  )
  if ( testQ$warning == '' ) testQ$warning = "Safe"
  
  if ( any( grepl( "lavaan WARNING: model has NOT converged!", as.character(testQ$warning) ) ) ) {
    
    ModelQ$ustart<-ifelse(ModelQ$free > 0, .01, ModelQ$ustart)
    
    testQ2 <- tryCatch.W.E(
      ModelQ_Results <- sem(
        model = ModelQ,
        sample.cov = S_LD,
        estimator = estimation,
        WLS.V = W_Reorder,
        sample.nobs = S_LD_nobs,
        start = ModelQ$ustart,
        optim.dx.tol = +Inf
      )
    )
  } else { testQ2<-testQ }
  
  if ( testQ2$warning == '' ) testQ2$warning = "Safe"
  
  if ( any( grepl( "lavaan WARNING: model has NOT converged!", as.character(testQ2$warning) ) ) ) {
    
    ModelQ$ustart<-ifelse(ModelQ$free > 0, .1, ModelQ$ustart)
    
    testQ3 <- tryCatch.W.E(
      ModelQ_Results <- sem(
        model = ModelQ,
        sample.cov = S_LD,
        estimator = estimation,
        WLS.V = W_Reorder,
        sample.nobs = S_LD_nobs,
        start = ModelQ$ustart,
        optim.dx.tol = +Inf
      )
    ) 
  } else { testQ3<-testQ2 }
  
  if ( testQ3$warning == '' ) testQ3$warning = "Safe"
  
  if ( !any( grepl( "lavaan WARNING: model has NOT converged!", as.character(testQ3$warning) ) ) ) {
    
    #pull the delta matrix (this doesn't depend on N)
    S2.delt_Q <- lavInspect(ModelQ_Results, "delta")
    
    ##weight matrix from stage 2
    S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V") 
    
    #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
    bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q) 
    
    #create the "lettuce" part of the sandwich
    lettuce_Q <- S2.W_Q%*%S2.delt_Q
    
    #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
    Ohtt_Q <- bread_Q %*% t(lettuce_Q)%*%V_Reorder%*%lettuce_Q%*%bread_Q  
    
    ##pull the sampling covariance matrix of the residual covariances and compute diagonal matrix of eigenvalues
    V_eta<- Ohtt_Q
    Eig2<-as.matrix(eigen(V_eta)$values)
    Eig<-diag(z)
    diag(Eig)<-Eig2
    
    #Pull P1 (the eigen vectors of V_eta)
    P1<-eigen(V_eta)$vectors
    
    ##Pull eta = vector of residual covariances
    eta_test<-parTable(ModelQ_Results)
    eta_test<-subset(eta_test, eta_test$free != 0)
    eta<-cbind(eta_test[,14])
    
    #Combining all the pieces from above:
    Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta

  } else { Q<-NA }
  
  ##now CFI
  ##run independence model
  testCFI <- tryCatch.W.E(
    fitCFI <- sem(
      modelCFI, sample.cov = S_LD,
      estimator = estimation,
      WLS.V = W_CFI,
      sample.nobs = S_LD_nobs,
      optim.dx.tol = +Inf
    )
  )
  if ( testCFI$warning == '' ) testCFI$warning = "Safe"
  
  if ( !any( grepl( "lavaan WARNING: model has NOT converged!", as.character(testCFI$warning) ) ) ) {
    
    ##code to estimate chi-square of independence model#
    #First pull the estimates from Step 2
    ModelQ_CFI <- parTable(fitCFI)
    p2<-length(ModelQ_CFI$free)-z
    
    ##fix variances and freely estimate covariances
    ModelQ_CFI$free <- c(rep(0, p2), 1:z)
    ModelQ_CFI$ustart <- ModelQ_CFI$est
    
    testCFI2 <- tryCatch.W.E(
      ModelQ_Results_CFI <- sem(
        model = ModelQ_CFI,
        sample.cov = S_LD,
        estimator = estimation,
        WLS.V = W_CFI,
        sample.nobs = S_LD_nobs,
        optim.dx.tol = +Inf
      )
    )
    if ( testCFI2$warning == '' ) testCFI2$warning = "Safe"
    
    if ( !any( grepl( "lavaan WARNING: model has NOT converged!", as.character(testCFI2$warning) ) ) ) {
      
      #pull the delta matrix (this doesn't depend on N)
      S2.delt_Q_CFI <- lavInspect(ModelQ_Results_CFI, "delta")
      
      ##weight matrix from stage 2
      S2.W_Q_CFI <- lavInspect(ModelQ_Results_CFI, "WLS.V") 
      
      #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
      bread_Q_CFI <- solve(t(S2.delt_Q_CFI)%*%S2.W_Q_CFI%*%S2.delt_Q_CFI) 
      
      #create the "lettuce" part of the sandwich
      lettuce_Q_CFI <- S2.W_Q_CFI%*%S2.delt_Q_CFI
      
      #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
      Ohtt_Q_CFI <- bread_Q_CFI %*% t(lettuce_Q_CFI)%*%V_Reorder2%*%lettuce_Q_CFI%*%bread_Q_CFI
      
      ##pull the sampling covariance matrix of the residual covariances and compute diagonal matrix of eigenvalues
      V_etaCFI<- Ohtt_Q_CFI
      Eig2_CFI<-as.matrix(eigen(V_etaCFI)$values)
      Eig_CFI<-diag(z)
      diag(Eig_CFI)<-Eig2_CFI
      
      #Pull P1 (the eigen vectors of V_eta)
      P1_CFI<-eigen(V_etaCFI)$vectors
      
      ##Pull eta = vector of residual covariances
      eta_test_CFI<-parTable(ModelQ_Results_CFI)
      eta_test_CFI<-subset(eta_test_CFI, eta_test_CFI$free != 0)
      eta_CFI<-cbind(eta_test_CFI[,14])
      
      #Combining all the pieces from above:
      Q_CFI<-t(eta_CFI)%*%P1_CFI%*%solve(Eig_CFI)%*%t(P1_CFI)%*%eta_CFI

    } else {Q_CFI<-"NA"}

  }
  
  ##to get standardized output
  ##transform the S covariance matrix to S correlation matrix
  D=sqrt(diag(diag(S_LD)))
  S_Stand=solve(D)%*%S_LD%*%solve(D)
  rownames(S_Stand)<-rownames(S_LD)
  colnames(S_Stand)<-colnames(S_Stand)
  
  #obtain diagonals of the original V matrix and take their sqrt to get SE's
  Dvcov<-sqrt(diag(V_LD))
  
  #calculate the ratio of the rescaled and original S matrices
  scaleO=as.vector(lowerTriangle((S_Stand/S_LD),diag=T))
  
  #rescale the SEs by the same multiples that the S matrix was rescaled by
  Dvcovl<-as.vector(Dvcov*t(scaleO))
  
  #obtain the sampling correlation matrix by standardizing the original V matrix
  Vcor<-cov2cor(V_LD)
  
  #rescale the sampling correlation matrix by the appropriate diagonals
  V_stand<-diag(Dvcovl)%*%Vcor%*%diag(Dvcovl)
  V_stand2<-diag(z)
  diag(V_stand2)<-diag(V_stand)
  W_stand<-solve(V_stand2[order,order])
  
  fit_stand <- sem(Model1, sample.cov = S_Stand, estimator = estimation, WLS.V = W_stand, sample.nobs = S_LD_nobs, optim.dx.tol = +Inf) 
  
  ##perform same procedures for sandwich correction as in the unstandardized case
  delt_stand <- lavInspect(fit_stand, "delta") 
  W_stand <- lavInspect(fit_stand, "WLS.V") 
  bread_stand <- solve(t(delt_stand)%*%W_stand %*%delt_stand) 
  lettuce_stand <- W_stand%*%delt_stand
  Vcov_stand<-as.matrix(V_stand[order,order])
  Ohtt_stand <- bread_stand %*% t(lettuce_stand)%*%Vcov_stand%*%lettuce_stand%*%bread_stand
  SE_stand <- as.matrix(sqrt(diag(Ohtt_stand)))
  
  if (Q_CFI != "NA") {
    CFI <- as.numeric( (
      (Q_CFI-lavInspect(fitCFI, "fit")["df"]) - (Q-lavInspect(Model1_Results, "fit")["df"])
    ) / (Q_CFI-lavInspect(fitCFI, "fit")["df"]) )
    CFI <- ifelse(CFI > 1, 1, CFI)
  } else { CFI <- NA }
  
  chisq<-Q
  df<-lavInspect(Model1_Results, "fit")["df"]
  AIC<-(Q + 2*lavInspect(Model1_Results, "fit")["npar"])
  SRMR<-lavInspect(Model1_Results, "fit")["srmr"]
  
  unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
  unstand<-subset(unstand, unstand$free != 0)                    
  unstand$free<-NULL
  
  stand<-data.frame(inspect(fit_stand,"list")[,c(8,14)])
  stand<-subset(stand, stand$free != 0)
  stand$free<-NULL
  
  modelfit<-cbind(chisq,df,AIC,CFI,SRMR)
  results<-cbind(unstand,SE,stand,SE_stand)
  
  ##name the columns of the results file
  colnames(results)=c("lhs","op","rhs","Unstandardized_Estimate","Unstandardized_SE","Standardized_Estimate","Standardized_SE")
  
  ##replace V1-VX general form in output with user provided trait names
  for(i in 1:nrow(results)){
    for(p in 1:length(traits)){
      results$lhs[[i]]<-ifelse(results$lhs[[i]] %in% S_names[[p]], gsub(results$lhs[[i]], traits[[p]], results$lhs[[i]]), results$lhs[[i]])
      results$rhs[[i]]<-ifelse(results$rhs[[i]] %in% S_names[[p]], gsub(results$rhs[[i]], traits[[p]], results$rhs[[i]]), results$rhs[[i]])
    }
  }
  
  ##name model fit columns
  colnames(modelfit)=c("chisq","df","AIC","CFI","SRMR")
  modelfit<-data.frame(modelfit)
  modelfit$p_chisq<-ifelse(modelfit$chisq != 'NA', modelfit$p_chisq<-pchisq(modelfit$chisq, modelfit$df,lower.tail=FALSE), modelfit$p_chisq<-NA)
  modelfit$chisq<-ifelse(modelfit$df == 0, modelfit$chisq == NA, modelfit$chisq)  
  modelfit$AIC<-ifelse(modelfit$df == 0, modelfit$AIC == NA, modelfit$AIC)  
  modelfit$p_chisq<-ifelse(modelfit$df == 0, modelfit$p_chisq == NA, modelfit$p_chisq)  
  
  order<-c(1,2,6,3,4,5)
  modelfit<-modelfit[,order]
  
  time_all<-proc.time()-time
  print(time_all[3])
  
  if(LD_sdiff > 0){
    print(paste("The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest
                difference in a cell between the smoothed and non-smoothed matrix was", LD_sdiff, sep = " "))
  }
  
  if(LD_sdiff2 > 0){
    print(paste("The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest
                difference in a cell between the smoothed and non-smoothed matrix was", LD_sdiff2, sep = " "))
  }
  
  if(!is.na( modelfit$CFI ) && modelfit$CFI < 0){
    warning(paste("CFI estimates below 0 should not be trusted, and indicate that the other model fit estimates should be interpreted with caution. A negative CFI estimates typically appears due to negative residual variances."))
  }
         
  results$p_value<-2*pnorm(abs(as.numeric(results$Unstandardized_Estimate)/as.numeric(results$Unstandardized_SE)),lower.tail=FALSE)
  results$p_value<-ifelse(results$p_value == 0, "< 5e-300", results$p_value)
  
  return(list(modelfit=modelfit,results=results))
  
}

