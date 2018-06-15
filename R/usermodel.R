usermodel <-function(covstruc,estimation="DWLS", model = "", CFIcalc=TRUE){ 
  time<-proc.time()
  
  #function to rearrange the sampling covariance matrix from original order to lavaan's order: 
  #'k' is the number of variables in the model
  #'fit' is the fit function of the regression model
  #'names' is a vector of variable names in the order you used
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
  
  ##read in the LD portion of the V (sampling covariance) matrix
  V_LD<-as.matrix(covstruc[[1]])
  
  ##read in the LD portion of the S (covariance) matrix
  S_LD<-as.matrix(covstruc[[2]])
  
  ##k = number of phenotypes in dataset (i.e., number of columns in LD portion of S matrix)
  k<-ncol(S_LD)
  
  ##size of V matrix used later in code to create diagonal V matrix
  z<-(k*(k+1))/2
  
  #function to creat row/column names for S_LD matrix
  write.names <- function(k, label = "V") {  
    varnames<-vector(mode="character",length=k)
    
    for (i in 1:k){
      varnames[i]<-paste(label,i,sep="")}
    
    return(varnames)
  }
  
  ##create the names
  S_names<-write.names(k=k)
  
  ##pull the column names specified in the munge function
  traits<-colnames(S_LD)
  
  ##replace trait names in user provided model with general form of V1-VX
  for(i in 1:length(traits)){
    model<-gsub(traits[[i]], S_names[[i]], model)
  }
  
  Model1<-model
  
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
  #transform V_LD matrix into a weight matrix: 
  W <- solve(V_LD)
  
  ##run the model
  ReorderModel <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W, sample.nobs = 2) 
  
  ##determine number of latent variables from writing extended model
  r<-nrow(lavInspect(ReorderModel, "cor.lv"))
  
  write.Model1 <- function(k, label = "V", label2 = "VF") {  
    
    ModelsatF<-""
    for (i in 1:(k-1)) {
      linestartc <- paste(" ", label2, i, "~~0*", label2, i+1,  sep = "")
      if (k-i >= 2) {
        linemidc <- ""
        for (j in (i+2):k) {
          linemidc <- paste(linemidc, "+0*", label2, j, sep = "")
        }
      } else {linemidc <- ""}
      ModelsatF <- paste(ModelsatF, linestartc, linemidc, " \n ", sep = "")
    } 
    
    if(r > 0){
      Model1b <- ""
      for (t in 1:r) {
        for (i in 1) {
          linestartb <- paste("F", t, " =~ 0*",label2, i, sep = "")  
          if ((k-1)-i > 0) {
            linemidb <- ""
            for (j in (i+1):k) {
              linemidb <- paste(linemidb, " + 0*", label2, j, sep = "")
            }
          } else {linemidb <- ""}
          
        }
        Model1b <- paste(Model1b, linestartb, linemidb, " \n ", sep = "")
      }
    }
    else {Model1b <- ""}
    
    Model2<-""
    for (p in 1:k) {
      linestart2 <- paste(label2, p, " =~ 1*", label, p, sep = "")
      Model2<-paste(Model2, linestart2, " \n ", sep = "")}
    
    Model3<-""
    for (p in 1:k) {
      linestart3 <- paste(label, p, "~~", label, p, sep = "")
      Model3<-paste(Model3, linestart3, " \n ", sep = "")}
    
    Model4<-""
    for (p in 1:k) {
      linestart4 <- paste(label2, p, "~~0*", label2, p, sep = "")
      Model4<-paste(Model4, linestart4, " \n ", sep = "")}
    
    Modelsat<-""
    for (i in 1:(k-1)) {
      linestartc <- paste("", label, i, "~~0*", label, i+1, sep = "")
      if (k-i >= 2) {
        linemidc <- ""
        for (j in (i+2):k) {
          linemidc <- paste("", linemidc, label, i, "~~0*", label, j, " \n ", sep="")
          
        }
      } else {linemidc <- ""}
      Modelsat <- paste(Modelsat, linestartc, " \n ", linemidc, sep = "")
    } 
    
    Model5<-paste(model, " \n ", ModelsatF, Model1b, Model2, Model3, Model4, Modelsat, sep = "")
    
    return(Model5)
  } 
  
  Model1<-write.Model1(k)
  
  ##function to remove duplicated elements between user/automatically specified Model Input
  tryCatch.W.E <- function(expr)
  {
    W <- NULL
    w.handler <- function(w){ # warning handler
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler),
         warning = W)
  }
  
  
  while(class(tryCatch.W.E(lavParseModelString(Model1))$value$message) != 'NULL'){
    u<-tryCatch.W.E(lavParseModelString(Model1))$value$message
    t<-paste(strsplit(u, ": ")[[1]][3], " \n ", sep = "")
    Model1<-str_replace(Model1, fixed(t), "")
  }
  
 if(CFIcalc==TRUE){
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
      linestartc <- paste(label, i, " ~~ 0*", label, i+1,  sep = "")
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
      linestartc <- paste(" ", label2, i, " ~~ 0*", label2, i+1,  sep = "")
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
  
  ##create inependence model for calculation of CFI
  modelCFI<-write.null(k)
  
  ##run CFI model so it knows the reordering for the independence model
  fitCFI <- sem(modelCFI, sample.cov = S_LD, estimator = "DWLS", WLS.V = W,sample.nobs=2)
  orderCFI <- rearrange(k = k, fit =  fitCFI, names =  rownames(S_LD))
  
  ##reorder matrix for independence (i.e., null) model for CFI calculation
  V_Reorder2 <- V_LD[orderCFI,orderCFI]
  V_Reorder2b<-diag(z)
  diag(V_Reorder2b)<-diag(V_Reorder2)
  W_CFI<-solve(V_Reorder2b)
  }
  
  ##code to write saturated model to check there are no redundancies
  ##with user provided model in later part of script
  write.test<-function(k, label = "V", label2 = "VF") {
    
    Modelsat<-"" 
    for (i in 1:(k)) {
      if (k-i >= 1) {
        linemidc <- ""
        for (j in (i+1):k) {
          linemidc <- paste(linemidc, label, i, "~~", label, j, " \n ", sep = "")
        }
      }else{linemidc<-""} 
      Modelsat <- paste(Modelsat, linemidc, sep = "")
    }
    
    Model4<-""
    for (p in 1:k) {
      linestart4 <- paste(label2, p, "~~", label2, p, sep = "")
      Model4<-paste(Model4, linestart4, " \n ", sep = "")}
    
    modelCFI<-paste(Modelsat, Model4)
    return(modelCFI)
  }
  
  modeltest<-data.frame(write.test(k))
  modeltest$write.test.k.<-as.character(modeltest$write.test.k.)
  modeltest2 <- cSplit(modeltest, "write.test.k.", sep = "\n", direction = "long") 
  modeltest2$write.test.k.<-as.character(modeltest2$write.test.k.)
  
  ReorderModel <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W, sample.nobs = 2) 
  
  ##save the ordering
  order <- rearrange(k = k, fit = ReorderModel, names = rownames(S_LD))
  
  ##reorder the weight (inverted V_LD) matrix
  V_Reorder<-V_LD[order,order]
  V_Reorderb<-diag(z)
  diag(V_Reorderb)<-diag(V_Reorder)
  W_Reorder<-solve(V_Reorderb)

  ##estimation for DWLS
  if(estimation=="DWLS"){
    
    print("Running primary model")
    ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
    Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2)
    
    #pull the delta matrix (this doesn't depend on N)
    ##note that while the delta matrix is reordered based on the ordering in the model specification
    ##that the lavaan output is also reordered so that this actually ensures that the results match up 
    S2.delt <- lavInspect(Model1_Results, "delta")
    
    ##weight matrix from stage 2. S2.W is not reordered by including something like model constraints
    S2.W <- lavInspect(Model1_Results, "WLS.V") 
    
    #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
    bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
    
    #create the "lettuce" part of the sandwich
    lettuce <- S2.W%*%S2.delt
    
    #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
    Ohtt <- bread %*% t(lettuce)%*%V_Reorder%*%lettuce%*%bread  
    
    #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
    SE <- as.matrix(sqrt(diag(Ohtt)))
    
    ModelQ_WLS <- parTable(Model1_Results)
    
    ##remove any parameter constraint labels
    ModelQ_WLS<-subset(ModelQ_WLS, ModelQ_WLS$plabel != "")
    
    for (i in 1:length(ModelQ_WLS)){
      if(((paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% modeltest2$write.test.k)) & ModelQ_WLS$est[i] == 0)
      {ModelQ_WLS$free[i] == 1} else{ModelQ_WLS$free[i] == 0} 
    }
    
    for (i in 1:nrow(ModelQ_WLS)){
      ModelQ_WLS$free[i]<-ifelse((paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% modeltest2$write.test.k) & ModelQ_WLS$est[i] == 0, 1, 0)
    }
    
    for (i in 1:nrow(ModelQ_WLS)){
      ModelQ_WLS$free[i]<-ifelse((paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% modeltest2$write.test.k) & ModelQ_WLS$est[i] != 0, 2, ModelQ_WLS$free[i])
    }
    
    test<-vector(mode="list",length=nrow(ModelQ_WLS))
    
    for(i in 1:nrow(ModelQ_WLS)){
      if(ModelQ_WLS$free[i] == 2) { 
        t<-paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "")
        t2<-gsub("V", "VF", t)
        test[[i]]<-t2}else{}
    }
    test2<-Filter(Negate(is.null), test)
    
    for (i in 1:nrow(ModelQ_WLS)){
      ModelQ_WLS$free[i]<-ifelse((paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% test2), 1, ModelQ_WLS$free[i])
    }
    
    ModelQ_WLS$free<-ifelse(ModelQ_WLS$free != 1, 0, ModelQ_WLS$free)
    
    #want to freely estimate the residual factor variances and the residual covariances
    z<-(k*(k+1))/2
    
    p<-length(ModelQ_WLS$free)-z
    
    ModelQ_WLS <- ModelQ_WLS[order(ModelQ_WLS$free),] 
    
    ModelQ_WLS$free <- c(rep(0, p),1:z)
    
    ModelQ_WLS$ustart <- ModelQ_WLS$est
    ModelQ_WLS$ustart<-ifelse(ModelQ_WLS$free > 0, .05, ModelQ_WLS$ustart)
    
    print("Calculating model chi-square")
    testQ<-tryCatch.W.E(ModelQ_Results_WLS <- sem(model = ModelQ_WLS, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs=2, start = ModelQ_WLS$ustart)) 
    testQ$warning$message[1]<-ifelse(is.null(testQ$warning$message), testQ$warning$message[1]<-"Safe", testQ$warning$message[1])
    testQ$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_WLS, "se")$theta[1,2]) == TRUE, testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ$warning$message[1])
    
    if(as.character(testQ$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
      
      ModelQ_WLS$ustart<-ifelse(ModelQ_WLS$free > 0, .01, ModelQ_WLS$ustart)
      
      testQ2<-tryCatch.W.E(ModelQ_Results_WLS <- sem(model = ModelQ_WLS, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2, start=ModelQ$ustart)) 
    }else{testQ2<-testQ}
    
    testQ2$warning$message[1]<-ifelse(is.null(testQ2$warning$message), testQ2$warning$message[1]<-"Safe", testQ2$warning$message[1])
    testQ2$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_WLS, "se")$theta[1,2]) == TRUE, testQ2$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ2$warning$message[1])
    
    if(as.character(testQ2$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
      
      ModelQ_WLS$ustart<-ifelse(ModelQ_WLS$free > 0, .1, ModelQ_WLS$ustart)
      
      testQ3<-tryCatch.W.E(ModelQ_Results_WLS <- sem(model = ModelQ_WLS, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2, start=ModelQ$ustart)) 
    }else{testQ3<-testQ2}
    
    testQ3$warning$message[1]<-ifelse(is.null(testQ3$warning$message), testQ3$warning$message[1]<-"Safe", testQ3$warning$message[1])
    testQ3$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_WLS, "se")$theta[1,2]) == TRUE, testQ3$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ3$warning$message[1])
    
    if(as.character(testQ3$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
      
      #pull the delta matrix (this doesn't depend on N)
      S2.delt_Q <- lavInspect(ModelQ_Results_WLS, "delta")
      
      ##weight matrix from stage 2
      S2.W_Q <- lavInspect(ModelQ_Results_WLS, "WLS.V") 
      
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
      eta_test<-parTable(ModelQ_Results_WLS)
      eta_test<-subset(eta_test, eta_test$free != 0)
      eta<-cbind(eta_test[,14])
      
      #Ronald's magic combining all the pieces from above:
      Q_WLS<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta}else{Q_WLS<-NA}
    
    if(CFIcalc == TRUE){
    print("Calculating CFI")
    ##now CFI
    ##run independence model
    testCFI<-tryCatch.W.E(fitCFI <- sem(modelCFI, sample.cov =  S_LD, estimator = "DWLS", WLS.V = W_CFI, sample.nobs=2))
    testCFI$warning$message[1]<-ifelse(is.null(testCFI$warning$message), testCFI$warning$message[1]<-"Safe", testCFI$warning$message[1])
    testCFI$warning$message[1]<-ifelse(is.na(inspect(fitCFI, "se")$theta[1,2]) == TRUE, testCFI$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testCFI$warning$message[1])
    
    if(as.character(testCFI$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
      
      ##code to estimate chi-square of independence model#
      #First pull the estimates from Step 2
      ModelQ_WLS_CFI <- parTable(fitCFI)
      p2<-length(ModelQ_WLS_CFI$free)-z
      
      ##fix variances and freely estimate covariances
      ModelQ_WLS_CFI$free <- c(rep(0, p2), 1:z)
      ModelQ_WLS_CFI$ustart <- ModelQ_WLS_CFI$est
      
      testCFI2<-tryCatch.W.E(ModelQ_Results_WLS_CFI <- sem(model = ModelQ_WLS_CFI, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_CFI, sample.nobs=2))
      testCFI2$warning$message[1]<-ifelse(is.null(testCFI2$warning$message), testCFI2$warning$message[1]<-"Safe", testCFI2$warning$message[1])
      testCFI2$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_WLS_CFI , "se")$theta[1,2]) == TRUE, testCFI2$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testCFI2$warning$message[1])
      
      if(as.character(testCFI2$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
        
        #pull the delta matrix (this doesn't depend on N)
        S2.delt_Q_CFI <- lavInspect(ModelQ_Results_WLS_CFI, "delta")
        
        ##weight matrix from stage 2
        S2.W_Q_CFI <- lavInspect(ModelQ_Results_WLS_CFI, "WLS.V") 
        
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
        eta_test_CFI<-parTable(ModelQ_Results_WLS_CFI)
        eta_test_CFI<-subset(eta_test_CFI, eta_test_CFI$free != 0)
        eta_CFI<-cbind(eta_test_CFI[,14])
        
        #Ronald's magic combining all the pieces from above:
        Q_CFI_WLS<-t(eta_CFI)%*%P1_CFI%*%solve(Eig_CFI)%*%t(P1_CFI)%*%eta_CFI}else{Q_CFI_WLS<-NA}}
    
    ##df of independence Model
          ## MICHEL ADDED THE df line as a TEMP fix 
      dfCFI <- (((k * (k + 1))/2) - k)
      df <- (k * (k + 1)/2) - max(parTable(Model1_Results)$free)
      
      
    if(!(is.na(Q_CFI_WLS)) & !(is.na(Q_WLS))){
      CFI<-as.numeric(((Q_CFI_WLS-dfCFI)-(Q_WLS-df))/(Q_CFI_WLS-dfCFI))
      CFI<-ifelse(CFI > 1, 1, CFI)
    }else{CFI<-NA}
    
     }
    
    print("Calculating Standardized Results")
    ##transform the S covariance matrix to S correlation matrix
    D=sqrt(diag(diag(S_LD)))
    S_Stand=solve(D)%*%S_LD%*%solve(D)
    rownames(S_Stand)<-rownames(S_LD)
    colnames(S_Stand)<-colnames(S_Stand)
    
    #obtain diagonals of the original V matrix and take their sqrt to get SE's
    Dvcov<-sqrt(diag(V_LD))
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO=as.vector(lowerTriangle((S_Stand/S_LD),diag=T))
    
    ## MAke sure that if ratio in NaN (devision by zero) we put the zero back in: ### TEMP STUPID MICHEL FIX!
    scaleO[is.nan(scaleO)] <- 0
    
    #rescale the SEs by the same multiples that the S matrix was rescaled by
    Dvcovl<-as.vector(Dvcov*t(scaleO))
    
    #obtain the sampling correlation matrix by standardizing the original V matrix
    Vcor<-cov2cor(V_LD)
    
    #rescale the sampling correlation matrix by the appropriate diagonals
    V_stand<-diag(Dvcovl)%*%Vcor%*%diag(Dvcovl)
    V_stand2<-diag(z)
    diag(V_stand2)<-diag(V_stand)
    
    ### make sure no value on the diagonal of V is 0 ### TEMP STUPID MICHEL FIX
    diag(V_stand2)[diag(V_stand2) == 0] <- 2e-9
    
    W_stand<-solve(V_stand2[order,order])
  
    DWLS.fit_stand <- sem(Model1, sample.cov = S_Stand, estimator = "DWLS", WLS.V = W_stand, sample.nobs = 2) 
    
    ##perform same procedures for sandwich correction as in the unstandardized case
    DWLS.delt_stand <- lavInspect(DWLS.fit_stand, "delta") 
    DWLS.W_stand <- lavInspect(DWLS.fit_stand, "WLS.V") 
    bread_stand <- solve(t(DWLS.delt_stand)%*%DWLS.W_stand %*%DWLS.delt_stand) 
    lettuce_stand <- DWLS.W_stand%*%DWLS.delt_stand
    Vcov_stand<-as.matrix(V_stand[order,order])
    Ohtt_stand <- bread_stand %*% t(lettuce_stand)%*%Vcov_stand%*%lettuce_stand%*%bread_stand
    SE_stand <- as.matrix(sqrt(diag(Ohtt_stand)))
    
    unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
    unstand<-subset(unstand, unstand$free != 0)                    
    unstand$free<-NULL
    
    stand<-data.frame(inspect(DWLS.fit_stand,"list")[,c(8,14)])
    stand<-subset(stand, stand$free != 0)
    stand$free<-NULL
    
    ##df of user Model
    df<-(k*(k+1)/2)-max(parTable(Model1_Results)$free)
    
    if(!(is.na(Q_WLS))){
      chisq<-Q_WLS
      AIC<-(Q_WLS + 2*max(parTable(Model1_Results)$free))}else{chisq<-NA
      AIC<-NA}
    
    print("Calculating SRMR")
    SRMR<-lavInspect(Model1_Results, "fit")["srmr"]
    
    if(CFIcalc == TRUE){
    modelfit<-cbind(chisq,df,AIC,CFI,SRMR)}else{modelfit<-cbind(chisq,df,AIC,SRMR)}
    
    results<-cbind(unstand,SE,stand,SE_stand)
    
  }
  
  ##ML estimation
  if(estimation=="ML"){
    
    print("Running primary model")
    ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
    Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "ML", sample.nobs = 200)
    
    #pull the delta matrix (this doesn't depend on N)
    S2.delt <- lavInspect(Model1_Results, "delta")
    
    ##weight matrix from stage 2. S2.W is not reordered by including something like model constraints
    S2.W <- lavInspect(Model1_Results, "WLS.V") 
    
    #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
    bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
    
    #create the "lettuce" part of the sandwich
    lettuce <- S2.W%*%S2.delt
    
    #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
    Ohtt <- bread %*% t(lettuce)%*%V_Reorder%*%lettuce%*%bread  
    
    #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
    SE <- as.matrix(sqrt(diag(Ohtt)))
    
    ModelQ_ML <- parTable(Model1_Results)
    
    ##remove any parameter constraint labels
    ModelQ_ML<-subset(ModelQ_ML, ModelQ_ML$plabel != "")
    
    for (i in 1:length(ModelQ_ML)){
      if(((paste(ModelQ_ML$lhs[i], ModelQ_ML$op[i], ModelQ_ML$rhs[i], sep = "") %in% modeltest2$write.test.k)) & ModelQ_ML$est[i] == 0)
      {ModelQ_ML$free[i] == 1} else{ModelQ_ML$free[i] == 0} 
    }
    
    for (i in 1:nrow(ModelQ_ML)){
      ModelQ_ML$free[i]<-ifelse((paste(ModelQ_ML$lhs[i], ModelQ_ML$op[i], ModelQ_ML$rhs[i], sep = "") %in% modeltest2$write.test.k) & ModelQ_ML$est[i] == 0, 1, 0)
    }
    
    for (i in 1:nrow(ModelQ_ML)){
      ModelQ_ML$free[i]<-ifelse((paste(ModelQ_ML$lhs[i], ModelQ_ML$op[i], ModelQ_ML$rhs[i], sep = "") %in% modeltest2$write.test.k) & ModelQ_ML$est[i] != 0, 2, ModelQ_ML$free[i])
    }
    
    test<-vector(mode="list",length=nrow(ModelQ_ML))
    
    for(i in 1:nrow(ModelQ_ML)){
      if(ModelQ_ML$free[i] == 2) { 
        t<-paste(ModelQ_ML$lhs[i], ModelQ_ML$op[i], ModelQ_ML$rhs[i], sep = "")
        t2<-gsub("V", "VF", t)
        test[[i]]<-t2}else{}
    }
    test2<-Filter(Negate(is.null), test)
    
    for (i in 1:nrow(ModelQ_ML)){
      ModelQ_ML$free[i]<-ifelse((paste(ModelQ_ML$lhs[i], ModelQ_ML$op[i], ModelQ_ML$rhs[i], sep = "") %in% test2), 1, ModelQ_ML$free[i])
    }
    
    ModelQ_ML$free<-ifelse(ModelQ_ML$free != 1, 0, ModelQ_ML$free)
    
    #want to freely estimate the residual factor variances and the residual covariances
    z<-(k*(k+1))/2
    
    p<-length(ModelQ_ML$free)-z
    
    ModelQ_ML <- ModelQ_ML[order(ModelQ_ML$free),] 
    
    ModelQ_ML$free <- c(rep(0, p),1:z)
    
    ModelQ_ML$ustart <- ModelQ_ML$est
    
    ModelQ_ML$ustart<-ifelse(ModelQ_ML$free > 0, .05, ModelQ_ML$ustart)
    
    print("Calculating model chi-square")
    testQ<-tryCatch.W.E(ModelQ_Results_ML <- sem(model = ModelQ_ML, sample.cov = S_LD, estimator = "ML", sample.nobs=200, start =  ModelQ_ML$ustart)) 
    testQ$warning$message[1]<-ifelse(is.null(testQ$warning$message), testQ$warning$message[1]<-"Safe", testQ$warning$message[1])
    
    testQ$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_ML, "se")$theta[1,2]) == TRUE, testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ$warning$message[1])
    
    if(as.character(testQ$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
      
      ModelQ_ML$ustart<-ifelse(ModelQ_ML$free > 0, .01, ModelQ_ML$ustart)
      
      testQ2<-tryCatch.W.E(ModelQ_Results_ML <- sem(model = ModelQ_ML, sample.cov = S_LD, estimator = "ML", sample.nobs=200, start =  ModelQ_ML$ustart)) 
    }else{testQ2<-testQ}
    
    
    testQ2$warning$message[1]<-ifelse(is.null(testQ2$warning$message), testQ2$warning$message[1]<-"Safe", testQ2$warning$message[1])
    testQ2$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_ML, "se")$theta[1,2]) == TRUE, testQ2$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ2$warning$message[1])
    
    if(as.character(testQ2$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
      
      ModelQ_ML$ustart<-ifelse(ModelQ_ML$free > 0, .1, ModelQ_ML$ustart)
      
      testQ3<-tryCatch.W.E(ModelQ_Results_ML <- sem(model = ModelQ_ML, sample.cov = S_LD, estimator = "ML", sample.nobs=200, start =  ModelQ_ML$ustart)) 
    }else{testQ3<-testQ2}
    
    testQ3$warning$message[1]<-ifelse(is.null(testQ3$warning$message), testQ3$warning$message[1]<-"Safe", testQ3$warning$message[1])
    testQ3$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_ML, "se")$theta[1,2]) == TRUE, testQ3$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ3$warning$message[1])
    
    if(as.character(testQ3$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
      
      #pull the delta matrix (this doesn't depend on N)
      S2.delt_Q <- lavInspect(ModelQ_Results_ML, "delta")
      
      ##weight matrix from stage 2
      S2.W_Q <- lavInspect(ModelQ_Results_ML, "WLS.V") 
      
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
      eta_test<-parTable(ModelQ_Results_ML)
      eta_test<-subset(eta_test, eta_test$free != 0)
      eta<-cbind(eta_test[,14])
      
      #Combining all the pieces from above:
      Q_ML<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta}else{Q_ML<-NA}
    
    if(CFIcalc == TRUE){
    print("Calculating CFI")
    ##now CFI
    ##run independence model
    testCFI<-tryCatch.W.E(fitCFI <- sem(modelCFI, sample.cov =  S_LD, estimator = "ML", sample.nobs=200))
    testCFI$warning$message[1]<-ifelse(is.null(testCFI$warning$message), testCFI$warning$message[1]<-"Safe", testCFI$warning$message[1])
    testCFI$warning$message[1]<-ifelse(is.na(inspect(fitCFI, "se")$theta[1,2]) == TRUE, testCFI$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testCFI$warning$message[1])
    
    if(as.character(testCFI$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
      
      ##code to estimate chi-square of independence model#
      #First pull the estimates from Step 2
      ModelQ_ML_CFI <- parTable(fitCFI)
      p2<-length(ModelQ_ML_CFI$free)-z
      
      ##fix variances and freely estimate covariances
      ModelQ_ML_CFI$free <- c(rep(0, p2), 1:z)
      ModelQ_ML_CFI$ustart <- ModelQ_ML_CFI$est
      
      testCFI2<-tryCatch.W.E(ModelQ_Results_ML_CFI <- sem(model = ModelQ_ML_CFI, sample.cov = S_LD, estimator = "ML", sample.nobs=200))
      
      testCFI2$warning$message[1]<-ifelse(is.null(testCFI2$warning$message), testCFI2$warning$message[1]<-"Safe", testCFI2$warning$message[1])
      testCFI2$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_ML_CFI , "se")$theta[1,2]) == TRUE, testCFI2$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testCFI2$warning$message[1])
      
      if(as.character(testCFI2$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
        
        #pull the delta matrix (this doesn't depend on N)
        S2.delt_Q_CFI <- lavInspect(ModelQ_Results_ML_CFI, "delta")
        
        ##weight matrix from stage 2
        S2.W_Q_CFI <- lavInspect(ModelQ_Results_ML_CFI, "WLS.V") 
        
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
        eta_test_CFI<-parTable(ModelQ_Results_ML_CFI)
        eta_test_CFI<-subset(eta_test_CFI, eta_test_CFI$free != 0)
        eta_CFI<-cbind(eta_test_CFI[,14])
        
        #Ronald's magic combining all the pieces from above:
        Q_CFI_ML<-t(eta_CFI)%*%P1_CFI%*%solve(Eig_CFI)%*%t(P1_CFI)%*%eta_CFI}else{Q_CFI_ML<-NA}}
    
    ##df of independence Model
    dfCFI<-(((k*(k+1))/2)-k)
    
    if(!(is.na(Q_CFI_ML)) & !(is.na(Q_ML))){
      CFI<-as.numeric(((Q_CFI_ML-dfCFI)-(Q_ML-df))/(Q_CFI_ML-dfCFI))
      CFI<-ifelse(CFI > 1, 1, CFI)
    }else{CFI<-NA}
    
    }
    
    ##transform the S covariance matrix to S correlation matrix
    D=sqrt(diag(diag(S_LD)))
    S_Stand=solve(D)%*%S_LD%*%solve(D)
    rownames(S_Stand)<-rownames(S_LD)
    colnames(S_Stand)<-colnames(S_Stand)
    
    #obtain diagonals of the original V matrix and take their sqrt to get SE's
    Dvcov<-sqrt(diag(V_LD))
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO=as.vector(lowerTriangle((S_Stand/S_LD),diag=T))
    
    ## MAke sure that if ratio in NaN (division by zero) we put the zero back in: ### TEMP STUPID MICHEL FIX!
    scaleO[is.nan(scaleO)] <- 0
    
    #rescale the SEs by the same multiples that the S matrix was rescaled by
    Dvcovl<-as.vector(Dvcov*t(scaleO))
    
    #obtain the sampling correlation matrix by standardizing the original V matrix
    Vcor<-cov2cor(V_LD)
    
    #rescale the sampling correlation matrix by the appropriate diagonals
    V_stand<-diag(Dvcovl)%*%Vcor%*%diag(Dvcovl)
    V_stand2<-diag(z)
    diag(V_stand2)<-diag(V_stand)
    
    ### make sure no value on the diagonal of V is 0 ### TEMP STUPID MICHEL FIX
    diag(V_stand2)[diag(V_stand2) == 0] <- 2e-9
    
    W_stand<-solve(V_stand2[order,order])
    
    print("Calculating Standardized Results")
    ML.fit_stand <- sem(Model1, sample.cov = S_Stand, estimator = "ML", sample.nobs = 200) 
    
    ##perform same procedures for sandwich correction as in the unstandardized case
    ML.delt_stand <- lavInspect(ML.fit_stand, "delta") 
    ML.W_stand <- lavInspect(ML.fit_stand, "WLS.V") 
    bread_stand <- solve(t(ML.delt_stand)%*%ML.W_stand %*%ML.delt_stand) 
    lettuce_stand <- ML.W_stand%*%ML.delt_stand
    Vcov_stand<-as.matrix(V_stand[order,order])
    Ohtt_stand <- bread_stand %*% t(lettuce_stand)%*%Vcov_stand%*%lettuce_stand%*%bread_stand
    SE_stand <- as.matrix(sqrt(diag(Ohtt_stand)))
    
    unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
    unstand<-subset(unstand, unstand$free != 0)                    
    unstand$free<-NULL
    
    stand<-data.frame(inspect(ML.fit_stand,"list")[,c(8,14)])
    stand<-subset(stand, stand$free != 0)
    stand$free<-NULL
    
    ##df of user model
    df<-(k*(k+1)/2)-max(parTable(Model1_Results)$free)
    
    if(!(is.na(Q_ML))){
      chisq<-Q_ML
      AIC<-(Q_ML + 2*max(parTable(Model1_Results)$free))}else{chisq<-NA
      AIC<-NA}
    
    print("Calculating SRMR")
    SRMR<-lavInspect(Model1_Results, "fit")["srmr"]
    
    if(CFIcalc == TRUE){
      modelfit<-cbind(chisq,df,AIC,CFI,SRMR)}else{modelfit<-cbind(chisq,df,AIC,SRMR)}
    
    results<-cbind(unstand,SE,stand,SE_stand)
    
  }
  
  ##name the columns of the results file
  colnames(results)=c("lhs","op","rhs","Unstandardized_Estimate","Unstandardized_SE","Standardized_Est","Standardized_SE")
  
  ##replace V1-VX general form in output with user provided trait names
  for(i in 1:nrow(results)){
    for(p in 1:length(traits)){
      results$lhs[[i]]<-ifelse(results$lhs[[i]] %in% S_names[[p]], gsub(results$lhs[[i]], traits[[p]], results$lhs[[i]]), results$lhs[[i]])
      results$rhs[[i]]<-ifelse(results$rhs[[i]] %in% S_names[[p]], gsub(results$rhs[[i]], traits[[p]], results$rhs[[i]]), results$rhs[[i]])
    }
  }
  
  ##name model fit columns
  if(CFIcalc == TRUE){
  colnames(modelfit)=c("chisq","df","AIC","CFI","SRMR")}else{colnames(modelfit)=c("chisq","df","AIC","SRMR")}
  
  modelfit<-data.frame(modelfit)
  
  modelfit$p_chisq<-ifelse(!(is.na(modelfit$chisq)), modelfit$p_chisq<-pchisq(modelfit$chisq, modelfit$df,lower.tail=FALSE), modelfit$p_chisq<-NA)
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
  
  
  return(list(modelfit=modelfit,results=results))
  
}
