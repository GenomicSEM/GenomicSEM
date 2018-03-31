
commonfactorGWAS <-function(Output,estimation="DWLS"){ 
  time<-proc.time()
  
  ##split the V and S matrices into as many (cores - 1) as are aviailable on the local computer
  V_Full<-(Output[[1]])
  S_Full<-(Output[[2]])
  
  #enter in k for number of phenotypes 
  k<-ncol(S_Full[[1]])-1
  
  f<-length(Output[[1]])
  
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
  
  #function to create lavaan syntax for a 1 factor model given k phenotypes
  write.Model1 <- function(k, label = "V") {  
    Model1 <- ""
    for (i in 1) {
      lineSNP <- paste(label, i, " ~ 0*SNP",sep = "")
      if (k-i > 0) {
        lineSNP2 <- " \n "
        for (j in (i+1):k) {
          lineSNP2 <- paste(lineSNP2, label, j, " ~ 0*SNP", " \n ", sep = "")
        }
      }
    } 
    
    for (i in 1) {
      linestart <- paste("F1"," =~ ",label, i, sep = "")  
      if (k-i > 0) {
        linemid <- ""
        for (j in (i+1):k) {
          linemid <- paste(linemid, " + ", label, j, sep = "")
        }
      } else {linemid <- ""}
    }
    
    Model1 <- paste(Model1, linestart, linemid, " \n ", "F1 ~ SNP", " \n ", lineSNP, lineSNP2, sep = "")
    return(Model1)
  } 
  
  ##create the model
  Model1 <- write.Model1(k)
  
  ##modification of trycatch that allows the results of a failed run to still be saved
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
  
  ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  for (i in 1) {
    
    #transform sampling covariance matrix into a weight matrix: 
    W <- solve(V_Full[[i]])
    
    S_Fullrun<-S_Full[[i]]
    
    ReorderModel <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2) 
    
    order <- rearrange(k = k+1, fit = ReorderModel, names = rownames(S_Full[[1]]))
  }
  
  ##name the columns of the results file
  results=as.data.frame(matrix(NA,ncol=10,nrow=f))
  colnames(results)=c("i","lhs","op","rhs","est","se", "se_c", "Q", "fail", "warning")
  
  ##estimation for 2S-DWLS-R
  if(estimation=="DWLS"){
  
    for (i in 1:f) { 
      
      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Full_Reorder <- V_Full[[i]][order,order]
      u<-nrow(V_Full_Reorder)
      V_Full_Reorderb<-diag(u)
      diag(V_Full_Reorderb)<-diag(V_Full_Reorder)
      
      ##invert the reordered sampling covariance matrix to create a weight matrix 
      W <- solve(V_Full_Reorderb) 
      
      #import the S_Full matrix for appropriate run
      S_Fullrun<-S_Full[[i]]
      
      ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
      test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2))
      
      #pull the delta matrix (this doesn't depend on N)
      S2.delt <- lavInspect(Model1_Results, "delta")
      
      ##weight matrix from stage 2
      S2.W <- lavInspect(Model1_Results, "WLS.V") 
      
      #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
      bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
      
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
      
      #fix the indicator loadings from Step 1, free the direct effects of the SNP on the indicators, and fix the factor residual variance
      ModelQ$free <- c(rep(0, k+1), 1:(k*2), 0, 0) 
      
      ##added##
      ModelQ$ustart <- ModelQ$est
      SNPresid<-resid(Model1_Results)$cov[k+1,1:k]
      
      for(t in 1:nrow(ModelQ)) {
        if(ModelQ$free[t] > 0 & ModelQ$free[t] <= k){
          ModelQ$ustart[t]<-SNPresid[ModelQ$free[t]]} else{}}

      #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
      testQ<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, start=ModelQ$ustart)) 
      testQ$warning$message[1]<-ifelse(is.null(testQ$warning$message), testQ$warning$message[1]<-"Safe", testQ$warning$message[1])
      
      if(as.character(testQ$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
        
        for(n in 1:nrow(ModelQ)) {
          if(ModelQ$free[n] > k & ModelQ$free[n] <= k*2){
            ModelQ$ustart[n]<-ModelQ$ustart[n]-.01} else{}}
        
        testQ2<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, start=ModelQ$ustart)) 
      }else{testQ2<-testQ}
      
      testQ2$warning$message[1]<-ifelse(is.null(testQ2$warning$message), testQ2$warning$message[1]<-"Safe", testQ2$warning$message[1])
      
      if(as.character(testQ2$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
        
        ModelQ$ustart<-ifelse(ModelQ$free > k, .05, ModelQ$ustart)
        
        testQ3<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, start=ModelQ$ustart)) 
      }else{testQ3<-testQ2}
      
      testQ3$warning$message[1]<-ifelse(is.null(testQ3$warning$message), testQ3$warning$message[1]<-"Safe", testQ3$warning$message[1])
      
      if(as.character(testQ3$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
        #pull the delta matrix for Q (this doesn't depend on N)
        S2.delt_Q <- lavInspect(ModelQ_Results, "delta")
        
        ##weight matrix from stage 2 for Q
        S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V") 
        
        #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
        bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q) 
        
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
        
        #Ronald's magic combining all the pieces from above:
        Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta} else{Q<-"Not Computed"}
      
      ##pull all the results into a single row
      results[i,]<-cbind(i,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
      
    }
  }

  ##2S-ML-R estimation
  if(estimation=="ML"){

    for (i in 1:f) { 
      
      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Full_Reorder <- V_Full[[i]][order,order]
      u<-nrow(V_Full_Reorder)
      V_Full_Reorderb<-diag(u)
      diag(V_Full_Reorderb)<-diag(V_Full_Reorder)
      
      ##invert the reordered sampling covariance matrix to create a weight matrix 
      W <- solve(V_Full_Reorderb) 
      
      #import the S_Full matrix for appropriate run
      S_Fullrun<-S_Full[[i]]
      
      ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
      test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML", sample.nobs = 200))
      
      #pull the delta matrix (this doesn't depend on N)
      S2.delt <- lavInspect(Model1_Results, "delta")
      
      ##weight matrix from stage 2
      S2.W <- lavInspect(Model1_Results, "WLS.V") 
      
      #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
      bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
      
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
      
      #fix the indicator loadings from Step 1, free the direct effects of the SNP on the indicators, and fix the factor residual variance
      ModelQ$free <- c(rep(0, k+1), 1:(k*2), 0, 0) 

      ##added##
      ModelQ$ustart <- ModelQ$est
      SNPresid<-resid(Model1_Results)$cov[k+1,1:k]

      for(t in 1:nrow(ModelQ)) {
        if(ModelQ$free[t] > 0 & ModelQ$free[t] <= k){
          ModelQ$ustart[t]<-SNPresid[ModelQ$free[t]]} else{}}
      
    #  for(n in 1:nrow(ModelQ)) {
     #   if(ModelQ$free[n] > k & ModelQ$free[n] <= k*2){
    #      for(y in 1:k){
    #      ModelQ$ustart[n]<-VARresid[y]}} else{}}
  
      
      #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
      testQ<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "ML", sample.nobs=200,start=ModelQ$ustart)) 
      testQ$warning$message[1]<-ifelse(is.null(testQ$warning$message), testQ$warning$message[1]<-"Safe", testQ$warning$message[1])
      
       if(as.character(testQ$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
      
      for(n in 1:nrow(ModelQ)) {
        if(ModelQ$free[n] > k & ModelQ$free[n] <= k*2){
          ModelQ$ustart[n]<-ModelQ$ustart[n]-.01} else{}}
      
      testQ2<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "ML", sample.nobs=200,start=ModelQ$ustart)) 
      }else{testQ2<-testQ}
     
      testQ2$warning$message[1]<-ifelse(is.null(testQ2$warning$message), testQ2$warning$message[1]<-"Safe", testQ2$warning$message[1])
      
      if(as.character(testQ2$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
        
        ModelQ$ustart<-ifelse(ModelQ$free > k, .05, ModelQ$ustart)
        
        testQ3<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "ML", sample.nobs=200,start=ModelQ$ustart)) 
      }else{testQ3<-testQ2}
      
      testQ3$warning$message[1]<-ifelse(is.null(testQ3$warning$message), testQ3$warning$message[1]<-"Safe", testQ3$warning$message[1])
      
      if(as.character(testQ3$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
      #pull the delta matrix for Q (this doesn't depend on N)
      S2.delt_Q <- lavInspect(ModelQ_Results, "delta")
      
      ##weight matrix from stage 2 for Q
      S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V") 
      
      #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
      bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q) 
      
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
      
      #Ronald's magic combining all the pieces from above:
      Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta} else{Q<-"Not Computed"}
      
      
      ##put the corrected standard error and Q in same dataset
      results[i,]<-cbind(i,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
      
      
    }
  }
  
  time_all<-proc.time()-time
  print(time_all[3])
  results$se <- NULL
  results2<-cbind(Output[[3]],results)
  
  return(list(results2))
  
}
