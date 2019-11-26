
commonfactorGWASpar <-function(covstruc,SNPs,estimation="DWLS",cores=NULL,toler=FALSE,SNPSE=FALSE,Output=NULL){ 
  time<-proc.time()
  
  print("Please note that an update was made to commonfactorGWASpar on 11/21/19 so that it combines addSNPs and commonfactorGWASpar.")
  
  
  if(class(covstruc[[1]]) != "matrix"){
    print("You are likely listing arguments in the order of a previous version of commonfactorGWAS. The current version of the function expects the 
          output from ldsc followed by the output from sumstats as the first two arguments.")    
    warning("You are likely listing arguments in the order of a previous version of commonfactorGWAS. The current version of the function expects the 
            output from ldsc followed by the output from sumstats as the first two arguments.")
  }
  
  
  if(is.null(cores)){
    ##if no default provided use 1 less than the total number of cores available so your computer will still function
    int <- detectCores() - 1
  }else{int<-cores}
  
  registerDoParallel(int)
  
  ##specify the cores should have access the local environment
  makeCluster(int, type="FORK")
  
  if(is.null(Output)){
  
  ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
  SNPs<-data.frame(SNPs)
  SNPs$A1<-as.character(SNPs$A1)
  SNPs$A2<-as.character(SNPs$A2)
  SNPs$SNP<-as.character(SNPs$SNP)
  
  #SNP variance
  varSNP=2*SNPs$MAF*(1-SNPs$MAF)  
  
  
  #small number because treating MAF as fixed
  if(SNPSE == FALSE){
    varSNPSE2=(.00000001)^2
  }
  
  if(SNPSE != FALSE){
    varSNPSE2 = SNPSE^2
  }
  
  V_LD<-as.matrix(covstruc[[1]])
  S_LD<-as.matrix(covstruc[[2]])
  I_LD<-as.matrix(covstruc[[3]])
  
  beta_SNP<-SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))] 
  SE_SNP<-SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))] 
  
  #enter in k for number of phenotypes
  k<-ncol(beta_SNP)
  
  #set univariate intercepts to 1 if estimated below 1
  diag(I_LD)<-ifelse(diag(I_LD)<= 1, 1, diag(I_LD))
  
  #f = number of SNPs in dataset
  f=nrow(beta_SNP) 
  
  
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
  
  ##pull the column names specified in the munge function
  traits<-colnames(S_LD)
  
  #function to create lavaan syntax for a 1 factor model given k phenotypes
  write.Model1 <- function(k, label = "V") {  
    Model1 <- ""
    for (i in 1) {
      lineSNP <- paste(traits[[i]], " ~ 0*SNP",sep = "")
      if (k-i > 0) {
        lineSNP2 <- " \n "
        for (j in (i+1):k) {
          lineSNP2 <- paste(lineSNP2, traits[[j]], " ~ 0*SNP", " \n ", sep = "")
        }
      }
    } 
    
    for (i in 1) {
      linestart <- paste("F1"," =~ ",traits[[i]], sep = "")  
      if (k-i > 0) {
        linemid <- ""
        for (j in (i+1):k) {
          linemid <- paste(linemid, " + ", traits[[j]], sep = "")
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
    #create empty shell of V_SNP matrix
    V_SNP<-diag(k)
    
    #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
    for (p in 1:nrow(coords)) { 
      x<-coords[p,1]
      y<-coords[p,2]
      if (x != y) { 
        V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i]^2)}
      if (x == y) {
        V_SNP[x,x]<-(SE_SNP[i,x]*I_LD[x,x]*varSNP[i])^2
      }
    }
    
    ##create shell of full sampling covariance matrix
    V_Full<-diag(((k+1)*(k+2))/2)
    
    ##input the ld-score regression region of sampling covariance from ld-score regression SEs
    V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD
    
    ##add in SE of SNP variance as first observation in sampling covariance matrix
    V_Full[1,1]<-varSNPSE2
    
    ##add in SNP region of sampling covariance matrix
    V_Full[2:(k+1),2:(k+1)]<-V_SNP
    
    k2<-nrow(V_Full)
    smooth2<-ifelse(eigen(V_Full)$values[k2] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)
    
    if(toler==FALSE){
      W <- solve(V_Full)
    }
    
    if(toler!=FALSE){
      W <- solve(V_Full,tol=toler)
    }
    
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
    S_Fullrun[1:(k+1),1]<-S_SNP
    S_Fullrun[1,1:(k+1)]<-t(S_SNP)
    
    ##pull in variables names specified in LDSC function and name first column as SNP
    colnames(S_Fullrun)<-c("SNP", colnames(S_LD))
    
    ##name rows like columns
    rownames(S_Fullrun)<-colnames(S_Fullrun)
    
    ##smooth to near positive definite if either V or S are non-positive definite
    ks<-nrow(S_Fullrun)
    smooth1<-ifelse(eigen(S_Fullrun)$values[ks] <= 0, S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat), S_Fullrun<-S_Fullrun)
    
    suppress<-tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf,optim.force.converged=TRUE)) 
    
    order <- rearrange(k = k+1, fit = ReorderModel, names = rownames(S_Fullrun))
  }
  
  SNPs2<-SNPs[,1:6]
  rm(SNPs)
  #SNPs2<-suppressWarnings(split(SNPs2,1:int))
  #split the V_SNP and S_SNP matrices into as many (cores - 1) as are aviailable on the local computer
  beta_SNP<-suppressWarnings(split(beta_SNP,1:int))
  SE_SNP<-suppressWarnings(split(SE_SNP,1:int))
  varSNP<-suppressWarnings(split(varSNP,1:int))
  
  ##estimation for 2S-DWLS-R
  if(estimation=="DWLS"){
    
    ##foreach parallel processing that rbinds results across cores
    results<-foreach(n = icount(int), .combine = 'rbind') %:% 
      
      foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = "lavaan") %dopar% { 
        
        #create empty shell of V_SNP matrix
        V_SNP<-diag(k)
        
        #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SE_SNP[[n]][i,y]*SE_SNP[[n]][i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[[n]][i]^2)}
          if (x == y) {
            V_SNP[x,x]<-(SE_SNP[[n]][i,x]*I_LD[x,x]*varSNP[[n]][i])^2
          }
        }
        
        ##create shell of full sampling covariance matrix
        V_Full<-diag(((k+1)*(k+2))/2)
        
        ##input the ld-score regression region of sampling covariance from ld-score regression SEs
        V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD
        
        ##add in SE of SNP variance as first observation in sampling covariance matrix
        V_Full[1,1]<-varSNPSE2
        
        ##add in SNP region of sampling covariance matrix
        V_Full[2:(k+1),2:(k+1)]<-V_SNP
        
        kv<-nrow(V_Full)
        smooth2<-ifelse(eigen(V_Full)$values[kv] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)
        
        #reorder sampling covariance matrix based on what lavaan expects given the specified model
        V_Full_Reorder <- V_Full[order,order]
        u<-nrow(V_Full_Reorder)
        V_Full_Reorderb<-diag(u)
        diag(V_Full_Reorderb)<-diag(V_Full_Reorder)
        
        ##invert the reordered sampling covariance matrix to create a weight matrix 
        if(toler==FALSE){
          W<- solve(V_Full_Reorderb)
        }
        
        if(toler!=FALSE){
          W <- solve(V_Full_Reorderb,tol=toler)
        }
        
        #create empty vector for S_SNP
        S_SNP<-vector(mode="numeric",length=k+1)
        
        #enter SNP variance from reference panel as first observation
        S_SNP[1]<-varSNP[[n]][i]
        
        #enter SNP covariances (standardized beta * SNP variance from refference panel)
        for (p in 1:k) {
          S_SNP[p+1]<-varSNP[[n]][i]*beta_SNP[[n]][i,p]
        }
        
        #create shell of the full S (observed covariance) matrix
        S_Fullrun<-diag(k+1)
        
        ##add the LD portion of the S matrix
        S_Fullrun[(2:(k+1)),(2:(k+1))]<-S_LD
        
        ##add in observed SNP variances as first row/column
        S_Fullrun[1:(k+1),1]<-S_SNP
        S_Fullrun[1,1:(k+1)]<-t(S_SNP)
        
        ##smooth to near positive definite if either V or S are non-positive definite
        ks<-nrow(S_Fullrun)
        smooth1<-ifelse(eigen(S_Fullrun)$values[ks] <= 0, S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat), S_Fullrun<-S_Fullrun)
        
        #name the columns
        colnames(S_Fullrun)<-c("SNP", colnames(S_LD))
        
        ##name rows like columns
        rownames(S_Fullrun)<-colnames(S_Fullrun)
        
        ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
        test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf))
        
        test$warning$message[1]<-ifelse(is.null(test$warning$message), test$warning$message[1]<-0, test$warning$message[1])
        
        if(class(test$value)[1] == "lavaan" & grepl("solution has NOT",  as.character(test$warning)) != TRUE){
          #pull the delta matrix (this doesn't depend on N)
          S2.delt <- lavInspect(Model1_Results, "delta")
          
          ##weight matrix from stage 2
          S2.W <- lavInspect(Model1_Results, "WLS.V") 
          
          #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
          if(toler==FALSE){
            bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
          }
          
          if(toler!=FALSE){
            bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler)
          }
          
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
          
          #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
          ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs=2, optim.dx.tol = +Inf) 
          
          #pull the delta matrix for Q (this doesn't depend on N)
          S2.delt_Q <- lavInspect(ModelQ_Results, "delta")
          
          ##weight matrix from stage 2 for Q
          S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V") 
          
          #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
          if(toler==FALSE){
            bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q) 
          }
          
          if(toler!=FALSE){
            bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q,tol=toler) 
          }
          
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
          Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta
          
          ##pull all the results into a single row
          cbind(i,n,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
        }else{
          cbind(i,n,rep(NA,7),ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
        } 
      }
  }
  
  
  ##2S-ML-R estimation
  if(estimation=="ML"){
    
    ##foreach parallel processing that rbinds results across cores  
    results<-foreach(n = icount(int), .combine = 'rbind') %:% 
      
      foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = "lavaan") %dopar% { 
        
        #create empty shell of V_SNP matrix
        V_SNP<-diag(k)
        
        #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SE_SNP[[n]][i,y]*SE_SNP[[n]][i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[[n]][i]^2)}
          if (x == y) {
            V_SNP[x,x]<-(SE_SNP[[n]][i,x]*I_LD[x,x]*varSNP[[n]][i])^2
          }
        }
        
        ##create shell of full sampling covariance matrix
        V_Full<-diag(((k+1)*(k+2))/2)
        
        ##input the ld-score regression region of sampling covariance from ld-score regression SEs
        V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD
        
        ##add in SE of SNP variance as first observation in sampling covariance matrix
        V_Full[1,1]<-varSNPSE2
        
        ##add in SNP region of sampling covariance matrix
        V_Full[2:(k+1),2:(k+1)]<-V_SNP
        
        kv<-nrow(V_Full)
        smooth2<-ifelse(eigen(V_Full)$values[kv] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)
        
        #reorder sampling covariance matrix based on what lavaan expects given the specified model
        V_Full_Reorder <- V_Full[order,order]
        u<-nrow(V_Full_Reorder)
        V_Full_Reorderb<-diag(u)
        diag(V_Full_Reorderb)<-diag(V_Full_Reorder)
        
        ##invert the reordered sampling covariance matrix to create a weight matrix 
        if(toler==FALSE){
          W<- solve(V_Full_Reorderb)
        }
        
        if(toler!=FALSE){
          W <- solve(V_Full_Reorderb,tol=toler)
        }
        
        #create empty vector for S_SNP
        S_SNP<-vector(mode="numeric",length=k+1)
        
        #enter SNP variance from reference panel as first observation
        S_SNP[1]<-varSNP[[n]][i]
        
        #enter SNP covariances (standardized beta * SNP variance from refference panel)
        for (p in 1:k) {
          S_SNP[p+1]<-varSNP[[n]][i]*beta_SNP[[n]][i,p]
        }
        
        #create shell of the full S (observed covariance) matrix
        S_Fullrun<-diag(k+1)
        
        ##add the LD portion of the S matrix
        S_Fullrun[(2:(k+1)),(2:(k+1))]<-S_LD
        
        ##add in observed SNP variances as first row/column
        S_Fullrun[1:(k+1),1]<-S_SNP
        S_Fullrun[1,1:(k+1)]<-t(S_SNP)
        
        ##smooth to near positive definite if either V or S are non-positive definite
        ks<-nrow(S_Fullrun)
        smooth1<-ifelse(eigen(S_Fullrun)$values[ks] <= 0, S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat), S_Fullrun<-S_Fullrun)
        
        #name the columns
        colnames(S_Fullrun)<-c("SNP", colnames(S_LD))
        
        ##name rows like columns
        rownames(S_Fullrun)<-colnames(S_Fullrun)
        
        ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
        test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML", sample.nobs = 200, optim.dx.tol = +Inf))
        
        test$warning$message[1]<-ifelse(is.null(test$warning$message), test$warning$message[1]<-0, test$warning$message[1])
        
        if(class(test$value)[1] == "lavaan" & grepl("solution has NOT",  as.character(test$warning)) != TRUE){
          #pull the delta matrix (this doesn't depend on N)
          S2.delt <- lavInspect(Model1_Results, "delta")
          
          ##weight matrix from stage 2
          S2.W <- lavInspect(Model1_Results, "WLS.V") 
          
          #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
          
          if(toler==FALSE){
            bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
          }
          
          if(toler!=FALSE){
            bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler) 
          }
          
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
          
          #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
          ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "ML", sample.nobs=200, optim.dx.tol = +Inf) 
          
          #pull the delta matrix for Q (this doesn't depend on N)
          S2.delt_Q <- lavInspect(ModelQ_Results, "delta")
          
          ##weight matrix from stage 2 for Q
          S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V") 
          
          #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specifie
          if(toler==FALSE){
            bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q) 
          }
          
          if(toler!=FALSE){
            bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q,tol=toler) 
          }
          
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
          Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta
          
          ##put the corrected standard error and Q in same dataset
          cbind(i,n,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
        }else{
          cbind(i,n,rep(NA,7),ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
        }
        
      }
  }
  
  ##name the columns of the results file
  colnames(results)=c("i","n","lhs","op","rhs","est","se", "se_c", "Q", "fail", "warning")
  
  ##sort results so it is in order of the output lists provided for the function
  results<- results[order(results$i, results$n),] 
  
  results$se <- NULL
  results2<-cbind(SNPs2,results)
  results2$Z_Estimate<-results2$est/results2$se_c
  results2$Pval_Estimate<-2*pnorm(abs(results2$Z_Estimate),lower.tail=FALSE)
  results2$Q_df<-k-1
  results2$Q_pval<-pchisq(results2$Q,results2$Q_df,lower.tail=FALSE)
  results2$i<-1:nrow(results2)
  results2$n<-NULL
  results2<-results2[,c(1:12,16,17,13,18,19,14,15)]
  time_all<-proc.time()-time
  print(time_all[3])
  return(results2)
  
  }
  
  if(!is.null(Output)){
    ##split the V and S matrices into as many (cores - 1) as are aviailable on the local computer
    V_Full<-suppressWarnings(split(Output[[1]],1:int))
    S_Full<-suppressWarnings(split(Output[[2]],1:int))
    
    #enter in k for number of phenotypes 
    k<-ncol(S_Full[[1]][[1]])-1
    
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
    ##pull the column names specified in the munge function
    traits<-colnames(S_Full[[1]][[1]])[2:(k+1)]
    
    #function to create lavaan syntax for a 1 factor model given k phenotypes
    write.Model1 <- function(k, label = "V") {  
      Model1 <- ""
      for (i in 1) {
        lineSNP <- paste(traits[[i]], " ~ 0*SNP",sep = "")
        if (k-i > 0) {
          lineSNP2 <- " \n "
          for (j in (i+1):k) {
            lineSNP2 <- paste(lineSNP2, traits[[j]], " ~ 0*SNP", " \n ", sep = "")
          }
        }
      } 
      
      for (i in 1) {
        linestart <- paste("F1"," =~ ",traits[[i]], sep = "")  
        if (k-i > 0) {
          linemid <- ""
          for (j in (i+1):k) {
            linemid <- paste(linemid, " + ", traits[[j]], sep = "")
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
      if(toler==FALSE){
        W<- solve(V_Full[[1]][[i]])
      }
      
      if(toler!=FALSE){
        W <- solve(V_Full[[1]][[i]],tol=toler)
      }
      
      S_Fullrun<-S_Full[[1]][[i]]
      
      test2<-tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf,optim.force.converged=TRUE)) 
      if(class(test2$value)[1]=="lavaan"){
        order <- rearrange(k = k+1, fit = ReorderModel, names = rownames(S_Full[[1]][[i]]))}else{
          i<-10
          
          if(toler==FALSE){
            W<- solve(V_Full[[1]][[i]])
          }
          
          if(toler!=FALSE){
            W <- solve(V_Full[[1]][[i]],tol=toler)
          }
          
          S_Fullrun<-S_Full[[1]][[i]]
          
          test2<-tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf,optim.force.converged=TRUE))
          
          order <- rearrange(k = k+1, fit = ReorderModel, names = rownames(S_Full[[1]][[i]]))
        }
    }
    
    ##estimation for 2S-DWLS-R
    if(estimation=="DWLS"){
      
      ##foreach parallel processing that rbinds results across cores
      results<-foreach(n = icount(int), .combine = 'rbind') %:% 
        
        foreach (i=1:length(V_Full[[n]]), .combine='rbind', .packages = "lavaan") %dopar% { 
          
          #reorder sampling covariance matrix based on what lavaan expects given the specified model
          V_Full_Reorder <- V_Full[[n]][[i]][order,order]
          u<-nrow(V_Full_Reorder)
          V_Full_Reorderb<-diag(u)
          diag(V_Full_Reorderb)<-diag(V_Full_Reorder)
          
          ##invert the reordered sampling covariance matrix to create a weight matrix 
          if(toler==FALSE){
            W<- solve(V_Full_Reorderb)
          }
          
          if(toler!=FALSE){
            W <- solve(V_Full_Reorderb,tol=toler)
          }
          
          #import the S_Full matrix for appropriate run
          S_Fullrun<-S_Full[[n]][[i]]
          
          ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
          test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf))
          
          test$warning$message[1]<-ifelse(is.null(test$warning$message), test$warning$message[1]<-0, test$warning$message[1])
          
          if(class(test$value)[1] == "lavaan" & grepl("solution has NOT",  as.character(test$warning)) != TRUE){
            #pull the delta matrix (this doesn't depend on N)
            S2.delt <- lavInspect(Model1_Results, "delta")
            
            ##weight matrix from stage 2
            S2.W <- lavInspect(Model1_Results, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
            if(toler==FALSE){
              bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
            }
            
            if(toler!=FALSE){
              bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler)
            }
            
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
            
            #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
            ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs=2, optim.dx.tol = +Inf) 
            
            #pull the delta matrix for Q (this doesn't depend on N)
            S2.delt_Q <- lavInspect(ModelQ_Results, "delta")
            
            ##weight matrix from stage 2 for Q
            S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
            if(toler==FALSE){
              bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q) 
            }
            
            if(toler!=FALSE){
              bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q,tol=toler) 
            }
            
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
            Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta
            
            ##pull all the results into a single row
            cbind(i,n,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
          }else{
            cbind(i,n,rep(NA,7),ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
          } 
        }
    }
    
    
    ##2S-ML-R estimation
    if(estimation=="ML"){
      
      ##foreach parallel processing that rbinds results across cores  
      results<-foreach(n = icount(int), .combine = 'rbind') %:% 
        
        foreach (i=1:length(V_Full[[n]]), .combine='rbind', .packages = "lavaan") %dopar% { 
          
          #reorder sampling covariance matrix based on what lavaan expects given the specified model
          V_Full_Reorder <- V_Full[[n]][[i]][order,order]
          
          #import the S_Full matrix for appropriate run
          S_Fullrun<-S_Full[[n]][[i]]
          
          ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
          test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML", sample.nobs = 200, optim.dx.tol = +Inf))
          
          test$warning$message[1]<-ifelse(is.null(test$warning$message), test$warning$message[1]<-0, test$warning$message[1])
          
          if(class(test$value)[1] == "lavaan" & grepl("solution has NOT",  as.character(test$warning)) != TRUE){
            #pull the delta matrix (this doesn't depend on N)
            S2.delt <- lavInspect(Model1_Results, "delta")
            
            ##weight matrix from stage 2
            S2.W <- lavInspect(Model1_Results, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
            
            if(toler==FALSE){
              bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
            }
            
            if(toler!=FALSE){
              bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler) 
            }
            
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
            
            #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
            ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "ML", sample.nobs=200, optim.dx.tol = +Inf) 
            
            #pull the delta matrix for Q (this doesn't depend on N)
            S2.delt_Q <- lavInspect(ModelQ_Results, "delta")
            
            ##weight matrix from stage 2 for Q
            S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specifie
            if(toler==FALSE){
              bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q) 
            }
            
            if(toler!=FALSE){
              bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q,tol=toler) 
            }
            
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
            Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta
            
            ##put the corrected standard error and Q in same dataset
            cbind(i,n,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
          }else{
            cbind(i,n,rep(NA,7),ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
          }
          
        }
    }
    
    ##name the columns of the results file
    colnames(results)=c("i","n","lhs","op","rhs","est","se", "se_c", "Q", "fail", "warning")
    
    ##sort results so it is in order of the output lists provided for the function
    results<- results[order(results$i, results$n),] 
    
    results$se <- NULL
    results2<-cbind(Output[[3]],results)
    results2$Z_Estimate<-results2$est/results2$se_c
    results2$Pval_Estimate<-2*pnorm(abs(results2$Z_Estimate),lower.tail=FALSE)
    results2$Q_df<-ncol(S_Full[[1]][[1]])-2
    results2$Q_pval<-pchisq(results2$Q,results2$Q_df,lower.tail=FALSE)
    results2$i<-1:nrow(results2)
    results2$n<-NULL
    results2<-results2[,c(1:12,16,17,13,18,19,14,15)]
    time_all<-proc.time()-time
    print(time_all[3])
    return(results2)
    
  }
  }
