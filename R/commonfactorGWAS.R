
commonfactorGWAS <-function(covstruc=NULL,SNPs=NULL,estimation="DWLS",cores=NULL,toler=FALSE,SNPSE=FALSE,parallel=TRUE,GC="standard",MPI=FALSE,TWAS=FALSE,smooth_check=FALSE){ 
  time<-proc.time()
  
  if(exists("Output")){
    stop("Please note that an update was made to commonfactorGWAS on 4/1/21 so that addSNPs output CANNOT be fed directly to the function. It now expects the 
            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments.")
  }
  
  
  if(class(SNPs)[1] == "character"){
    print("You are likely listing arguments (e.g. Output = ...) in the order of a previous version of commonfactorGWAS. The function no longer accepts output from addSNPs. The current version of the function is faster and saves memory. It expects the 
          output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?commonfactorGWAS for help on proper usage.")
    warning("You are likely listing arguments (e.g. Output = ...) in the order of a previous version of commonfactorGWAS.  The function no longer accepts output from addSNPs. The current version of the function is faster and saves memory. It expects the 
            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?commonfactorGWAS for help on proper usage.")
  }else{
    if(is.null(SNPs) | is.null(covstruc)){
      print("You are likely listing arguments (e.g. Output = ...) in the order of a previous version of commonfactorGWAS. The function no longer accepts output from addSNPs. The current version of the function is faster and saves memory. It expects the 
          output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?commonfactorGWAS for help on proper usage.")
      warning("You are likely listing arguments (e.g. Output = ...) in the order of a previous version of commonfactorGWAS.  The function no longer accepts output from addSNPs. The current version of the function is faster and saves memory. It expects the 
            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?commonfactorGWAS for help on proper usage.") 
    }
  }
  
  
  Operating<-Sys.info()[['sysname']]
  
  ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
  SNPs<-data.frame(SNPs)
  
  if(TWAS == FALSE){
    SNPs$A1<-as.character(SNPs$A1)
    SNPs$A2<-as.character(SNPs$A2)
    SNPs$SNP<-as.character(SNPs$SNP)
    
    #SNP variance
    varSNP=2*SNPs$MAF*(1-SNPs$MAF)
  }
  
  if(TWAS == TRUE){
    SNPs$Gene<-as.character(SNPs$Gene)
    SNPs$Panel<-as.character(SNPs$Panel)
    varSNP=SNPs$HSQ
  }
  
  #small number because treating MAF as fixed
  if(SNPSE == FALSE){
    varSNPSE2=(.0005)^2
  }
  
  if(SNPSE != FALSE){
    varSNPSE2 = SNPSE^2
  }
  
  V_LD<-as.matrix(covstruc[[1]])
  S_LD<-as.matrix(covstruc[[2]])
  I_LD<-as.matrix(covstruc[[3]])
  
  check_names<-str_detect(colnames(S_LD), "-")
  if(any(check_names==TRUE)){warning("Your trait names specified when running the ldsc function include mathematical arguments (e.g., + or -) that will be misread by lavaan. Please rename the traits.")}
  
  beta_SNP<-SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))] 
  SE_SNP<-SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))] 
  
  #enter in k for number of phenotypes
  k<-ncol(beta_SNP)
  
  #print warning if number of traits are unequal across SNPs and LDSC output
  if(ncol(beta_SNP) != ncol(S_LD)){
    stop("There are different numbers of traits in the sumstats and ldsc output. Please verify that the same summary statistics have been provided to both functions before running commonfactorGWAS.") 
  }
  
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
  
  ##pull the coordinates of the I_LD matrix to loop making the V_SNP matrix
  coords<-which(I_LD != 'NA', arr.ind= T)
  
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
    
    suppress<-tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf,optim.force.converged=TRUE,control=list(iter.max=1))) 
    
    order <- rearrange(k = k+1, fit = ReorderModel, names = rownames(S_Fullrun))
  }
  
  
  if(TWAS == FALSE){
    SNPs2<-SNPs[,1:6]}
  if(TWAS == TRUE){
    SNPs2<-SNPs[,1:3]
  }
  rm(SNPs)
  
  if(parallel==FALSE){
    
    ##name the columns of the results file
    if(smooth_check==FALSE){
    results=as.data.frame(matrix(NA,ncol=10,nrow=f))
    colnames(results)=c("i","lhs","op","rhs","est","se", "se_c", "Q", "fail", "warning")
    }
    if(smooth_check==TRUE){
      results=as.data.frame(matrix(NA,ncol=11,nrow=f))
      colnames(results)=c("i","lhs","op","rhs","est","se", "se_c", "Q", "fail", "warning","Z_smooth")
    }
    
    for (i in 1:f) { 
      
      #create empty shell of V_SNP matrix
      V_SNP<-diag(k)
      
      #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
      
      #double GC correctiong using univariate LDSC intercepts
      if(GC == "conserv"){
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i]^2)}
          if (x == y) {
            V_SNP[x,x]<-(SE_SNP[i,x]*I_LD[x,x]*varSNP[i])^2
          }
        }
      }
      
      #single GC correction using sqrt of univariate LDSC intercepts
      if(GC == "standard"){
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*varSNP[i]^2)}
          if (x == y) {
            V_SNP[x,x]<-(SE_SNP[i,x]*sqrt(I_LD[x,x])*varSNP[i])^2
          }
        }
      }
      
      #no GC correction
      if(GC == "none"){
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*varSNP[i]^2)}
          if (x == y) {
            V_SNP[x,x]<-(SE_SNP[i,x]*varSNP[i])^2
          }
        }
      }
      
      if(smooth_check == TRUE){
        if(GC == "conserv"){
          Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*diag(I_LD))
        }
        if(GC=="standard"){
          Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*sqrt(diag(I_LD))) 
        }
        if(GC=="none"){
          Z_pre<-beta_SNP[i,]/SE_SNP[i,] 
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
      if(eigen(V_Full)$values[kv] <= 0){
        V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat)
        V_smooth<-1
      }
    
      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Full_Reorder <- V_Full[order,order]
      u<-nrow(V_Full_Reorder)
      W<-diag(u)
      diag(W)<-diag(V_Full_Reorder)
      
      ##invert the reordered sampling covariance matrix to create a weight matrix 
      if(toler==FALSE){
        W<- solve(W)
      }
      
      if(toler!=FALSE){
        W <- solve(W,tol=toler)
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
      
      colnames(S_Fullrun)<-c("SNP", colnames(S_LD))
      
      ##name rows like columns
      rownames(S_Fullrun)<-colnames(S_Fullrun)
      
      ##smooth to near positive definite if either V or S are non-positive definite
      ks<-nrow(S_Fullrun)
      
      if(eigen(S_Fullrun)$values[ks] <= 0){
        S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat)
        S_smooth<-1
      }
      
      if(smooth_check == TRUE){
        if(exists("S_smooth") | exists("V_smooth")){
          SE_smooth<-matrix(0, ks, ks)
          SE_smooth[lower.tri(SE_smooth,diag=TRUE)] <-sqrt(diag(V_Full))
          Z_smooth<-(S_Fullrun/SE_smooth)[2:ks,1]
          Z_smooth<-max(abs(Z_smooth-Z_pre))
        }else{Z_smooth<-0}
      }
      
      ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
      if(estimation == "DWLS"){
        test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf))
      }
      
      if(estimation == "ML"){
        test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML",sample.nobs = 200, optim.dx.tol = +Inf, sample.cov.rescale=FALSE))
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
          testQ<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2,  optim.dx.tol = +Inf))
        }
        
        if(estimation == "ML"){
          testQ<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "ML", sample.nobs = 200, optim.dx.tol = +Inf, sample.cov.rescale=FALSE))
        }
        
        testQ$warning$message[1]<-ifelse(is.null(testQ$warning$message), testQ$warning$message[1]<-"Safe", testQ$warning$message[1])
        testQ$warning$message[1]<-ifelse(exists("ModelQ_Results") == FALSE, testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!", ifelse(is.na(inspect(ModelQ_Results, "se")$theta[1,2]) == TRUE ,testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!" , testQ$warning$message[1]))
        
        if(as.character(testQ$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
          #pull the delta matrix for Q (this doesn't depend on N)
          S2.delt_Q <- lavInspect(ModelQ_Results, "delta")
          
          ##weight matrix from stage 2 for Q
          S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V") 
          
          #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
          Q_catch<-tryCatch.W.E(bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q,tol=toler)) 
          
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
        if(smooth_check == FALSE){
        results[i,]<-data.frame(i,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value)[1] == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning)[1] == 'NULL', 0, as.character(test$warning$message[1])),stringsAsFactors = FALSE)
        }
        
       if(smooth_check == TRUE){
        results[i,]<-data.frame(i,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value)[1] == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning)[1] == 'NULL', 0, as.character(test$warning$message[1])),Z_smooth,stringsAsFactors = FALSE)
       }
        
      }else{
        if(smooth_check == FALSE){
        results[i,]<-data.frame(i,inspect(Model1_Results,"list")[k+1,-c(1,5:13,15)],t(rep(NA,3)),ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning)[1] == 'NULL', 0, as.character(test$warning$message[1])),stringsAsFactors = FALSE)
        }
        if(smooth_check == TRUE){
        results[i,]<-data.frame(i,inspect(Model1_Results,"list")[k+1,-c(1,5:13,15)],t(rep(NA,3)),ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning)[1] == 'NULL', 0, as.character(test$warning$message[1])),Z_smooth,stringsAsFactors = FALSE)
        }
      } 
    }

    time_all<-proc.time()-time
    print(time_all[3])
    results$se <- NULL
    results<-cbind(SNPs2,results)
    results$Z_Estimate<-results$est/results$se_c
    results$Pval_Estimate<-2*pnorm(abs(results$Z_Estimate),lower.tail=FALSE)
    results$Q_df<-k-1
    results$Q_pval<-pchisq(results$Q,results$Q_df,lower.tail=FALSE)
    
    if(TWAS == FALSE & smooth_check == FALSE){
      results<-results[,c(1:12,16,17,13,18,19,14,15)]
      }
    if(TWAS == FALSE & smooth_check == TRUE){
      results<-results[,c(1:12,17,18,13,19,20,14,15,16)]
      }
    if(TWAS == TRUE & smooth_check == FALSE){
      results<-results[,c(1:9,13,14,10,15,16,11,12)]
      results$rhs<-rep("Gene",nrow(results))
    }
    if(TWAS == TRUE & smooth_check == TRUE){
      results<-results[,c(1:9,14,15,10,16,17,11,12,13)]
      results$rhs<-rep("Gene",nrow(results))
    }
    return(results)
    
  }
  
  
  if(parallel == TRUE & Operating != "Windows"){
    
    if(is.null(cores)){
      ##if no default provided use 1 less than the total number of cores available so your computer will still function
      int <- detectCores() - 1
    }else{int<-cores}
    
    if(MPI == FALSE){
      
      registerDoParallel(int)
      
      ##specify the cores should have access to the local environment
      makeCluster(int, type="FORK")
    }
    
    if(MPI == TRUE){
      #register MPI
      cluster <- getMPIcluster()
      
      #register cluster; no makecluster as ibrun already starts the MPI process. 
      registerDoParallel(cluster)
      
    }
    
    #split the V_SNP and S_SNP matrices into as many (cores - 1) as are aviailable on the local computer
    beta_SNP<-suppressWarnings(split(beta_SNP,1:int))
    SE_SNP<-suppressWarnings(split(SE_SNP,1:int))
    varSNP<-suppressWarnings(split(varSNP,1:int))
    
    ##foreach parallel processing that rbinds results across cores
    results<-foreach(n = icount(int), .combine = 'rbind') %:% 
      
      foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = "lavaan") %dopar% { 
        
        #create empty shell of V_SNP matrix
        V_SNP<-diag(k)
        
        #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
        if(GC == "conserv"){
          for (p in 1:nrow(coords)) { 
            x<-coords[p,1]
            y<-coords[p,2]
            if (x != y) { 
              V_SNP[x,y]<-(SE_SNP[[n]][i,y]*SE_SNP[[n]][i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[[n]][i]^2)}
            if (x == y) {
              V_SNP[x,x]<-(SE_SNP[[n]][i,x]*I_LD[x,x]*varSNP[[n]][i])^2
            }
          }
        }
        
        if(GC == "standard"){
          for (p in 1:nrow(coords)) { 
            x<-coords[p,1]
            y<-coords[p,2]
            if (x != y) { 
              V_SNP[x,y]<-(SE_SNP[[n]][i,y]*SE_SNP[[n]][i,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*varSNP[[n]][i]^2)}
            if (x == y) {
              V_SNP[x,x]<-(SE_SNP[[n]][i,x]*sqrt(I_LD[x,x])*varSNP[[n]][i])^2
            }
          }
        }
        
        if(GC == "none"){
          for (p in 1:nrow(coords)) { 
            x<-coords[p,1]
            y<-coords[p,2]
            if (x != y) { 
              V_SNP[x,y]<-(SE_SNP[[n]][i,y]*SE_SNP[[n]][i,x]*I_LD[x,y]*varSNP[[n]][i]^2)}
            if (x == y) {
              V_SNP[x,x]<-(SE_SNP[[n]][i,x]*varSNP[[n]][i])^2
            }
          }
        }
        
        if(smooth_check == TRUE){
          if(GC == "conserv"){
            Z_pre<-beta_SNP[[n]][i,]/(SE_SNP[[n]][i,]*diag(I_LD))
          }
          if(GC=="standard"){
            Z_pre<-beta_SNP[[n]][i,]/(SE_SNP[[n]][i,]*sqrt(diag(I_LD))) 
          }
          if(GC=="none"){
            Z_pre<-beta_SNP[[n]][i,]/SE_SNP[[n]][i,] 
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
        if(eigen(V_Full)$values[kv] <= 0){
          V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat)
          V_smooth<-1
        }
        
        #reorder sampling covariance matrix based on what lavaan expects given the specified model
        V_Full_Reorder <- V_Full[order,order]
        u<-nrow(V_Full_Reorder)
        W<-diag(u)
        diag(W)<-diag(V_Full_Reorder)
        
        ##invert the reordered sampling covariance matrix to create a weight matrix 
        if(toler==FALSE){
          W<- solve(W)
        }
        
        if(toler!=FALSE){
          W <- solve(W,tol=toler)
        }
        
        #create empty vector for S_SNP
        S_SNP<-vector(mode="numeric",length=k+1)
        
        #enter SNP variance from reference panel as first observation
        S_SNP[1]<-varSNP[[n]][i]
        
        #enter SNP covariances (standardized beta * SNP variance from reference panel)
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
        if(eigen(S_Fullrun)$values[ks] <= 0){
          S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat)
          S_smooth<-1
        }
        
        if(smooth_check == TRUE){
          if(exists("S_smooth") | exists("V_smooth")){
            SE_smooth<-matrix(0, ks, ks)
            SE_smooth[lower.tri(SE_smooth,diag=TRUE)] <-sqrt(diag(V_Full))
            Z_smooth<-(S_Fullrun/SE_smooth)[2:ks,1]
            Z_smooth<-max(abs(Z_smooth-Z_pre))
          }else{Z_smooth<-0}
        }
        
        #name the columns
        colnames(S_Fullrun)<-c("SNP", colnames(S_LD))
        
        ##name rows like columns
        rownames(S_Fullrun)<-colnames(S_Fullrun)
        
        ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
        if(estimation == "DWLS"){
          test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf))
        }
        
        if(estimation == "ML"){
          test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML",sample.nobs = 200, optim.dx.tol = +Inf, sample.cov.rescale=FALSE))
        }
        
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
          if(estimation == "DWLS"){
            testQ<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs=2, optim.dx.tol = +Inf)) 
          }
          if(estimation == "ML"){
            testQ<-tryCatch.W.E(ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Fullrun, estimator = "ML", sample.nobs=200, optim.dx.tol = +Inf, sample.cov.rescale=FALSE))
          }
          
          testQ$warning$message[1]<-ifelse(is.null(testQ$warning$message), testQ$warning$message[1]<-"Safe", testQ$warning$message[1])
          testQ$warning$message[1]<-ifelse(exists("ModelQ_Results") == FALSE, testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!", ifelse(is.na(inspect(ModelQ_Results, "se")$theta[1,2]) == TRUE ,testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!" , testQ$warning$message[1]))
          
          if(as.character(testQ$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
            #pull the delta matrix for Q (this doesn't depend on N)
            S2.delt_Q <- lavInspect(ModelQ_Results, "delta")
            
            ##weight matrix from stage 2 for Q
            S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
            Q_catch<-tryCatch.W.E(bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q,tol=toler)) 
            
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
          
          if(smooth_check==TRUE){
            cbind(i,n,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1],Z_smooth)
          }else{cbind(i,n,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q, ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])}
          
        }else{
          se<-NA
          se_c<-NA
          Q<-NA
          ##pull all the results into a single row
          if(smooth_check == FALSE){
          cbind(i,n,inspect(Model1_Results,"list")[k+1,-c(1,5:13,15)],se,se_c,Q,ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1])
          }else{cbind(i,n,inspect(Model1_Results,"list")[k+1,-c(1,5:13,15)],se,se_c,Q,ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1],  ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1],Z_smooth)}
            }
          
      }
    
  
    ##name the columns of the results file
    if(smooth_check == FALSE){
    colnames(results)=c("i","n","lhs","op","rhs","est","se", "se_c", "Q", "fail", "warning")
    }
    if(smooth_check == TRUE){
      colnames(results)=c("i","n","lhs","op","rhs","est","se", "se_c", "Q", "fail", "warning","Z_smooth")
    }
    
    ##sort results so it is in order of the output lists provided for the function
    results<- results[order(results$i, results$n),] 
    
    results$se <- NULL
    results<-cbind(SNPs2,results)
    results$Z_Estimate<-results$est/results$se_c
    results$Pval_Estimate<-2*pnorm(abs(results$Z_Estimate),lower.tail=FALSE)
    results$Q_df<-k-1
    results$Q_pval<-pchisq(results$Q,results$Q_df,lower.tail=FALSE)
    results$i<-1:nrow(results)
    results$n<-NULL
    
    if(TWAS == FALSE & smooth_check == FALSE){
      results<-results[,c(1:12,16,17,13,18,19,14,15)]
      }
    if(TWAS == FALSE & smooth_check == TRUE){
      results<-results[,c(1:12,17,18,13,19,20,14,15,16)]
      }
    if(TWAS == TRUE & smooth_check == FALSE){
      results<-results[,c(1:9,13,14,10,15,16,11,12)]
      results$rhs<-rep("Gene",nrow(results))
    }
    if(TWAS == TRUE & smooth_check == TRUE){
      results<-results[,c(1:9,14,15,10,16,17,11,12,13)]
      results$rhs<-rep("Gene",nrow(results))
    }

    time_all<-proc.time()-time
    print(time_all[3])
    return(results)
    
  }
  
  if(parallel == TRUE & Operating == "Windows"){
    stop("Parallel processing is not currently available for Windows operating systems. Please set the parallel argument to FALSE, or switch to a Linux or Mac operating system.")
  }
  
}
