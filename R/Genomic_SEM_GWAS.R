require(lavaan)

#requires these items as input from Michels LDSC software:
#I_LD = matrix of univariate and bivariate ldsc intercepts
#S_LD = genetic covariancae matrix of ldsc
#V_LD = genetic covairance matrix of ldsc

#requires these items from cleaned univariate GWAS summary stats along with MAFs from reference panel:
#beta_SNP = matrix of univariate GWAS betas (to be standardized)
#SE_SNP = matrix of univariate GWAS SEs
#MAF = vector of MAFs from reference panel

##estimation options = "DWLS" or "ML", with default for DWLS

Genomic_SEM_GWAS<-function(S_LD, V_LD, I_LD, beta_SNP, SE_SNP, MAF, estimation = "DWLS"){
  
  time<-proc.time()
  
  S_LD<-as.matrix(S_LD)
  V_LD<-as.matrix(V_LD)
  I_LD<-as.matrix(I_LD)
  
  #enter in k for number of phenotypes 
  k<-ncol(beta_SNP)
  
  #f = number of SNPs in dataset
  f=nrow(beta_SNP) 
  
  ##add in f column to match up lavaan results in case only a subset of SNPs are being run
  beta_SNP$f<-(1:f)
  SE_SNP$f<-(1:f)
  
  #SNP variance (updated with 1KG phase 3 MAFs)
  varSNP=2*MAF*(1-MAF)  
  
  #small number because treating MAF as fixed
  varSNPSE2=(.00000001)^2
  
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
  
  Model1 <- write.Model1(k=ncol(I_LD))
  
  #function to creat row/column names for S_full matrix
  write.names <- function(k, label = "V") {  
    varnames<-vector(mode="character",length=k+1)
    
    for (i in 1){
      varnames[1]<-c("SNP")}
    
    for (j in i:k) { 
      varnames[j+1]<-paste(label,j,sep="")}
    
    return(varnames)
  }
  
  S_names<-write.names(k=ncol(I_LD))
  
  ##run one model that specifies the factor structure so that lavaan knows how to rearrange the sampling covariance matrix
  for (i in 1) {
    
    #create empty vector for S_SNP
    S_SNP<-vector(mode="numeric",length=k+1)
    
    #enter SNP variance from reference panel as first observation
    S_SNP[1]<-varSNP[i,1]
    
    #enter SNP covariances (standardized beta * SNP variance from refference panel)
    for (p in 1:k) {
      S_SNP[p+1]<-varSNP[i,1]*beta_SNP[i,p]
    }
    
    #create shell of the full S (observed covariance) matrix
    S_Full<-diag(k+1)
    
    ##add the LD portion of the S matrix
    S_Full[(2:(k+1)),(2:(k+1))]<-S_LD
    
    ##add in observed SNP variances as first row/column
    S_Full[1:(k+1),1]<-S_SNP
    S_Full[1,1:(k+1)]<-t(S_SNP)
    
    ##name the columns/rows using the naming function defined outside of the loop
    rownames(S_Full) <- S_names
    colnames(S_Full) <- S_names
    
    #create empty shell of V_SNP matrix
    V_SNP<-diag(k)
    
    ##pull the coordinates of the I_LD matrix to loop making the V_SNP matrix
    coords<-which(I_LD != 'NA', arr.ind= T)
    
    #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
    for (p in 1:nrow(coords)) { 
      x<-coords[p,1]
      y<-coords[p,2]
      if (x != y) { 
        V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i,1]^2)}
      if (x == y) {
        V_SNP[x,x]<-(SE_SNP[i,x]*I_LD[x,x]*varSNP[i,1])^2
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
    
    #transform sampling covariance matrix into a weight matrix: 
    W <- solve(V_Full) 
    
    ReorderModel <- sem(Model1, sample.cov = S_Full, estimator = "DWLS", WLS.V = W, 
                        sample.nobs = 2) 
    
    order <- rearrange(k = ncol(beta_SNP), fit = ReorderModel, names = rownames(S_Full))
  }
  
  ##create results dataset
  RESULTS=as.data.frame(matrix(NA,ncol=8,nrow=f))
  colnames(RESULTS)=c("f","lhs","op","rhs","est","se", "se_c", "Q")
  
  ##create dataframe of runs that throw errors due to non-convergence or negative residuals
  fail<-as.data.frame(matrix(NA,ncol=2,nrow=f))
  colnames(fail)=c("fail", "run")
  
  if(estimation=="DWLS"){
  
  for (i in 1:f) {
    
    #create empty vector for S_SNP
    S_SNP<-vector(mode="numeric",length=k+1)
    
    #enter SNP variance from reference panel as first observation
    S_SNP[1]<-varSNP[i,1]
    
    #enter SNP covariances (standardized beta * SNP variance from refference panel)
    for (p in 1:k) {
      S_SNP[p+1]<-varSNP[i,1]*beta_SNP[i,p]
    }
    
    #create shell of the full S (observed covariance) matrix
    S_Full<-diag(k+1)
    
    ##add the LD portion of the S matrix
    S_Full[(2:(k+1)),(2:(k+1))]<-S_LD
    
    ##add in observed SNP variances as first row/column
    S_Full[1:(k+1),1]<-S_SNP
    S_Full[1,1:(k+1)]<-t(S_SNP)
    
    ##name the columns/rows using the naming function defined outside of the loop
    rownames(S_Full) <- S_names
    colnames(S_Full) <- S_names
    
    #create empty shell of V_SNP matrix
    V_SNP<-diag(k)
    
    ##pull the coordinates of the I_LD matrix to loop making the V_SNP matrix
    coords<-which(I_LD != 'NA', arr.ind= T)
    
    #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
    for (p in 1:nrow(coords)) { 
      x<-coords[p,1]
      y<-coords[p,2]
      if (x != y) { 
        V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i,1]^2)}
      if (x == y) {
        V_SNP[x,x]<-(SE_SNP[i,x]*I_LD[x,x]*varSNP[i,1])^2
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
    
    #reorder sampling covariance matrix based on what lavaan expects given the specified model
    V_Full_Reorder <- V_Full[order,order]
    
    ##invert the reordered sampling covariance matrix to create a weight matrix 
    W <- solve(V_Full_Reorder) 
    
    ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
    tryCatch(Model1_Results <- sem(Model1, sample.cov = S_Full, estimator = "DWLS", WLS.V = W, sample.nobs = 2),
             warning = function(w) {
               print(w)
               fail[i,]<<-cbind(1, i)},
             error = function(e) { 
               print(e)
               fail[i,]<<-cbind(2, i)})
    
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
    
    ##code to estimate Ronald's Q##
    #First pull the estimates from Step 1
    ModelQ <- parTable(Model1_Results)
    
    #fix the indicator loadings from Step 1, free the direct effects of the SNP on the indicators, and fix the factor residual variance
    ModelQ$free <- c(rep(0, k+1), 1:(k*2), 0, 0) 
    
    #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
    ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Full, estimator = "DWLS", WLS.V = W, sample.nobs=2) 
    
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
    Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta
    
    ##put the corrected standard error and Q in same dataset
    RESULTS[i,]=cbind(i,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q)
    
    print(i)
  }
  }
  
  if(estimation=="ML"){
    
    for (i in 1:f) {
      
      #create empty vector for S_SNP
      S_SNP<-vector(mode="numeric",length=k+1)
      
      #enter SNP variance from reference panel as first observation
      S_SNP[1]<-varSNP[i,1]
      
      #enter SNP covariances (standardized beta * SNP variance from refference panel)
      for (p in 1:k) {
        S_SNP[p+1]<-varSNP[i,1]*beta_SNP[i,p]
      }
      
      #create shell of the full S (observed covariance) matrix
      S_Full<-diag(k+1)
      
      ##add the LD portion of the S matrix
      S_Full[(2:(k+1)),(2:(k+1))]<-S_LD
      
      ##add in observed SNP variances as first row/column
      S_Full[1:(k+1),1]<-S_SNP
      S_Full[1,1:(k+1)]<-t(S_SNP)
      
      ##name the columns/rows using the naming function defined outside of the loop
      rownames(S_Full) <- S_names
      colnames(S_Full) <- S_names
      
      #create empty shell of V_SNP matrix
      V_SNP<-diag(k)
      
      ##pull the coordinates of the I_LD matrix to loop making the V_SNP matrix
      coords<-which(I_LD != 'NA', arr.ind= T)
      
      #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
      for (p in 1:nrow(coords)) { 
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) { 
          V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i,1]^2)}
        if (x == y) {
          V_SNP[x,x]<-(SE_SNP[i,x]*I_LD[x,x]*varSNP[i,1])^2
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
      
      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Full_Reorder <- V_Full[order,order]
      
      ##invert the reordered sampling covariance matrix to create a weight matrix 
      W <- solve(V_Full_Reorder) 
      
      ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
      tryCatch(Model1_Results <- sem(Model1, sample.cov = S_Full, estimator = "ML", sample.nobs = 200),
               warning = function(w) {
                 print(w)
                 fail[i,]<<-cbind(1, i)},
               error = function(e) { 
                 print(e)
                 fail[i,]<<-cbind(2, i)})
      
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
      
      ##code to estimate Ronald's Q##
      #First pull the estimates from Step 1
      ModelQ <- parTable(Model1_Results)
      
      #fix the indicator loadings from Step 1, free the direct effects of the SNP on the indicators, and fix the factor residual variance
      ModelQ$free <- c(rep(0, k+1), 1:(k*2), 0, 0) 
      
      #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
      ModelQ_Results <- sem(model = ModelQ, sample.cov = S_Full, estimator = "ML", sample.nobs=200) 
      
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
      Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta
      
      ##put the corrected standard error and Q in same dataset
      RESULTS[i,]=cbind(i,inspect(Model1_Results,"list")[k+1,-c(1,5:13)],se_c,Q)
      
      print(i)
    }
  }
  
  time_all<-proc.time()-time
  print(time_all[3])
  
  ##combine list of runs that produced errors/warnings
  RESULTS<-cbind(RESULTS, fail)
  
  return(list(RESULTS))
  
  
}
