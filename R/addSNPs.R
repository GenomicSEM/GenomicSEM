
addSNPs <-function(covstruc, SNPs, SNPSE = FALSE,parallel=FALSE,cores=NULL,GC="standard"){
  
  print("Please note that an update was made on 11/21/19 that combine addSNPs and the multivariate GWAS functions into a single step. Therefore, addSNPs is no longer a necessary function.")    
  warning("Please note that an update was made on 11/21/19 that combine addSNPs and the multivariate GWAS functions into a single step. Therefore, addSNPs is no longer a necessary function.")
  
  time<-proc.time()
  
  V_LD<-as.matrix(covstruc[[1]])
  S_LD<-as.matrix(covstruc[[2]])
  I_LD<-as.matrix(covstruc[[3]])
  
  SNPs<-data.frame(SNPs)
  beta_SNP<-SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))] 
  SE_SNP<-SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))] 
  
  #set univariate intercepts to 1 if estimated below 1
  diag(I_LD)<-ifelse(diag(I_LD)<= 1, 1, diag(I_LD))
  
  #enter in k for number of phenotypes 
  k<-ncol(beta_SNP)
  
  #f = number of SNPs in dataset
  f=nrow(beta_SNP) 
  
  #SNP variance (updated with 1KG phase 3 MAFs)
  varSNP=2*SNPs$MAF*(1-SNPs$MAF)  
  
  #small number because treating MAF as fixed
  if(SNPSE == FALSE){
    varSNPSE2=(.0005)^2
  }
  
  if(SNPSE != FALSE){
    varSNPSE2 = SNPSE^2
  }
  
  if(parallel == TRUE){
    if(is.null(cores)){
      ##if no default provided use 1 less than the total number of cores available so your computer will still function
      int <- parallel::detectCores() - 1
    }else{int<-cores}
    
    print("Creating Sampling Covariance [V] matrices")
    V_Full_List<-mclapply(X=1:f,FUN=function(X){
      
      i<-X
      
      #create empty shell of V_SNP matrix
      V_SNP<-diag(k)
      
      ##pull the coordinates of the I_LD matrix to loop making the V_SNP matrix
      coords<-which(I_LD != 'NA', arr.ind= T)
      
      #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
      
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
      
      ##create shell of full sampling covariance matrix
      V_Full<-.get_V_full(k, V_LD, varSNPSE2, V_SNP)
      
      k2<-nrow(V_Full)
      smooth2<-ifelse(eigen(V_Full)$values[k2] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)
      
      ##store the full V to a list of V_full matrices
      V_Full
      
    },mc.cores=int)
    
    print("Creating Genetic Covariance [S] matrices")
    S_Full_List<-mclapply(X=1:f,FUN=function(X){
      
      i<-X
      
      #create empty vector for S_SNP
      S_SNP<-vector(mode="numeric",length=k+1)
      
      #enter SNP variance from reference panel as first observation
      S_SNP[1]<-varSNP[i]
      
      #enter SNP covariances (standardized beta * SNP variance from refference panel)
      for (p in 1:k) {
        S_SNP[p+1]<-varSNP[i]*beta_SNP[i,p]
      }
      
      #create shell of the full S (observed covariance) matrix
      S_Full<-diag(k+1)
      
      ##add the LD portion of the S matrix
      S_Full[(2:(k+1)),(2:(k+1))]<-S_LD
      
      ##add in observed SNP variances as first row/column
      S_Full[1:(k+1),1]<-S_SNP
      S_Full[1,1:(k+1)]<-t(S_SNP)
      
      ##pull in variables names specified in LDSC function and name first column as SNP
      colnames(S_Full)<-c("SNP", colnames(S_LD))
      
      ##name rows like columns
      rownames(S_Full)<-colnames(S_Full)
      
      ##smooth to near positive definite if either V or S are non-positive definite
      ks<-nrow(S_Full)
      smooth1<-ifelse(eigen(S_Full)$values[ks] <= 0, S_Full<-as.matrix((nearPD(S_Full, corr = FALSE))$mat), S_Full<-S_Full)
      
      ##store the full S to a list of S_full matrices
      S_Full
    },mc.cores=int)}
  
  if(parallel == FALSE){
    #make empty matrices for S_full
    S_Full_List<-vector(mode="list",length=f)
    V_Full_List<-vector(mode="list",length=f)
    
    
    for (i in 1:f) {
      
      #create empty vector for S_SNP
      S_SNP<-vector(mode="numeric",length=k+1)
      
      #enter SNP variance from reference panel as first observation
      S_SNP[1]<-varSNP[i]
      
      #enter SNP covariances (standardized beta * SNP variance from refference panel)
      for (p in 1:k) {
        S_SNP[p+1]<-varSNP[i]*beta_SNP[i,p]
      }
      
      #create shell of the full S (observed covariance) matrix
      S_Full<-diag(k+1)
      
      ##add the LD portion of the S matrix
      S_Full[(2:(k+1)),(2:(k+1))]<-S_LD
      
      ##add in observed SNP variances as first row/column
      S_Full[1:(k+1),1]<-S_SNP
      S_Full[1,1:(k+1)]<-t(S_SNP)
      
      ##pull in variables names specified in LDSC function and name first column as SNP
      colnames(S_Full)<-c("SNP", colnames(S_LD))
      
      ##name rows like columns
      rownames(S_Full)<-colnames(S_Full)
      
      ##smooth to near positive definite if either V or S are non-positive definite
      ks<-nrow(S_Full)
      smooth1<-ifelse(eigen(S_Full)$values[ks] <= 0, S_Full<-as.matrix((nearPD(S_Full, corr = FALSE))$mat), S_Full<-S_Full)
      
      ##store the full S to a list of S_full matrices
      S_Full_List[[i]]<-S_Full
      
      #create empty shell of V_SNP matrix
      V_SNP<-diag(k)
      
      ##pull the coordinates of the I_LD matrix to loop making the V_SNP matrix
      coords<-which(I_LD != 'NA', arr.ind= T)
      
      #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
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
      
      ##store the full V to a list of V_full matrices
      V_Full_List[[i]]<-V_Full
      
      if(i == 1){
        print(i)
      }else{
        if(i %% 10000==0) {
          print(i)
        }}
    }}
  
  ##save the rsnumbers, MAF, A1/A2, and BP
  SNPs2<-SNPs[,1:6]
  
  time_all<-proc.time()-time
  print(time_all[3])
  
  return(Output <- list(V_Full=V_Full_List,S_Full=S_Full_List,RS=SNPs2))
  
}
