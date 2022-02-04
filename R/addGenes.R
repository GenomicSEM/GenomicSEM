

addGenes <-function(covstruc, Genes, GC="standard"){
  time<-proc.time()
  
  V_LD<-as.matrix(covstruc[[1]])
  S_LD<-as.matrix(covstruc[[2]])
  I_LD<-as.matrix(covstruc[[3]])
  
  Genes<-data.frame(Genes)
  beta_Gene<-Genes[,grep("beta.",fixed=TRUE,colnames(Genes))] 
  SE_Gene<-Genes[,grep("se.",fixed=TRUE,colnames(Genes))] 
  
  #set univariate intercepts to 1 if estimated below 1
  diag(I_LD)<-ifelse(diag(I_LD)<= 1, 1, diag(I_LD))
  
  #enter in k for number of phenotypes 
  k<-ncol(beta_Gene)
  
  #f = number of Genes in dataset
  f=nrow(beta_Gene) 
  
  #make empty matrices for S_full
  S_Full_List<-vector(mode="list",length=f)
  V_Full_List<-vector(mode="list",length=f)
  
  #Gene variance (assume het HSQ is correct)
  varGene= Genes$HSQ
  
  #small number because treating MAF as fixed
  varGeneSE2=(.0005)^2
  
  #function to creat row/column names for S_full matrix
  write.names <- function(k, label = "V") {  
    varnames<-vector(mode="character",length=k+1)
    
    for (i in 1){
      varnames[1]<-c("Gene")}
    
    for (j in i:k) { 
      varnames[j+1]<-paste(label,j,sep="")}
    
    return(varnames)
  }
  
  S_names<-write.names(k=ncol(I_LD))
  
  for (i in 1:f) {
    
    #create empty vector for S_Gene
    S_Gene<-vector(mode="numeric",length=k+1)
    
    #enter Gene variance from reference panel as first observation
    S_Gene[1]<-varGene[i]
    
    #enter Gene covariances (standardized beta * Gene variance from refference panel)
    for (p in 1:k) {
      S_Gene[p+1]<-varGene[i]*beta_Gene[i,p]
    }
    
    #create shell of the full S (observed covariance) matrix
    S_Full<-diag(k+1)
    
    ##add the LD portion of the S matrix
    S_Full[(2:(k+1)),(2:(k+1))]<-S_LD
    
    ##add in observed Gene variances as first row/column
    S_Full[1:(k+1),1]<-S_Gene
    S_Full[1,1:(k+1)]<-t(S_Gene)
    
    ##name the columns/rows using the naming function defined outside of the loop
    rownames(S_Full) <- S_names
    colnames(S_Full) <- S_names
    
    ##smooth to near positive definite if either V or S are non-positive definite
    ks<-nrow(S_Full)
    smooth1<-ifelse(eigen(S_Full)$values[ks] <= 0, S_Full<-as.matrix((nearPD(S_Full, corr = FALSE))$mat), S_Full<-S_Full)
    
    ##store the full S to a list of S_full matrices
    S_Full_List[[i]]<-S_Full
    
    #create empty shell of V_Gene matrix
    V_Gene<-diag(k)
    
    ##pull the coordinates of the I_LD matrix to loop making the V_Gene matrix
    coords<-which(I_LD != 'NA', arr.ind= T)
    
    #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by Gene variance from reference panel
    if(GC == "conserv"){
    for (p in 1:nrow(coords)) { 
      x<-coords[p,1]
      y<-coords[p,2]
      if (x != y) { 
        V_Gene[x,y]<-(SE_Gene[i,y]*SE_Gene[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varGene[i]^2)}
      if (x == y) {
        V_Gene[x,x]<-(SE_Gene[i,x]*I_LD[x,x]*varGene[i])^2
      }
    }
    }
    if(GC == "standard"){
      for (p in 1:nrow(coords)) { 
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) { 
          V_Gene[x,y]<-(SE_Gene[i,y]*SE_Gene[i,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*varGene[i]^2)}
        if (x == y) {
          V_Gene[x,x]<-(SE_Gene[i,x]*sqrt(I_LD[x,x])*varGene[i])^2
        }
      }
    }
    
 
    
    if(GC == "none"){
      for (p in 1:nrow(coords)) { 
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) { 
          V_Gene[x,y]<-(SE_Gene[i,y]*SE_Gene[i,x]*I_LD[x,y]*varGene[i]^2)}
        if (x == y) {
          V_Gene[x,x]<-(SE_Gene[i,x]*varGene[i])^2
        }
      }
    }
    
    V_Full <- .get_V_full(k, V_LD, varGeneSE2, V_Gene)

    k2<-nrow(V_Full)
    smooth2<-ifelse(eigen(V_Full)$values[k2] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)
    
    
    ##store the full V to a list of V_full matrices
    V_Full_List[[i]]<-V_Full
    
    print(i)
  }

  ##save the GeneID and panel
  Genes2<-Genes[,c(1,2)]
  
  return(Output <- list(V_Full=V_Full_List,S_Full=S_Full_List,ID=Genes2))
  
}
