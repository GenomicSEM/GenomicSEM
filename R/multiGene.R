##multiGene takes as arguments:
#1.covstruc from ldsc
#2. Genes from sumstats
#3. LD matrix from package of users choice
#4. [optional] varGene_SE in case 1e-9 is too small 
#5. [optional] Genelist if the user wants to create all sumstats at once, but only subset certain rs numbers for multiGene



multiGene <-function(covstruc, Genes, LD, GeneSE = "F", Genelist = "F"){
  time<-proc.time()
  i = 1
  V_LD<-as.matrix(covstruc[[1]])
  S_LD<-as.matrix(covstruc[[2]])
  I_LD<-as.matrix(covstruc[[3]])
  
  Genes<-data.frame(Genes)
  
  LD<-as.matrix(LD)
  
  ##take Genes and go from long to wide? or just assume that the full dataset is all the Genes wanted...
  ##option to only pull Genes on a list, otherwise assume using full dataset
  if(Genelist != "F"){
    Genes<-subset(Genes, Genes$ID %in% Genelist)
  }
  
  beta_Genes<-Genes[,grep("beta.",fixed=TRUE,colnames(Genes))] 
  SE_Genes<-Genes[,grep("se.",fixed=TRUE,colnames(Genes))] 
  
  if(nrow(beta_Genes) != nrow(LD)){
    print("ERROR: The number of Genes in your dataset is not equal to the number of rows in your LD matrix. 
          Please check your input.")
  }
  
  #set univariate intercepts to 1 if estimated below 1
  diag(I_LD)<-ifelse(diag(I_LD)<= 1, 1, diag(I_LD))
  
  #enter in k for number of phenotypes 
  k<-ncol(beta_Genes)
  
  #f = number of Genes in dataset
  f=nrow(beta_Genes) 
  
  #Gene variance (heritability)
  varGenes =Genes$HSQ 
  
  #small number because treating MAF as fixed
  if(GeneSE == "F"){
    varGeneSE2=(.00000001)^2
  }
  
  ##if user provides own GeneSE use that instead 
  if(GeneSE != "F"){
    varGeneSE2 = GeneSE^2
  }
  
  #create shell of the full S (observed covariance) matrix = k phenotypes + f Genes
  S_Full<-diag(k+f)
  
  #enter Gene covariance with phenotypes (i.e., univariates sumstats) for Gene 1
  for (p in 1:k) {
    S_Full[p+f,1]<-varGenes[i]*beta_Genes[i,p]
  }
  
  u<-1
  e<-1
  #enter Gene covariance with phenotypes (i.e., univariates sumstats) for remaining Genes
  for(u in (i+1):f){
    for(e in 1:k){
      S_Full[e+f,u]<-varGenes[u]*beta_Genes[u,e]}}
  
  ##add the LDSC portion of the S matrix
  S_Full[((f+1):nrow(S_Full)),((f+1):nrow(S_Full))]<-S_LD
  
  #turn LD (Gene-Gene correlations) into LD_cov (Gene-Gene covariances)
  ##pull the coordinates of the I_LD matrix to loop making the V_Gene matrix
  LD_coords<-which(LD != 'NA', arr.ind= T)
  
  ##create empty matrix for LD between Genes
  LD_Cov<-diag(f)
  
  p<-1
  #loop to create LD_Cov with Gene variances on diagonal and Gene covariances on off-diagonal
  for (p in 1:nrow(LD_coords)) { 
    x<-LD_coords[p,1]
    y<-LD_coords[p,2]
    if (x != y) { 
      LD_Cov[x,y]<-(LD[x,y]*sqrt(varGenes[x])*sqrt(varGenes[y]))}
    if (x == y) {
      LD_Cov[x,x]<-varGenes[x]
    }
  }
  
  ##add in LD info across Genes 
  S_Full[1:f,1:f]<-LD_Cov
  
  #make S_Full symmetric
  makeSymm <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  
  S_Full<-makeSymm(S_Full)
  
  ##name Genes by rs#
  S_names<-Genes$Gene
  
  ##name columns of S_Full
  colnames(S_Full)<-c(as.character(S_names), colnames(S_LD))
  
  ##name rows like columns
  rownames(S_Full)<-colnames(S_Full)
  
  ##smooth to near positive definite if either V or S are non-positive definite
  ks<-nrow(S_Full)
  smooth1<-ifelse(eigen(S_Full)$values[ks] <= 0, S_Full<-as.matrix((nearPD(S_Full, corr = FALSE))$mat), S_Full<-S_Full)
  print("S matrix created")
  
  ##create shell of full sampling covariance matrix
  V_Full<-diag(((k+f)*(k+f+1))/2)
  
  ##create empty shell of V_Gene = length of V_Full - length of V_LD
  V_Gene<-diag(length(diag(V_Full))-length(diag(V_LD)))
  
  ##add in sampling variance of first Gene as first observation
  V_Gene[1,1]<-varGeneSE2
  
  ##fill in the remaining Gene SEs and Gene-Gene SEs (very small values = treating as fixed)
  if(f > 1){t<-f
  for(b in 2:f){
    ##fills in Gene1-Gene2, Gene1-Gene3, Gene1-GeneX
    V_Gene[b,b]<-varGeneSE2
    for(a in 1:(f+1-b)){
      u<-t+k
      V_Gene[u+1,u+1]<-varGeneSE2 
      if(a == (f+1-b)){t<-u+a}else{V_Gene[u+a+1,u+a+1]<-varGeneSE2}}}}
  
  ##pull coordinatoes of Gene-phenotype covariances (all remaining diagonals in V_Gene that = 1)     
  coords<-which(V_Gene == 1, arr.ind= T)
  
  print("Filling in diagonals of V matrix")
  ##fill in diagonals of V_Gene with Gene-phenotype sampling covariances
  r<-1
  for(l in 1:nrow(SE_Genes)){
    for(c in 1:nrow(I_LD)){
      x<-coords[r]
      V_Gene[x,x]<-(SE_Genes[l,c]*I_LD[c,c]*varGenes[l])^2
      r<-r+1
    }}
  
  ##fill in cross-trait within-Gene sampling covariance
  ##pull coordinates with length =
  ##number of unique off-diagonal elements within a Gene * the number of Genes = ((k*(k+1))/2-k)*f
  r<-1
  u<-1
  coords2=as.data.frame(matrix(NA,ncol=2,nrow=(((k*(k+1))/2-k)*f)))
  for(t in 1:(length(diag(V_Gene))-1)){
    for(u in 1:k){
      m<-t+u
      if(m <= length(diag(V_Gene))){
        if(varGeneSE2 %in% (diag(V_Gene)[t:m])){}else{
          coords2[r,]<-c(t,t+u)
          r<-r+1
        }
      }}}
  
  #correct the SEs ahead of time
  SE_Genes2<-as.data.frame(matrix(NA,ncol=k,nrow=f))
  b<-1
  for(b in 1:nrow(I_LD)){
    SE_Genes2[,b]<-SE_Genes[,b]*I_LD[b,b]
  }
  
  ##now multiply by Gene variance ahead of time
  c<-1
  for(c in 1:length(varGenes)){
    SE_Genes2[c,]<-SE_Genes2[c,]*varGenes[c]}
  
  ##pull coordinates in I_LD
  coords3<-data.frame(which(I_LD != 'NA', arr.ind= T))
  
  #take only off-diagonal coordinates of I_LD
  coords3<-subset(coords3, coords3$row < coords3$col)
  
  #turn off-diagonals of I_LD into a vector
  I_LD2<-I_LD[lower.tri(I_LD,diag=FALSE)]
  
  ##fill in actual pieces of V_Gene with cross-trait within-Gene sampling covariances
  print("Filling in cross-trait within-Gene sampling covariances of V matrix")
  r<-1
  u<-1
  p<-1
  for(p in 1:f){
    for(u in 1:((k*(k+1))/2-k)){
       x<-coords2[r,1]
      y<-coords2[r,2]
      x2<-coords3[u,1]
      y2<-coords3[u,2]
      V_SNP[y,x]<-(SE_Genes2[p,x2]*SE_Genes2[p,y2]*I_LD2[u])
      r<-r+1
      
    }
  }  
  
  
  
  if(sum(abs(LD[lower.tri(LD)])) > 0){
    ##get coordinates for cross-Gene within-trait sampling covariances
    #f*(f+1)/2-f = number of smaller matrices within V_Gene
    #*k = number of diagonal elements within the smaller matrices
    b<-(f*(f+1)/2-f)*k
    coords4=as.data.frame(matrix(NA,ncol=2,nrow=b))
    m<-1
    d<-k+1
    s<-1
    j<-1
    n<-k+1
    o<-1
    w<-1
    e<-1
    h<-1
    
    for(j in 1:(f-1)){
      if(j == 1){
        for(h in 1:k){
          ##takes k+1 and first observation
          coords4[m,]<-c(coords[d,1],coords[m,1])   
          d<-d+1
          m<-m+1
        }
      }
      s<-1
      if(j > 1){
        for(o in 1:k){
          w<-(k*j)+o
          coords4[m,]<-c(coords[w,1], coords[s,1])
          m<-m+1
          n<-s+k
          s<-s+1
          for(t in 1:(j-1)){
            coords4[m,]<-c(coords[w,1], coords[n,1])
            m<-m+1
            n<-(s-1)+k+t*k
          }}}}
    
    ##create vector of cross-Gene within-trait sampling covariances
    samp<-as.data.frame(matrix(NA,ncol=1,nrow=b))
    
    h<-1
    p<-1
    y<-1
    u<-1
    for(h in 1:f){
      for(p in 1:k){
        if(h == 1){
          samp[u,]<-SE_Genes2[h,p]*SE_Genes2[h+1,p]*LD[h,h+1]
          u<-u+1}
        if(h >= 3){
          for(y in 1:(h-1)){
            samp[u,]<-SE_Genes2[h,p]*SE_Genes2[y,p]*LD[h,y]
            u<-u+1
          }
        }else{}
      }
    }
    
    print("Filling in cross-Gene within-trait sampling covariances of V matrix")
    ##add in cross-Gene within-trait sampling covariances
    for(w in 1:nrow(coords4)){
      x<-coords4[w,1]
      y<-coords4[w,2]
      V_Gene[x,y]<-samp[w,]}
    
    print("Pulling coordinates of cross-Gene cross-trait sampling covariances in V matrix. 
          For large numbers of Genes or traits (e.g., > 50) this may take > 20 minutes.")
    ##pull coordinates of cross-Gene cross-trait sampling covariances
    coords6=as.data.frame(matrix(NA,ncol=2,nrow=(((k*k)-k)*(((f*f)-f)/2))))
    t<-1
    u<-1
    p<-1
    n<-nrow(V_Gene)
    for(u in 1:nrow(V_Gene)){
      print(c(u, "of", n))
      for(p in 1:nrow(V_Gene)){
        if(diag(V_Gene)[u] != varGeneSE2 & diag(V_Gene)[p] != varGeneSE2){
          if (u != p){
            if (u > p) {
              coords6[t,]<-c(u,p)
              t<-t+1
            }}}else{}
      }
    }
    
    coords7<-data.frame(coords4$V2,coords4$V1)
    coords8<-data.frame(coords2$V2,coords2$V1)
    colnames(coords7)<-c("V1", "V2")
    colnames(coords8)<-c("V1", "V2")
    
    ##remove redundant vectors from cross-Gene within-trait
    #and reordered versions of same position e.g., x,y -> y,x
    t<-anti_join(coords6,coords4, by = c("V1", "V2"))
    t2<-anti_join(t,coords7,by = c("V1", "V2"))
    t3<-anti_join(t2,coords2,by = c("V1", "V2"))
    t4<-anti_join(t3, coords8,by = c("V1", "V2"))
    
    #turn LD into a vector so its easier
    LD2<-LD[lower.tri(LD,diag=FALSE)]
    
    #turn CTI into a vector to match order of t4
    CTI<-as.vector(I_LD)
    
    #remove diagonal elements
    CTI<-subset(CTI, CTI < 1)
    
    u<-((f*f)-f)/2
    w<-1
    
    print("Filling in cross-Gene cross-trait sampling covariances of V matrix. 
          For large numbers of Genes or traits (e.g., > 50) this may take > 10 minutes.")
    ##add cross-Gene cross-trait sampling covariances to V_Gene
    for(l in 1:u){
      print(c(l, "of", u))
      for(n in 1:length(CTI)){
        x<-t4[w,1]
        y<-t4[w,2]
        V_Gene[x,y]<-sqrt(diag(V_Gene)[x])*sqrt(diag(V_Gene)[y])*CTI[n]*LD2[u]  
        w<-w+1
      }
    }
  }else{print("The LD among all Genes was 0. Please check this was intended")}
  
  ##input the ld-score regression region of sampling covariance from ld-score regression SEs
  V_Full[(nrow(V_Gene)+1):nrow(V_Full),(nrow(V_Gene)+1):nrow(V_Full)]<-V_LD
  
  ##add in Gene region of sampling covariance matrix
  V_Full[1:nrow(V_Gene),1:nrow(V_Gene)]<-V_Gene
  
  V_Full<-makeSymm(V_Full)
  
  k2<-nrow(V_Full)
  smooth2<-ifelse(eigen(V_Full)$values[k2] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)
  
  ##save the rsnumbers, MAF, A1/A2, and BP
  Genes2<-Genes[,c(1,4)]
  
  ##save in order that can be used with both usermodel and usergGWAS
  return(Output <- list(V_Full=V_Full,S_Full=S_Full,Info=Genes2))
  
  }
