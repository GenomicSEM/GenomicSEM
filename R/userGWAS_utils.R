# Title     : Minor supporting functions of userGWAS.R
# Objective : Primarily to prevent clutter in userGWAS.R
# Created by: mze210
# Created on: 20/08/2021

.get_V_SNP <- function(SE_SNP, I_LD, varSNP, GC, coords, k, i) {
     V_SNP<-diag(k)
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
    return(V_SNP)
}

.get_Z_pre <- function(i, beta_SNP, SE_SNP, I_LD, GC) {
    if(GC == "conserv"){
        Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*diag(I_LD))
    }
    if(GC=="standard"){
        Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*sqrt(diag(I_LD)))
    }
    if(GC=="none"){
        Z_pre<-beta_SNP[i,]/SE_SNP[i,]
    }
    return(Z_pre)
}

.get_V_full <- function(k, V_LD, varSNPSE2, V_SNP) {
    ##create shell of full sampling covariance matrix
    V_Full<-diag(((k+1)*(k+2))/2)

    ##input the ld-score regression region of sampling covariance from ld-score regression SEs
    V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD

    ##add in SE of SNP variance as first observation in sampling covariance matrix
    V_Full[1,1]<-varSNPSE2

    ##add in SNP region of sampling covariance matrix
    V_Full[2:(k+1),2:(k+1)]<-V_SNP
    return(V_Full)
}
