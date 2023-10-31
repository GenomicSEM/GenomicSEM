#--------------------------------------------------------------------------------------------------------------------------
# Author: Javier de la Fuente
# Date: 10-31-2023
#
# Filename: paLDSC.R
#
# Purpose: Defining a function to perform Parallel Analysis (PA) on LDSC genetic correlation matrices. 
# The method compares the eigenvalues generated from the eigen decomposition of the LDSC genetic 
# correlation matrix to the eigenvalues of a Monte-Carlo simulated null correlation matrix with random noise drawn from 
# the multivariate LDSC sampling distribution V. The suggested number of factors to be extracted
# corresponds with the last component with a larger eigenvalue than the same component in the null
# correlation matrix.
# The mandatory arguments of the function are S and V, corresponding either
# with the genetic correlation and standardized sampling distribution matrices (i.e., S_Stand and V_Stand matrices
# from the LDSC output), or the genetic covariance and unstandardized sampling distribution matrices from the LDSC output (S and V matrices, 
# respectively, which can be obtained by setting stand = TRUE in the LDSC function). 
# The following optional arguments are also implemented:
#  - r: defines the number of replications for the Monte-Carlo simulations of null correlation matrices (500 by default).
#  - p: defines the percentile for the simulated eigen-values.
#  - diag: defines whether diagonallized PA should be conducted (FALSE by default). Diagonallized PA assumes uncorrelated sampling
#          errors for the simulated null correlation matrices, as would happen if the data were pure noise. 
#          Thus, if diag = TRUE, the function will only use the diagonal of the sampling distribution V (i.e., the sampling variances
#          or the phenotypes in the original genetic correlation matrix) to introduce variation in the simulated null correlation matrices.
#          If diag = FALSE (by default), the whole multivariate sampling distribution LDSC V matrix will be used to simulate the null 
#          correlation matrices, thus assuming correlated sampling errors.
#  - fa: defines whether the eigenvalues should also be computed from a common factor solution from a explotory factor analysis of 
#        n factors (by default 1) using the factor method fm (by default minimum residual). By default this argument is set to FALSE, since
#        there is some degree of variation on the factor analysis eigenvalues depending on the number of factors extracted.
#  - fm: factor method for the exploratory factor analysis if fa = T (by default minimum residual).
#  - nfactors: number of factors to be extracted in the exploratory factor analysis if fa = T (by default = 1).
#  - save.pdf: whether the scree-plots derived from the PA function should be saved into a .pdf file.
#--------------------------------------------------------------------------------------------------------------------------

paLDSC <- function(S = S, V = V, r = NULL, p = NULL, save.pdf = F, diag = F, fa = F,
                   fm = NULL, nfactors = NULL) {
  list.of.packages <- c("ggplot2", "MASS","matrixStats","gdata","psych","matrixStats","egg","ggpubr","Matrix")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  invisible(lapply(list.of.packages, library,character.only = TRUE))
  if (is.null(r)){
    r <- 500
  }
  if (is.null(p)){
    p <- .95
  }
  if (is.null(fm)){
    fm <- "minres"
  }
  if (is.null(nfactors)){
    nfactors <- 1
  }
  
  #---- Get Dimensions of Matrices ----#
  k=dim(S)[1] #k phenotypes
  kstar=k*(k+1)/2 #kstar unique variances/covariances
  Svec=lowerTriangle(S,diag=T) #vectorize S
  SNULL=(0*S)
  diag(SNULL)=diag(S)
  SNULLvec=lowerTriangle(SNULL,diag=T) #vectorize S null
  #---- Parallel Analysis ----#
  EIG=as.data.frame(matrix(NA,nrow=k,ncol=r))
  for (i in 1:r) {
    Sample_null=mvrnorm(n=1,mu=SNULLvec,Sigma=V) #Simulate a null vectorized correlation matrix with noise drawn from the multivariate sample distribution
    Sample_null_M=matrix(0,ncol=k,nrow=k) #Turn it back into a matrix
    lowerTriangle(Sample_null_M,diag=T)=Sample_null
    upperTriangle(Sample_null_M,diag=F)=upperTriangle(t(Sample_null_M))
    EIG[,i]=eigen(Sample_null_M)$values #Store the eigenvalues
    cat("Running parallel analysis. Replication number: ",i,"\n")
  }
  #---- Find the p percentile for permuted and observed eigenvalues ----#
  Parallel_values=rowQuantiles(as.matrix(EIG),probs=p)
  Observed_PCA_values=eigen(S)$values
  #Create data frame from observed eigenvalue data
  obs = data.frame(Observed_PCA_values)
  obs$type = c('Observed Data')
  obs$num = c(row.names(obs))
  obs$num = as.numeric(obs$num)
  colnames(obs) = c('eigenvalue', 'type', 'num')
  #Create data frame from permuted data
  simPA = data.frame(Parallel_values)
  simPA$type = paste("Simulated data (",(p*100),"th %ile)",sep = "")
  simPA$num = c(row.names(obs))
  simPA$num = as.numeric(simPA$num)
  colnames(simPA) = c('eigenvalue', 'type', 'num')
  eigendatPA = rbind(obs,simPA)
  nfactPA <- min(which((eigendatPA[1:k,1] < eigendatPA[(k+1):(k*2),1] ) == TRUE))-1
  #Create vector for observed minus permuted data and # of factors
  obsPA <- data.frame(obs[1]-simPA[1])
  obsPA$type = paste("Observed minus simulated data (",(p*100),"th %ile)",sep = "")
  obsPA$num = c(row.names(obs))
  obsPA$num = as.numeric(obsPA$num)
  colnames(obsPA) = c('eigenvalue', 'type', 'num')
  nfactobsPA <- which(obsPA < 0)[1]-1
  
  if (diag == T) { 
    Vd_stand=0*V
    diag(Vd_stand)=diag(V)
    #---- Diagonalized Parallel Analysis ----#
    EIGd=as.data.frame(matrix(NA,nrow=k,ncol=r))
    for (i in 1:r) {
      Sample_null=mvrnorm(n=1,mu=SNULLvec,Sigma=Vd_stand) #Simulate a null vectorized correlation matrtix with noise drawn from the multivariate sample distribution of S
      Sample_null_M=matrix(0,ncol=k,nrow=k) #Turn it back into a matrix
      lowerTriangle(Sample_null_M,diag=T)=Sample_null
      upperTriangle(Sample_null_M,diag=F)=upperTriangle(t(Sample_null_M))
      EIGd[,i]=eigen(Sample_null_M)$values #Store the eigen values
      cat("Running diagonalized parallel analysis. Replication number: ",i,"\n")
    }
    
    #---- Find the p percentil for diagonalized permuted data----#
    Paralleld_values=rowQuantiles(as.matrix(EIGd),probs=p)
    #Create data frame from diagonalized permuted data
    simPAd = data.frame(Paralleld_values)
    simPAd$type = paste("Simulated data (",(p*100),"th %ile)",sep = "")
    simPAd$num = c(row.names(obs))
    simPAd$num = as.numeric(simPAd$num)
    colnames(simPAd) = c('eigenvalue', 'type', 'num')
    eigendatPAd = rbind(obs,simPAd)
    nfactPAd <- min(which((eigendatPAd[1:k,1] < eigendatPAd[(k+1):(k*2),1] ) == TRUE))-1
    #Create vector for observed minus diagonalized permuted data and suggested number of components to be extracted
    obsPAd <- data.frame(obs[1]-simPAd[1])
    obsPAd$type = paste("Observed minus simulated data (",(p*100),"th %ile)",sep = "")
    obsPAd$num = c(row.names(obs))
    obsPAd$num = as.numeric(obsPAd$num)
    colnames(obsPAd) = c('eigenvalue', 'type', 'num')
    nfactobsPAd <- which(obsPAd < 0)[1]-1
  }  
  
  #--- If eigenvalues from FA are required ---#
  if (fa == T) {
    Ssmooth<-as.matrix((nearPD(S, corr = T))$mat) #Smooth S matrix
    Observed_FA_values=fa(Ssmooth, fm = fm, nfactors = nfactors,SMC = FALSE, warnings = FALSE, rotate = "none")$values
    #Create data frame for observed FA eigenvalue data
    obsFA = data.frame(Observed_FA_values)
    obsFA$type = c('Observed Data')
    obsFA$num = c(row.names(obsFA))
    obsFA$num = as.numeric(obsFA$num)
    colnames(obsFA) = c('eigenvalue', 'type', 'num')
    #---- FA Parallel Analysis ----#
    EIGfa=as.data.frame(matrix(NA,nrow=k,ncol=r))
    for (i in 1:r) {
      Sample_null=mvrnorm(n=1,mu=SNULLvec,Sigma=V) #Simulate a null vectorized correlation matrtix with noise drawn from the multivariate sample distribution of S
      Sample_null_M=matrix(0,ncol=k,nrow=k) #Turn it back into a matrix
      lowerTriangle(Sample_null_M,diag=T)=Sample_null
      upperTriangle(Sample_null_M,diag=F)=upperTriangle(t(Sample_null_M))
      Ssmooth<-as.matrix((nearPD(Sample_null_M, corr = T))$mat)
      EIGfa[,i]=fa(Ssmooth, fm = fm, nfactors = nfactors,SMC = FALSE, warnings = FALSE, rotate = "none")$values #store the eigen values
      cat("Running FA parallel analysis. Replication number: ",i,"\n")
    } 
    
    #---- Find the quantile for FA PA, diagonalized FA PA, and observed eigen values ----#
    Parallel_fa_values=rowQuantiles(as.matrix(EIGfa),probs=p)
    #Create data frame from simulated PA FA data
    simPAfa = data.frame(Parallel_fa_values)
    simPAfa$type = paste("Simulated data (",(p*100),"th %ile)",sep = "")
    simPAfa$num = c(row.names(obs))
    simPAfa$num = as.numeric(simPAfa$num)
    colnames(simPAfa) = c('eigenvalue', 'type', 'num')
    eigendatPAfa = rbind(obsFA,simPAfa)
    nfactPAfa <- min(which((eigendatPAfa[1:k,1] < eigendatPAfa[(k+1):(k*2),1]) == TRUE))-1
    if(nfactPAfa==Inf){
      nfactPAfa <- 1
    }
    #Create vector for observed minus fa PA simulated data and # of factors
    obsPAfa <- data.frame(obsFA[1]-simPAfa[1])
    obsPAfa$type = paste("Observed minus simulated data (",(p*100),"th %ile)",sep = "")
    obsPAfa$num = c(row.names(obsFA))
    obsPAfa$num = as.numeric(obsPAfa$num)
    colnames(obsPAfa) = c('eigenvalue', 'type', 'num')
    nfactobsPAfa <- which(obsPAfa < 0)[1]-1
    
    if (diag == T){ 
      #---- FA Diagonalized Parallel Analysis ----#
      EIGdfa=as.data.frame(matrix(NA,nrow=k,ncol=r))
      for (i in 1:r) {
        Sample_null=mvrnorm(n=1,mu=SNULLvec,Sigma=Vd_stand) #simulate a null vectorized correlation matrtix with noise drawn from the multivariate sample distribution of S
        Sample_null_M=matrix(0,ncol=k,nrow=k) #turn it back into a matrix
        lowerTriangle(Sample_null_M,diag=T)=Sample_null
        upperTriangle(Sample_null_M,diag=F)=upperTriangle(t(Sample_null_M))
        Ssmooth<-as.matrix((nearPD(Sample_null_M, corr = T))$mat)
        EIGdfa[,i]=fa(Ssmooth, fm = fm, nfactors =  nfactors,SMC = FALSE, warnings = FALSE,rotate = "none")$values #store the eigen values    cat("Running diagonalized parallel analysis. Replication number: ",i,"\n")
        cat("Running FA diagonalized parallel analysis. Replication number: ",i,"\n")
      }
      #---- Find the quantile for FA PA, diagonalized FA PA, and observed eigen values ----#
      Parallel_dfa_values=rowQuantiles(as.matrix(EIGdfa),probs=p)
      #Create data frame from diagonalized PA FA simulated data
      simPAdfa = data.frame(Parallel_dfa_values)
      simPAdfa$type = paste("Simulated data (",(p*100),"th %ile)",sep = "")
      simPAdfa$num = c(row.names(obs))
      simPAdfa$num = as.numeric(simPAdfa$num)
      colnames(simPAdfa) = c('eigenvalue', 'type', 'num')
      eigendatPAdfa = rbind(obsFA,simPAdfa)
      nfactPAdfa <- min(which((eigendatPAdfa[1:k,1] < eigendatPAdfa[(k+1):(k*2),1] ) == TRUE))-1
      if(nfactPAdfa==Inf){
        nfactPAdfa <- 1
      }
      #Create vector for observed minus diagonalized permuted data and # of factors
      obsPAdfa <- data.frame(obsFA[1]-simPAdfa[1])
      obsPAdfa$type = paste("Observed minus simulated data (",(p*100),"th %ile)",sep = "")
      obsPAdfa$num = c(row.names(obsFA))
      obsPAdfa$num = as.numeric(obsPAdfa$num)
      colnames(obsPAdfa) = c('eigenvalue', 'type', 'num')
      nfactobsPAdfa <- which(obsPAdfa < 0)[1]-1
    }
  }
  #---- Scree plots ----#
  apatheme=theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          text=element_text(family='serif'),
          legend.title=element_blank(),
          legend.position=c(.7,.8),
          axis.line.x = element_line(color='black'),
          axis.line.y = element_line(color='black'))
  #---- Scree plot ----#
  if (isTRUE(all(diag(S) != rep(1, nrow(S))))){ # Verifica si la condición es verdadera
    pPA = ggplot(eigendatPA, aes(x = num, y = eigenvalue, shape = type)) +
      geom_line() +
      geom_point(size = 3) +
      scale_y_continuous(name = 'Eigenvalue') +
      scale_x_continuous(name = 'Component Number', breaks = min(1:k):max(1:k)) +
      scale_shape_manual(values = c(2, 17)) +
      geom_vline(xintercept = nfactPA, linetype = 'dashed') +
      apatheme
  } else {
    pPA = ggplot(eigendatPA, aes(x = num, y = eigenvalue, shape = type)) +
      geom_line() +
      geom_point(size = 3) +
      scale_y_continuous(name = 'Eigenvalue') +
      scale_x_continuous(name = 'Component Number', breaks = min(1:k):max(1:k)) +
      scale_shape_manual(values = c(2, 17)) +
      geom_hline(yintercept = 1) +
      apatheme
  }
  #---- Scree plot Observed minus permuted ----#
  pobsPA = ggplot(obsPA, aes(x=num, y=eigenvalue, shape=type)) +
    geom_line()+
    geom_point(size=3)+
    scale_y_continuous(name='Difference in Eigenvalues') +
    scale_x_continuous(name='Component Number', breaks=min(1:k):max(1:k))+
    scale_shape_manual(values=c(17)) +
    geom_vline(xintercept = nfactobsPA, linetype = 'dashed')+
    geom_hline(yintercept = 0)+
    apatheme
  
  if (diag == T){
    #---- Scree plot diagonalized PA ----#
    pPAd = ggplot(eigendatPAd, aes(x=num, y=eigenvalue, shape=type)) +
      geom_line()+
      geom_point(size=3)+
      scale_y_continuous(name='Diagonalized Eigenvalue')+
      scale_x_continuous(name='Component Number', breaks=min(1:k):max(1:k))+
      scale_shape_manual(values=c(2,17)) +
      geom_vline(xintercept = nfactPAd, linetype = 'dashed')+
      geom_hline(yintercept = 1)+
      apatheme
    #---- Scree plot Observed minus diagonalized permuted data ----#
    pobsPAd = ggplot(obsPAd, aes(x=num, y=eigenvalue, shape=type)) +
      geom_line()+
      geom_point(size=3)+
      scale_y_continuous(name='Difference in diagonalized Eigenvalues')+
      scale_x_continuous(name='Component Number', breaks=min(1:k):max(1:k))+
      scale_shape_manual(values=c(17)) +
      geom_vline(xintercept = nfactobsPAd, linetype = 'dashed')+
      geom_hline(yintercept = 0)+
      apatheme
  }
  
  if (isTRUE(fa)){
    #---- Scree plot FA PA ----#
    pPAfa = ggplot(eigendatPAfa, aes(x=num, y=eigenvalue, shape=type)) +
      geom_line()+
      geom_point(size=3)+
      scale_y_continuous(name='Eigenvalue from FA Solution')+
      scale_x_continuous(name='Factor Number', breaks=min(1:k):max(1:k))+
      scale_shape_manual(values=c(2,17)) +
      geom_vline(xintercept = nfactPAfa, linetype = 'dashed')+
      geom_hline(yintercept = 1)+
      apatheme
    #---- Scree plot Observed minus FA PA ----#
    pobsPAfa = ggplot(obsPAfa, aes(x=num, y=eigenvalue, shape=type)) +
      geom_line()+
      geom_point(size=3)+
      scale_y_continuous(name='Difference in FA Eigenvalues') +
      scale_x_continuous(name='Factor Number', breaks=min(1:k):max(1:k))+
      scale_shape_manual(values=c(17)) +
      geom_vline(xintercept = nfactobsPAfa, linetype = 'dashed')+
      geom_hline(yintercept = 0)+
      apatheme
    
    if (diag == T){
      #---- Scree plot FA diagonalized PA ----#
      pPAdfa = ggplot(eigendatPAdfa, aes(x=num, y=eigenvalue, shape=type)) +
        geom_line()+
        geom_point(size=3)+
        scale_y_continuous(name='Eigenvalue from diagonalized FA Solution')+
        scale_x_continuous(name='Factor Number', breaks=min(1:k):max(1:k))+
        scale_shape_manual(values=c(2,17)) +
        geom_vline(xintercept = nfactPAdfa, linetype = 'dashed')+
        geom_hline(yintercept = 1)+
        apatheme
      #---- Scree plot Observed minus FA diagonalized PA ----#
      pobsPAdfa = ggplot(obsPAdfa, aes(x=num, y=eigenvalue, shape=type)) +
        geom_line()+
        geom_point(size=3)+
        scale_y_continuous(name='Difference in diagonalized FA Eigenvalues')+
        scale_x_continuous(name='Factor Number', breaks=min(1:k):max(1:k))+
        scale_shape_manual(values=c(17)) +
        geom_vline(xintercept = nfactobsPAdfa, linetype = 'dashed')+
        geom_hline(yintercept = 0)+
        apatheme
    }
  }
  
  print(pPA)
  print(pobsPA)
  cat("------------------------------------------------------------------------","\n")
  cat("Parallel Analysis suggests extracting",nfactPA,"components","\n")
  cat("------------------------------------------------------------------------","\n")
  
  if (diag == T){
    print(pPAd)
    print(pobsPAd)
    cat("Diagonalized Parallel Analysis suggests extracting",nfactPAd,"components","\n")
    cat("------------------------------------------------------------------------","\n")
  }
  
  if (fa == T){
    print(pPAfa)
    print(pobsPAfa)
    cat("------------------------------------------------------------------------","\n")
    cat("Parallel Analysis suggests extracting",nfactPAfa,"factors","\n")
    cat("------------------------------------------------------------------------","\n")
    if (diag == T){
      print(pPAdfa)
      print(pobsPAdfa)
      cat("------------------------------------------------------------------------","\n")
      cat("Diagonalized Parallel Analysis suggests extracting",nfactPAdfa,"factors","\n")
      cat("------------------------------------------------------------------------","\n")
    }
  }
  
  
  if (isTRUE(save.pdf)){
    figurePCA <- ggarrange(pPA, pobsPA, 
                           labels = c(NA,NA),
                           ncol = 1, nrow = 2)
    ggexport(figurePCA, filename = "PA_LDSC.pdf",res = 300, verbose = F)
    if (diag == T){
      figurePCAd <- ggarrange(pPAd, pobsPAd, 
                              labels = c(NA,NA),
                              ncol = 1, nrow = 2)
      ggexport(figurePCAd, filename = "Diagonalized_PA_LDSC.pdf",res = 300, verbose = F)
    }
    if (fa == T){
      figureFA <- ggarrange(pPAfa, pobsPAfa, 
                            labels = c(NA,NA),
                            ncol = 1, nrow = 2)
      ggexport(figureFA, filename = "FA_PA_LDSC.pdf",res = 300, verbose = F)
      if (diag == T){
        figureFAd <- ggarrange(pPAdfa, pobsPAdfa, 
                               labels = c(NA,NA),
                               ncol = 1, nrow = 2)
        ggexport(figureFAd, filename = "FA_Diagonalized_PA_LDSC.pdf",res = 300, verbose = F)
      }
    }
  }
}
