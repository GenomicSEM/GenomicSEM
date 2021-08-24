userGWAS <- function(covstruc=NULL, SNPs=NULL, estimation="DWLS", model="", modelchi=TRUE, printwarn=TRUE,
sub=FALSE,cores=NULL, toler=FALSE, SNPSE=FALSE, parallel=TRUE, GC="standard", MPI=FALSE,
smooth_check=FALSE, TWAS=FALSE, std.lv=FALSE){
  time<-proc.time()
  Operating<-Sys.info()[['sysname']]
  if(parallel == TRUE & Operating == "Windows"){
    stop("Parallel processing is not currently available for Windows operating systems. Please set the parallel argument to FALSE, or switch to a Linux or Mac operating system.")
  }

  if(exists("Output")){
    stop("Please note that an update was made to commonfactorGWAS on 4/1/21 so that addSNPs output CANNOT be fed directly to the function. It now expects the 
            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments.")
  }

  ##determine if the model is likely being listed in quotes and print warning if so
  test<-c(str_detect(model, "~"),str_detect(model, "="),str_detect(model, "\\+"))
  if(all(test) == FALSE){
    warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned stopped running after returning an error.")
  }

  print("Please note that an update was made to userGWAS on 11/21/19 so that it combines addSNPs and userGWAS.")

  if(class(SNPs)[1] == "character"){
    print("You are likely listing arguments in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS. The current version of the function is faster and saves memory. It expects the 
          output from ldsc followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usag")
    warning("You are likely listing arguments (e.g. Output = ...) in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS. The current version of the function is faster and saves memory. It expects the 
            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usage")
  }else{
    if(is.null(SNPs) | is.null(covstruc)){
      print("You may be listing arguments in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS;if you already did this and the function ran then you can disregard this warning. The current version of the function is faster and saves memory. It expects the 
            output from ldsc followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usag")
      warning("You may be listing arguments (e.g. Output = ...) in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS; ; if you already did this and the function ran then you can disregard this warning. The current version of the function is faster and saves memory. It expects the 
              output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usage")
    }
  }

  #remove white spacing on subset argument so will exact match lavaan representation of parameter
  if(!(sub[[1]])==FALSE){
    sub<-str_replace_all(sub, fixed(" "), "")
  }

  ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
  SNPs<-data.frame(SNPs)

  if (TWAS) {
    SNPs$Gene<-as.character(SNPs$Gene)
    SNPs$Panel<-as.character(SNPs$Panel)
    varSNP=SNPs$HSQ
  } else {
    SNPs$A1<-as.character(SNPs$A1)
    SNPs$A2<-as.character(SNPs$A2)
    SNPs$SNP<-as.character(SNPs$SNP)

    #SNP variance
    varSNP=2*SNPs$MAF*(1-SNPs$MAF)
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

  beta_SNP<-SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))]
  SE_SNP<-SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))]

  #enter in k for number of phenotypes
  k<-ncol(beta_SNP)

  #set univariate intercepts to 1 if estimated below 1
  diag(I_LD)<-ifelse(diag(I_LD)<= 1, 1, diag(I_LD))

  Model1<-model


  coords<-which(I_LD != 'NA', arr.ind= T)
  i<-1
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

  kv<-nrow(V_Full)
  smooth2<-ifelse(eigen(V_Full)$values[kv] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)

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
  if(TWAS == FALSE){
    colnames(S_Full)<-c("SNP", colnames(S_LD))
  }

  if(TWAS == TRUE){
    colnames(S_Full)<-c("Gene", colnames(S_LD))
  }

  ##name rows like columns
  rownames(S_Full)<-colnames(S_Full)

  ##smooth to near positive definite if either V or S are non-positive definite
  ks<-nrow(S_Full)
  smooth1<-ifelse(eigen(S_Full)$values[ks] <= 0, S_Full<-as.matrix((nearPD(S_Full, corr = FALSE))$mat), S_Full<-S_Full)

  k2<-ncol(S_Full)

  ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  for (i in 1) {

    if(toler==FALSE){
      W<- solve(V_Full)
    }

    if(toler!=FALSE){
      W <- solve(V_Full,tol=toler)
    }

    if(std.lv == FALSE){
      test2 <- .tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Full, estimator = "DWLS",
                             WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf,optim.force.converged=TRUE,
                             control=list(iter.max=1)))
    }
    if(std.lv == TRUE){
      test2 <- .tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Full, estimator = "DWLS",
                                                 WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf,
                                                 optim.force.converged=TRUE,control=list(iter.max=1),std.lv=TRUE))
    }

    order <- .rearrange(k = k2, fit = ReorderModel, names = rownames(S_Full))

    suppressWarnings(df<-lavInspect(ReorderModel, "fit")["df"])
    suppressWarnings(npar<-lavInspect(ReorderModel, "fit")["npar"])

  }

  if(!(TWAS)){
    SNPs2<-SNPs[,1:6]
  } else {
    SNPs2<-SNPs[,1:3]
  }

  rm(SNPs)

  if(parallel==FALSE){
    #make empty list object for model results if not saving specific model parameter
    if(sub[[1]]==FALSE){
      Results_List<-vector(mode="list",length=nrow(beta_SNP))
    }

    if(TWAS){
      print("Starting TWAS Estimation")
    } else {
      print("Starting GWAS Estimation")
    }

    for (i in 1:nrow(beta_SNP)) {
      if(i == 1){
        cat(paste0("Running Model: ", i, "\n"))
      }else{
        if(i %% 1000==0) {
          cat(paste0("Running Model: ", i, "\n"))
        }}
      final2 <- .userGWAS_analysis(i, k, SNPs2, beta_SNP, SE_SNP, varSNP, GC, coords, smooth_check, TWAS, printwarn)
      if(!(sub[[1]])==FALSE){
        final3<-as.data.frame(matrix(NA,ncol=ncol(final2),nrow=length(sub)))
        final3[1:length(sub),]<-final2[1,]
        if(i == 1){
          Results_List<-vector(mode="list",length=length(sub))
          for(y in 1:length(sub)){
            Results_List[[y]]<-as.data.frame(matrix(NA,ncol=ncol(final3),nrow=f))
            colnames(Results_List[[y]])<-colnames(final3)
            Results_List[[y]][1,]<-final3[y,]
          }
        }else{
          for(y in 1:nrow(final3)){
            Results_List[[y]][i,]<-final3[y,]
          }
        }
      }else{
        Results_List[[i]]<-final2
      }
    }

    time_all<-proc.time()-time
    print(time_all[3])

    return(Results_List)

  }

  if(parallel == TRUE & Operating != "Windows"){

    if(is.null(cores)){
      ##if no default provided use 1 less than the total number of cores available so your computer will still function
      int <- detectCores() - 1
    }else{
      int<-cores
    }

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

    SNPs2<-suppressWarnings(split(SNPs2,1:int))
    #split the V_SNP and S_SNP matrices into as many (cores - 1) as are aviailable on the local computer
    beta_SNP<-suppressWarnings(split(beta_SNP,1:int))
    SE_SNP<-suppressWarnings(split(SE_SNP,1:int))
    varSNP<-suppressWarnings(split(varSNP,1:int))

    if(TWAS){
      print("Starting TWAS Estimation")
    } else {
      print("Starting TWAS Estimation")
    }


    results<-foreach(n = icount(int), .combine = 'rbind') %:%

    foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = "lavaan") %dopar% {

      .userGWAS_analysis(i, k, SNPs2[[n]], beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]],
                         GC, coords, smooth_check, TWAS, printwarn)

    }

    ##sort results so it is in order of the output lists provided for the function
    results<- results[order(results$i, results$n),]
    results$n<-NULL

    if(!(sub[[1]])==FALSE){
      results$i<-NULL
      Results_List<-vector(mode="list",length=length(sub))
      for(y in 1:length(sub)){
        Results_List[[y]]<-as.data.frame(matrix(NA,ncol=ncol(results),nrow=nrow(results)/length(sub)))
        colnames(Results_List[[y]])<-colnames(results)
        Results_List[[y]]<-subset(results, paste0(results$lhs, results$op, results$rhs, sep = "") %in% sub[[y]] | is.na(results$lhs))
      }
      rm(results)
    }

    if(TWAS == FALSE){
      if(sub[[1]]==FALSE){
        names<-unique(results$SNP)
        Results_List<-vector(mode="list",length=length(names))
        for(y in 1:length(names)){
          Results_List[[y]]<-subset(results, results$SNP == names[[y]])
          Results_List[[y]]$Model_Number<-NULL
        }
        rm(results)
        rm(names)
      }
    }


    if(TWAS == TRUE){
      if(sub[[1]]==FALSE){
        names<-unique(results$Panel)
        Results_List<-vector(mode="list",length=length(names))
        for(y in 1:length(names)){
          Results_List[[y]]<-subset(results, results$Panel == names[[y]])
          Results_List[[y]]$Model_Number<-NULL
        }
        rm(results)
        rm(names)
      }
    }

    time_all<-proc.time()-time
    print(time_all[3])
    return(Results_List)

  }
}
