userGWAS <- function(covstruc=NULL, SNPs=NULL, estimation="DWLS", model="", printwarn=TRUE,
                     sub=FALSE,cores=NULL, toler=FALSE, SNPSE=FALSE, parallel=TRUE, GC="standard", MPI=FALSE,
                     smooth_check=FALSE, TWAS=FALSE, std.lv=FALSE,fix_measurement=TRUE,Q_SNP=FALSE){
  
  # Set toler to machine precision to enable passing this to solve() directly
  if (!toler) toler <- .Machine$double.eps
  .check_one_of(estimation, c("DWLS", "ML"))
  .check_boolean(printwarn)
  if (!is.null(cores)) .check_range(cores, min=0, max=Inf)
  .check_range(toler, min=0, max=Inf)
  .check_boolean(parallel)
  .check_one_of(GC, c("standard", "conserv", "none"))
  .check_boolean(MPI)
  .check_boolean(smooth_check)
  .check_boolean(TWAS)
  .check_boolean(std.lv)
  # Sanity checks finished
  time <- proc.time()
  Operating <- Sys.info()[['sysname']]
  
  if(MPI == TRUE & Operating == "Windows"){
    stop("MPI is not currently available for Windows operating systems. Please set the MPI argument to FALSE, or switch to a Linux or Mac operating system.")
  }
  
  if(exists("Output")){
    stop("Please note that an update was made to commonfactorGWAS on 4/1/21 so that addSNPs output CANNOT be fed directly to the function. It now expects the
            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments.")
  }
  
  ##determine if the model is likely being listed in quotes and print warning if so
  test <- c(str_detect(model, "~"),str_detect(model, "="),str_detect(model, "\\+"))
  if (!all(test)){
    warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned stopped running after returning an error.")
  }
  
  print("Please note that an update was made to userGWAS on Sept 1 2023  so that the default behavior is to fix the measurement model using the fix_measurement argument.")
  
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
  if(sub[[1]] != FALSE){
    sub <- str_replace_all(sub, fixed(" "), "")
  }
  
  ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
  SNPs <- data.frame(SNPs)
  
  if (TWAS) {
    SNPs$Gene <- as.character(SNPs$Gene)
    SNPs$Panel <- as.character(SNPs$Panel)
    varSNP <- SNPs$HSQ
  } else {
    SNPs$A1 <- as.character(SNPs$A1)
    SNPs$A2 <- as.character(SNPs$A2)
    SNPs$SNP <- as.character(SNPs$SNP)
    
    #SNP variance
    varSNP <- 2*SNPs$MAF*(1-SNPs$MAF)
  }
  
  #small number because treating MAF as fixed
  if(SNPSE == FALSE){
    varSNPSE2 <- (.0005)^2
  }else{varSNPSE2 <- SNPSE^2}

  V_LD <- as.matrix(covstruc[[1]])
  S_LD <- as.matrix(covstruc[[2]])
  I_LD <- as.matrix(covstruc[[3]])
  Model1 <- model
  
  if(fix_measurement){
      
      #name rownames as column names
      rownames(S_LD)<-colnames(S_LD)
      
      #split by lines of code
      lines <- strsplit(model, "\n")[[1]]
      
      # Use grep to find lines containing "SNP" or "Gene" and exclude them
      if(TWAS){
        filtered_lines <- lines[!grepl(c("Gene"), lines)]
      }else{filtered_lines <- lines[!grepl("SNP", lines)]}
      
      #remove ghost parameters as this could include SNP/Gene effects 
      filtered_lines<-filtered_lines[!grepl(c(":="), filtered_lines)]
      
      # Join the filtered lines back into a single text string
      noSNPmodel <- paste(filtered_lines, collapse = "\n")
      
      #smooth S and V matrices if necessary
      smoothS<-ifelse(eigen(S_LD)$values[nrow(S_LD)] <= 0, S_LD<-as.matrix((nearPD(S_LD, corr = FALSE))$mat), S_LD<-S_LD)
      smoothV<-ifelse(eigen(V_LD)$values[nrow(V_LD)] <= 0, V_LD<-as.matrix((nearPD(V_LD, corr = FALSE))$mat), V_LD<-V_LD)
      
      #estimate a no SNP model to get the estimates for the measurement model
      W <- solve(V_LD,tol=toler)
      
      #estimate with incorrectly ordered V to get internal representation of V
      testnoSNP <- .tryCatch.W.E(ReorderModelnoSNP <- sem(noSNPmodel, sample.cov = S_LD, estimator = "DWLS",
                                                          WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf, optim.force.converged=TRUE
                                                          ,control=list(iter.max=1),std.lv=std.lv))
      
      #obtain V ordering for this model
      order <- .rearrange(k = ncol(S_LD), fit = ReorderModelnoSNP, names = colnames(S_LD))
      
      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Reorder <- V_LD[order, order]
      u <- nrow(V_Reorder)
      W_Reorder<- diag(u)
      diag(W_Reorder) <- diag(V_Reorder)
      
      ##invert the reordered sampling covariance matrix to create a weight matrix
      W_Reorder <- solve(W_Reorder, tol=toler)
      
      #estimate the measurement model
      if(estimation == "DWLS"){
        ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
        emptynoSNP<-.tryCatch.W.E(Model1_Results <- sem(noSNPmodel, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2,optim.dx.tol = +Inf,std.lv=std.lv))
      }
      
      if(estimation == "ML"){
          emptynoSNP<-.tryCatch.W.E(Model1_Results <- sem(noSNPmodel, sample.cov = S_LD, estimator = "ML",  sample.nobs = 200,optim.dx.tol = +Inf,sample.cov.rescale=FALSE,std.lv=std.lv))
      }
      
      #pull the model results
      Model1<- parTable(Model1_Results)
      
      #fix any parameter that is not a (residual) variance
      for(p in 1:nrow(Model1)){
        Model1$free[p]<-ifelse(Model1$lhs[p] != Model1$rhs[p], 0, Model1$free[p]) 
      }
      
    }
  
  beta_SNP <- SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))]
  SE_SNP <- SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))]
  
  #enter in k for number of phenotypes
  n_phenotypes <- ncol(beta_SNP)
  
  #set univariate intercepts to 1 if estimated below 1
  diag(I_LD) <- ifelse(diag(I_LD)<= 1, 1, diag(I_LD))
  
  coords <- which(I_LD != 'NA', arr.ind= T)
  
  V_SNP <- .get_V_SNP(SE_SNP, I_LD, varSNP, GC, coords, n_phenotypes, 1)
  
  V_full <- .get_V_full(n_phenotypes, V_LD, varSNPSE2, V_SNP)
  
  kv <- nrow(V_full)
  smooth2 <- ifelse(eigen(V_full)$values[kv] <= 0, V_full <- as.matrix((nearPD(V_full, corr = FALSE))$mat), V_full <- V_full)
  
  S_Full<-.get_S_Full(n_phenotypes,S_LD,varSNP,beta_SNP,TWAS,1)
  
  ##smooth to near positive definite if either V or S are non-positive definite
  ks <- nrow(S_Full)
  smooth1 <- ifelse(eigen(S_Full)$values[ks] <= 0, S_Full <- as.matrix((nearPD(S_Full, corr = FALSE))$mat), S_Full <- S_Full)
  
  k2 <- ncol(S_Full)
  
  ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  for (i in 1) {
    
    W <- solve(V_full,tol=toler)
    
    test2 <- .tryCatch.W.E(ReorderModel <- sem(model, sample.cov = S_Full, estimator = "DWLS",
                                               WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf, optim.force.converged=TRUE
                                               ,control=list(iter.max=1),std.lv=std.lv))
    
    if(fix_measurement){
      #pull the model with SNP effects
      withSNP<-parTable(ReorderModel)
      
      #rbind back in Gene/SNP effects
      if(TWAS){
        for(p in 1:nrow(withSNP)){
          if(withSNP$rhs[p] == "Gene" | withSNP$lhs[p] == "Gene"){
            Model1<-rbind(Model1,withSNP[p,])
          }
        }
      }else{ 
        for(p in 1:nrow(withSNP)){
          if(withSNP$rhs[p] == "SNP" | withSNP$lhs[p] == "SNP"){
            Model1<-rbind(Model1,withSNP[p,])
          }
        }
      }
      
      #rbind back in ghost parameters if relevant
      for(p in 1:nrow(withSNP)){
        if(withSNP$op[p] == ":="){
          Model1<-rbind(Model1,withSNP[p,])
        }
      }
      
      #estimate model with SNP effects and fixed measurement model to get ordering of V
      test3 <- .tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Full, estimator = "DWLS",
                                                 WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf, optim.force.converged=TRUE,
                                                 control=list(iter.max=1),std.lv=std.lv))
    }
    
    #final ordering to use for multivariate GWAS model
    order <- .rearrange(k = k2, fit = ReorderModel, names = rownames(S_Full))
    
    suppressWarnings(df <- lavInspect(ReorderModel, "fit")["df"])
    suppressWarnings(npar <- lavInspect(ReorderModel, "fit")["npar"])
    
  }
  
  if(TWAS){
    SNPs <- SNPs[,1:3]
  } else {
    SNPs <- SNPs[,1:6]
  }
  
  f <- nrow(beta_SNP)
  # Run a single SNP to obtain base Lavaan model object
  LavModel1 <- .userGWAS_main(i=1, cores=1, n_phenotypes, n=1, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order, SNPs, beta_SNP, SE_SNP, varSNP, GC,
                              coords, smooth_check, TWAS, printwarn, toler, estimation, sub, Model1, df, npar, returnlavmodel=TRUE,Q_SNP=Q_SNP,model=model)
  if(!parallel){
    #make empty list object for model results if not saving specific model parameter
    if(sub[[1]]==FALSE){
      Results_List <- vector(mode="list",length=nrow(beta_SNP))
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
        }
      }
      final2 <- .userGWAS_main(i, cores=1, n_phenotypes, 1, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order, SNPs, beta_SNP, SE_SNP, varSNP, GC,
                               coords, smooth_check, TWAS, printwarn, toler, estimation, sub, Model1, df, npar, basemodel=LavModel1,Q_SNP=Q_SNP,model=model)
      final2$i <- NULL
      if(sub[[1]] != FALSE){
        final3 <- as.data.frame(matrix(NA,ncol=ncol(final2),nrow=length(sub)))
        final3[1:length(sub),] <- final2[1:length(sub),]
        if(i == 1){
          Results_List <- vector(mode="list",length=length(sub))
          for(y in 1:length(sub)){
            Results_List[[y]] <- as.data.frame(matrix(NA,ncol=ncol(final3),nrow=f))
            colnames(Results_List[[y]]) <- colnames(final2)
            Results_List[[y]][1,] <- final3[y,]
          }
        }else{
          for(y in 1:nrow(final3)){
            Results_List[[y]][i,] <- final3[y,]
          }
        }
      }else{
        Results_List[[i]] <- final2
      }
    }
    
    time_all <- proc.time()-time
    print(time_all[3])
    return(Results_List)
    
  } else {
    if(is.null(cores)){
      ##if no default provided use 1 less than the total number of cores available so your computer will still function
      int <- min(c(nrow(SNPs2), detectCores() - 1))
    }else{
      if (cores > nrow(SNPs))
        warning(paste0("Provided number of cores was greater than number of SNPs, reverting to cores=",nrow(SNPs2)))
      int <- min(c(cores, nrow(SNPs)))
    }
    if(MPI){
      #register MPI
      cl <- getMPIcluster()
      #register cluster; no makecluster as ibrun already starts the MPI process.
      registerDoParallel(cl)
    } else {
      ##specify the cores should have access to the local environment
      if (Operating != "Windows") {
        cl <- makeCluster(int, type="FORK")
      } else {
        cl <- makeCluster(int, type="PSOCK")
      }
      registerDoParallel(cl)
      on.exit(stopCluster(cl))
    }
    
    #split the V_SNP and S_SNP matrices into as many (cores - 1) as are aviailable on the local computer
    SNPs <- suppressWarnings(split(SNPs,1:int))
    beta_SNP <- suppressWarnings(split(beta_SNP,1:int))
    SE_SNP <- suppressWarnings(split(SE_SNP,1:int))
    varSNP <- suppressWarnings(split(varSNP,1:int))
    
    if(TWAS){
      print("Starting TWAS Estimation")
    } else {
      print("Starting GWAS Estimation")
    }
    if (Operating != "Windows") {
      results <- foreach(n = icount(int), .combine = 'rbind') %:%
        foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = "lavaan") %dopar% 
        .userGWAS_main(i, int, n_phenotypes, n, I_LD, V_LD, S_LD,
                       std.lv, varSNPSE2, order, SNPs[[n]], beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]], GC,
                       coords, smooth_check, TWAS, printwarn, toler, estimation, sub, Model1, df, npar, basemodel=LavModel1,Q_SNP=Q_SNP,model=model)
    } else {
      #Util-functions have to be explicitly passed to the analysis function in PSOCK cluster
      utilfuncs <- list()
      utilfuncs[[".tryCatch.W.E"]] <- .tryCatch.W.E
      utilfuncs[[".get_V_SNP"]] <- .get_V_SNP
      utilfuncs[[".get_Z_pre"]] <- .get_Z_pre
      utilfuncs[[".get_V_full"]] <- .get_V_full
      results <- foreach(n = icount(int), .combine = 'rbind') %:%
        foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = c("lavaan", "gdata"),
                 .export=c(".userGWAS_main")) %dopar% {
                   .userGWAS_main(i, int, n_phenotypes, n, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order,
                                  SNPs[[n]], beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]], GC, coords,
                                  smooth_check, TWAS, printwarn, toler, estimation, sub, Model1,
                                  df, npar, utilfuncs, basemodel=LavModel1,Q_SNP=Q_SNP,model=model)
                 }
    }
    
    ##sort results so it is in order of the output lists provided for the function
    results <-  results[order(results$i),]
    results$i <- NULL
    if(sub[[1]] != FALSE){
      Results_List <- vector(mode="list", length=length(sub))
      for(y in 1:length(sub)){
        Results_List[[y]] <- as.data.frame(matrix(NA,ncol=ncol(results),nrow=nrow(results)/length(sub)))
        colnames(Results_List[[y]]) <- colnames(results)
        Results_List[[y]] <- subset(results, paste0(results$lhs, results$op, results$rhs, sep = "") %in% sub[[y]] | is.na(results$lhs))
      }
      rm(results)
    }
    
    if(sub[[1]] == FALSE){
      if(TWAS){
        names <- unique(results$Panel)
        Results_List <- vector(mode="list", length=length(names))
        for(y in 1:length(names)){
          Results_List[[y]] <- subset(results, results$Panel == names[[y]])
          Results_List[[y]]$Model_Number <- NULL
        }
      }else{
        names <- unique(results$SNP)
        Results_List <- vector(mode="list", length=length(names))
        for(y in 1:length(names)){
          Results_List[[y]] <- subset(results, results$SNP == names[[y]])
          Results_List[[y]]$Model_Number <- NULL
        }
      }
      rm(results)
      rm(names)
    }
    
    time_all <- proc.time()-time
    print(time_all[3])
    return(Results_List)
    
  }
}
