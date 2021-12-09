userGWAS <- function(covstruc=NULL, SNPs=NULL, estimation="DWLS", model="", printwarn=TRUE,
                     sub=FALSE,cores=NULL, toler=FALSE, SNPSE=FALSE, parallel=TRUE, GC="standard", MPI=FALSE,
                     smooth_check=FALSE, TWAS=FALSE, std.lv=FALSE){
  # Set toler to machine precision to enable passing this to solve() directly
  if (!toler) toler <- .Machine$double.eps
  # Sanity checks
  #TODO: add check for covstruc
  #TODO: add check for SNPs
  #TODO: add check for sub
  .check_one_of(estimation, c("DWLS", "ML"))
  .check_boolean(printwarn)
  if (!is.null(cores)) .check_range(cores, min=0, max=Inf, allowNA=FALSE)
  .check_range(toler, min=0, max=Inf)
  .check_boolean(SNPSE)
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
  }

  if(SNPSE != FALSE){
    varSNPSE2 <- SNPSE^2
  }

  V_LD <- as.matrix(covstruc[[1]])
  S_LD <- as.matrix(covstruc[[2]])
  I_LD <- as.matrix(covstruc[[3]])

  beta_SNP <- SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))]
  SE_SNP <- SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))]

  #enter in k for number of phenotypes
  n_phenotypes <- ncol(beta_SNP)

  #set univariate intercepts to 1 if estimated below 1
  diag(I_LD) <- ifelse(diag(I_LD)<= 1, 1, diag(I_LD))

  Model1 <- model


  coords <- which(I_LD != 'NA', arr.ind= T)
  i <- 1  # TODO: reason for this?
  #create empty shell of V_SNP matrix
  V_SNP <- diag(n_phenotypes)

  #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
  for (p in 1:nrow(coords)) {
    x <- coords[p,1]
    y <- coords[p,2]
    if (x != y) {
      V_SNP[x,y] <- (SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i]^2)}
    if (x == y) {
      V_SNP[x,x] <- (SE_SNP[i,x]*I_LD[x,x]*varSNP[i])^2
    }
  }

  V_full <- .get_V_full(n_phenotypes, V_LD, varSNPSE2, V_SNP)

  kv <- nrow(V_full)
  smooth2 <- ifelse(eigen(V_full)$values[kv] <= 0, V_full <- as.matrix((nearPD(V_full, corr = FALSE))$mat), V_full <- V_full)

  #create empty vector for S_SNP
  S_SNP <- vector(mode="numeric",length=n_phenotypes+1)

  #enter SNP variance from reference panel as first observation
  S_SNP[1] <- varSNP[i]

  #enter SNP covariances (standardized beta * SNP variance from refference panel)
  for (p in 1:n_phenotypes) {
    S_SNP[p+1] <- varSNP[i]*beta_SNP[i,p]
  }

  #create shell of the full S (observed covariance) matrix
  S_Full <- diag(n_phenotypes+1)

  ##add the LD portion of the S matrix
  S_Full[(2:(n_phenotypes+1)),(2:(n_phenotypes+1))] <- S_LD

  ##add in observed SNP variances as first row/column
  S_Full[1:(n_phenotypes+1),1] <- S_SNP
  S_Full[1,1:(n_phenotypes+1)] <- t(S_SNP)

  ##pull in variables names specified in LDSC function and name first column as SNP
  if(TWAS){
    colnames(S_Full) <- c("Gene", colnames(S_LD))
  } else {
    colnames(S_Full) <- c("SNP", colnames(S_LD))
  }

  ##name rows like columns
  rownames(S_Full) <- colnames(S_Full)

  ##smooth to near positive definite if either V or S are non-positive definite
  ks <- nrow(S_Full)
  smooth1 <- ifelse(eigen(S_Full)$values[ks] <= 0, S_Full <- as.matrix((nearPD(S_Full, corr = FALSE))$mat), S_Full <- S_Full)

  k2 <- ncol(S_Full)

  ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  for (i in 1) {

    W <- solve(V_full,tol=toler)

    test2 <- .tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Full, estimator = "DWLS",
    WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf, optim.force.converged=TRUE
      ,control=list(iter.max=1),std.lv=std.lv))

    order <- .rearrange(k = k2, fit = ReorderModel, names = rownames(S_Full))

    suppressWarnings(df <- lavInspect(ReorderModel, "fit")["df"])
    suppressWarnings(npar <- lavInspect(ReorderModel, "fit")["npar"])

  }

  if(!(TWAS)){
    SNPs2 <- SNPs[,1:6]
  } else {
    SNPs2 <- SNPs[,1:3]
  }

  rm(SNPs)
  f <- nrow(beta_SNP)
  if(parallel==FALSE){
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
      n <- 1
      final2 <- .userGWAS_main(i, n_phenotypes, n, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order, SNPs2, beta_SNP, SE_SNP, varSNP, GC,
      coords, smooth_check, TWAS, printwarn, toler, estimation, sub, Model1, df, npar)
      if(sub[[1]] != FALSE){
        final3 <- as.data.frame(matrix(NA,ncol=ncol(final2),nrow=length(sub)))
        final3[1:length(sub),] <- final2[1,]
        if(i == 1){
          Results_List <- vector(mode="list",length=length(sub))
          for(y in 1:length(sub)){
            Results_List[[y]] <- as.data.frame(matrix(NA,ncol=ncol(final3),nrow=f))
            colnames(Results_List[[y]]) <- colnames(final3)
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

  }

  if(parallel == TRUE){

    if(is.null(cores)){
      ##if no default provided use 1 less than the total number of cores available so your computer will still function
      int <- detectCores() - 1
    }else{
      int <- cores
    }
    if(MPI){
      #register MPI
      cl <- getMPIcluster()
      #register cluster; no makecluster as ibrun already starts the MPI process.
      registerDoParallel(cl)
    } else {
      registerDoParallel(int)
      ##specify the cores should have access to the local environment
      if (Operating != "Windows") {
        cl <- makeCluster(int, type="FORK")
      } else {
        cl <- makeCluster(int, type="PSOCK")
      }
    }

    SNPs2 <- suppressWarnings(split(SNPs2,1:int))
    #split the V_SNP and S_SNP matrices into as many (cores - 1) as are aviailable on the local computer
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
      foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = "lavaan") %dopar% {
        .userGWAS_main(i, n_phenotypes, n, I_LD, V_LD, S_LD,
        std.lv, varSNPSE2, order, SNPs2[[n]], beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]], GC,
        coords, smooth_check, TWAS, printwarn, toler, estimation, sub, Model1, df, npar)
      }
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
        .userGWAS_main(i, n_phenotypes, n, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order,
        SNPs2[[n]], beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]], GC, coords,
        smooth_check, TWAS, printwarn, toler, estimation, sub, Model1,
        df, npar, utilfuncs)
      }
    }
    if (!MPI) {
      stopCluster(cl)
    }

    ##sort results so it is in order of the output lists provided for the function
    results <-  results[order(results$i, results$n),]
    results$n <- NULL

    if(sub[[1]] != FALSE){
      results$i <- NULL
      Results_List <- vector(mode="list", length=length(sub))
      for(y in 1:length(sub)){
        Results_List[[y]] <- as.data.frame(matrix(NA,ncol=ncol(results),nrow=nrow(results)/length(sub)))
        colnames(Results_List[[y]]) <- colnames(results)
        Results_List[[y]] <- subset(results, paste0(results$lhs, results$op, results$rhs, sep = "") %in% sub[[y]] | is.na(results$lhs))
      }
      rm(results)
    }
    if (TWAS) {
      if(!sub[[1]]){
        names <- unique(results$Panel)
        Results_List <- vector(mode="list", length=length(names))
        for(y in 1:length(names)){
          Results_List[[y]] <- subset(results, results$Panel == names[[y]])
          Results_List[[y]]$Model_Number <- NULL
        }
        rm(results)
        rm(names)
      }
    } else {
      if(sub[[1]]==FALSE){
        names <- unique(results$SNP)
        Results_List <- vector(mode="list", length=length(names))
        for(y in 1:length(names)){
          Results_List[[y]] <- subset(results, results$SNP == names[[y]])
          Results_List[[y]]$Model_Number <- NULL
        }
        rm(results)
        rm(names)
      }
    }

    time_all <- proc.time()-time
    print(time_all[3])
    return(Results_List)

  }
}
