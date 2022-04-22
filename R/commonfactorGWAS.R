commonfactorGWAS <-function(covstruc=NULL,SNPs=NULL,estimation="DWLS",cores=NULL,toler=FALSE,SNPSE=FALSE,parallel=TRUE,
                            GC="standard",MPI=FALSE,TWAS=FALSE,smooth_check=FALSE){
  # Set toler to machine precision to enable passing this to solve() directly
  if (!toler) toler <- .Machine$double.eps
  # Sanity checks
  #TODO: add check for covstruc
  #TODO: add check for SNPs
  .check_one_of(estimation, c("DWLS", "ML"))
  if (!is.null(cores)) .check_range(cores, min=0, max=Inf, allowNA=FALSE)
  .check_range(toler, min=0, max=Inf)
  # .check_boolean(SNPSE)
  .check_boolean(parallel)
  .check_one_of(GC, c("standard", "conserv", "none"))
  .check_boolean(MPI)
  .check_boolean(TWAS)
  .check_boolean(smooth_check)
  # Sanity checks finished

  time <- proc.time()
  Operating <- Sys.info()[['sysname']]
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

  ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
  SNPs <- data.frame(SNPs)

  if (TWAS) {
    SNPs$Gene  <- as.character(SNPs$Gene)
    SNPs$Panel <- as.character(SNPs$Panel)
    varSNP <- SNPs$HSQ
  } else {
    SNPs$A1  <- as.character(SNPs$A1)
    SNPs$A2  <- as.character(SNPs$A2)
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

  check_names <- str_detect(colnames(S_LD), "-")
  if(any(check_names)){warning("Your trait names specified when running the ldsc function include mathematical arguments (e.g., + or -) that will be misread by lavaan. Please rename the traits.")}

  beta_SNP <- SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))]
  SE_SNP <- SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))]

  #enter in k for number of phenotypes
  k <- ncol(beta_SNP)

  #print warning if number of traits are unequal across SNPs and LDSC output
  if(ncol(beta_SNP) != ncol(S_LD)){
    stop("There are different numbers of traits in the sumstats and ldsc output. Please verify that the same summary statistics have been provided to both functions before running commonfactorGWAS.")
  }

  #set univariate intercepts to 1 if estimated below 1
  diag(I_LD) <- ifelse(diag(I_LD)<= 1, 1, diag(I_LD))

  #f = number of SNPs in dataset
  f <- nrow(beta_SNP)

  ##pull the column names specified in the munge function
  traits <- colnames(S_LD)

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

  ##pull the coordinates of the I_LD matrix to loop making the V_SNP matrix
  coords <- which(I_LD != 'NA', arr.ind= T)

  ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  for (i in 1) {
    V_SNP <- .get_V_SNP(SE_SNP, I_LD, varSNP, "conserv", coords, k, i)

    ##create shell of full sampling covariance matrix
    V_Full <- .get_V_full(k, V_LD, varSNPSE2, V_SNP)

    k2 <- nrow(V_Full)
    smooth2 <- ifelse(eigen(V_Full)$values[k2] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)

    W <- solve(V_Full,tol=toler)

    #create empty vector for S_SNP
    S_SNP <- vector(mode="numeric",length=k+1)

    #enter SNP variance from reference panel as first observation
    S_SNP[1] <- varSNP[i]

    #enter SNP covariances (standardized beta * SNP variance from refference panel)
    for (p in 1:k) {
      S_SNP[p+1] <- varSNP[i]*beta_SNP[i,p]
    }

    #create shell of the full S (observed covariance) matrix
    S_Fullrun <- diag(k+1)

    ##add the LD portion of the S matrix
    S_Fullrun[(2:(k+1)),(2:(k+1))] <- S_LD

    ##add in observed SNP variances as first row/column
    S_Fullrun[1:(k+1),1] <- S_SNP
    S_Fullrun[1,1:(k+1)] <- t(S_SNP)

    ##pull in variables names specified in LDSC function and name first column as SNP
    colnames(S_Fullrun) <- c("SNP", colnames(S_LD))

    ##name rows like columns
    rownames(S_Fullrun) <- colnames(S_Fullrun)

    ##smooth to near positive definite if either V or S are non-positive definite
    ks <- nrow(S_Fullrun)
    smooth1 <- ifelse(eigen(S_Fullrun)$values[ks] <= 0, S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat), S_Fullrun<-S_Fullrun)

    suppress <- .tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf,optim.force.converged=TRUE,control=list(iter.max=1)))

    order <- .rearrange(k = k+1, fit = ReorderModel, names = rownames(S_Fullrun))
  }

  if(TWAS){
    SNPs2 <- SNPs[,1:3]
  } else {
    SNPs2 <- SNPs[,1:6]
  }

  rm(SNPs)
  ##name the columns of the results file
  if(smooth_check){
    colnamesresults <- c("i", "lhs", "op", "rhs", "est", "se", "se_c", "Q", "fail", "warning", "Z_smooth")
  } else {
    colnamesresults <- c("i", "lhs", "op", "rhs", "est", "se", "se_c", "Q", "fail", "warning")
  }
  LavModel1 <- .commonfactorGWAS_main(1, cores=1, 1, S_LD, V_LD, I_LD, beta_SNP, SE_SNP, varSNP, varSNPSE2, GC, coords, k, smooth_check,Model1, toler, estimation, order,returnlavmodel=TRUE)
  if(!parallel){
    if(smooth_check){
      results <- as.data.frame(matrix(NA, ncol=11, nrow=f))
    } else {
      results <- as.data.frame(matrix(NA, ncol=10, nrow=f))
    }
    colnames(results) <- colnamesresults
    for (i in 1:f) {
      results[i, ] <- .commonfactorGWAS_main(i,cores=1, 1, S_LD, V_LD, I_LD, beta_SNP, SE_SNP, varSNP, varSNPSE2, GC, coords, k, smooth_check,Model1, toler, estimation, order, basemodel=LavModel1)
    }
  }
  if(parallel){
    if(is.null(cores)){
      ##if no default provided use 1 less than the total number of cores available so your computer will still function
      int <- min(c(nrow(SNPs2), detectCores() - 1))
    }else{
      if (cores > nrow(SNPs2))
        warning(paste0("Provided number of cores was greater than number of SNPs, reverting to cores=",nrow(SNPs2)))
      int <- min(c(nrow(SNPs2), cores))
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

    #split the V_SNP and S_SNP matrices into n_cores
    beta_SNP <- suppressWarnings(split(beta_SNP,1:int))
    SE_SNP   <- suppressWarnings(split(SE_SNP,1:int))
    varSNP   <- suppressWarnings(split(varSNP,1:int))
    ##foreach parallel processing that rbinds results across cores
    if (Operating != "Windows") {
      results <- foreach(n = icount(int), .combine = 'rbind') %:%
        foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = "lavaan") %dopar% {
        .commonfactorGWAS_main(i,cores=int, n, S_LD, V_LD, I_LD, beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]], varSNPSE2, GC, coords, k, smooth_check,Model1, toler, estimation, order, basemodel=LavModel1)
      }
    } else {
      utilfuncs <- list()
      utilfuncs[[".tryCatch.W.E"]] <- .tryCatch.W.E
      utilfuncs[[".get_V_SNP"]] <- .get_V_SNP
      utilfuncs[[".get_V_full"]] <- .get_V_full
      results <- foreach(n = icount(int), .combine = 'rbind') %:%
        foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = c("lavaan", "gdata"),
                 .export=c(".commonfactorGWAS_main")) %dopar% {
        .commonfactorGWAS_main(i,cores=int, n, S_LD, V_LD, I_LD, beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]], varSNPSE2, GC, coords, k, smooth_check,Model1, toler, estimation, order, utilfuncs, basemodel=LavModel1)
      }
    }
    colnames(results) <- colnamesresults
    results <- results[order(results$i),]
  }
  results$se <- NULL
  results <- cbind(SNPs2,results)
  results$Z_Estimate <- results$est/results$se_c
  results$Pval_Estimate <- 2*pnorm(abs(results$Z_Estimate),lower.tail=FALSE)
  results$Q_df <- k-1
  results$Q_pval <- pchisq(results$Q,results$Q_df,lower.tail=FALSE)

  if (!TWAS & !smooth_check){
    results <- results[,c(1:12,16,17,13,18,19,14,15)]
  } else if (!TWAS & smooth_check){
    results <- results[,c(1:12,17,18,13,19,20,14,15,16)]
  } else if (TWAS & !smooth_check){
    results <- results[,c(1:9,13,14,10,15,16,11,12)]
    results$rhs <- rep("Gene",nrow(results))
  } else if (TWAS & smooth_check){
    results <- results[,c(1:9,14,15,10,16,17,11,12,13)]
    results$rhs <- rep("Gene",nrow(results))
  }
  time_all <- proc.time()-time
  print(time_all[3])
  # Fix last two column names, these are incorrectly labelled in parallel operation
  if (smooth_check) {
    colnames(results)[(length(colnames(results))-2):(length(colnames(results))-1)] <- c("fail", "warning")
  } else {
    colnames(results)[(length(colnames(results))-1):length(colnames(results))] <- c("fail", "warning")
  }
  return(results)
}
