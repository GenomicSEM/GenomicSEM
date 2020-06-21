commonfactorGWAS <- function(covstruc=NULL,SNPs=NULL,estimation="DWLS",cores=NULL,toler=.Machine$double.eps,SNPSE=NA,parallel=TRUE,output=NULL,GC="standard",MPI=FALSE,nearPD.toler=1.E-7,logfile="") {

  start.time<-proc.time()

  cat("Please note that as of 11/21/2019 commonfactorGWAS was updated to combine addSNPs and commonfactorGWAS.\n\n")

  if (class(SNPs) == "character") {
    cat(
      "You are likely listing arguments (e.g. output = ...) in the order of a previous version of commonfactorGWAS.",
      "The output of addSNPs() may still be used provided the argument is explicitely spelled out (i.e. output = ...)",
      "The current version of the function is faster and saves memory. It expects the output from ldsc()",
      "(i.e. covstruc = ...) followed by the output from sumstats() (i.e., SNPs = ...) as first two arguments.",
      "See ?commonfactorGWAS for help on proper usage.", sep='\n'
    )
    warning( paste(
      "You are likely listing arguments (e.g. output = ...) in the order of a previous version of commonfactorGWAS.",
      "The output of addSNPs() may still be used provided the argument is explicitely spelled out (i.e. output = ...)",
      "The current version of the function is faster and saves memory. It expects the output from ldsc()",
      "(i.e. covstruc = ...) followed by the output from sumstats() (i.e., SNPs = ...) as first two arguments.",
      "See ?commonfactorGWAS for help on proper usage.", sep='\n'
    ) )
  } else if (is.null(SNPs) || is.null(covstruc)) {
    cat(
      "You may be listing arguments (e.g. output = ...) in the order of a previous version of commonfactorGWAS.",
      "The output of addSNPs() may still be used provided the argument is explicitely spelled out (i.e. output = ...)",
      "The current version of the function is faster and saves memory. It expects the output from ldsc()",
      "(i.e. covstruc = ...) followed by the output from sumstats() (i.e., SNPs = ...) as first two arguments.",
      "See ?commonfactorGWAS for help on proper usage.", sep='\n'
    )
    warning( paste(
      "You may be listing arguments (e.g. output = ...) in the order of a previous version of commonfactorGWAS.",
      "The output of addSNPs() may still be used provided the argument is explicitely spelled out (i.e. output = ...)",
      "The current version of the function is faster and saves memory. It expects the output from ldsc()",
      "(i.e. covstruc = ...) followed by the output from sumstats() (i.e., SNPs = ...) as first two arguments.",
      "See ?commonfactorGWAS for help on proper usage.", sep='\n'
    ) )
  }

  ## detect operating system
  Operating<-Sys.info()[['sysname']]

  foreach.eh = 'remove'

  ## detect number of cores for parallel execution.
  ## if no number is provided use 1 less than the total number of cores available
  if ( is.null(cores) ) {
    int <- detectCores() - 1
  } else { int <- cores }

  .cls=NULL
  if ( parallel ) {
    ## register parallel backend
    if ( Operating == "Windows" ) {
      .cls = makeCluster(int, outfile=logfile)
      warning( paste(
        "Parallel processing is under development for Windows operating systems. Should you run into any issues,",
        "please, set the parallel argument to FALSE, or switch to a Linux or Mac operating system."
      ) )
    } else {
      if ( MPI ) {
        .cls <- getMPIcluster()  # no makecluster as ibrun already starts the MPI process
      } else {
        .cls <- makeCluster( int, type='FORK', outfile=logfile )  # makecluster with fork
      }
    }
    ## register cluster
    registerDoParallel(.cls)
    `%do.operator%` = `%dopar%`
    cat( 'processing variants in parallel..\n' )
  } else {
    `%do.operator%` = `%do%`
    cat( 'processing variants sequentially..\n' )
  }

  ##modification of trycatch that allows the results of a failed run to still be saved
  tryCatch.W.E <- function(expr)
  {
    .warning <- NULL
    w.handler <- function(w) { # warning handler
      .warning <<- c( .warning, w$message )
      invokeRestart("muffleWarning")
    }
    list(
      value=withCallingHandlers(
        tryCatch( expr, error=function(e) e ),
        warning=w.handler
      ),
      warning=paste( .warning, collapse='\n' )
    )
  }

  loadNamespace( 'lavaan' )
  results.colnames = c(names(SNPs)[1:6],"lhs","op","rhs","est","se","se_c","Q","fail","warning")

  if ( is.null(output) ) {

    ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
    SNPs<-data.frame(SNPs)
    SNPs$A1<-as.character(SNPs$A1)
    SNPs$A2<-as.character(SNPs$A2)
    SNPs$SNP<-as.character(SNPs$SNP)

    #SNP variance
    SNPvar=2*SNPs$MAF*(1-SNPs$MAF)

    #small number because treating MAF as fixed
    SNPvar.se2 = ifelse( is.na(SNPSE), (.0005)^2, SNPSE^2 )

    V_LD<-as.matrix(covstruc[[1]])
    S_LD<-as.matrix(covstruc[[2]])
    I_LD<-as.matrix(covstruc[[3]])

    #set the number of observations (LD-score regression blocks) behind the covariance matrix
    S_LD_nobs<-covstruc$n  #[typically n=200] this causes long fitting time in 'DWLS' -> __NOBS_OVERRIDE__

    check_names<-str_detect(colnames(S_LD), "-")
    if (any(check_names==TRUE)) { warning(
      "The trait names you specified when running the ldsc function include mathematical symbols (e.g., + or -) that",
      "will be misread by lavaan. Please rename them to avoid unexpected behaviours."
    ) }

    SNPbeta<-SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))]
    SNPse<-SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))]

    #enter in k for number of phenotypes
    k<-ncol(SNPbeta)

    #set univariate intercepts to 1 if estimated below 1
    diag(I_LD)<-ifelse(diag(I_LD)<=1, 1, diag(I_LD))

    #f = number of SNPs in dataset
    f=nrow(SNPbeta)

    #function to rearrange the sampling covariance matrix from original order to lavaan's order:
    # 'k' is the number of variables in the model
    # 'fit' is the fit function of the regression model
    # 'names' is a vector of variable names in the order you used
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

    ##pull the column names specified in the munge function
    traits<-colnames(S_LD)

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
    coords<-which(I_LD != 'NA', arr.ind= T)

    cat( "calibrating lavaan's factor structure.. " )
    ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
    for (i in 1) {
      #create empty shell of V_SNP matrix
      V_SNP<-diag(k)
      
      #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
      for (p in 1:nrow(coords)) { 
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) { 
          V_SNP[x,y]<-(SNPse[i,y]*SNPse[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*SNPvar[i]^2)}
        if (x == y) {
          V_SNP[x,x]<-(SNPse[i,x]*I_LD[x,x]*SNPvar[i])^2
        }
      }

      ##create shell of full sampling covariance matrix
      V_Full<-diag(((k+1)*(k+2))/2)

      ##input the ld-score regression region of sampling covariance from ld-score regression SEs
      V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD

      ##add in SE of SNP variance as first observation in sampling covariance matrix
      V_Full[1,1]<-SNPvar.se2

      ##add in SNP region of sampling covariance matrix
      V_Full[2:(k+1),2:(k+1)]<-V_SNP

      k2<-nrow(V_Full)
      if ( eigen(V_Full)$values[k2] <= 0 )
        V_Full <- as.matrix( nearPD(V_Full, corr = FALSE)$mat )

      W <- solve(V_Full,tol=toler)

      #create empty vector for S_SNP
      S_SNP<-vector(mode="numeric",length=k+1)

      #enter SNP variance from reference panel as first observation
      S_SNP[1]<-SNPvar[i]

      #enter SNP covariances (standardized beta * SNP variance from refference panel)
      for (p in 1:k) {
        S_SNP[p+1]<-SNPvar[i]*SNPbeta[i,p]
      }

      #create shell of the full S (observed covariance) matrix
      S_Fullrun<-diag(k+1)

      ##add the LD portion of the S matrix
      S_Fullrun[(2:(k+1)),(2:(k+1))]<-S_LD

      ##add in observed SNP variances as first row/column
      S_Fullrun[1:(k+1),1]<-S_SNP
      S_Fullrun[1,1:(k+1)]<-t(S_SNP)

      ##pull in variables names specified in LDSC function and name first column as SNP
      colnames(S_Fullrun)<-c("SNP", colnames(S_LD))

      ##name rows like columns
      rownames(S_Fullrun)<-colnames(S_Fullrun)

      ##smooth to near positive definite if either V or S are non-positive definite
      ks<-nrow(S_Fullrun)
      if ( eigen(S_Fullrun)$values[ks] <= 0 )
        S_Fullrun <- as.matrix( nearPD(S_Fullrun, corr = FALSE)$mat )

      suppress <- tryCatch.W.E(
        ReorderModel <- sem(
          Model1,
          sample.cov = S_Fullrun,
          estimator = "DWLS",
          WLS.V = W,
          sample.nobs = 2,
          optim.dx.tol = +Inf,
          optim.force.converged = TRUE,
          control = list(iter.max=1)
        )
      )

      order <- rearrange(k=k+1, fit=ReorderModel, names=rownames(S_Fullrun))
    }
    cat( "done.\n" )

    workChunks = rep_len( 1:int, dim(SNPs)[1] )
    foreach.opts = list( chunkSize=as.integer( dim(SNPs)[1]/int ) )

    ## foreach (parallel) processing that rbinds results across cores
    cat( "running lavaan common factor models including individual variants..\n" )
    results <- foreach (
      SNPs.it=iterators::isplit(SNPs[,1:6],workChunks),
      SNPbeta.it=iterators::isplit(SNPbeta,workChunks),
      SNPse.it=iterators::isplit(SNPse,workChunks),
      SNPvar.it=iterators::isplit(SNPvar,workChunks),
      .combine=rbind, .errorhandling=foreach.eh, .inorder=F,
      .noexport=c('SNPs','SNPbeta','SNPse','SNPvar'),
      .options.nws=foreach.opts, .packages='lavaan'
      ) %:% foreach (
      SNPs.itw=iterators::iter(SNPs.it$value, by='row'),
      SNPbeta.itw=iterators::iter(SNPbeta.it$value, by='row'),
      SNPse.itw=iterators::iter(SNPse.it$value, by='row'),
      SNPvar.itw=iterators::iter(SNPvar.it$value, by='row'),
      .combine=rbind, .multicombine=T, .errorhandling=foreach.eh, .inorder=F,
      .noexport=c('SNPs','SNPbeta','SNPse','SNPvar'), .packages='lavaan'
    ) %do.operator% {
    
      #create empty shell of V_SNP matrix
      V_SNP<-diag(k)
      
      cat( 'processing', SNPs.itw[,1], '\n' )
      #double GC correction using univariate LDSC intercepts
      if(GC == "conserv"){
      #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
      for (p in 1:nrow(coords)) { 
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) { 
          V_SNP[x,y]<-(SNPse.itw[,y]*SNPse.itw[,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*SNPvar.itw^2)
        }
        if (x == y) {
          V_SNP[x,x]<-(SNPse.itw[,x]*I_LD[x,x]*SNPvar.itw)^2
        }
      }
      }
      
      #single GC correction using sqrt of univariate LDSC intercepts
      if(GC == "standard"){
      #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SNPse.itw[,y]*SNPse.itw[,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*SNPvar.itw^2)
          }
          if (x == y) {
            V_SNP[x,x]<-(SNPse.itw[,x]*sqrt(I_LD[x,x])*SNPvar.itw)^2
          }
        }
      }
      
      #no GC correction
      if(GC == "none"){
      #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SNPse.itw[,y]*SNPse.itw[,x]*I_LD[x,y]*SNPvar.itw^2)
          }
          if (x == y) {
            V_SNP[x,x]<-(SNPse.itw[,x]*SNPvar.itw)^2
          }
        }
      }

      ##create shell of full sampling covariance matrix
      V_Full<-diag(((k+1)*(k+2))/2)

      ##input the ld-score regression region of sampling covariance from ld-score regression SEs
      V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD

      ##add in SE of SNP variance as first observation in sampling covariance matrix
      V_Full[1,1]<-SNPvar.se2

      ##add in SNP region of sampling covariance matrix
      V_Full[2:(k+1),2:(k+1)]<-V_SNP

      kv<-nrow(V_Full)
      if ( eigen(V_Full)$values[kv] <= 0 ) {
        cat( 'smoothing non positive definite matrix\n' )
        V_Full <- as.matrix( nearPD(V_Full, corr=FALSE)$mat )
      }

      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Full_Reorder <- V_Full[order,order]

      #create empty vector for S_SNP
      S_SNP<-vector(mode="numeric",length=k+1)

      #enter SNP variance from reference panel as first observation
      S_SNP[1]<-SNPvar.itw

      #enter SNP covariances (standardized beta * SNP variance from refference panel)
      for (p in 1:k) {
        S_SNP[p+1]<-SNPvar.itw*SNPbeta.itw[,p]
      }

      #create shell of the full S (observed covariance) matrix
      S_Fullrun<-diag(k+1)

      ##add the LD portion of the S matrix
      S_Fullrun[(2:(k+1)),(2:(k+1))]<-S_LD

      ##add in observed SNP variances as first row/column
      S_Fullrun[1:(k+1),1]<-S_SNP
      S_Fullrun[1,1:(k+1)]<-t(S_SNP)

      #name the columns
      colnames(S_Fullrun)<-c("SNP", colnames(S_LD))

      ##name rows like columns
      rownames(S_Fullrun)<-colnames(S_Fullrun)

      ##smooth to near positive definite if either V or S are non-positive definite
      ks<-nrow(S_Fullrun)
      if ( eigen(S_Fullrun)$values[ks] <= 0 )
          S_Fullrun <- as.matrix( nearPD(S_Fullrun, corr = FALSE)$mat )

      ##set estimation-specific parameters
      S_LD_nobs = 200  # see __NOBS_OVERRIDE__
      W = NULL
      if ( estimation == "DWLS" ) {
        S_LD_nobs = 2  # see __NOBS_OVERRIDE__
        u <- nrow(V_Full_Reorder)
        V_Full_Reorderb <- diag(u)
        diag(V_Full_Reorderb) <- diag(V_Full_Reorder)
        ##invert the reordered sampling covariance matrix to create a weight matrix
        W <- solve(V_Full_Reorderb,tol=toler)
      }

      # initialize output values
      se<-NA
      se_c<-NA
      Q<-NA

      ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error.
      test <- tryCatch.W.E(
        Model1_Results <- sem(
          Model1,
          sample.cov = S_Fullrun,
          estimator = estimation,
          WLS.V = W,
          sample.nobs = S_LD_nobs,
          optim.dx.tol = +Inf
        )
      )

      if ( test$warning == '' ) test$warning = 0

      if ( class(test$value)[1] == "lavaan" && !any( grepl("solution has NOT",  as.character(test$warning)) ) ) {

        #pull the delta matrix (this doesn't depend on N)
        S2.delt <- lavInspect(Model1_Results, "delta")

        ##weight matrix from stage 2
        S2.W <- lavInspect(Model1_Results, "WLS.V")

        #the "bread" part of the sandwich is the naive parameter estimates covariance matrix that would only be correct
        #if the fit function were correctly specified
        ### NdF: toler is now by default set to .Machine$double.eps; in the original implementation the default FALSE
        ###      value was read as zero, possibly making this calculation too tolerant to singularities
        bread <- solve(t(S2.delt) %*% S2.W %*% S2.delt, tol=toler)

        #create the "lettuce" part of the sandwich
        lettuce <- S2.W %*% S2.delt

        #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
        Ohtt <- bread %*% t(lettuce) %*% V_Full_Reorder %*% lettuce %*% bread

        #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
        SE <- as.matrix(sqrt(diag(Ohtt)))

        ##pull the model se estimate
        se <- inspect(Model1_Results,"list")[k+1,15]

        ##pull the corrected SE for SNP effect on P-factor
        se_c <- SE[k,1]

        ##code to estimate Q_SNP##
        #First pull the estimates from Step 1
        ModelQ <- parTable(Model1_Results)

        #fix the indicator loadings from Step 1, free the direct effects of the SNP on the indicators, and fix the factor residual variance
        ModelQ$free <- c(rep(0, k+1), 1:(k*2), 0, 0)

        ##added##
        ModelQ$ustart <- ModelQ$est
        SNPresid <- resid(Model1_Results)$cov[k+1,1:k]

        for (t in 1:nrow(ModelQ)) {
          if (ModelQ$free[t] > 0 && ModelQ$free[t] <= k) {
            ModelQ$ustart[t]<-SNPresid[ModelQ$free[t]]} else{}
          }

        #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
        testQ <- tryCatch.W.E(
          ModelQ_Results <- sem(
            model = ModelQ,
            sample.cov = S_Fullrun,
            estimator = estimation,
            WLS.V = W,
            sample.nobs = S_LD_nobs,
            start = ModelQ$ustart,
            optim.dx.tol = +Inf
          )
        )

        if ( testQ$warning == '' ) testQ$warning = "Safe"

        if ( any( grepl("lavaan WARNING: model has NOT converged!", as.character(testQ$warning)) ) ) {

          Qvec = ModelQ$free > k && ModelQ$free <= k*2
          ModelQ$ustart[ Qvec ] <- ModelQ$ustart[ Qvec ] - .01

          testQ2 <- tryCatch.W.E(
            ModelQ_Results <- sem(
              model = ModelQ,
              sample.cov = S_Fullrun,
              estimator = estimation,
              WLS.V = W,
              sample.nobs = S_LD_nobs,
              start = ModelQ$ustart,
              optim.dx.tol = +Inf
            )
          )
        } else { testQ2<-testQ }

        if ( testQ2$warning == '' ) testQ2$warning = "Safe"

        if ( any( grepl("lavaan WARNING: model has NOT converged!", as.character(testQ2$warning)) ) ) {

          Qvec = ModelQ$free > k
          ModelQ$ustart[ Qvec ] = .05

          testQ3 <- tryCatch.W.E(
            ModelQ_Results <- sem(
              model = ModelQ,
              sample.cov = S_Fullrun,
              estimator = estimation,
              WLS.V = W,
              sample.nobs = S_LD_nobs,
              start = ModelQ$ustart,
              optim.dx.tol = +Inf
            )
          )
        } else { testQ3<-testQ2 }

        if ( testQ3$warning == '' ) testQ3$warning = "Safe"
        
        if ( !any( grepl("lavaan WARNING: model has NOT converged!", as.character(testQ3$warning)) ) ) {

          #pull the delta matrix for Q (this doesn't depend on N)
          S2.delt_Q <- lavInspect(ModelQ_Results, "delta")

          ##weight matrix from stage 2 for Q
          S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V")

          #the "bread" part of the sandwich is the naive parameter estimates covariance matrix that would only be correct
          #if the fit function were correctly specified
          ### NdF: toler is now by default set to .Machine$double.eps; in the original implementation the default FALSE
          ###      value was read as zero, possibly making this calculation too tolerant to singularities
          bread_Q <- solve(t(S2.delt_Q) %*% S2.W_Q %*% S2.delt_Q, tol=toler)

          #create the "lettuce" part of the sandwich
          lettuce_Q <- S2.W_Q %*% S2.delt_Q

          #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
          Ohtt_Q <- bread_Q %*% t(lettuce_Q) %*% V_Full_Reorder %*% lettuce_Q %*% bread_Q

          ##compute diagonal matrix (Ron calls this lambda, we call it Eig) of eigenvalues of the sampling covariance matrix of the model residuals (V_eta)
          V_eta<-Ohtt_Q[1:k,1:k]
          Eig2<-as.matrix(eigen(V_eta)$values)
          Eig<-diag(k)
          diag(Eig)<-Eig2

          #Pull P1 (the eigen vectors of V_eta)
          P1<-eigen(V_eta)$vectors

          ##Pull eta = vector of direct effects of the SNP (Model Residuals)
          eta<-cbind(inspect(ModelQ_Results,"list")[(k+2):(2*k+1),14])

          #Ronald's magic combining all the pieces from above:
          Q <- t(eta) %*% P1 %*% solve(Eig,tol=toler) %*% t(P1) %*% eta

        } else { Q<-"Not Computed" }

      }

      ##pull all the results into a single row
      ##put the corrected standard error and Q in one line
      df <- cbind(
        SNPs.itw,
        inspect(Model1_Results,"list")[k+1,-c(1,5:13,15)],
        se,
        se_c,
        Q,
        ifelse( class(test$value)[1] == "lavaan", 0, as.character(test$value) ),
        as.character(test$warning),
        stringsAsFactors=FALSE
      )

      return(df)

    }

    cat( "done.\n" )

  }

  if (!is.null(output)) {

    #enter in k for number of phenotypes
    k<-ncol(output[[2]][[1]])-1

    check_names<-str_detect(colnames(output[[2]][[1]]), "-")
    if (any(check_names==TRUE)) { warning(
      "The trait names you specified when running the ldsc function include mathematical symbols (e.g., + or -) that",
      "will be misread by lavaan. Please rename them to avoid unexpected behaviours."
    ) }

    #function to rearrange the sampling covariance matrix from original order to lavaan's order:
    # 'k' is the number of variables in the model
    # 'fit' is the fit function of the regression model
    # 'names' is a vector of variable names in the order you used
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

    ##pull the column names specified in the munge function
    traits<-colnames(output[[2]][[1]])[2:(k+1)]

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

    ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
    i = 1
    calibrate.lavaan = 0
    calibrate.lavaan.test = NULL
    while ( calibrate.lavaan == 0 || class( calibrate.lavaan.test$value )[1] != 'lavaan' ) {

      calibrate.lavaan = 1

      i = i+1
      if ( i > length( output[[2]] ) ) {
        warning( 'could not calibrate lavaan.\n' )
        break()
      }

      S_Fullrun<-output[[2]][[i]]

      #transform sampling covariance matrix into a weight matrix:
      W <- solve(output[[1]][[i]],tol=toler)

      calibrate.lavaan.test <- tryCatch.W.E(
        ReorderModel <- sem(
          Model1,
          sample.cov = S_Fullrun,
          estimator = "DWLS",
          WLS.V = W,
          sample.nobs = 2,
          optim.dx.tol = +Inf,
          optim.force.converged = TRUE,
          control = list(iter.max=1)
        )
      )
      order <- rearrange(k=k+1, fit=ReorderModel, names=rownames(output[[2]][[i]]))

    }

    ##set the number of observations (LD-score regression blocks) behind the covariance matrix
    S_LD_nobs<-covstruc$n  #[typically n=200] this causes long fitting time in 'DWLS' -> __NOBS_OVERRIDE__

    workChunks = rep_len( 1:int, dim(SNPs)[1] )

    ## foreach (parallel) processing that rbinds results across cores
    cat( "running lavaan common factor models including individual variants..\n" )
    results <- foreach (
      V_Full.it=iterators::isplit(output[[1]],workChunks),
      S_Full.it=iterators::isplit(output[[2]],workChunks),
      SNPs.it=iterators::isplit(output[[3]],workChunks),
      .combine=rbind, .errorhandling=foreach.eh,
      .noexport='output', .packages='lavaan'
      ) %:% foreach (
      V_Full.itw=iterators::iter(V_Full.it$value, by='row'),
      S_Full.itw=iterators::iter(S_Full.it$value, by='row'),
      SNPs.itw=iterators::iter(SNPs.it$value, by='row'),
      .combine=rbind, .multicombine=T, .errorhandling=foreach.eh,
      .noexport='output', .packages='lavaan'
    ) %do.operator% {

      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Full_Reorder <- V_Full.itw[order,order]

      ##set estimation-specific parameters
      S_LD_nobs = 200  # see __NOBS_OVERRIDE__
      W = NULL
      if ( estimation == "DWLS" ) {
        S_LD_nobs = 2  # see __NOBS_OVERRIDE__
        u<-nrow(V_Full_Reorder)
        V_Full_Reorderb<-diag(u)
        diag(V_Full_Reorderb)<-diag(V_Full_Reorder)
        ##invert the reordered sampling covariance matrix to create a weight matrix
        W <- solve(V_Full_Reorderb,tol=toler)
      }

      #import the S_Full matrix for appropriate run
      S_Fullrun<-S_Full.itw

      # initialize output values
      se<-NA
      se_c<-NA
      Q<-NA

      ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error.
      test <- tryCatch.W.E(
        Model1_Results <- sem(
          Model1,
          sample.cov = S_Fullrun,
          estimator = estimation,
          WLS.V = W,
          sample.nobs = S_LD_nobs,
          optim.dx.tol = +Inf
        )
      )

      if ( test$warning == '' ) test$warning = 0

      if ( class(test$value)[1] == "lavaan" && !any( grepl("solution has NOT",  as.character(test$warning)) ) ) {

        #pull the delta matrix (this doesn't depend on N)
        S2.delt <- lavInspect(Model1_Results, "delta")

        ##weight matrix from stage 2
        S2.W <- lavInspect(Model1_Results, "WLS.V")

        #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
        bread <- solve(t(S2.delt) %*% S2.W %*% S2.delt, tol=toler)

        #create the "lettuce" part of the sandwich
        lettuce <- S2.W %*% S2.delt

        #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
        Ohtt <- bread %*% t(lettuce) %*% V_Full_Reorder %*% lettuce %*% bread

        #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
        SE <- as.matrix(sqrt(diag(Ohtt)))

        ##pull the model se estimate
        se <- inspect(Model1_Results,"list")[k+1,15]

        ##pull the corrected SE for SNP effect on P-factor
        se_c <- SE[k,1]

        ##code to estimate Q_SNP##
        #First pull the estimates from Step 1
        ModelQ <- parTable(Model1_Results)

        #fix the indicator loadings from Step 1, free the direct effects of the SNP on the indicators, and fix the factor residual variance
        ModelQ$free <- c(rep(0, k+1), 1:(k*2), 0, 0)

        ##added##
        ModelQ$ustart <- ModelQ$est
        SNPresid <- resid(Model1_Results)$cov[k+1,1:k]

        for (t in 1:nrow(ModelQ)) {
          if (ModelQ$free[t] > 0 && ModelQ$free[t] <= k) {
            ModelQ$ustart[t]<-SNPresid[ModelQ$free[t]]} else{}
          }

        #run the updated common and independent pathways model with fixed indicator loadings and free direct effects. these direct effects are the model residuals
        testQ <- tryCatch.W.E(
          ModelQ_Results <- sem(
            model = ModelQ,
            sample.cov = S_Fullrun,
            estimator = estimation,
            WLS.V = W,
            sample.nobs = S_LD_nobs,
            start = ModelQ$ustart,
            optim.dx.tol = +Inf
          )
        )

        if ( testQ$warning == '' ) testQ$warning = "Safe"

        if ( any( grepl("lavaan WARNING: model has NOT converged!", as.character(testQ$warning)) ) ) {

          for(n in 1:nrow(ModelQ)) {
            if(ModelQ$free[n] > k & ModelQ$free[n] <= k*2){
              ModelQ$ustart[n]<-ModelQ$ustart[n]-.01} else{}}

          testQ2 <- tryCatch.W.E(
            ModelQ_Results <- sem(
              model = ModelQ,
              sample.cov = S_Fullrun,
              estimator = estimation,
              WLS.V = W,
              sample.nobs = S_LD_nobs,
              start = ModelQ$ustart,
              optim.dx.tol = +Inf
            )
          )
        } else { testQ2<-testQ }

        if ( testQ2$warning == '' )  testQ2$warning = "Safe"

        if ( any( grepl("lavaan WARNING: model has NOT converged!", as.character(testQ2$warning)) ) ) {

          ModelQ$ustart <- ifelse( ModelQ$free > k, .05, ModelQ$ustart )

          testQ3 <- tryCatch.W.E(
            ModelQ_Results <- sem(
              model = ModelQ,
              sample.cov = S_Fullrun,
              estimator = estimation,
              WLS.V = W,
              sample.nobs = S_LD_nobs,
              start = ModelQ$ustart,
              optim.dx.tol = +Inf
            )
          )
        } else { testQ3<-testQ2 }

        if ( testQ3$warning == '' ) testQ3$warning = "Safe"

        if ( !any( grepl("lavaan WARNING: model has NOT converged!", as.character(testQ3$warning)) ) ) {

          #pull the delta matrix for Q (this doesn't depend on N)
          S2.delt_Q <- lavInspect(ModelQ_Results, "delta")

          ##weight matrix from stage 2 for Q
          S2.W_Q <- lavInspect(ModelQ_Results, "WLS.V")

          #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
          bread_Q <- solve(t(S2.delt_Q) %*% S2.W_Q %*% S2.delt_Q, tol=toler)

          #create the "lettuce" part of the sandwich
          lettuce_Q <- S2.W_Q %*% S2.delt_Q

          #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
          Ohtt_Q <- bread_Q %*% t(lettuce_Q) %*% V_Full_Reorder %*% lettuce_Q %*% bread_Q

          ##compute diagonal matrix (Ron calls this lambda, we call it Eig) of eigenvalues of the sampling covariance matrix of the model residuals (V_eta)
          V_eta<-Ohtt_Q[1:k,1:k]
          Eig2<-as.matrix(eigen(V_eta)$values)
          Eig<-diag(k)
          diag(Eig)<-Eig2

          #Pull P1 (the eigen vectors of V_eta)
          P1<-eigen(V_eta)$vectors

          ##Pull eta = vector of direct effects of the SNP (Model Residuals)
          eta<-cbind(inspect(ModelQ_Results,"list")[(k+2):(2*k+1),14])

          #Ronald's magic combining all the pieces from above:
          Q <- t(eta) %*% P1 %*% solve(Eig,tol=toler) %*% t(P1) %*% eta

        } else { Q<-"Not Computed" }

      }

      ##pull all the results into a single row
      ##put the corrected standard error and Q in one line
      df <- cbind(
        SNPs.itw,
        inspect(Model1_Results,"list")[k+1,-c(1,5:13,15)],
        se,
        se_c,
        Q,
        ifelse( class(test$value)[1] == "lavaan", 0, as.character(test$value) ),
        as.character(test$warning),
        stringsAsFactors=FALSE
      )

      return(df)

    }

    cat( "done.\n" )

  }

  if ( !is.null( results ) ) {
    results = data.frame( results )
    colnames(results) = results.colnames
    results$Z_Estimate <- results$est/results$se_c
    results$Pval_Estimate <- 2*pnorm( abs(results$Z_Estimate), lower.tail=FALSE )
    results$Q_df <- k-1
    results$Q_pval <- ifelse( is.finite( results$Q ), pchisq( results$Q, results$Q_df, lower.tail=FALSE ), NA )
#   [1]  "SNP"           "CHR"           "BP"            "MAF"          
#   [5]  "A1"            "A2"            "lhs"           "op"           
#   [9]  "rhs"           "est"           "se"            "se_c"         
#   [13] "Q"             "fail"          "warning"       "Z_Estimate"   
#   [17] "Pval_Estimate" "Q_df"          "Q_pval"       
    results <- results[, c(1:10,12,16,17,13,18,19,14,15) ]
  }

  if ( !is.null( .cls ) ) stopCluster( .cls )
  time_all<-proc.time()-start.time
  print(time_all[3])

  return(results)

}

