usermodelGWAS <- function(covstruc=NULL,SNPs=NULL,estimation="DWLS",model="",modelchi=FALSE,printwarn=TRUE,sub=FALSE,cores=NULL,toler=FALSE,SNPSE=FALSE,parallel=TRUE,output=NULL,GC="standard",MPI=FALSE) {

  time<-proc.time()

  ##determine if the model is likely being listed in quotes and print warning if so
  test<-c(str_detect(model, "~"),str_detect(model, "="),str_detect(model, "\\+"))
  if(all(test) == FALSE){
    warning( paste(
      "Your model name may be listed in quotes; please remove the quotes and try re-running if the function stopped",
      "running after returning an error.",
      sep='\n'
    ) )
  }

  cat("Please note that as of 11/21/2019 usermodelGWAS was updated to combine addSNPs and usermodelGWAS.\n\n")

  if (class(SNPs) == "character") {
    cat(
      "You are likely listing arguments (e.g. output = ...) in the order of a previous version of usermodelGWAS.",
      "The output of addSNPs() may still be used provided the argument is explicitely spelled out (i.e. output = ...)",
      "The current version of the function is faster and saves memory. It expects the output from ldsc()",
      "(i.e. covstruc = ...) followed by the output from sumstats() (i.e., SNPs = ...) as first two arguments.",
      "See ?usermodelGWAS for help on proper usage.", sep='\n'
    )
    warning( paste(
      "You are likely listing arguments (e.g. output = ...) in the order of a previous version of usermodelGWAS.",
      "The output of addSNPs() may still be used provided the argument is explicitely spelled out (i.e. output = ...)",
      "The current version of the function is faster and saves memory. It expects the output from ldsc()",
      "(i.e. covstruc = ...) followed by the output from sumstats() (i.e., SNPs = ...) as first two arguments.",
      "See ?usermodelGWAS for help on proper usage.", sep='\n'
    ) )
  } else if (is.null(SNPs) || is.null(covstruc)) {
    cat(
      "You may be listing arguments (e.g. output = ...) in the order of a previous version of usermodelGWAS.",
      "The output of addSNPs() may still be used provided the argument is explicitely spelled out (i.e. output = ...)",
      "The current version of the function is faster and saves memory. It expects the output from ldsc()",
      "(i.e. covstruc = ...) followed by the output from sumstats() (i.e., SNPs = ...) as first two arguments.",
      "See ?usermodelGWAS for help on proper usage.", sep='\n'
    )
    warning( paste(
      "You may be listing arguments (e.g. output = ...) in the order of a previous version of usermodelGWAS.",
      "The output of addSNPs() may still be used provided the argument is explicitely spelled out (i.e. output = ...)",
      "The current version of the function is faster and saves memory. It expects the output from ldsc()",
      "(i.e. covstruc = ...) followed by the output from sumstats() (i.e., SNPs = ...) as first two arguments.",
      "See ?usermodelGWAS for help on proper usage.", sep='\n'
    ) )
  }

  ## detect operating system
  Operating<-Sys.info()[['sysname']]

  ## detect number of cores for parallel execution.
  ## if no number is provided use 1 less than the total number of cores available
  if ( is.null(cores) ) {
    int <- detectCores() - 1
  } else { int <- cores }

  .cls=NULL
  if ( parallel ) {
    ## register parallel backend
    if ( Operating == "Windows" ) {
      .cls = makeCluster(int)
      warning( paste(
        "Parallel processing is under development for Windows operating systems. Should you run into any issues,",
        "please, set the parallel argument to FALSE, or switch to a Linux or Mac operating system."
      ) )
    } else {
      if ( MPI ) {
        .cls <- getMPIcluster()                  # no makecluster as ibrun already starts the MPI process
      } else {
        .cls <- makeCluster( int, type='FORK' )  # makecluster with fork
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

  if(is.null(output)){

    ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
    SNPs<-data.frame(SNPs)
    SNPs$A1<-as.character(SNPs$A1)
    SNPs$A2<-as.character(SNPs$A2)
    SNPs$SNP<-as.character(SNPs$SNP)

    ##SNP variance
    SNPvar=2*SNPs$MAF*(1-SNPs$MAF)

    ##small number because treating MAF as fixed
    if(SNPSE == FALSE){
      SNPvar.se2=(.0005)^2
    }

    if(SNPSE != FALSE){
      SNPvar.se2 = SNPSE^2
    }

    V_LD<-as.matrix(covstruc[[1]])
    S_LD<-as.matrix(covstruc[[2]])
    I_LD<-as.matrix(covstruc[[3]])

    SNPbeta<-SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))]
    SNPse<-SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))]

    ##enter in k for number of phenotypes
    k<-ncol(SNPbeta)

    ##set univariate intercepts to 1 if estimated below 1
    diag(I_LD)<-ifelse(diag(I_LD)<= 1, 1, diag(I_LD))

    ##set the number of observations (typically LD-score bins) behind the covariance matrix
    S_LD_nobs<-covstruc$n  # this is what it should be (typically n-200), but causes long fitting time in 'DWLS' -> __NOBS_OVERRIDE__

    ##number of SNPs in dataset
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

    Model1<-model

    if(modelchi == TRUE){
      #function to creat row/column names for S_LD matrix
      write.names <- function(k, label = "V") {
        varnames<-vector(mode="character",length=k)

        for (i in 1:k){
          varnames[i]<-paste(label,i,sep="")}

        return(varnames)
      }

      ##pull the coordinates of the I_LD matrix to loop making the V_SNP matrix
      coords<-which(I_LD != 'NA', arr.ind= T)

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

        kv<-nrow(V_Full)
        smooth2<-ifelse(eigen(V_Full)$values[kv] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)

        #create empty vector for S_SNP
        S_SNP<-vector(mode="numeric",length=k+1)

        #enter SNP variance from reference panel as first observation
        S_SNP[1]<-SNPvar[i]

        #enter SNP covariances (standardized beta * SNP variance from refference panel)
        for (p in 1:k) {
          S_SNP[p+1]<-SNPvar[i]*SNPbeta[i,p]
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

        ##pull the column names specified in the munge function
        traits<-colnames(S_Full)

        if(toler==FALSE){
          W_test <- solve(V_Full)
        }

        if(toler!=FALSE){
          W_test <- solve(V_Full,tol=toler)
        }

        ##run the model to determine number of latent variables
        suppress <- tryCatch.W.E(
          ReorderModel1 <- sem(
            model,
            sample.cov = S_Full,
            estimator = "DWLS",
            WLS.V = W_test,
            sample.nobs = 2,
            warn = FALSE,
            optim.dx.tol = +Inf,
            optim.force.converged = TRUE,
            control = list(iter.max=1)
          )
        )

        if ( class(suppress$value)[1] != "lavaan" ) {
          i<-20
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
          smooth2<-ifelse(eigen(V_Full)$values[k2] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)

          #create empty vector for S_SNP
          S_SNP<-vector(mode="numeric",length=k+1)

          #enter SNP variance from reference panel as first observation
          S_SNP[1]<-SNPvar[i]

          #enter SNP covariances (standardized beta * SNP variance from refference panel)
          for (p in 1:k) {
            S_SNP[p+1]<-SNPvar[i]*SNPbeta[i,p]
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

          if(toler==FALSE){
            W_test <- solve(V_Full)
          }

          if(toler!=FALSE){
            W_test <- solve(V_Full,tol=toler)
          }

          ##run the model to determine number of latent variables
          suppress <- tryCatch.W.E(
            ReorderModel1 <- sem(
              model,
              sample.cov = S_Full,
              estimator = "DWLS",
              WLS.V = W_test,
              sample.nobs = 2,
              warn = FALSE,
              optim.dx.tol = +Inf,
              optim.force.converged = TRUE,
              control = list(iter.max=1)
            )
          )
        }

      }

      k2<-ncol(S_Full)

      ##create the names
      S_names<-write.names(k=k2)

      ##add bracketing so gsub knows to replace exact cases
      traits2<-traits
      for(i in 1:length(traits)){
        traits2[[i]]<-paste0("\\<", traits[[i]],"\\>",sep="")
      }

      ##replace trait names in user provided model with general form of V1-VX
      for(i in 1:length(traits)){
        model<-gsub(traits2[[i]], S_names[[i]], model)
      }

      ##determine number of latent variables from writing extended model.
      r<-nrow(lavInspect(ReorderModel1, "cor.lv"))
      lat_labs<-colnames(lavInspect(ReorderModel1, "cor.lv"))

      write.Model1 <- function(k2, label = "V", label2 = "VF") {

        ModelsatF<-""
        for (i in 1:(k2-1)) {
          linestartc <- paste(" ", label2, i, "~~0*", label2, i+1,  sep = "")
          if (k2-i >= 2) {
            linemidc <- ""
            for (j in (i+2):k2) {
              linemidc <- paste(linemidc, "+0*", label2, j, sep = "")
            }
          } else {linemidc <- ""}
          ModelsatF <- paste(ModelsatF, linestartc, linemidc, " \n ", sep = "")
        }

        if(r > 0){
          Model1b <- ""
          for (t in 1:r) {
            for (i in 1) {
              linestartb <- paste(lat_labs[t], " =~ 0*",label2, i, sep = "")
              if ((k2-1)-i > 0 | k2 == 2) {
                linemidb <- ""
                for (j in (i+1):k2) {
                  linemidb <- paste(linemidb, " + 0*", label2, j, sep = "")
                }
              } else {linemidb <- ""}

            }
            Model1b <- paste(Model1b, linestartb, linemidb, " \n ", sep = "")
          }
        }
        else {Model1b <- ""}

        if(r > 0){
          Model1c <- ""
          for (t in 1:r) {
            for (i in 1) {
              linestartc <- paste(lat_labs[t], " ~~ 0*",label2, i, sep = "")
              if ((k2-1)-i > 0 | k2 == 2) {
                linemidc <- ""
                for (j in (i+1):k2) {
                  linemidc <- paste(linemidc, " + 0*", label2, j, sep = "")
                }
              } else {linemidc <- ""}

            }
            Model1c <- paste(Model1c, linestartc, linemidc, " \n ", sep = "")
          }
        }
        else {Model1c <- ""}

        Model2<-""
        for (p in 1:k2) {
          linestart2 <- paste(label2, p, " =~ 1*", label, p, sep = "")
          Model2<-paste(Model2, linestart2, " \n ", sep = "")}

        Model3<-""
        for (p in 1:k2) {
          linestart3 <- paste(label, p, "~~", label, p, sep = "")
          Model3<-paste(Model3, linestart3, " \n ", sep = "")}

        Model4<-""
        for (p in 1:k2) {
          linestart4 <- paste(label2, p, "~~0*", label2, p, sep = "")
          Model4<-paste(Model4, linestart4, " \n ", sep = "")}

        Modelsat<-""
        for (i in 1:(k2-1)) {
          linestartc <- paste("", label, i, "~~0*", label, i+1, sep = "")
          if (k2-i >= 2) {
            linemidc <- ""
            for (j in (i+2):k2) {
              linemidc <- paste("", linemidc, label, i, "~~0*", label, j, " \n ", sep="")

            }
          } else {linemidc <- ""}
          Modelsat <- paste(Modelsat, linestartc, " \n ", linemidc, sep = "")
        }

        Model5<-paste(model, " \n ", ModelsatF, Model1b, Model1c, Model2, Model3, Model4, Modelsat, sep = "")

        return(Model5)
      }

      Model1<-write.Model1(k2)

      while ( class(u<-tryCatch.W.E(lavParseModelString(Model1))$value$message) != 'NULL' ) {
        t<-paste(strsplit(u, ": ")[[1]][3], " \n ", sep = "")
        Model1<-str_replace(Model1, fixed(t), "")
      }

      ##code to write fake model to make sure elements estimated in user model
      ##do not overlap with saturated model, in which case an alternative
      ##specification is used that uses, for example, covariances among residual factors
      write.test<-function(k2, label = "V", label2 = "VF") {

        Modelsat<-""
        for (i in 1:(k2)) {
          if (k2-i >= 1) {
            linemidc <- ""
            for (j in (i+1):k2) {
              linemidc <- paste(linemidc, label, i, "~~", label, j, " \n ", sep = "")
            }
          }else{linemidc<-""}
          Modelsat <- paste(Modelsat, linemidc, sep = "")
        }

        Modelsat2<-""
        for (i in 1:(k2)) {
          if (k2-i >= 1) {
            linemidc2 <- ""
            for (j in (i+1):k2) {
              linemidc2 <- paste(linemidc2, label, j, "~~", label, i, " \n ", sep = "")
            }
          }else{linemidc2<-""}
          Modelsat2 <- paste(Modelsat2, linemidc2, sep = "")
        }

        Model4<-""
        for (p in 1:k2) {
          linestart4 <- paste(label2, p, "~~", label2, p, sep = "")
          Model4<-paste(Model4, linestart4, " \n ", sep = "")}

        modelCFI<-paste(Modelsat, Modelsat2, Model4)
        return(modelCFI)
      }

      modeltest<-data.frame(write.test(k2))
      modeltest$write.test.k2.<-as.character(modeltest$write.test.k2.)
      modeltest2 <- cSplit(modeltest, "write.test.k2.", sep = "\n", direction = "long")
      modeltest2$write.test.k2.<-as.character(modeltest2$write.test.k2.)

      z<-(k2*(k2+1))/2
    }

    if(modelchi==FALSE){
      coords<-which(I_LD != 'NA', arr.ind= T)
      i<-1
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

      kv<-nrow(V_Full)
      smooth2<-ifelse(eigen(V_Full)$values[kv] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)

      #create empty vector for S_SNP
      S_SNP<-vector(mode="numeric",length=k+1)

      #enter SNP variance from reference panel as first observation
      S_SNP[1]<-SNPvar[i]

      #enter SNP covariances (standardized beta * SNP variance from refference panel)
      for (p in 1:k) {
        S_SNP[p+1]<-SNPvar[i]*SNPbeta[i,p]
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

      k2<-ncol(S_Full)
    }

    ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
    for (i in 1) {

      if(toler==FALSE){
        W <- solve(V_Full)
      }

      if(toler!=FALSE){
        W <- solve(V_Full,tol=toler)
      }

      if(modelchi == TRUE){
        ##name the columns and rows of the S matrix in general format V1-VX
        rownames(S_Full) <- S_names
        colnames(S_Full) <- S_names
      }

      suppress <- tryCatch.W.E(
        ReorderModel <- sem(
          Model1,
          sample.cov = S_Full,
          estimator = "DWLS",
          WLS.V = W,
          sample.nobs = 2,
          optim.dx.tol = +Inf,
          optim.force.converged = TRUE,
          control = list(iter.max=1)
        )
      )

      suppressWarnings(df<-lavInspect(ReorderModel, "fit")["df"])
      suppressWarnings(npar<-lavInspect(ReorderModel, "fit")["npar"])
      order <- rearrange(k = k2, fit = ReorderModel, names = rownames(S_Full))

    }

    workChunks = rep_len( 1:int, dim(SNPs)[1] )

    ## foreach (parallel) processing that rbinds results across cores
    cat( "running lavaan common factor models including individual variants..\n" )
    results <- foreach (
      SNPs.it=iterators::isplit(SNPs[,1:6],workChunks),
      SNPbeta.it=iterators::isplit(SNPbeta,workChunks),
      SNPse.it=iterators::isplit(SNPse,workChunks),
      SNPvar.it=iterators::isplit(SNPvar,workChunks),
      .combine=rbind, .errorhandling='remove',
      .noexport=c('SNPs','SNPbeta','SNPse','SNPvar'), .packages='lavaan'
      ) %:% foreach (
      SNPs.itw=iterators::iter(SNPs.it$value, by='row'),
      SNPbeta.itw=iterators::iter(SNPbeta.it$value, by='row'),
      SNPse.itw=iterators::iter(SNPse.it$value, by='row'),
      SNPvar.itw=iterators::iter(SNPvar.it$value, by='row'),
      .combine=rbind, .multicombine=T, .errorhandling='remove',
      .noexport=c('SNPs','SNPbeta','SNPse','SNPvar'), .packages='lavaan'
    ) %do.operator% {

        #create empty shell of V_SNP matrix
        V_SNP<-diag(k)
        
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
        smooth2<-ifelse(eigen(V_Full)$values[kv] <= 0, V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat), V_Full<-V_Full)

        #reorder sampling covariance matrix based on what lavaan expects given the specified model
        V_Full_Reorder <- V_Full[order,order]

        W = NULL
        S_LD_nobs = 200  # <- __NOBS_OVERRIDE__
        if ( estimation == 'DWLS' ) {
          S_LD_nobs = 2  # <- __NOBS_OVERRIDE__
          u<-nrow(V_Full_Reorder)
          V_Full_Reorderb<-diag(u)
          diag(V_Full_Reorderb)<-diag(V_Full_Reorder)

          ##invert the reordered sampling covariance matrix to create a weight matrix
          if(toler==FALSE){
            W <- solve(V_Full_Reorderb)
          }

          if(toler!=FALSE){
            W <- solve(V_Full_Reorderb,tol=toler)
          }
        }

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

        ##smooth to near positive definite if either V or S are non-positive definite
        ks<-nrow(S_Fullrun)
        smooth1<-ifelse(eigen(S_Fullrun)$values[ks] <= 0, S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat), S_Fullrun<-S_Fullrun)

        if(modelchi == TRUE){
          ##name the columns and rows of the S matrix in general format V1-VX
          rownames(S_Fullrun) <- S_names
          colnames(S_Fullrun) <- S_names
        }

        if(modelchi == FALSE){
          #name the columns
          colnames(S_Fullrun)<-c("SNP", colnames(S_LD))

          ##name rows like columns
          rownames(S_Fullrun)<-colnames(S_Fullrun)
        }

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

        Model_table <- parTable(Model1_Results)

        resid_var <- 0
        if ( estimation=='DWLS' && class(test$value)[1] == "lavaan" && !any( grepl("solution has NOT", as.character(test$warning)) ) ) {
          resid_var_set <- subset( Model_table,
            Model_table$op == "~~" & Model_table$free != 0 & Model_table$lhs == Model_table$rhs
          )
          resid_var <- min( resid_var_set$est )
        } else { resid_var <- -9 }

        if ( estimation != 'DWLS' || resid_var > 0 ) {

          #pull the delta matrix (this doesn't depend on N)
          S2.delt <- lavInspect(Model1_Results, "delta")

          ##weight matrix from stage 2
          S2.W <- lavInspect(Model1_Results, "WLS.V")

          #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
          bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler)

          #create the "lettuce" part of the sandwich
          lettuce <- S2.W%*%S2.delt

          #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
          Ohtt <- bread %*% t(lettuce)%*%V_Full_Reorder%*%lettuce%*%bread

          #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
          SE <- as.matrix(sqrt(diag(Ohtt)))

          #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
          if ( estimation == 'ML' ) {
            if(":=" %in% Model_table$op){
              print("SEs of ghost parameters are not available for ML estimation")
              #pull the ghost parameter point estiamte
              ghost<-subset(Model_table, Model_table$op == ":=")[,c(2:4,8,11,14)]
              se.ghost<-rep("SE of ghost parameters not available for ML estimation", count(":=" %in% Model_table$op)$freq)
              ##combine with delta method SE
              ghost2<-cbind(ghost,se.ghost)
              colnames(ghost2)[7]<-"SE"
            }
          }
          if ( estimation == 'DWLS' ) {
            if(":=" %in% Model_table$op & !(NA %in% Model_table$se)){
              #variance-covariance matrix of parameter estimates, q-by-q (this is the naive one)
              vcov <- lavInspect(Model1_Results, "vcov")

              #internal lavaan representation of the model
              lavmodel <- Model1_Results@Model

              #lavaan representation of the indirect effect
              func <- lavmodel@def.function

              #vector of parameter estimates
              x <- lav_model_get_parameters(lavmodel, type = "free")

              #vector of indirect effect derivatives evaluated @ parameter estimates
              Jac <- lav_func_jacobian_complex(func = func, x = x)

              #replace vcov here with our corrected one. this gives parameter variance
              var.ind <- Jac %*% vcov %*% t(Jac)

              #square root of parameter variance = parameter SE.
              se.ghost <- sqrt(diag(var.ind))

              #pull the ghost parameter point estiamte
              ghost<-subset(Model_table, Model_table$op == ":=")[,c(2:4,8,11,14)]

              #combine with delta method SE
              ghost2<-cbind(ghost,se.ghost)
              colnames(ghost2)[7]<-"SE"
            }else{
              if(":=" %in% Model_table$op & (NA %in% Model_table$se)){
                se.ghost<-rep("SE could not be computed", count(":=" %in% Model_table$op)$freq)
                ghost<-subset(Model_table, Model_table$op == ":=")[,c(2:4,8,11,14)]
                ghost2<-cbind(ghost,se.ghost)
                colnames(ghost2)[7]<-"SE"}else{}}
          }

          if (modelchi == TRUE) {

            ModelQ_table <- parTable(Model1_Results)

            ##remove any parameter constraint labels or ghost parameters
            ModelQ_table<-subset(ModelQ_table, ModelQ_table$plabel != "")

            ##identify components of saturated model not already estimated and label as 1
            for (g in 1:nrow(ModelQ_table)){
              ModelQ_table$free[g]<-ifelse((paste(ModelQ_table$lhs[g], ModelQ_table$op[g], ModelQ_table$rhs[g], sep = "") %in% modeltest2$write.test.k) & ModelQ_table$est[g] == 0, 1, 0)
            }

            ##identify components of saturated model that were already estimated in user model and label as 2
            for (g in 1:nrow(ModelQ_table)){
              ModelQ_table$free[g]<-ifelse((paste(ModelQ_table$lhs[g], ModelQ_table$op[g], ModelQ_table$rhs[g], sep = "") %in% modeltest2$write.test.k) & ModelQ_table$est[g] != 0, 2, ModelQ_table$free[g])
            }

            tester<-vector(mode="list",length=nrow(ModelQ_table))

            ##replace components of saturdated model already estimated with residaul factor estimates
            for(g in 1:nrow(ModelQ_table)){
              if(ModelQ_table$free[g] == 2) {
                t<-paste(ModelQ_table$lhs[g], ModelQ_table$op[g], ModelQ_table$rhs[g], sep = "")
                t2<-gsub("V", "VF", t)
                tester[[g]]<-t2}else{}
            }

            test2<-Filter(Negate(is.null), tester)

            for (g in 1:nrow(ModelQ_table)){
              ModelQ_table$free[g]<-ifelse((paste(ModelQ_table$lhs[g], ModelQ_table$op[g], ModelQ_table$rhs[g], sep = "") %in% test2), 1, ModelQ_table$free[g])
            }

            ModelQ_table$free<-ifelse(ModelQ_table$free != 1, 0, ModelQ_table$free)

            #want to freely estimate the residual factor variances and the residual covariances
            p<-length(ModelQ_table$free)-z

            ModelQ_table <- ModelQ_table[order(ModelQ_table$free),]

            ModelQ_table$free <- c(rep(0, p),1:z)

            ModelQ_table$ustart <- ModelQ_table$est
            ModelQ_table$ustart<-ifelse(ModelQ_table$free > 0, .05, ModelQ_table$ustart)

            testQ <- tryCatch.W.E(
              ModelQ_Results_table <- sem(
                model = ModelQ_table,
                sample.cov = S_Fullrun,
                estimator = estimation,
                WLS.V = W,
                sample.nobs = S_LD_nobs,
                start = ModelQ_table$ustart,
                optim.dx.tol = +Inf
              )
            )
            if ( testQ$warning == '' ) testQ$warning = "Safe"
            testQ$warning <- ifelse( is.na(inspect(ModelQ_Results_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testQ$warning )

            if (as.character(testQ$warning) == "lavaan WARNING: model has NOT converged!") {

              ModelQ_table$ustart<-ifelse(ModelQ_table$free > 0, .01, ModelQ_table$ustart)

#                 NdF: W_Reorder is undefined; TEMPORARY FIX: replace with W
#                 testQ2<-tryCatch.W.E(ModelQ_Results_table <- sem(model = ModelQ_table, sample.cov = S_LD, estimator = estimation, WLS.V = W_Reorder, sample.nobs = S_LD_nobs, start=ModelQ$ustart, optim.dx.tol = +Inf))
              testQ2<-tryCatch.W.E(ModelQ_Results_table <- sem(model = ModelQ_table, sample.cov = S_LD, estimator = estimation, WLS.V = W, sample.nobs = S_LD_nobs, start=ModelQ$ustart, optim.dx.tol = +Inf))
            } else { testQ2<-testQ }

            if ( testQ2$warning == '' ) testQ2$warning = "Safe"
            testQ2$warning <- ifelse( is.na(inspect(ModelQ_Results_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testQ2$warning )

            if (as.character(testQ2$warning) == "lavaan WARNING: model has NOT converged!") {

              ModelQ_table$ustart<-ifelse(ModelQ_table$free > 0, .1, ModelQ_table$ustart)

#                 NdF: W_Reorder is undefined; TEMPORARY FIX: replace with W
#                 testQ3<-tryCatch.W.E(ModelQ_Results_table <- sem(model = ModelQ_table, sample.cov = S_LD, estimator = estimation, WLS.V = W_Reorder, sample.nobs = S_LD_nobs, start=ModelQ$ustart, optim.dx.tol = +Inf))
              testQ3<-tryCatch.W.E(ModelQ_Results_table <- sem(model = ModelQ_table, sample.cov = S_LD, estimator = estimation, WLS.V = W, sample.nobs = S_LD_nobs, start=ModelQ$ustart, optim.dx.tol = +Inf))
            } else { testQ3<-testQ2 }

            if ( testQ3$warning == '' ) testQ3$warning = "Safe"
            testQ3$warning <- ifelse( is.na(inspect(ModelQ_Results_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testQ3$warning )

            ModelQ_table2 <- parTable(ModelQ_Results_table)

            if (as.character(testQ3$warning) != "lavaan WARNING: model has NOT converged!" && !(NA %in% ModelQ_table2$se)) {

              #pull the delta matrix (this doesn't depend on N)
              S2.delt_Q <- lavInspect(ModelQ_Results_table, "delta")

              ##weight matrix from stage 2
              S2.W_Q <- lavInspect(ModelQ_Results_table, "WLS.V")

              #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
              bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q,tol=toler)

              #create the "lettuce" part of the sandwich
              lettuce_Q <- S2.W_Q%*%S2.delt_Q

              #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
              Ohtt_Q <- bread_Q %*% t(lettuce_Q)%*%V_Full_Reorder%*%lettuce_Q%*%bread_Q

              ##pull the sampling covariance matrix of the residual covariances and compute diagonal matrix of eigenvalues
              V_eta<- Ohtt_Q
              Eig2<-as.matrix(eigen(V_eta)$values)
              Eig<-diag(z)
              diag(Eig)<-Eig2

              #Pull P1 (the eigen vectors of V_eta)
              P1<-eigen(V_eta)$vectors

              ##Pull eta = vector of residual covariances
              eta_test<-parTable(ModelQ_Results_table)
              eta_test<-subset(eta_test, eta_test$free != 0)
              eta<-cbind(eta_test[,14])

              #Ronald's magic combining all the pieces from above:
              Qstat<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta
            } else { Qstat<-NA }
            
          }

          ##remove parameter constraints, ghost parameters, and fixed effects from output to merge with SEs
          unstand<-subset(Model_table, Model_table$plabel != "" & Model_table$free > 0)[,c(2:4,8,11,14)]

          ##combine ghost parameters with rest of output
          if(exists("ghost2") == "TRUE"){
            unstand2<-rbind(cbind(unstand,SE),ghost2)
          }else{unstand2<-cbind(unstand,SE)}

          ##add in fixed effects and parameter constraints to output
          other<-subset(Model_table, (Model_table$plabel == "" & Model_table$op != ":=") | (Model_table$free == 0 & Model_table$plabel != ""))[,c(2:4,8,11,14)]
          other$SE<-rep(NA, nrow(other))

          ##combine fixed effects and parameter constraints with output if there are any
          if(nrow(other) > 0){
            final<-rbind(unstand2,other)
          }else{final<-unstand2}

          #reorder based on row numbers so it is in order the user provided
          final$index <- as.numeric(row.names(final))
          final<-final[order(final$index), ]
          final$index<-NULL

          ##add in p-values
          if(class(final$SE) != "factor"){
            final$Z_Estimate<-final$est/final$SE
            final$Pval_Estimate<-2*pnorm(abs(final$Z_Estimate),lower.tail=FALSE)
          }else{
            final$SE<-as.character(final$SE)
            final$Z_Estimate<-NA
            final$Pval_Estimate<-NA
          }

          if(modelchi == TRUE){
            ##replace V1-VX general form in output with user provided trait names
            for(g in 1:nrow(final)){
              for(p in 1:length(traits)){
                final$lhs[[g]]<-ifelse(final$lhs[[g]] %in% S_names[[p]], gsub(final$lhs[[g]], traits[[p]], final$lhs[[g]]), final$lhs[[g]])
                final$rhs[[g]]<-ifelse(final$rhs[[g]] %in% S_names[[p]], gsub(final$rhs[[g]], traits[[p]], final$rhs[[g]]), final$rhs[[g]])
              }}
            ##subest to only the pieces the user specified
            ##this is overly long and confusing with the model fit components added in
            Modeltest <- parTable(ReorderModel1)
            for(t in 1:nrow(final)){
              final$test[t]<-ifelse(paste(final$lhs[t], final$op[t], final$rhs[t], sep = "") %in% paste(Modeltest$lhs, Modeltest$op, Modeltest$rhs, sep = ""), 1, 0)
              final$test[t]<-ifelse(paste(final$lhs[t], final$op[t], final$rhs[t], sep = "") %in% paste(Modeltest$rhs, Modeltest$op, Modeltest$lhs, sep = ""), 1, final$test[t])
            }
            final<-subset(final, final$test==1)
            final$test<-NULL

            ##add in model fit components to each row
            if(!(is.na(Qstat))){
              final$chisq<-rep(Qstat,nrow(final))
              final$chisq_df<-df
              final$chisq_pval<-pchisq(final$chisq,final$chisq_df,lower.tail=FALSE)
              final$AIC<-rep(Qstat + 2*npar,nrow(final))}else{final$chisq<-rep(NA, nrow(final))
              final$chisq_df<-rep(NA,nrow(final))
              final$chisq_pval<-rep(NA,nrow(final))
              final$AIC<-rep(NA, nrow(final))}
          }

          ##add in error and warning messages
          if(printwarn == TRUE){
            final$error<-ifelse( "lavaan" %in% class(test$value), 0, as.character(test$value$message) )
            final$warning<-ifelse( test$warning == '', 0, as.character(test$warning) )
          }

          ##combine with rs-id, BP, CHR, etc.
          final2<-cbind(SNPs.itw,final,row.names=NULL)

          if(!(sub[[1]])==FALSE){
            final2<-subset(final2, paste0(final2$lhs, final2$op, final2$rhs, sep = "") %in% sub)
          }else{##pull results -- NdF: final2$est == NA equals TRUE or FALSE
            final2$est<-ifelse(final2$op == "<" | final2$op == ">" | final2$op == ">=" | final2$op == "<=", final2$est == NA, final2$est)
          }

          ##results to be put into the output
          final2

        } else {
          if(modelchi == TRUE){
            final<-data.frame(t(rep(NA, 13)))
            if(printwarn == TRUE){
              final$error<-ifelse( "lavaan" %in% class(test$value), 0, as.character(test$value$message) )
              if(resid_var != -9){
                final$error<-paste0(
                  "This lavaan run produced negative (residual) variances for either latent or observed variables.",
                  "You may discard the run for this SNP, re-run the model with constraints to keep variances above 0,",
                  "or specify an alternative model.",
                  sep='\n'
                )
              }
              final$warning<-ifelse( test$warning == '', 0, as.character(test$warning) )
            }

            ##combine results with SNP, CHR, BP, A1, A2 for particular model
            final2<-cbind(SNPs.itw,final,row.names=NULL)
            colnames(final2)<-c("SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate", "chisq", "chisq_df", "chisq_pval", "AIC", "error", "warning")
          }

          if(modelchi == FALSE){
            final<-data.frame(t(rep(NA, 9)))
            if(printwarn == TRUE){
              final$error<-ifelse( "lavaan" %in% class(test$value), 0, as.character(test$value$message) )
              if(resid_var != -9){
                final$error<-paste0(
                  "This lavaan run produced negative (residual) variances for either latent or observed variables.",
                  "You may discard the run for this SNP, re-run the model with constraints to keep variances above 0,",
                  "or specify an alternative model.",
                  sep='\n'
                )
              }
              final$warning<-ifelse( test$warning == '', 0, as.character(test$warning) )
            }

            ##combine results with SNP, CHR, BP, A1, A2 for particular model
            final2<-cbind(SNPs.itw,final,row.names=NULL)
            colnames(final2)<-c("SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate", "error", "warning")
          }
          final2
        }

      }


    if(!(sub[[1]])==FALSE){
      Results_List<-vector(mode="list",length=length(sub))
      for(y in 1:length(sub)){
        Results_List[[y]]<-as.data.frame(matrix(NA,ncol=ncol(results),nrow=nrow(results)/length(sub)))
        colnames(Results_List[[y]])<-colnames(results)
        Results_List[[y]]<-subset(results, paste0(results$lhs, results$op, results$rhs, sep = "") %in% sub[[y]] | is.na(results$lhs))
      }
      rm(results)
    }


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

    time_all<-proc.time()-time
    print(time_all[3])
    return(Results_List)

  }

  if(!is.null(output)){

    ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
    output$RS$SNP<-as.character(output$RS$SNP)
    output$RS$A1<-as.character(output$RS$A1)
    output$RS$A2<-as.character(output$RS$A2)

    ##split the V and S matrices into as many (cores - 1) as are aviailable on the local computer
    V_Full<-suppressWarnings(split(output[[1]],1:int))
    S_Full<-suppressWarnings(split(output[[2]],1:int))
    SNPs<-suppressWarnings(split(output[[3]],1:int))

    ##enter in k for number of columns in S matrix
    k<-ncol(S_Full[[1]][[1]])

    ##number of models to run = number of distinct S/V matrices
    f=length(output[[1]])

    ##set the number of observations (typically LD-score bins) behind the covariance matrix
    S_LD_nobs<-covstruc$n  # this is what it should be (typically n-200), but causes long fitting time in 'DWLS' -> __NOBS_OVERRIDE__

    rm(output)

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

    Model1<-model

    if(modelchi == TRUE){
      #function to creat row/column names for S_LD matrix
      write.names <- function(k, label = "V") {
        varnames<-vector(mode="character",length=k)

        for (i in 1:k){
          varnames[i]<-paste(label,i,sep="")}

        return(varnames)
      }

      if(toler==FALSE){
        W_test <- solve(V_Full[[1]][[1]])
      }

      if(toler!=FALSE){
        W_test <- solve(V_Full[[1]][[1]],tol=toler)
      }

      S_Fulltest<-S_Full[[1]][[1]]

      ##run the model to determine number of latent variables
      suppress <- tryCatch.W.E(
        ReorderModel1 <- sem(
          model,
          sample.cov = S_Fulltest,
          estimator = "DWLS",
          WLS.V = W_test,
          sample.nobs = 2,
          warn = FALSE,
          optim.dx.tol = +Inf,
          optim.force.converged = TRUE,
          control = list(iter.max=1)
        )
      )

      ##pull the column names specified in the munge function
      traits<-colnames(S_Full[[1]][[1]])

      ##create the names
      S_names<-write.names(k=k)

      ##add bracketing so gsub knows to replace exact cases
      traits2<-traits
      for(i in 1:length(traits)){
        traits2[[i]]<-paste0("\\<", traits[[i]],"\\>",sep="")
      }

      ##replace trait names in user provided model with general form of V1-VX
      for(i in 1:length(traits)){
        model<-gsub(traits2[[i]], S_names[[i]], model)
      }

      ##determine number of latent variables from writing extended model.
      r<-nrow(lavInspect(ReorderModel1, "cor.lv"))
      lat_labs<-colnames(lavInspect(ReorderModel1, "cor.lv"))

      write.Model1 <- function(k, label = "V", label2 = "VF") {

        ModelsatF<-""
        for (i in 1:(k-1)) {
          linestartc <- paste(" ", label2, i, "~~0*", label2, i+1,  sep = "")
          if (k-i >= 2) {
            linemidc <- ""
            for (j in (i+2):k) {
              linemidc <- paste(linemidc, "+0*", label2, j, sep = "")
            }
          } else {linemidc <- ""}
          ModelsatF <- paste(ModelsatF, linestartc, linemidc, " \n ", sep = "")
        }

        if(r > 0){
          Model1b <- ""
          for (t in 1:r) {
            for (i in 1) {
              linestartb <- paste(lat_labs[t], " =~ 0*",label2, i, sep = "")
              if ((k-1)-i > 0 | k == 2) {
                linemidb <- ""
                for (j in (i+1):k) {
                  linemidb <- paste(linemidb, " + 0*", label2, j, sep = "")
                }
              } else {linemidb <- ""}

            }
            Model1b <- paste(Model1b, linestartb, linemidb, " \n ", sep = "")
          }
        }
        else {Model1b <- ""}

        if(r > 0){
          Model1c <- ""
          for (t in 1:r) {
            for (i in 1) {
              linestartc <- paste(lat_labs[t], " ~~ 0*",label2, i, sep = "")
              if ((k-1)-i > 0 | k == 2) {
                linemidc <- ""
                for (j in (i+1):k) {
                  linemidc <- paste(linemidc, " + 0*", label2, j, sep = "")
                }
              } else {linemidc <- ""}

            }
            Model1c <- paste(Model1c, linestartc, linemidc, " \n ", sep = "")
          }
        }
        else {Model1c <- ""}

        Model2<-""
        for (p in 1:k) {
          linestart2 <- paste(label2, p, " =~ 1*", label, p, sep = "")
          Model2<-paste(Model2, linestart2, " \n ", sep = "")}

        Model3<-""
        for (p in 1:k) {
          linestart3 <- paste(label, p, "~~", label, p, sep = "")
          Model3<-paste(Model3, linestart3, " \n ", sep = "")}

        Model4<-""
        for (p in 1:k) {
          linestart4 <- paste(label2, p, "~~0*", label2, p, sep = "")
          Model4<-paste(Model4, linestart4, " \n ", sep = "")}

        Modelsat<-""
        for (i in 1:(k-1)) {
          linestartc <- paste("", label, i, "~~0*", label, i+1, sep = "")
          if (k-i >= 2) {
            linemidc <- ""
            for (j in (i+2):k) {
              linemidc <- paste("", linemidc, label, i, "~~0*", label, j, " \n ", sep="")

            }
          } else {linemidc <- ""}
          Modelsat <- paste(Modelsat, linestartc, " \n ", linemidc, sep = "")
        }

        Model5<-paste(model, " \n ", ModelsatF, Model1b, Model1c, Model2, Model3, Model4, Modelsat, sep = "")

        return(Model5)
      }

      Model1<-write.Model1(k)

      while ( class(u<-tryCatch.W.E(lavParseModelString(Model1))$value$message ) != 'NULL') {
        t<-paste(strsplit(u, ": ")[[1]][3], " \n ", sep = "")
        Model1<-str_replace(Model1, fixed(t), "")
      }

      ##code to write fake model to make sure elements estimated in user model
      ##do not overlap with saturated model, in which case an alternative
      ##specification is used that uses, for example, covariances among residual factors
      write.test<-function(k, label = "V", label2 = "VF") {

        Modelsat<-""
        for (i in 1:(k)) {
          if (k-i >= 1) {
            linemidc <- ""
            for (j in (i+1):k) {
              linemidc <- paste(linemidc, label, i, "~~", label, j, " \n ", sep = "")
            }
          }else{linemidc<-""}
          Modelsat <- paste(Modelsat, linemidc, sep = "")
        }

        Modelsat2<-""
        for (i in 1:(k)) {
          if (k-i >= 1) {
            linemidc2 <- ""
            for (j in (i+1):k) {
              linemidc2 <- paste(linemidc2, label, j, "~~", label, i, " \n ", sep = "")
            }
          }else{linemidc2<-""}
          Modelsat2 <- paste(Modelsat2, linemidc2, sep = "")
        }

        Model4<-""
        for (p in 1:k) {
          linestart4 <- paste(label2, p, "~~", label2, p, sep = "")
          Model4<-paste(Model4, linestart4, " \n ", sep = "")}

        modelCFI<-paste(Modelsat, Modelsat2, Model4)
        return(modelCFI)
      }


      modeltest<-data.frame(write.test(k))
      modeltest$write.test.k.<-as.character(modeltest$write.test.k.)
      modeltest2 <- cSplit(modeltest, "write.test.k.", sep = "\n", direction = "long")
      modeltest2$write.test.k.<-as.character(modeltest2$write.test.k.)

      z<-(k*(k+1))/2
    }

    ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
    for (i in 1) {

      if(toler==FALSE){
        W <- solve(V_Full[[1]][[i]])
      }

      if(toler!=FALSE){
        W <- solve(V_Full[[1]][[i]],tol=toler)
      }

      if(modelchi == TRUE){
        ##name the columns and rows of the S matrix in general format V1-VX
        rownames(S_Full[[1]][[i]]) <- S_names
        colnames(S_Full[[1]][[i]]) <- S_names
      }

      S_Fullrun<-S_Full[[1]][[i]]

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

      if ( class(suppress$value)[1] == "lavaan" ) {
        suppressWarnings(df<-lavInspect(ReorderModel, "fit")["df"])
        suppressWarnings(npar<-lavInspect(ReorderModel, "fit")["npar"])
        order <- rearrange(k = k, fit = ReorderModel, names = rownames(S_Full[[1]][[i]]))
      } else {

        i<-10

        if(toler==FALSE){
          W <- solve(V_Full[[1]][[i]])
        }

        if(toler!=FALSE){
          W <- solve(V_Full[[1]][[i]],tol=toler)
        }

        if(modelchi == TRUE){
          ##name the columns and rows of the S matrix in general format V1-VX
          rownames(S_Full[[1]][[i]]) <- S_names
          colnames(S_Full[[1]][[i]]) <- S_names
        }

        S_Fullrun<-S_Full[[1]][[i]]

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

        suppressWarnings(df<-lavInspect(ReorderModel, "fit")["df"])
        suppressWarnings(npar<-lavInspect(ReorderModel, "fit")["npar"])
        order <- rearrange(k = k, fit = ReorderModel, names = rownames(S_Full[[1]][[i]]))

      }

    }

    workChunks = rep_len( 1:int, dim(SNPs)[1] )

    ## foreach (parallel) processing that rbinds results across cores
    cat( "running lavaan common factor models including individual variants..\n" )
    results <- foreach (
      V_Full.it=iterators::isplit(output[[1]],workChunks),
      S_Full.it=iterators::isplit(output[[2]],workChunks),
      SNPs.it=iterators::isplit(output[[3]],workChunks),
      .combine=rbind, .errorhandling='remove',
      .noexport='output', .packages='lavaan'
      ) %:% foreach (
      V_Full.itw=iterators::iter(V_Full.it$value, by='row'),
      S_Full.itw=iterators::iter(S_Full.it$value, by='row'),
      SNPs.itw=iterators::iter(SNPs.it$value, by='row'),
      .combine=rbind, .multicombine=T, .errorhandling='remove',
      .noexport='output', .packages='lavaan'
    ) %do.operator% {

        #reorder sampling covariance matrix based on what lavaan expects given the specified model
        V_Full_Reorder <- V_Full.itw[order,order]

        W = NULL
        S_LD_nobs = 200  # <- __NOBS_OVERRIDE__
        if ( estimation == 'DWLS' ) {
          S_LD_nobs = 2  # <- __NOBS_OVERRIDE__
          u<-nrow(V_Full_Reorder)
          V_Full_Reorderb<-diag(u)
          diag(V_Full_Reorderb)<-diag(V_Full_Reorder)

          ##invert the reordered sampling covariance matrix to create a weight matrix
          if(toler==FALSE){
            W <- solve(V_Full_Reorderb)
          }

          if(toler!=FALSE){
            W <- solve(V_Full_Reorderb,tol=toler)
          }
        }

        #import the S_Full matrix for appropriate run
        S_Fullrun<-S_Full.itw

        if(modelchi == TRUE){
          ##name the columns and rows of the S matrix in general format V1-VX
          rownames(S_Fullrun) <- S_names
          colnames(S_Fullrun) <- S_names
        }

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

        Model_table <- parTable(Model1_Results)

        resid_var <- 0
        if ( estimation=='DWLS' && class(test$value)[1] == "lavaan" && !any( grepl("solution has NOT", as.character(test$warning)) ) ) {
          resid_var_set <- subset( Model_table,
            Model_table$op == "~~" & Model_table$free != 0 & Model_table$lhs == Model_table$rhs
          )
          resid_var <- min( resid_var_set$est )
        } else { resid_var <- -9 }

        if ( estimation != 'DWLS' || resid_var > 0 ) {

          #pull the delta matrix (this doesn't depend on N)
          S2.delt <- lavInspect(Model1_Results, "delta")

          ##weight matrix from stage 2
          S2.W <- lavInspect(Model1_Results, "WLS.V")

          #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
          bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler)

          #create the "lettuce" part of the sandwich
          lettuce <- S2.W%*%S2.delt

          #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
          Ohtt <- bread %*% t(lettuce)%*%V_Full_Reorder%*%lettuce%*%bread

          #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
          SE <- as.matrix(sqrt(diag(Ohtt)))

          #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
          if ( estimation == 'ML' ) {
            if(":=" %in% Model_table$op){
              print("SEs of ghost parameters are not available for ML estimation")
              #pull the ghost parameter point estiamte
              ghost<-subset(Model_table, Model_table$op == ":=")[,c(2:4,8,11,14)]
              se.ghost<-rep("SE of ghost parameters not available for ML estimation", count(":=" %in% Model_table$op)$freq)
              ##combine with delta method SE
              ghost2<-cbind(ghost,se.ghost)
              colnames(ghost2)[7]<-"SE"
            }
          }
          if ( estimation == 'DWLS' ) {
            if(":=" %in% Model_table$op & !(NA %in% Model_table$se)){
              #variance-covariance matrix of parameter estimates, q-by-q (this is the naive one)
              vcov <- lavInspect(Model1_Results, "vcov")

              #internal lavaan representation of the model
              lavmodel <- Model1_Results@Model

              #lavaan representation of the indirect effect
              func <- lavmodel@def.function

              #vector of parameter estimates
              x <- lav_model_get_parameters(lavmodel, type = "free")

              #vector of indirect effect derivatives evaluated @ parameter estimates
              Jac <- lav_func_jacobian_complex(func = func, x = x)

              #replace vcov here with our corrected one. this gives parameter variance
              var.ind <- Jac %*% vcov %*% t(Jac)

              #square root of parameter variance = parameter SE.
              se.ghost <- sqrt(diag(var.ind))

              #pull the ghost parameter point estiamte
              ghost<-subset(Model_table, Model_table$op == ":=")[,c(2:4,8,11,14)]

              ##combine with delta method SE
              ghost2<-cbind(ghost,se.ghost)
              colnames(ghost2)[7]<-"SE"
            }else{
              if(":=" %in% Model_table$op & (NA %in% Model_table$se)){
                se.ghost<-rep("SE could not be computed", count(":=" %in% Model_table$op)$freq)
                ghost<-subset(Model_table, Model_table$op == ":=")[,c(2:4,8,11,14)]
                ghost2<-cbind(ghost,se.ghost)
                colnames(ghost2)[7]<-"SE"
              }
            }
          }

          if (modelchi == TRUE) {

            ModelQ_table <- parTable(Model1_Results)

            ##remove any parameter constraint labels or ghost parameters
            ModelQ_table<-subset(ModelQ_table, ModelQ_table$plabel != "")

            ##identify components of saturated model not already estimated and label as 1
            for (g in 1:nrow(ModelQ_table)){
              ModelQ_table$free[g]<-ifelse((paste(ModelQ_table$lhs[g], ModelQ_table$op[g], ModelQ_table$rhs[g], sep = "") %in% modeltest2$write.test.k) & ModelQ_table$est[g] == 0, 1, 0)
            }

            ##identify components of saturated model that were already estimated in user model and label as 2
            for (g in 1:nrow(ModelQ_table)){
              ModelQ_table$free[g]<-ifelse((paste(ModelQ_table$lhs[g], ModelQ_table$op[g], ModelQ_table$rhs[g], sep = "") %in% modeltest2$write.test.k) & ModelQ_table$est[g] != 0, 2, ModelQ_table$free[g])
            }

            tester<-vector(mode="list",length=nrow(ModelQ_table))

            ##replace components of saturdated model already estimated with residaul factor estimates
            for(g in 1:nrow(ModelQ_table)){
              if(ModelQ_table$free[g] == 2) {
                t<-paste(ModelQ_table$lhs[g], ModelQ_table$op[g], ModelQ_table$rhs[g], sep = "")
                t2<-gsub("V", "VF", t)
                tester[[g]]<-t2}else{}
            }

            test2<-Filter(Negate(is.null), tester)

            for (g in 1:nrow(ModelQ_table)){
              ModelQ_table$free[g]<-ifelse((paste(ModelQ_table$lhs[g], ModelQ_table$op[g], ModelQ_table$rhs[g], sep = "") %in% test2), 1, ModelQ_table$free[g])
            }

            ModelQ_table$free<-ifelse(ModelQ_table$free != 1, 0, ModelQ_table$free)

            #want to freely estimate the residual factor variances and the residual covariances
            p<-length(ModelQ_table$free)-z

            ModelQ_table <- ModelQ_table[order(ModelQ_table$free),]

            ModelQ_table$free <- c(rep(0, p),1:z)

            ModelQ_table$ustart <- ModelQ_table$est
            ModelQ_table$ustart<-ifelse(ModelQ_table$free > 0, .05, ModelQ_table$ustart)

            testQ <- tryCatch.W.E(
              ModelQ_Results_table <- sem(
                model = ModelQ_table,
                sample.cov = S_Fullrun,
                estimator = estimation,
                WLS.V = W,
                sample.nobs = S_LD_nobs,
                start = ModelQ_table$ustart,
                optim.dx.tol = +Inf
              )
            )
            if ( testQ$warning == '' ) testQ$warning = "Safe"
            testQ$warning <- ifelse( is.na(inspect(ModelQ_Results_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testQ$warning )

            if (as.character(testQ$warning) == "lavaan WARNING: model has NOT converged!") {

              ModelQ_table$ustart<-ifelse(ModelQ_table$free > 0, .01, ModelQ_table$ustart)

#                 NdF: W_Reorder is undefined; TEMPORARY FIX: replace with W
#                 testQ2<-tryCatch.W.E(ModelQ_Results_table <- sem(model = ModelQ_table, sample.cov = S_LD, estimator = estimation, WLS.V = W_Reorder, sample.nobs = S_LD_nobs, start=ModelQ$ustart, optim.dx.tol = +Inf))
              testQ2<-tryCatch.W.E(ModelQ_Results_table <- sem(model = ModelQ_table, sample.cov = S_LD, estimator = estimation, WLS.V = W, sample.nobs = S_LD_nobs, start=ModelQ$ustart, optim.dx.tol = +Inf))
            } else { testQ2<-testQ }

            if ( testQ2$warning == '' ) testQ2$warning = "Safe"
            testQ2$warning <- ifelse( is.na(inspect(ModelQ_Results_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testQ2$warning )

            if (as.character(testQ2$warning) == "lavaan WARNING: model has NOT converged!") {

              ModelQ_table$ustart<-ifelse(ModelQ_table$free > 0, .1, ModelQ_table$ustart)

#                 NdF: W_Reorder is undefined; TEMPORARY FIX: replace with W
#                 testQ3<-tryCatch.W.E(ModelQ_Results_table <- sem(model = ModelQ_table, sample.cov = S_LD, estimator = estimation, WLS.V = W_Reorder, sample.nobs = S_LD_nobs, start=ModelQ$ustart, optim.dx.tol = +Inf))
              testQ3<-tryCatch.W.E(ModelQ_Results_table <- sem(model = ModelQ_table, sample.cov = S_LD, estimator = estimation, WLS.V = W, sample.nobs = S_LD_nobs, start=ModelQ$ustart, optim.dx.tol = +Inf))
            } else { testQ3<-testQ2 }

            if ( testQ3$warning == '' ) testQ3$warning = "Safe"
            testQ3$warning <- ifelse( is.na(inspect(ModelQ_Results_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testQ3$warning )

            ModelQ_table2 <- parTable(ModelQ_Results_table)

            if (as.character(testQ3$warning) != "lavaan WARNING: model has NOT converged!" && !(NA %in% ModelQ_table2$se)) {

              #pull the delta matrix (this doesn't depend on N)
              S2.delt_Q <- lavInspect(ModelQ_Results_table, "delta")

              ##weight matrix from stage 2
              S2.W_Q <- lavInspect(ModelQ_Results_table, "WLS.V")

              #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
              bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q,tol=toler)

              #create the "lettuce" part of the sandwich
              lettuce_Q <- S2.W_Q%*%S2.delt_Q

              #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
              Ohtt_Q <- bread_Q %*% t(lettuce_Q)%*%V_Full_Reorder%*%lettuce_Q%*%bread_Q

              ##pull the sampling covariance matrix of the residual covariances and compute diagonal matrix of eigenvalues
              V_eta<- Ohtt_Q
              Eig2<-as.matrix(eigen(V_eta)$values)
              Eig<-diag(z)
              diag(Eig)<-Eig2

              #Pull P1 (the eigen vectors of V_eta)
              P1<-eigen(V_eta)$vectors

              ##Pull eta = vector of residual covariances
              eta_test<-parTable(ModelQ_Results_table)
              eta_test<-subset(eta_test, eta_test$free != 0)
              eta<-cbind(eta_test[,14])

              #Ronald's magic combining all the pieces from above:
              Qstat<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta

            } else { Qstat<-NA }

          }

          ##remove parameter constraints, ghost parameters, and fixed effects from output to merge with SEs
          unstand<-subset(Model_table, Model_table$plabel != "" & Model_table$free > 0)[,c(2:4,8,11,14)]

          ##combine ghost parameters with rest of output
          if(exists("ghost2") == "TRUE"){
            unstand2<-rbind(cbind(unstand,SE),ghost2)
          }else{unstand2<-cbind(unstand,SE)}

          ##add in fixed effects and parameter constraints to output
          other<-subset(Model_table, (Model_table$plabel == "" & Model_table$op != ":=") | (Model_table$free == 0 & Model_table$plabel != ""))[,c(2:4,8,11,14)]
          other$SE<-rep(NA, nrow(other))

          ##combine fixed effects and parameter constraints with output if there are any
          if(nrow(other) > 0){
            final<-rbind(unstand2,other)
          }else{final<-unstand2}

          #reorder based on row numbers so it is in order the user provided
          final$index <- as.numeric(row.names(final))
          final<-final[order(final$index), ]
          final$index<-NULL

          ##add in p-values
          if(class(final$SE) != "factor"){
            final$Z_Estimate<-final$est/final$SE
            final$Pval_Estimate<-2*pnorm(abs(final$Z_Estimate),lower.tail=FALSE)
          }else{
            final$SE<-as.character(final$SE)
            final$Z_Estimate<-NA
            final$Pval_Estimate<-NA}

          if(modelchi == TRUE){
            ##replace V1-VX general form in output with user provided trait names
            for(g in 1:nrow(final)){
              for(p in 1:length(traits)){
                final$lhs[[g]]<-ifelse(final$lhs[[g]] %in% S_names[[p]], gsub(final$lhs[[g]], traits[[p]], final$lhs[[g]]), final$lhs[[g]])
                final$rhs[[g]]<-ifelse(final$rhs[[g]] %in% S_names[[p]], gsub(final$rhs[[g]], traits[[p]], final$rhs[[g]]), final$rhs[[g]])
              }}
            ##subest to only the pieces the user specified
            ##this is overly long and confusing with the model fit components added in
            Modeltest <- parTable(ReorderModel1)
            for(t in 1:nrow(final)){
              final$test[t]<-ifelse(paste(final$lhs[t], final$op[t], final$rhs[t], sep = "") %in% paste(Modeltest$lhs, Modeltest$op, Modeltest$rhs, sep = ""), 1, 0)
              final$test[t]<-ifelse(paste(final$lhs[t], final$op[t], final$rhs[t], sep = "") %in% paste(Modeltest$rhs, Modeltest$op, Modeltest$lhs, sep = ""), 1, final$test[t])
            }
            final<-subset(final, final$test==1)
            final$test<-NULL

            ##add in model fit components to each row
            if(!(is.na(Qstat))){
              final$chisq<-rep(Qstat,nrow(final))
              final$chisq_df<-df
              final$chisq_pval<-pchisq(final$chisq,final$chisq_df,lower.tail=FALSE)
              final$AIC<-rep(Qstat + 2*npar,nrow(final))}else{final$chisq<-rep(NA, nrow(final))
              final$chisq_df<-rep(NA,nrow(final))
              final$chisq_pval<-rep(NA,nrow(final))
              final$AIC<-rep(NA, nrow(final))}
          }

          ##add in error and warning messages
          if(printwarn == TRUE){
            final$error<-ifelse( "lavaan" %in% class(test$value), 0, as.character(test$value$message) )
            final$warning<-ifelse( test$warning == '', 0, as.character(test$warning) )
          }

          ##combine with rs-id, BP, CHR, etc.
          final2<-cbind(SNPs.itw,final,row.names=NULL)

          if(!(sub[[1]])==FALSE){
            final2<-subset(final2, paste0(final2$lhs, final2$op, final2$rhs, sep = "") %in% sub)
          }else{##pull results -- NdF: final2$est == NA equals TRUE or FALSE
            final2$est<-ifelse(final2$op == "<" | final2$op == ">" | final2$op == ">=" | final2$op == "<=", final2$est == NA, final2$est)
          }

          ##results to be put into the output
          final2

        }else{
          if(modelchi == TRUE){
            final<-data.frame(t(rep(NA, 13)))
            if(printwarn == TRUE){
              final$error<-ifelse( "lavaan" %in% class(test$value), 0, as.character(test$value$message) )
              if(resid_var != -9){
                final$error<-paste0(
                  "This lavaan run produced negative (residual) variances for either latent or observed variables.",
                  "You may discard the run for this SNP, re-run the model with constraints to keep variances above 0,",
                  "or specify an alternative model.",
                  sep='\n'
                )
              }
              final$warning<-ifelse( test$warning == '', 0, as.character(test$warning) )
            }

            ##combine results with SNP, CHR, BP, A1, A2 for particular model
            final2<-cbind(SNPs.itw,final,row.names=NULL)
            colnames(final2)<-c("SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate", "chisq", "chisq_df", "chisq_pval", "AIC", "error", "warning")
          }

          if(modelchi == FALSE){
            final<-data.frame(t(rep(NA, 9)))
            if(printwarn == TRUE){
              final$error<-ifelse( "lavaan" %in% class(test$value), 0, as.character(test$value$message) )
              if(resid_var != -9){
                final$error<-paste0(
                  "This lavaan run produced negative (residual) variances for either latent or observed variables.",
                  "You may discard the run for this SNP, re-run the model with constraints to keep variances above 0,",
                  "or specify an alternative model.",
                  sep='\n'
                )
              }
              final$warning<-ifelse( test$warning == '', 0, as.character(test$warning) )
            }

            ##combine results with SNP, CHR, BP, A1, A2 for particular model
            final2<-cbind(SNPs.itw,final,row.names=NULL)
            colnames(final2)<-c("SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate", "error", "warning")
          }
          final2
        }

      }


    if(!(sub[[1]])==FALSE){
      Results_List<-vector(mode="list",length=length(sub))
      for(y in 1:length(sub)){
        Results_List[[y]]<-as.data.frame(matrix(NA,ncol=ncol(results),nrow=nrow(results)/length(sub)))
        colnames(Results_List[[y]])<-colnames(results)
        Results_List[[y]]<-subset(results, paste0(results$lhs, results$op, results$rhs, sep = "") %in% sub[[y]] | is.na(results$lhs))
      }
      rm(results)
    }


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

    time_all<-proc.time()-time
    print(time_all[3])

    return(Results_List)

  }

}

userGWAS <- usermodelGWAS

