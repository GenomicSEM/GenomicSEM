usermodel <- function(covstruc, estimation="DWLS", model = "", CFIcalc=TRUE, std.lv=FALSE, imp_cov=FALSE, fix_resid=TRUE) { 
  time<-proc.time()
  ##determine if the model is likely being listed in quotes and print warning if so
  test<-c(str_detect(model, "~"),str_detect(model, "="),str_detect(model, "\\+"))
  if(any(test) != TRUE){
    warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned an error about not locating the ReorderModel.")
  }
  
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
  
  ##read in the LD portion of the V (sampling covariance) matrix
  V_LD<-as.matrix(covstruc$V)
  
  ##read in the LD portion of the S (covariance) matrix
  S_LD<-as.matrix(covstruc$S)
  
  ##set the number of observations (typically LD-score bins) behind the covariance matrix
  S_LD_nobs<-covstruc$n  # this is what it should be, but see -> __NOBS_OVERRIDE__
  
  ##k = number of phenotypes in dataset (i.e., number of columns in LD portion of S matrix)
  k<-ncol(S_LD)
  
  ##size of V matrix used later in code to create diagonal V matrix
  z<-(k*(k+1))/2
  
  #function to creat row/column names for S_LD matrix
  write.names <- function(k, label = "V") {
    varnames<-vector(mode="character",length=k)
    
    for (i in 1:k){
      varnames[i]<-paste(label,i,sep="")
    }
    
    return(varnames)
  }
  
  ##create the names
  S_names<-write.names(k=k)
  
  ##pull the column names specified in the munge function
  traits<-colnames(S_LD)
  
  ##rename in general form in case using previous version of Genomic SEM
  if(is.null(traits)){
    traits<-S_names}
  
  ##add bracketing so gsub knows to replace exact cases
  traits2<-traits
  for(i in 1:length(traits)){
    traits2[[i]]<-paste0("\\<", traits[[i]],"\\>",sep="")
  }
  
  ##replace trait names in user provided model with general form of V1-VX
  for(i in 1:length(traits)){
    model<-gsub(traits2[[i]], S_names[[i]], model)
  }
  
  Model1<-model
  
  ##name the columns and rows of the S matrix
  rownames(S_LD) <- S_names
  colnames(S_LD) <- S_names
  
  cat( 'renaming V matrix entries..\n' )
  ##name columns of V to remove any variables not used in the current analysis
  y<-expand.grid(S_names,S_names)
  y<-y[!duplicated(apply(y,1,function(x) paste(sort(x),collapse=''))),]
  V_Names<-paste(y$Var1,y$Var2,sep=" ")
  colnames(V_LD)<-V_Names
  rownames(V_LD)<-V_Names
  
  ##determine whether all variables in S are in the model
  ##if not, remove them from S_LD and V_LD for this particular run
  remove2<-c()
  w<-1
  
  ##also for exact cases
  for(i in 1:length(S_names)){
    S_names[[i]]<-paste0("\\b", S_names[[i]],"\\b",sep="")
  }
  i<-1
  for(i in 1:length(S_names)){
    b<-grepl(S_names[i], model)
    if(b == FALSE){
      remove<-paste0("\\b", colnames(S_LD)[i],"\\b",sep="")
      remove2[w]<-i
      V_LD <- V_LD[-grep(pattern=remove[1],row.names(V_LD)),-grep(pattern=remove[1],colnames(V_LD))]
      w<-w+1
    }
  }
  
  if(is.null(remove2) == FALSE){
    S_LD<-S_LD[-remove2,-remove2]
    traits<-traits[-remove2]
  }
  
  ##redefine k and z and model names after removing non-used variables
  k<-ncol(S_LD)
  z<-(k*(k+1))/2
  before<-colnames(S_LD)
  S_names<-write.names(k=k)
  colnames(S_LD)<-S_names
  rownames(S_LD)<-S_names
  
  for(i in 1:length(before)){
    before[[i]]<-paste0("\\<", before[[i]],"\\>",sep="")
  }
  
  for(i in 1:length(traits)){
    model<-gsub(before[[i]], S_names[[i]], model)
  }
  Model1<-model
  
  cat( 'smoothing V and S matrices..\n' )
  ##smooth to near positive definite if either V or S are non-positive definite
  ks<-nrow(S_LD)
  S_LDb<-S_LD
  smooth1<-ifelse(eigen(S_LD)$values[ks] <= 0, S_LD<-as.matrix((nearPD(S_LD, corr = FALSE))$mat), S_LD<-S_LD)
  LD_sdiff<-max(abs(S_LD-S_LDb))

  kv<-nrow(V_LD)
  V_LDb<-V_LD
  smooth2<-ifelse(eigen(V_LD)$values[kv] <= 0, V_LD<-as.matrix((nearPD(V_LD, corr = FALSE))$mat), V_LD<-V_LD)
  LD_sdiff2<-max(abs(V_LD-V_LDb))

  SE_pre<-matrix(0, k, k)
  SE_pre[lower.tri(SE_pre,diag=TRUE)] <-sqrt(diag(V_LDb))
  
  SE_post<-matrix(0, k, k)
  SE_post[lower.tri(SE_post,diag=TRUE)] <-sqrt(diag(V_LD))
  
  Z_pre<-S_LDb/SE_pre
  Z_post<-S_LD/SE_post
  Z_diff<-max(abs(Z_pre-Z_post),na.rm=T)

  cat( 'initializing lavaan model factor structure..\n' )
  ##run model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  #transform V_LD matrix into a weight matrix: 
  W <- solve(V_LD)
  
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

  ##run the model
  
  empty2 <- tryCatch.W.E(
    ReorderModel1 <- sem(
      Model1,
      sample.cov = S_LD,
      estimator = "DWLS",  #should this not depend on the estimation method?
      WLS.V = W,           #should this?
      sample.nobs = 2,
      warn = FALSE,
      std.lv = std.lv,
      optim.dx.tol = +Inf,
      optim.force.converged = TRUE,
      control = list(iter.max=1)
    )
  )
  
  if ( !( "lavaan" %in% class(empty2$value) ) ) {
    latentcorr <- any( grepl("not defined:", empty2$value$message[1][1]) )
    latentcorr2 <- any( grepl("unknown", empty2$value$message[1][1]) )
    if ( latentcorr || latentcorr2 ) {
      warning( paste(
        "The function may have stopped either because a variable has been misnamed in the model or because you have",
        "tried to estimate a correlation between an observed and latent variable. In the latter case, one workaround",
        "is to define a latent variable solely by the observed variable.", sep='\n'
      ) )
    }
  }
  
  ##determine number of latent variables from writing extended model
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
    } else {Model1b <- ""}
    
    
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
    } else {Model1c <- ""}
    
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
  
  cat( 'pruning lavaan model..\n' )
  ##code to remove duplicated elements between user/automatically specified Model Input
  while(class(tryCatch.W.E(lavParseModelString(Model1))$value$message) != 'NULL'){
    u<-tryCatch.W.E(lavParseModelString(Model1))$value$message
    t<-paste(strsplit(u, ": ")[[1]][3], " \n ", sep = "")
    Model1<-str_replace(Model1, fixed(t), "")
  }
  
  if (CFIcalc==TRUE) {

    print("Calculating CFI..")

    ##code to write null model for calculation of CFI
    write.null<-function(k, label = "V", label2 = "VF") {
      Model3<-""
      for (p in 1:k) {
        linestart3 <- paste(label, p, " ~~ ", label, p, sep = "")
        Model3<-paste(Model3, linestart3, " \n ", sep = "")}
      
      Model2<-""
      for (p in 1:k) {
        linestart2 <- paste(label2, p, " =~ 1*", label, p, sep = "")
        Model2<-paste(Model2, linestart2, " \n ", sep = "")}
      
      Modelsat<-""
      for (i in 1:(k-1)) {
        linestartc <- paste(label, i, " ~~ 0*", label, i+1,  sep = "")
        if (k-i >= 2) {
          linemidc <- ""
          for (j in (i+2):k) {
            linemidc <- paste(linemidc, " + 0*", label, j, sep = "")
          }
        } else {linemidc <- ""}
        Modelsat <- paste(Modelsat, linestartc, linemidc, " \n ", sep = "")
      }
      
      ModelsatF<-""
      for (i in 1:(k-1)) {
        linestartc <- paste(" ", label2, i, " ~~ 0*", label2, i+1,  sep = "")
        if (k-i >= 2) {
          linemidc <- ""
          for (j in (i+2):k) {
            linemidc <- paste(linemidc, " + 0*", label2, j, sep = "")
          }
        } else {linemidc <- ""}
        ModelsatF <- paste(ModelsatF, linestartc, linemidc, " \n ", sep = "")
      }
      
      Model4<-""
      for (p in 1:k) {
        linestart4 <- paste(label2, p, " ~~ 0*", label2, p, sep = "")
        Model4<-paste(Model4, linestart4, " \n ", sep = "")}
      
      modelCFI<-paste(Model3, Model2, ModelsatF, Modelsat, Model4)

      return(modelCFI)
    }
    
    ##create independence model for calculation of CFI
    modelCFI<-write.null(k)
    
    ##run CFI model so it knows the reordering for the independence model
    empty <- tryCatch.W.E(
      fitCFI <- sem(
        modelCFI,
        sample.cov = S_LD,
        estimator = "DWLS",  #should this not depend on the estimation method?
        WLS.V = W,           #should this?
        sample.nobs = 2,
        optim.dx.tol = +Inf
      )
    )
    
    orderCFI <- rearrange(k = k, fit = fitCFI, names = rownames(S_LD))
    
    ##reorder matrix for independence (i.e., null) model for CFI calculation
    V_Reorder2 <- V_LD[orderCFI,orderCFI]
    V_Reorder2b<-diag(z)
    diag(V_Reorder2b)<-diag(V_Reorder2)
    W_CFI<-solve(V_Reorder2b)

  }
  
  ##code to write saturated model to check there are no redundancies
  ##with user provided model in later part of script
  write.test<-function(k, label = "V", label2 = "VF") {
    
    Modelsat<-""
    for (i in 1:(k)) {
      if (k-i >= 1) {
        linemidc <- ""
        for (j in (i+1):k) {
          linemidc <- paste(linemidc, label, i, "~~", label, j, " \n ", sep = "")
        }
      } else {linemidc<-""}
      Modelsat <- paste(Modelsat, linemidc, sep = "")
    }
    
    Modelsat2<-""
    for (i in 1:(k)) {
      if (k-i >= 1) {
        linemidc2 <- ""
        for (j in (i+1):k) {
          linemidc2 <- paste(linemidc2, label, j, "~~", label, i, " \n ", sep = "")
        }
      } else {linemidc2<-""}
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
  
  empty3 <- tryCatch.W.E(
    ReorderModel <- sem(
      Model1,
      sample.cov = S_LD,
      estimator = "DWLS",  #should this not depend on the estimation method?
      WLS.V = W,           #should this?
      sample.nobs = 2,
      warn = FALSE,
      std.lv = std.lv,
      optim.dx.tol = +Inf,
      optim.force.converged = TRUE,
      control = list(iter.max=1)
    )
  )
  
  if ( !( "lavaan" %in% class(empty3$value) ) ) {
    warning( paste(
        "The function has stopped due to convergence issues for your primary model.",
        "Please, contact us with your specific model and variables used or try specifying an",
        "alternative model",
        sep='\n'
    ) )
  }
  
  ##save the ordering
  order <- rearrange(k = k, fit = ReorderModel, names = rownames(S_LD))
  
  ##reorder the weight (inverted V_LD) matrix
  V_Reorder<-V_LD[order,order]

  ##set estimation-specific parameters
  S_LD_nobs = 200  # see __NOBS_OVERRIDE__
  W_Reorder = NULL
  letter.off = 1 # letter offset
  if ( estimation == "DWLS" ) {
    S_LD_nobs = 2  # see __NOBS_OVERRIDE__
    V_Reorderb<-diag(z)
    diag(V_Reorderb)<-diag(V_Reorder)
    W_Reorder<-solve(V_Reorderb)
    letter.off = 2
  }

  cat("Running primary model..\n")
  ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error.
  empty4 <- tryCatch.W.E(
    Model1_Results <- sem(
      Model1,
      sample.cov = S_LD,
      estimator = estimation,
      WLS.V = W_Reorder,
      sample.nobs = S_LD_nobs,
      std.lv = std.lv,
      optim.dx.tol = +Inf
    )
  )

  if ( empty4$warning == '' ) empty4$warning = 0
  
  if ( fix_resid ) {
    if ( class(empty4$value)[1] == "simpleError" || any( grepl("solution has NOT", as.character(empty4$warning)) ) ) {
      cat(
        "The model as initially specified failed to converge. A lower bound of 0 on residual",
        "variances has been automatically added to try and troubleshoot this.",
        sep='\n'
      )
      
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
        } else {Model1b <- ""}
        
        
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
        } else {Model1c <- ""}
        
        #create unique combination of letters for residual variance parameter labels
        n<-combn(letters,4)[,sample(1:14000, k, replace=FALSE)]

        Model2<-""
        for (p in 1:k) {
          linestart2 <- paste(label2, p, " =~ 1*", label, p, sep = "")
          Model2<-paste(Model2, linestart2, " \n ", sep = "")}
        
        Model3<-""
        for (p in 1:k) {
          linestart3a <- paste(label, p, " ~~ ",  paste(n[,p],collapse=""), "*", label, p, sep = "")
          linestart3b <- paste(paste(n[,p],collapse=""), " > .001", sep = "")
          Model3<-paste(Model3, linestart3a, " \n ", linestart3b, " \n ", sep = "")
        }
        
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
      
      empty4 <- tryCatch.W.E(
        Model1_Results <- sem(
          Model1,
          sample.cov = S_LD,
          estimator = estimation,
          WLS.V = W_Reorder,
          sample.nobs = S_LD_nobs,
          std.lv = std.lv,
          optim.dx.tol = +Inf
        )
      )
    }
  }
  
  if ( class(empty4$value)[1] == "simpleError" ) {
    print("The model failed to converge on a solution. Please try specifying an alternative model")
  }
  
  if ( any( grepl("solution has NOT", as.character(empty4$warning)) ) ) {
    print("The model failed to converge on a solution. Please try specifying an alternative model.")
  }
  
  ##save model implied matrix and difference between observed and model implied S_LD matrix
  if (imp_cov) {
    implied<-as.matrix(fitted(Model1_Results))[1]
    implied_order<-colnames(S_LD)
    implied[[1]]<-implied[[1]][implied_order,implied_order]
    f<-S_LD
    implied2<-f-implied[[1]]
  }
  
  ##pull the delta matrix (this doesn't depend on N)
  ##note that while the delta matrix is reordered based on the ordering in the model
  ##specification the lavaan output is also reordered. this actually ensures that the results
  ##match up
  S2.delta <- lavInspect(Model1_Results, "delta")
  
  ##weight matrix from stage 2. S2.W is not reordered by including something like model constraints
  S2.W <- lavInspect(Model1_Results, "WLS.V")
  
  #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that
  #would only be correct if the fit function were correctly specified
  bread2<-tryCatch.W.E(bread <- solve(t(S2.delta)%*%S2.W%*%S2.delta))
  
  if ( !( "matrix" %in% class(bread2$value) ) ) {
    print( "Error: The primary model did not converge!" )
    print( bread2$value )
    cat(
      "Additional warnings or errors are likely being printed by lavaan.",
      "The model output is also printed below (without standard errors) in case this is helpful",
      "for troubleshooting. Please note that these results should not be interpreted.\n",
      sep='\n'
    )
    unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
    unstand<-subset(unstand, unstand$free != 0)
    unstand$free<-NULL
    results<-unstand
    colnames(results)=c("lhs","op","rhs","Unstandardized_Estimate")
    
    ##replace V1-VX general form in output with user provided trait names
    for(i in 1:nrow(results)){
      for(p in 1:length(traits)){
        results$lhs[[i]]<-ifelse(results$lhs[[i]] %in% S_names[[p]], gsub(results$lhs[[i]], traits[[p]], results$lhs[[i]]), results$lhs[[i]])
        results$rhs[[i]]<-ifelse(results$rhs[[i]] %in% S_names[[p]], gsub(results$rhs[[i]], traits[[p]], results$rhs[[i]]), results$rhs[[i]])
      }
    }
    
    print(results)

    return( invisible(
      list( lavaan.out=empty4 )
    ) )
  }
  
  if ( "matrix" %in% class(bread2$value) ) {
    #create the "lettuce" part of the sandwich
    lettuce <- S2.W%*%S2.delta
    
    #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
    Ohtt <- bread %*% t(lettuce)%*%V_Reorder%*%lettuce%*%bread
    
    #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
    SE <- as.matrix(sqrt(diag(Ohtt)))
    
    Model_table <- parTable(Model1_Results)
    
    constraints<-subset(Model_table$label, Model_table$label != "")
    constraints2<-duplicated(constraints)
    
    #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
    if(estimation == "ML"){
      if(":=" %in% Model_table$op){
        print("SEs of ghost parameters are not available for ML estimation")
      }
    }
    if(estimation == "DWLS"){
      if(":=" %in% Model_table$op && !is.na(SE[1])){
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
        
      } else {
        se.ghost<-NA
        if(":=" %in% Model_table$op & is.na(se.ghost[1])){
          se.ghost<-rep("SE could not be computed", count(":=" %in% Model_table$op)$freq)
          ghost<-subset(Model_table, Model_table$op == ":=")[,c(2:4,8,11,14)]
          ghost2<-cbind(ghost,se.ghost)
          colnames(ghost2)[7]<-"SE"
        }
      }
    
      ##check whether the correlations among latent variables are positive definite
      if (r > 1) {
        empty<-tryCatch.W.E(corr.lowtri<-lowerTriangle(lavInspect(Model1_Results,"cor.lv")[1:r,1:r]))
        t<-max(corr.lowtri)
        t2<-min(corr.lowtri)
      } else {
        t<-1
        t2<--1
      }
    
      if(t > 1 || t2 < -1 || t == "NaN"){
        cat(
          "Error: The primary model produced correlations among your latent variables that are",
          "either greater than 1 or less than -1, or to have negative variances. Consequently,",
          "model fit estimates could not be computed and the results should likely not be",
          "interpreted. Results are provided below to enable troubleshooting.",
          "A model constraint that forces the latent correlations to be above -1, below 1, or to",
          "have positive variances is suggested.",
          sep='\n'
        )
        
        unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
        unstand<-subset(unstand, unstand$free != 0)
        unstand$free<-NULL
        results<-unstand
        colnames(results)=c("lhs","op","rhs","Unstandardized_Estimate")
        
        ##replace V1-VX general form in output with user provided trait names
        for(i in 1:nrow(results)){
          for(p in 1:length(traits)){
            results$lhs[[i]]<-ifelse(results$lhs[[i]] %in% S_names[[p]], gsub(results$lhs[[i]], traits[[p]], results$lhs[[i]]), results$lhs[[i]])
            results$rhs[[i]]<-ifelse(results$rhs[[i]] %in% S_names[[p]], gsub(results$rhs[[i]], traits[[p]], results$rhs[[i]]), results$rhs[[i]])
          }
        }
        if(exists("ghost2") == "TRUE"){
          ghost2$free<-NULL
          ghost2$label<-NULL
          unstand2<-rbind(cbind(results,SE),ghost2)
        } else {unstand2<-cbind(results,SE)}
        
        print(unstand2)

        return( invisible(
          list( lavaan.out=empty4 )
        ) )
      }
    }
    ModelQ_table <- parTable(Model1_Results)
    
    ##remove any parameter constraint labels
    ModelQ_table<-subset(ModelQ_table, ModelQ_table$plabel != "")
    
    for (i in 1:length(ModelQ_table)){
      if(((paste(ModelQ_table$lhs[i], ModelQ_table$op[i], ModelQ_table$rhs[i], sep = "") %in% modeltest2$write.test.k)) & ModelQ_table$est[i] == 0)
      {ModelQ_table$free[i] == 1} else {ModelQ_table$free[i] == 0}
    }
    
    for (i in 1:nrow(ModelQ_table)){
      ModelQ_table$free[i]<-ifelse((paste(ModelQ_table$lhs[i], ModelQ_table$op[i], ModelQ_table$rhs[i], sep = "") %in% modeltest2$write.test.k) & ModelQ_table$est[i] == 0, 1, 0)
    }
    
    for (i in 1:nrow(ModelQ_table)){
      ModelQ_table$free[i]<-ifelse((paste(ModelQ_table$lhs[i], ModelQ_table$op[i], ModelQ_table$rhs[i], sep = "") %in% modeltest2$write.test.k) & ModelQ_table$est[i] != 0, 2, ModelQ_table$free[i])
    }
    
    test<-vector(mode="list",length=nrow(ModelQ_table))
    
    for(i in 1:nrow(ModelQ_table)){
      if(ModelQ_table$free[i] == 2) {
        t<-paste(ModelQ_table$lhs[i], ModelQ_table$op[i], ModelQ_table$rhs[i], sep = "")
        t2<-gsub("V", "VF", t)
        test[[i]]<-t2
      }
    }
    test2<-Filter(Negate(is.null), test)
    
    for (i in 1:nrow(ModelQ_table)){
      ModelQ_table$free[i]<-ifelse((paste(ModelQ_table$lhs[i], ModelQ_table$op[i], ModelQ_table$rhs[i], sep = "") %in% test2), 1, ModelQ_table$free[i])
    }
    
    ModelQ_table$free<-ifelse(ModelQ_table$free != 1, 0, ModelQ_table$free)
    
    #want to freely estimate the residual factor variances and the residual covariances
    z<-(k*(k+1))/2
    
    p<-length(ModelQ_table$free)-z
    
    ModelQ_table <- ModelQ_table[order(ModelQ_table$free),]
    
    ModelQ_table$free <- c(rep(0, p),1:z)
    
    ModelQ_table$ustart <- ModelQ_table$est
    ModelQ_table$ustart<-ifelse(ModelQ_table$free > 0, .05, ModelQ_table$ustart)
    
    print("Calculating model chi-square..")

    testQ <- tryCatch.W.E(
      ModelQ_Results_table <- sem(
        model = ModelQ_table,
        sample.cov = S_LD,
        estimator = estimation,
        WLS.V = W_Reorder,
        sample.nobs = S_LD_nobs,
        start = ModelQ_table$ustart,
        std.lv = std.lv,
        optim.dx.tol = +Inf
      )
    )
    
    if ( testQ$warning == '' ) testQ$warning = "Safe"
    testQ$warning <- ifelse( is.na(inspect(ModelQ_Results_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testQ$warning )
    
    if ( any( grepl("lavaan WARNING: model has NOT converged!", as.character(testQ$warning)) ) ) {
      
      ModelQ_table$ustart<-ifelse(ModelQ_table$free > 0, .01, ModelQ_table$ustart)
      testQ2 <- tryCatch.W.E(
        ModelQ_Results_table <- sem(
          model = ModelQ_table,
          sample.cov = S_LD,
          estimator = estimation,
          WLS.V = W_Reorder,
          sample.nobs = S_LD_nobs,
          start = ModelQ_table$ustart, # this was 'start = ModelQ$ustart' for DWLS
          std.lv = std.lv,
          optim.dx.tol = +Inf
        )
      )

    } else {testQ2<-testQ}
    
    if ( testQ2$warning == '' ) testQ2$warning = "Safe"
    testQ2$warning <- ifelse( is.na(inspect(ModelQ_Results_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testQ2$warning )
    
    if ( any( grepl("lavaan WARNING: model has NOT converged!", as.character(testQ2$warning)) ) ) {
      
      ModelQ_table$ustart<-ifelse(ModelQ_table$free > 0, .1, ModelQ_table$ustart)
      testQ3 <- tryCatch.W.E(
        ModelQ_Results_table <- sem(
          model = ModelQ_table,
          sample.cov = S_LD,
          estimator = estimation,
          WLS.V = W_Reorder,
          sample.nobs = S_LD_nobs,
          start = ModelQ_table$ustart, # this was 'start = ModelQ$ustart' for DWLS
          std.lv = std.lv,
          optim.dx.tol = +Inf
        )
      )

    } else {testQ3<-testQ2}
    
    if ( testQ3$warning == '' ) testQ3$warning = "Safe"
    testQ3$warning <- ifelse( is.na(inspect(ModelQ_Results_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testQ3$warning )
    
    if ( !any( grepl("lavaan WARNING: model has NOT converged!", as.character(testQ3$warning)) ) ) {
      
      #pull the delta matrix (this doesn't depend on N)
      S2.delta_Q <- lavInspect(ModelQ_Results_table, "delta")
      
      ##weight matrix from stage 2
      S2.W_Q <- lavInspect(ModelQ_Results_table, "WLS.V")
      
      #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that
      #would only be correct if the fit function were correctly specified
      bread_Q <- solve(t(S2.delta_Q)%*%S2.W_Q%*%S2.delta_Q)
      
      #create the "lettuce" part of the sandwich
      lettuce_Q <- S2.W_Q%*%S2.delta_Q
      
      #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
      Ohtt_Q <- bread_Q %*% t(lettuce_Q)%*%V_Reorder%*%lettuce_Q%*%bread_Q
      
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
      Q<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta
    } else {
      Q<-"The follow-up chi-square model did not converge"
    }
    
    if (CFIcalc == TRUE) {

      print("Calculating CFI..")

      ##set estimation-specific parameters
      S_LD_nobs = 200  # see __NOBS_OVERRIDE__
      W_CFI = NULL
      letter.off = 1 # letter offset
      if ( estimation == "DWLS" ) {
        S_LD_nobs = 2  # see __NOBS_OVERRIDE__
        V_Reorder2b<-diag(z)
        diag(V_Reorder2b)<-diag(V_Reorder2)
        W_CFI<-solve(V_Reorder2b)
        letter.off = 2
      }

      ##now CFI
      ##run independence model
      testCFI <- tryCatch.W.E(
        fitCFI <- sem(
          modelCFI,
          sample.cov = S_LD,
          estimator = estimation,
          WLS.V = W_CFI,
          sample.nobs = S_LD_nobs,
          optim.dx.tol = +Inf
        )
      )

      if ( testCFI$warning == '' ) testCFI$warning = "Safe"
      testCFI$warning <- ifelse( is.na(inspect(fitCFI, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testCFI$warning )
      
      if ( !any( grepl("lavaan WARNING: model has NOT converged!", as.character(testCFI$warning)) ) ) {
        
        ##code to estimate chi-square of independence model#
        #First pull the estimates from Step 2
        ModelQ_table_CFI <- parTable(fitCFI)
        p2<-length(ModelQ_table_CFI$free)-z
        
        ##fix variances and freely estimate covariances
        ModelQ_table_CFI$free <- c(rep(0, p2), 1:z)
        ModelQ_table_CFI$ustart <- ModelQ_table_CFI$est
        
        testCFI2 <- tryCatch.W.E(
          ModelQ_Results_CFI_table <- sem(
            model = ModelQ_table_CFI,
            sample.cov = S_LD,
            estimator = estimation,
            WLS.V = W_CFI,
            sample.nobs = S_LD_nobs,
            optim.dx.tol = +Inf
          )
        )

        if ( testCFI2$warning == '' ) testCFI2$warning = "Safe"
        testCFI2$warning <- ifelse( is.na(inspect(ModelQ_Results_CFI_table, "se")$theta[1,2]), "lavaan WARNING: model has NOT converged!", testCFI2$warning )
        
        if ( !any( grepl("lavaan WARNING: model has NOT converged!", as.character(testCFI2$warning)) ) ) {
          
          #pull the delta matrix (this doesn't depend on N)
          S2.delta_Q_CFI <- lavInspect(ModelQ_Results_CFI_table, "delta")
          
          ##weight matrix from stage 2
          S2.W_Q_CFI <- lavInspect(ModelQ_Results_CFI_table, "WLS.V")
          
          #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates
          #that would only be correct if the fit function were correctly specified
          bread_Q_CFI <- solve(t(S2.delta_Q_CFI)%*%S2.W_Q_CFI%*%S2.delta_Q_CFI)
          
          #create the "lettuce" part of the sandwich
          lettuce_Q_CFI <- S2.W_Q_CFI%*%S2.delta_Q_CFI
          
          #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
          Ohtt_Q_CFI <- bread_Q_CFI %*% t(lettuce_Q_CFI)%*%V_Reorder2%*%lettuce_Q_CFI%*%bread_Q_CFI
          
          ##pull the sampling covariance matrix of the residual covariances and compute diagonal matrix of eigenvalues
          V_etaCFI<- Ohtt_Q_CFI
          Eig2_CFI<-as.matrix(eigen(V_etaCFI)$values)
          Eig_CFI<-diag(z)
          diag(Eig_CFI)<-Eig2_CFI
          
          #Pull P1 (the eigen vectors of V_eta)
          P1_CFI<-eigen(V_etaCFI)$vectors
          
          ##Pull eta = vector of residual covariances
          eta_test_CFI<-parTable(ModelQ_Results_CFI_table)
          eta_test_CFI<-subset(eta_test_CFI, eta_test_CFI$free != 0)
          eta_CFI<-cbind(eta_test_CFI[,14])
          
          #Ronald's magic combining all the pieces from above:
          Q_CFI<-t(eta_CFI)%*%P1_CFI%*%solve(Eig_CFI)%*%t(P1_CFI)%*%eta_CFI
        } else {
          Q_CFI<-"The null (i.e. independence) model did not converge"
        }
      }
      
      ##df of independence Model
      dfCFI <- (((k * (k + 1))/2) - k)
      
      ##df of user model
      df <- (k * (k + 1)/2) - max(parTable(Model1_Results)$free)
      
      if(!(is.character(Q_CFI)) & !(is.character(Q))){
        CFI<-as.numeric(((Q_CFI-dfCFI)-(Q-df))/(Q_CFI-dfCFI))
        CFI<-ifelse(CFI > 1, 1, CFI)
      } else {
        CFI<-"Either the chi-square or null (i.e. independence) model did not converge"
      }
      
    }
    
    print("Calculating standardized results..")

    ##transform the S covariance matrix to S correlation matrix
    D=sqrt(diag(diag(S_LD)))
    S_stand=solve(D)%*%S_LD%*%solve(D)
    rownames(S_stand)<-rownames(S_LD)
    colnames(S_stand)<-colnames(S_stand)
    
    #obtain diagonals of the original V matrix and take their sqrt to get SE's
    Dvcov<-sqrt(diag(V_LD))
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO=as.vector(lowerTriangle((S_stand/S_LD),diag=T))
    
    ## Make sure that if ratio in NaN (division by zero) we put the zero back in: ### TEMP STUPID MICHEL FIX!
    scaleO[is.nan(scaleO)] <- 0
    
    #rescale the SEs by the same multiples that the S matrix was rescaled by
    Dvcovl<-as.vector(Dvcov*t(scaleO))
    
    #obtain the sampling correlation matrix by standardizing the original V matrix
    Vcor<-cov2cor(V_LD)
    
    #rescale the sampling correlation matrix by the appropriate diagonals
    V_stand<-diag(Dvcovl)%*%Vcor%*%diag(Dvcovl)
    
    ##set estimation-specific parameters
    S_LD_nobs = 200  # see __NOBS_OVERRIDE__
    W_stand = NULL
    letter.off = 1 # letter offset
    if ( estimation == "DWLS" ) {
      S_LD_nobs = 2  # see __NOBS_OVERRIDE__
      V_stand2<-diag(z)
      diag(V_stand2)<-diag(V_stand)
      ### make sure no value on the diagonal of V is 0 ### TEMP STUPID MICHEL FIX
      diag(V_stand2)[diag(V_stand2) == 0] <- 2e-9
      W_stand<-solve(V_stand2[order,order])
      letter.off = 2
    }

    sem.fit_stand <- sem(
      Model1,
      sample.cov = S_stand,
      estimator = estimation,
      WLS.V = W_stand,
      sample.nobs = S_LD_nobs,
      std.lv = std.lv,
      optim.dx.tol = +Inf
    )
    
    ##perform same procedures for sandwich correction as in the unstandardized case
    sem.delta_stand <- lavInspect(sem.fit_stand, "delta")
    sem.W_stand <- lavInspect(sem.fit_stand, "WLS.V")
    bread_stand <- solve(t(sem.delta_stand)%*%sem.W_stand %*%sem.delta_stand)
    lettuce_stand <- sem.W_stand%*%sem.delta_stand
    Vcov_stand<-as.matrix(V_stand[order,order])
    Ohtt_stand <- bread_stand %*% t(lettuce_stand)%*%Vcov_stand%*%lettuce_stand%*%bread_stand
    SE_stand <- as.matrix(sqrt(diag(Ohtt_stand)))
    
    unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
    unstand<-subset(unstand, unstand$free != 0)
    unstand$free<-NULL
    
    stand<-data.frame(inspect(sem.fit_stand,"list")[,c(8,14)])
    stand<-subset(stand, stand$free != 0)
    stand$free<-NULL
    
    #df of user model
    df<-lavInspect(Model1_Results, "fit")["df"]
    
    Model_table_stand <- parTable(sem.fit_stand)

    #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
    if(estimation == "DWLS"){
      if(":=" %in% Model_table_stand$op & !(NA %in% Model_table_stand$se)){
        #variance-covariance matrix of parameter estimates, q-by-q (this is the naive one)
        vcov <- lavInspect(sem.fit_stand, "vcov")
        
        #internal lavaan representation of the model
        lavmodel <- sem.fit_stand@Model
        
        #lavaan representation of the indirect effect
        func <- lavmodel@def.function
        
        #vector of parameter estimates
        x <- lav_model_get_parameters(lavmodel, type = "free")
        
        #vector of indirect effect derivatives evaluated @ parameter estimates
        Jac <- lav_func_jacobian_complex(func = func, x = x)
        
        #replace vcov here with our corrected one. this gives parameter variance
        var.ind <- Jac %*% vcov %*% t(Jac)
        
        #square root of parameter variance = parameter SE.
        se.ghost_stand <- sqrt(diag(var.ind))
        
        #pull the ghost parameter point estiamte
        ghost_stand<-subset(Model_table_stand,  Model_table_stand$op == ":=")[,c(2:4,8,11,14)]
        
        ##combine with delta method SE
        ghost2_stand<-cbind(ghost_stand,se.ghost_stand)
        colnames(ghost2_stand)[7]<-"SE_stand"
      } else {
        if(":=" %in% Model_table_stand$op & (NA %in% Model_table_stand$se)){
          se.ghost_stand<-rep("SE could not be computed", count(":=" %in% Model_table_stand$op)$freq)
          ghost_stand<-subset(Model_table_stand, Model_table_stand$op == ":=")[,c(2:4,8,11,14)]
          ghost2_stand<-cbind(ghost_stand,se.ghost_stand)
          colnames(ghost2_stand)[7]<-"SE_stand"
        }
      }
      
      ##combine ghost parameters with rest of output
      if(exists("ghost2") == "TRUE"){
        ghost2$free<-NULL
        ghost2$label<-NULL
        unstand2<-rbind(cbind(unstand,SE),ghost2)
      } else {unstand2<-cbind(unstand,SE)}
      
      ##combine ghost parameters with rest of output
      if(exists("ghost2_stand") == "TRUE"){
        ghost2_stand[,1:5]<-NULL
        stand2<-rbind(cbind(stand,SE_stand),ghost2_stand)
      } else {stand2<-cbind(stand,SE_stand)}
    }  
    if(!(is.character(Q))){
      chisq<-as.numeric(Q)
      AIC<-(Q + 2*lavInspect(Model1_Results, "fit")["npar"])
    } else {
      chisq<-Q
      AIC<-NA
    }
    
    print("Calculating SRMR..")
    
    SRMR<-lavInspect(Model1_Results, "fit")["srmr"]
    
    if(CFIcalc == TRUE){
      modelfit<-cbind(chisq,df,AIC,CFI,SRMR)
    } else {modelfit<-cbind(chisq,df,AIC,SRMR)}
    
    std_all<-standardizedSolution(sem.fit_stand)
    std_all<-subset(std_all, std_all$est.std != "NaN" & std_all$est.std != 0)
    
    results<-cbind(unstand,SE,stand,SE_stand)

    ##add in fixed effects
    base_model<-data.frame(inspect(ReorderModel1, "list")[,c(2:4,8,14)])
    base_model<-subset(base_model, !(paste0(base_model$lhs, base_model$op,base_model$rhs) %in% paste0(unstand$lhs, unstand$op, unstand$rhs)))
    base_model<-subset(base_model, base_model$op == "=~" | base_model$op == "~~" | base_model$op == "~")
    if(nrow(base_model) > 0){
      base_model$free<-NULL
      base_model$SE<-""
      base_model[6]<-base_model$est
      base_model$SE_stand<-""
      #base_model[8]<-base_model$est
      colnames(base_model)<-colnames(results)
      results<-rbind(results,base_model)
    }
    std_all<-subset(std_all,  paste0(std_all$lhs, std_all$op, std_all$rhs) %in% paste0(results$lhs, results$op, results$rhs))
    std_all$order<-paste0(std_all$lhs, std_all$op, std_all$rhs)
    std_all<-data.frame(std_all$est.std,std_all$order)
    colnames(std_all)<-c("est.std","order")
    results$order<-paste0(results$lhs,results$op,results$rhs)
    results<-suppressWarnings(merge(results,std_all,by="order"))
    results$order<-NULL
  }
  
  if ( "matrix" %in% class(bread2$value) ) {
    ##name the columns of the results file
    colnames(results)=c("lhs","op","rhs","Unstand_Est","Unstand_SE","STD_Genotype","STD_Genotype_SE", "STD_All")
    
    ##replace V1-VX general form in output with user provided trait names
    for(i in 1:nrow(results)){
      for(p in 1:length(traits)){
        results$lhs[[i]]<-ifelse(results$lhs[[i]] %in% S_names[[p]], gsub(results$lhs[[i]], traits[[p]], results$lhs[[i]]), results$lhs[[i]])
        results$rhs[[i]]<-ifelse(results$rhs[[i]] %in% S_names[[p]], gsub(results$rhs[[i]], traits[[p]], results$rhs[[i]]), results$rhs[[i]])
      }
    }
    
    ##name model fit columns
    if(CFIcalc == TRUE){
      colnames(modelfit)=c("chisq","df","AIC","CFI","SRMR")} else {colnames(modelfit)=c("chisq","df","AIC","SRMR")}
    
    modelfit<-data.frame(modelfit)
    
    
    if(!(is.character(modelfit$chisq)) & !(is.factor(modelfit$chisq))){
      modelfit$chisq<-as.numeric(as.character(modelfit$chisq))
      modelfit$df<-as.numeric(as.character(modelfit$df))
      modelfit$p_chisq<-ifelse(!(is.character(modelfit$chisq)), modelfit$p_chisq<-pchisq(modelfit$chisq, modelfit$df,lower.tail=FALSE), modelfit$p_chisq<-NA)
      modelfit$chisq<-ifelse(modelfit$df == 0, modelfit$chisq == NA, modelfit$chisq)
      modelfit$AIC<-ifelse(modelfit$df == 0, modelfit$AIC == NA, modelfit$AIC)
      modelfit$p_chisq<-ifelse(modelfit$df == 0, modelfit$p_chisq == NA, modelfit$p_chisq)
      modelfit$SRMR<-ifelse(modelfit$df == 0, modelfit$SRMR == NA, modelfit$SRMR)
      if(CFIcalc == TRUE){
        order<-c(1,2,6,3,4,5)
        modelfit<-modelfit[,order]
        if(!(is.factor(modelfit$CFI))){
          if(modelfit$CFI < 0){
            warning( paste(
              "CFI estimates below 0 should not be trusted, and indicate that the other model fit estimates should be",
              "interpreted with caution. Negative CFI estimates typically appear due to negative residual variances.",
              sep='\n'
            ) )
          }
        }
        modelfit$CFI<-ifelse(modelfit$df == 0, modelfit$CFI == NA, modelfit$CFI)
      } else {
        order<-c(1,2,5,3,4)
        modelfit<-modelfit[,order]
      }
    }
    
    time_all<-proc.time()-time
    print(time_all[3])
    
    if(modelfit$df == 0){
      cat( "Model fit statistics are all printed as NA since you specified a fully saturated model (i.e., df = 0).\n" )
    }
    
    if(LD_sdiff > 0){
      cat( paste0(
        "The S matrix was smoothed prior to model estimation since it was not positive definite. The largest\n",
        "absolute difference in a cell between the smoothed and non-smoothed matrices was ", LD_sdiff, ". As a\n",
        "result of the smoothing, the largest Z-statistic change for the genetic covariances was ", Z_diff, ".\n",
        "We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS.\n",
      ) )
    }
    
    if(LD_sdiff > .025){
      warning( paste(
        "A difference greater than .025 was observed pre- and post-smoothing in the genetic covariance matrix. This",
        "reflects a large difference and results should be interpreted with caution!! Such differences can often",
        "result from including low-powered GWASes. You might consider removing those from the model. If you are going",
        "to run a multivariate GWAS, we strongly recommend setting the smooth_check argument to true to check",
        "smoothing for each SNP.",
        sep='\n'
      ) )
    }
    
    if(Z_diff > .025){
      warning( paste(
        "A difference greater than .025 was observed pre- and post-smoothing in the standardized genetic covariance",
        "matrix. This, reflects a large difference and results should be interpreted with caution!! Such differences",
        "can often result from including low-powered GWASes. You might consider removing those from the model. If you",
        "are going to run a multivariate GWAS, we strongly recommend setting the smooth_check argument to true to",
        "check smoothing for each SNP.",
        sep='\n'
      ) )
    }
    
    if(LD_sdiff2 > 0){
      cat( paste0(
        "The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest\n",
        "absolute difference in a cell between the smoothed and non-smoothed matrices was ", LD_sdiff2, ". As a\n",
        "result of the smoothing, the largest Z-statistic change for the genetic covariances was ", Z_diff,  ".\n",
        "We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS.\n"
      ) )
    }
    
    if(any(constraints2 == TRUE)){
      cat( paste(
        "Please note that when equality constraints are used in the current version of Genomic SEM, the standardized",
        "output will also impose the same constraints.",
        sep='\n'
      ) )
    }

    results$p_value<-2*pnorm(abs(as.numeric(results$Unstand_Est)/as.numeric(results$Unstand_SE)),lower.tail=FALSE)
    results$p_value<-ifelse(results$p_value == 0, "<5e-300", results$p_value)
    
    if (imp_cov) {
      ##replace general form of V1-VX with trait names in model implied and residual covariance matrix
      colnames(implied[[1]])<-traits
      rownames(implied[[1]])<-traits
      colnames(implied2)<-traits
      rownames(implied2)<-traits
      resid_cov<-list()
      resid_cov[[1]]<-implied[[1]]
      resid_cov[[2]]<-implied2
      names(resid_cov) <- c("Model Implied Covariance Matrix", "Residual Covariance Matrix: Calculated as Observed Cov - Model Implied Cov")
      return(list(modelfit=modelfit,results=results,resid_cov=resid_cov,lavaan.out=empty4))
    } else {
      return(list(modelfit=modelfit,results=results,lavaan.out=empty4))
    }
    
  }
  
}

