enrich <-function(s_covstruc, model = "",params,fix= "regressions",std.lv=FALSE,rm_flank=TRUE,tau=FALSE,base=TRUE,toler=NULL,fixparam=NULL){ 
  time<-proc.time()
  ##determine if the model is likely being listed in quotes and print warning if so
  test<-c(str_detect(model, "~"),str_detect(model, "="),str_detect(model, "\\+"))
  if(any(test) != TRUE){
    warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned an error about not locating the ReorderModel.")
  }
  
  if(tau == FALSE){
    ##read in the LD portion of the V (sampling covariance) matrix for the baseline annotation
    V_LD<-as.matrix(s_covstruc$V[[1]])
    
    ##read in the LD portion of the S (covariance) matrix for the baseline annotation
    S_LD<-as.matrix((s_covstruc$S[[1]]))
  }
  
  if(tau == TRUE){
    ##read in the LD portion of the V (sampling covariance) matrix for the baseline annotation
    V_LD<-as.matrix(s_covstruc$V_Tau[[1]])
    
    ##read in the LD portion of the S (covariance) matrix for the baseline annotation
    S_LD<-as.matrix((s_covstruc$S_Tau[[1]]))
    
  }
  
  rownames(S_LD)<-colnames(S_LD)
  base1<-names(s_covstruc[[1]][1])
  
  ##name the columns and rows of the S matrix
  S_names<-colnames(S_LD)
  
  ##name columns of V to remove any variables not used in the current analysis
  y<-expand.grid(S_names,S_names)
  y<-y[!duplicated(apply(y,1,function(x) paste(sort(x),collapse=''))),]
  V_Names<-paste(y$Var1,y$Var2,sep=" ")
  colnames(V_LD)<-V_Names
  rownames(V_LD)<-V_Names
  
  ##determine whether all variables in S are in the model
  ##if not, remove them from S_LD and V_LD for this particular run
  ##also for exact cases
  for(i in 1:length(S_names)){
    S_names[[i]]<-paste0("\\b", S_names[[i]],"\\b",sep="")
  }
  
  w<-1
  remove2<-c()
  for(i in 1:length(S_names)){
    b<-grepl(S_names[i], model)
    if(b == FALSE){
      remove<-paste0("\\b", colnames(S_LD)[i],"\\b",sep="")
      remove2[w]<-i
      V_LD <- V_LD[-grep(pattern=remove[1],row.names(V_LD)),-grep(pattern=remove[1],colnames(V_LD))]
      w<-w+1
    }else{}
  }
  
  if(is.null(remove2) == FALSE){
    S_LD<-S_LD[-remove2,-remove2]
  }
  
  Model1<-model
  
  print(paste0(base1, " is assumed to be the baseline annotation that includes all SNPs."))
  
  ##k = number of phenotypes in dataset (i.e., number of columns in LD portion of S matrix)
  k<-ncol(S_LD)
  
  ##size of V matrix used later in code to create diagonal V matrix
  z<-(k*(k+1))/2
  
  ##smooth to near positive definite if either V or S are non-positive definite
  ks<-nrow(S_LD)
  S_LDb<-S_LD
  smooth1<-ifelse(eigen(S_LD)$values[ks] <= 0, S_LD<-as.matrix((nearPD(S_LD, corr = FALSE))$mat), S_LD<-S_LD)
  diff<-(abs(S_LD-S_LDb))
  LD_sdiff<-max(diff)
  
  kv<-nrow(V_LD)
  V_LDb<-V_LD
  smooth2<-ifelse(eigen(V_LD)$values[kv] <= 0, V_LD<-as.matrix((nearPD(V_LD, corr = FALSE))$mat), V_LD<-V_LD)
  diff2<-(abs(V_LD-V_LDb))
  LD_sdiff2<-max(diff2)
  
  SE_pre<-matrix(0, k, k)
  SE_pre[lower.tri(SE_pre,diag=TRUE)] <-sqrt(diag(V_LDb))
  
  SE_post<-matrix(0, k, k)
  SE_post[lower.tri(SE_post,diag=TRUE)] <-sqrt(diag(V_LD))
  
  Z_pre<-S_LDb/SE_pre
  Z_post<-S_LD/SE_post
  Z_diff<-max(abs(Z_pre-Z_post),na.rm=T)
  
  ##run model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  #transform V_LD matrix into a weight matrix: 
  W <- solve(V_LD)
  
  ##run the model
  if(std.lv == FALSE){
    empty2<-.tryCatch.W.E(ReorderModel1 <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W, sample.nobs = 2,warn=FALSE, optim.dx.tol = +Inf,optim.force.converged=TRUE))
  }
  
  if(std.lv == TRUE){
    empty2<-.tryCatch.W.E(ReorderModel1 <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W, sample.nobs = 2,warn=FALSE,std.lv=TRUE, optim.dx.tol = +Inf,optim.force.converged=TRUE))
  }
  
  if(class(empty2$value)[1] != "lavaan"){
    latentcorr<-grepl("not defined:", empty2$value$message[1][1])
    latentcorr2<-grepl("unknown", empty2$value$message[1][1])
    if(latentcorr == TRUE | latentcorr2 == TRUE){
      warning(paste("The function may have stopped either because a variable has been misnamed in the model or because you have tried to estimate a correlation between an observed and latent variable. In the latter case, one workaround
                    is to define a latent variable solely by the observed variable."))
    }}
  
  ##save the ordering
  order <- .rearrange(k = k, fit = ReorderModel1, names = rownames(S_LD))
  
  ##reorder the weight (inverted V_LD) matrix
  V_Reorder<-V_LD[order,order]
  V_Reorderb<-diag(z)
  diag(V_Reorderb)<-diag(V_Reorder)
  W_Reorder<-solve(V_Reorderb)
  
  print("Running model for baseline annotation")
  check<-1
  ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
  if(std.lv == FALSE){
    empty4<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2,optim.dx.tol = +Inf))
  }
  
  if(std.lv == TRUE){
    empty4<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2,std.lv=TRUE, optim.dx.tol = +Inf))
  }
  
  empty4$warning$message[1]<-ifelse(is.null(empty4$warning$message), empty4$warning$message[1]<-0, empty4$warning$message[1])
  
  if(class(empty4$value)[1] == "simpleError" | grepl("solution has NOT",  as.character(empty4$warning)) == TRUE){
    print("The model as initially specified failed to converge for the baseline annotation. A lower bound of 0 on residual variances has been automatically added to try and troubleshoot this.")
    
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
            if ((k-1)-i > 0 | k ==2) {
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
        linestart3a <- paste(label, p, " ~~ ", letters[p], letters[p],letters[p+2], "*", label, p, sep = "")
        linestart3b <- paste(letters[p], letters[p],letters[p+2], " > .001", sep = "")
        Model3<-paste(Model3, linestart3a, " \n ", linestart3b, " \n ", sep = "")}
      
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
    
    if(std.lv == FALSE){
      empty4<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2, optim.dx.tol = +Inf))
    }
    
    if(std.lv == TRUE){
      empty4<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2,std.lv=TRUE, optim.dx.tol = +Inf))
    }
  }
  
  if(class(empty4$value)[1] == "simpleError"){
    print("The model failed to converge on a solution for the baseline annotation. Please try specifying an alternative model")
    check<-2}
  
  if(!(is.null(empty4$warning))){
    if(grepl("solution has NOT",  as.character(empty4$warning)) == TRUE){
      print("The model failed to converge on a solution for the baseline annotation. Please try specifying an alternative model.")
      check<-2
    }}
  
  
  if(check == 1){
    
    if(base==TRUE){
      
      base_results<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
      base_results<-subset(base_results, base_results$free != 0)                    
      base_results$free<-NULL
      
      ##fixed effects
      base_model<-data.frame(inspect(ReorderModel1, "list")[,c(2:4,8,14)])
      base_model<-subset(base_model,  !(paste0(base_model$lhs, base_model$op,base_model$rhs) %in% paste0(base_results$lhs, base_results$op, base_results$rhs)))
      base_model<-subset(base_model, base_model$op == "=~" | base_model$op == "~~" | base_model$op == "~")
      
      #calculate SEs
      S2.delt <- lavInspect(Model1_Results, "delta")
      
      ##weight matrix from stage 2. S2.W is not reordered by including something like model constraints
      S2.W <- lavInspect(Model1_Results, "WLS.V") 
      
      #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
      if(is.null(toler)){
        bread_check<-.tryCatch.W.E(bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt))
      }else{ bread_check<-.tryCatch.W.E(bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler))}
      
      if(class(bread_check$value)[1] != "matrix"){
        warning("SEs could not be computed for baseline model. Results may not be trustworthy")
      }else{
        
        lettuce <- S2.W%*%S2.delt
        
        #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
        Ohtt <- bread %*% t(lettuce)%*%V_Reorder%*%lettuce%*%bread  
        
        #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
        SE <- as.matrix(sqrt(diag(Ohtt)))
        
        base_results<-cbind(base_results,SE)
        
        if(nrow(base_model) > 0){
          base_model$SE<-""
        }
      }
      
      if(nrow(base_model) > 0){
        base_model$free<-NULL
        colnames(base_model)<-colnames(base_results)
        base_results<-rbind(base_results,base_model)
      }
    }
    
    ModelQ_WLS <- parTable(Model1_Results)
    
    #fix all values from baseline model
    ModelQ_WLS$free<-0
    
    #remove white space from parameters for easier matching
    params<-str_replace_all(params, fixed(" "), "") 
    
    x<-1:nrow(ModelQ_WLS)
    u<-1
    ##freely estimate specified parameters, (residual) variances, and correlations
    #fix regression estimates including factor loadings
    #likelymost useful for estimating enrichment of factor variances
    if(fix == "regressions"){
      for(i in 1:nrow(ModelQ_WLS)){
        if((paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% params | ModelQ_WLS$op[i] == "~~") & !(paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% fixparam)) {
          ModelQ_WLS$free[i]<-x[u]
          u<-u+1
        }else{ModelQ_WLS$free[i]<-0}
      }
    }

    ##freely estimate specified parameters, regressions (including factor loadings), and (residual) variances
    #fix covariances
    if(fix == "covariances"){
      for(i in 1:nrow(ModelQ_WLS)){
        if((paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% params | ModelQ_WLS$op[i] == "~" |  ModelQ_WLS$op[i] == "=~" | ModelQ_WLS$lhs[i] == ModelQ_WLS$rhs[i]) & !(paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% fixparam)) {
          ModelQ_WLS$free[i]<-x[u]
          u<-u+1
        }else{ModelQ_WLS$free[i]<-0}
      }
    }
    
    ##freely estimate specified parameters, regressions (including factor loadings), and covariances
    #fix (residual) variances
    if(fix == "variances"){
      for(i in 1:nrow(ModelQ_WLS)){
        if((paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% params | ModelQ_WLS$op[i] == "~" |  ModelQ_WLS$op[i] == "=~" | (ModelQ_WLS$op[i] == "~~" & (ModelQ_WLS$lhs[i] != ModelQ_WLS$rhs[i])))  & !(paste(ModelQ_WLS$lhs[i], ModelQ_WLS$op[i], ModelQ_WLS$rhs[i], sep = "") %in% fixparam)) {
          ModelQ_WLS$free[i]<-x[u]
          u<-u+1
        }else{ModelQ_WLS$free[i]<-0}
      }
    }
    
    if(max(ModelQ_WLS$free) == 0){
      stop("All parameters are being fixed from the baseline model. Please specify arguments so at least one free parameter is estimated.")
    }
    
    if(max(ModelQ_WLS$free) == nrow(ModelQ_WLS)){
      warning("All parameters are being freely estimated from the baseline model. Enrichment results should likely not be interpreted.")
    }
    
    if(base==TRUE){
      Merge_base<-data.frame(paste0(ModelQ_WLS$lhs,ModelQ_WLS$op,ModelQ_WLS$rhs,sep=""),ModelQ_WLS$free)
      colnames(Merge_base)<-c("Merge","Fixed_Enrich")
      base_results$Merge<-paste0(base_results$lhs,base_results$op,base_results$rhs,sep="")
      base_results<-merge(base_results,Merge_base,by="Merge")
      base_results$Fixed_Enrich<-ifelse(base_results$Fixed_Enrich == 0, "Yes", "No")
      base_results<-base_results[order(base_results$Fixed_Enrich),]
      base_results$Merge<-NULL
      rm(Merge_base)
    }
    
    ModelQ_WLS$ustart <- ModelQ_WLS$est
    ModelQ_WLS$ustart<-ifelse(ModelQ_WLS$free > 0, 1, ModelQ_WLS$ustart)
    ModelQ_WLS<-subset(ModelQ_WLS,ModelQ_WLS$op != "==")
    
    print("Confirming fixed model reproduces estimate from freely estimated model for baseline annotation.")
    
    if(std.lv == FALSE){
      testQ<-.tryCatch.W.E(ModelQ_Results_WLS <- sem(model = ModelQ_WLS, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs=2, optim.dx.tol = +Inf))
    }
    
    if(std.lv == TRUE){
      testQ<-.tryCatch.W.E(ModelQ_Results_WLS <- sem(model = ModelQ_WLS, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs=2, std.lv=TRUE, optim.dx.tol = +Inf))
    }
    
    test1<-subset(ModelQ_WLS, ModelQ_WLS$free != 0)
    
    ModelQ_WLS2 <- parTable(ModelQ_Results_WLS)
    test2<-subset(ModelQ_WLS2,ModelQ_WLS2$free != 0)
    
    if(max(round(test1$est-test2$est)) != 0){
      warning("The model re-estimated in the baseline annotation is producing different values when fixed and freely estimated. This may be an indication that the model is poorly specified.")
    }
    
    #subset to only parameters being estimated for enrichment for scaling in later code
    test1<-subset(test1, paste(test1$lhs, test1$op, test1$rhs, sep = "") %in% params)  
    
    print(paste0("Beginning estimation of enrichment for ", length(s_covstruc$S), " functional annotations."))
    
    #add in baseline annotation for selection variable that excludes flanking and continuous annots
    Select<-data.frame(c("Base",s_covstruc$Select$V1),c(1,s_covstruc$Select$V2))
    colnames(Select)<-c("Annot","Select")
    
    for(n in 1:length(s_covstruc$S)){
      
      #only run for binary and non-flanking annotations
      if(Select$Select[n] == 1){
        
        if(tau == FALSE){
          ##read in the LD portion of the V (sampling covariance) matrix for the baseline annotation
          V_LD<-as.matrix(s_covstruc$V[[n]])
          
          ##read in the LD portion of the S (covariance) matrix for the baseline annotation
          S_LD<-as.matrix((s_covstruc$S[[n]]))
        }
        
        if(tau == TRUE){
          ##read in the LD portion of the V (sampling covariance) matrix for the baseline annotation
          V_LD<-as.matrix(s_covstruc$V_Tau[[n]])
          
          ##read in the LD portion of the S (covariance) matrix for the baseline annotation
          S_LD<-as.matrix((s_covstruc$S_Tau[[n]]))
          
        }
        
        ##remove variables not used in the model from S and V
        colnames(V_LD)<-V_Names
        rownames(V_LD)<-V_Names
        
        w<-1
        remove2<-c()
        for(i in 1:length(S_names)){
          b<-grepl(S_names[i], model)
          if(b == FALSE){
            remove<-paste0("\\b", colnames(S_LD)[i],"\\b",sep="")
            remove2[w]<-i
            V_LD <- V_LD[-grep(pattern=remove[1],row.names(V_LD)),-grep(pattern=remove[1],colnames(V_LD))]
            w<-w+1
          }else{}
        }
        
        if(is.null(remove2) == FALSE){
          S_LD<-S_LD[-remove2,-remove2]
        }
        
        
        #check smoothing for S and V
        if(all(diag(S_LD) < 0) == FALSE) {
          S_LDb<-S_LD
          
          smooth1<-ifelse(eigen(S_LD)$values[ks] <= 0, S_LD<-as.matrix((nearPD(S_LD, corr = FALSE))$mat), S_LD<-S_LD)
          diff<-(abs(S_LD-S_LDb))
          LD_sdiff<-max(diff)
          
          V_LDb<-V_LD
          smooth2<-ifelse(eigen(V_LD)$values[kv] <= 0, V_LD<-as.matrix((nearPD(V_LD, corr = FALSE))$mat), V_LD<-V_LD)
          diff2<-(abs(V_LD-V_LDb))
          LD_sdiff2<-max(diff2)
          
          SE_pre<-matrix(0, k, k)
          SE_pre[lower.tri(SE_pre,diag=TRUE)] <-sqrt(diag(V_LDb))
          
          SE_post<-matrix(0, k, k)
          SE_post[lower.tri(SE_post,diag=TRUE)] <-sqrt(diag(V_LD))
          
          Z_pre<-S_LDb/SE_pre
          Z_post<-S_LD/SE_post
          Z_diff<-max(abs(Z_pre[lower.tri(Z_pre,diag=TRUE)]-Z_post[lower.tri(Z_post,diag=TRUE)]),na.rm=T)
          
          V_Reorder<-V_LD[order,order]
          V_Reorderb<-diag(z)
          diag(V_Reorderb)<-diag(V_Reorder)
          W_Reorder<-solve(V_Reorderb)
          
          part_warn<-.tryCatch.W.E(ModelPart_Results <- sem(ModelQ_WLS, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2, optim.dx.tol = +Inf))
          
          part_warn$warning$message[1]<-ifelse(is.null(part_warn$warning$message), part_warn$warning$message[1]<-0, part_warn$warning$message[1])
          
          if(class(part_warn$value)[1] != "simpleError" & grepl("solution has NOT",  as.character(part_warn$warning)) != TRUE){
            ##that the lavaan output is also reordered so that this actually ensures that the results match up 
            S2.delt <- lavInspect(ModelPart_Results, "delta")
            
            ##weight matrix from stage 2. S2.W is not reordered by including something like model constraints
            S2.W <- lavInspect(ModelPart_Results, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
            if(is.null(toler)){
              bread_check<-.tryCatch.W.E(bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt))
            }else{bread_check<-.tryCatch.W.E(bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler))}
            
            if(class(bread_check$value)[1] == "matrix"){
              
              lettuce <- S2.W%*%S2.delt
              
              #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
              Ohtt <- bread %*% t(lettuce)%*%V_Reorder%*%lettuce%*%bread  
              
              #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
              SE <- as.matrix(sqrt(diag(Ohtt)))
              
              ##replace any labels with the actual parameter name
              for(e in 1:nrow(SE)){
                for(w in 1:nrow(ModelQ_WLS)){
                  if(rownames(SE)[e] == ModelQ_WLS$label[w]){
                    rownames(SE)[e]<-paste(ModelQ_WLS$lhs[w], ModelQ_WLS$op[w], ModelQ_WLS$rhs[w], sep = "")
                  }
                }
              }
              
              unstand<-data.frame(inspect(ModelPart_Results, "list")[,c(2:4,8,14)])
              unstand<-subset(unstand, unstand$free != 0)                    
              unstand$free<-NULL
              unstand<-subset(unstand, paste(unstand$lhs, unstand$op, unstand$rhs, sep = "") %in% params)
              SE<-subset(SE, rownames(SE) %in% params)
              
              results<-cbind(as.character(names(s_covstruc$S[n])),unstand,SE,LD_sdiff,Z_diff)
              results[,1]<-as.character(results[,1])
              for(y in 1:nrow(results)){
                #calculate enrichment with null = 1 divided by proportional size of annotation
                results$enrichment[y]<-(results$est[y]/test1$est[y])/s_covstruc$Prop$Prop[n]
                
                #caculate enrichment SE
                results$enrichment_se[y]<-(results$SE[y]/abs(test1$est[y]))/s_covstruc$Prop$Prop[n]
                
              }
              
              #compute 1-trailed p-value subtracting null of 1 from enrichment estimate
              results$enrichment_p<-pnorm(((results$enrichment-1)/results$enrichment_se),lower.tail=FALSE)
              
              results$est<-NULL
              results$SE<-NULL
              results$error<-ifelse(class(part_warn$value) == "lavaan", 0, as.character(part_warn$value$message))[1]
              results$warning<-ifelse(class(part_warn$warning) == 'NULL', 0, as.character(part_warn$warning$message))[1]
              
              if(n == 1){
                Results_List<-vector(mode="list",length=nrow(results))
                for(y in 1:nrow(results)){
                  Results_List[[y]]<-as.data.frame(matrix(NA,ncol=ncol(results),nrow=length(s_covstruc$S)))
                  colnames(Results_List[[y]])<-c("Annotation", "lhs", "op", "rhs", "Cov_Smooth", "Z_smooth", "Enrichment", "Enrichment_SE", "Enrichment_p_value", "Error", "Warning")
                  Results_List[[y]][1,]<-results[y,]
                }
              }else{
                for(y in 1:nrow(results)){
                  Results_List[[y]][n,]<-results[y,]
                }
              }
            }else{
              for(y in 1:length(params)){
                final<-data.frame(as.character(names(s_covstruc$S[n])), test1$lhs[y], test1$op[y], test1$rhs[y],LD_sdiff, Z_diff, NA,NA,NA)
                final$error<-ifelse(class(part_warn$value) == "lavaan", 0, as.character(part_warn$value$message))[1]
                final$warning<-ifelse(class(part_warn$warning) == 'NULL', 0, as.character(part_warn$warning$message))[1]
                colnames(final)<-c("Annotation", "lhs", "op", "rhs", "Cov_Smooth", "Z_smooth", "Enrichment", "Enrichment_SE", "Enrichment_p_value", "Error", "Warning")
                Results_List[[y]][n,]<-final
              }
            }
          }else{
            for(y in 1:length(params)){
              final<-data.frame(as.character(names(s_covstruc$S[n])), test1$lhs[y], test1$op[y], test1$rhs[y],LD_sdiff, Z_diff, NA,NA,NA)
              final$error<-ifelse(class(part_warn$value) == "lavaan", 0, as.character(part_warn$value$message))[1]
              final$warning<-ifelse(class(part_warn$warning) == 'NULL', 0, as.character(part_warn$warning$message))[1]
              colnames(final)<-c("Annotation", "lhs", "op", "rhs", "Cov_Smooth", "Z_smooth", "Enrichment", "Enrichment_SE", "Enrichment_p_value", "Error", "Warning")
              Results_List[[y]][n,]<-final
            }
          }
        }else{
          for(y in 1:length(params)){
            final<-data.frame(as.character(names(s_covstruc$S[n])),  test1$lhs[y], test1$op[y], test1$rhs[y],NA,NA,NA,NA,NA)
            final$error<-0
            final$warning<-c("This annotation was not analyzed as all heritability estimates were below 0.")
            Results_List[[y]][n,]<-final
          }
        }
      }else{
        for(y in 1:length(params)){
          final<-data.frame(as.character(names(s_covstruc$S[n])),  test1$lhs[y], test1$op[y], test1$rhs[y],NA,NA,NA,NA,NA)
          final$error<-0
          final$warning<-c("This annotation was not analyzed as it is either a continuous or flanking annotation.")
          Results_List[[y]][n,]<-final
        }
      }
    }
    
    if(rm_flank==TRUE){
      flank1<-nrow(Results_List[[1]])
      for(i in 1:length(Results_List)){
        Results_List[[i]]<-subset(Results_List[[i]], Results_List[[i]]$Warning != "This annotation was not analyzed as it is either a continuous or flanking annotation.")
      }
      flank2<-flank1-nrow(Results_List[[1]])
      print(paste0(flank2, " annotations were removed from the output because they were either continuous or flanking window annotatoins."))
    }
    
    if(base == TRUE){
      Results_List[[((length(Results_List)+1))]]<-base_results
    }
    
    #turn Results_List into single dataframe if only analyzing single model parameter
    if(length(Results_List) == 1 & base == FALSE){
      Results_List<-data.frame(Results_List[[y]])
    }
    
    
    return(Results_List)
  }
  
  time_all<-proc.time()-time
  print(time_all[3])
  
}
