
userGWASpar<-function(Output,estimation="DWLS",model="",modelchi=FALSE,printwarn=TRUE,sub=FALSE,cores=NULL){ 
  time<-proc.time()
  
  ##determine if the model is likely being listed in quotes and print warning if so
  test<-c(str_detect(model, "~"),str_detect(model, "="),str_detect(model, "\\+"))
  if(all(test) == FALSE){
    warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned stopped running after returning an error.")
  }
  
  if(is.null(cores)){
    ##if no default provided use 1 less than the total number of cores available so your computer will still function
    int <- detectCores() - 1
  }else{int<-cores}
  
  registerDoParallel(int)
  
  ##specify the cores should have access to the local environment
  makeCluster(int, type="FORK")
  
  ##split the V and S matrices into as many (cores - 1) as are aviailable on the local computer
  V_Full<-suppressWarnings(split(Output[[1]],1:int))
  S_Full<-suppressWarnings(split(Output[[2]],1:int))
  SNPs<-suppressWarnings(split(Output[[3]],1:int))
  
  #enter in k for number of columns in S matrix
  k<-ncol(S_Full[[1]][[1]])
  
  ##number of models to run = number of distinct S/V matrices
  f<-length(Output[[1]])
  
  ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
  Output$RS$SNP<-as.character(Output$RS$SNP)
  Output$RS$A1<-as.character(Output$RS$A1)
  Output$RS$A2<-as.character(Output$RS$A2)
  
  #function to rearrange the sampling covariance matrix from original order to lavaan's order: 
  #'k' is the number of variables in the model
  #'fit' is the fit function of the regression model
  #'names' is a vector of variable names in the order you used
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
  
  ##modification of trycatch that allows the results of a failed run to still be saved
  tryCatch.W.E <- function(expr)
  {
    W <- NULL
    w.handler <- function(w){ # warning handler
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler),
         warning = W)
  }
  
  
  if(modelchi == TRUE){
    #function to creat row/column names for S_LD matrix
    write.names <- function(k, label = "V") {  
      varnames<-vector(mode="character",length=k)
      
      for (i in 1:k){
        varnames[i]<-paste(label,i,sep="")}
      
      return(varnames)
    }
    
    W_test <- solve(V_Full[[1]][[1]])
    
    S_Fulltest<-S_Full[[1]][[1]]
    
    ##run the model to determine number of latent variables
    suppress<-tryCatch.W.E(ReorderModel1 <- sem(model, sample.cov = S_Fulltest, estimator = "DWLS", WLS.V = W_test, sample.nobs = 2,warn=FALSE, optim.dx.tol = +Inf)) 
    
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
    
    while(class(tryCatch.W.E(lavParseModelString(Model1))$value$message) != 'NULL'){
      u<-tryCatch.W.E(lavParseModelString(Model1))$value$message
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
    
    #transform sampling covariance matrix into a weight matrix: 
    W <- solve(V_Full[[1]][[i]])
    
    if(modelchi == TRUE){
      ##name the columns and rows of the S matrix in general format V1-VX
      rownames(S_Full[[1]][[i]]) <- S_names
      colnames(S_Full[[1]][[i]]) <- S_names
    }
    
    S_Fullrun<-S_Full[[1]][[i]]
    
    test2<-tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf))
    
    order <- rearrange(k = k, fit = ReorderModel, names = rownames(S_Full[[1]][[i]]))
  }
  
  
  ##estimation for WLS
  if(estimation=="DWLS"){
    
    results<-foreach(n = icount(int), .combine = 'rbind') %:% 
      
      foreach (i=1:length(V_Full[[n]]), .combine='rbind', .packages = "lavaan") %dopar% { 
        
        #reorder sampling covariance matrix based on what lavaan expects given the specified model
        V_Full_Reorder <- V_Full[[n]][[i]][order,order]
        u<-nrow(V_Full_Reorder)
        V_Full_Reorderb<-diag(u)
        diag(V_Full_Reorderb)<-diag(V_Full_Reorder)
        
        ##invert the reordered sampling covariance matrix to create a weight matrix 
        W <- solve(V_Full_Reorderb) 
        
        #import the S_Full matrix for appropriate run
        S_Fullrun<-S_Full[[n]][[i]]
        
        if(modelchi == TRUE){
          ##name the columns and rows of the S matrix in general format V1-VX
          rownames(S_Fullrun) <- S_names
          colnames(S_Fullrun) <- S_names
        }
        
        ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
        test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf))
        
        Model_WLS <- parTable(Model1_Results)
        
        if(NA %in% Model_WLS$se){
          SE<-rep("SE could not be computed", max(Model_WLS$free))}else{
            #pull the delta matrix (this doesn't depend on N)
            S2.delt <- lavInspect(Model1_Results, "delta")
            
            ##weight matrix from stage 2
            S2.W <- lavInspect(Model1_Results, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
            bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
            
            #create the "lettuce" part of the sandwich
            lettuce <- S2.W%*%S2.delt
            
            #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
            Ohtt <- bread %*% t(lettuce)%*%V_Full_Reorder%*%lettuce%*%bread  
            
            #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
            SE <- as.matrix(sqrt(diag(Ohtt)))}
        
        #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
        if(":=" %in% Model_WLS$op & !(NA %in% Model_WLS$se)){
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
          ghost<-subset(Model_WLS, Model_WLS$op == ":=")[,c(2:4,8,11,14)]
          
          ##combine with delta method SE
          ghost2<-cbind(ghost,se.ghost)
          colnames(ghost2)[7]<-"SE"
        }else{
          if(":=" %in% Model_WLS$op & (NA %in% Model_WLS$se)){
            se.ghost<-rep("SE could not be computed", count(":=" %in% Model_WLS$op)$freq)
            ghost<-subset(Model_WLS, Model_WLS$op == ":=")[,c(2:4,8,11,14)]
            ghost2<-cbind(ghost,se.ghost)
            colnames(ghost2)[7]<-"SE"}else{}} 
        
        if(modelchi == TRUE){
          
          ModelQ_WLS <- parTable(Model1_Results)
          
          ##remove any parameter constraint labels or ghost parameters
          ModelQ_WLS<-subset(ModelQ_WLS, ModelQ_WLS$plabel != "")
          
          ##identify components of saturated model not already estimated and label as 1
          for (g in 1:nrow(ModelQ_WLS)){
            ModelQ_WLS$free[g]<-ifelse((paste(ModelQ_WLS$lhs[g], ModelQ_WLS$op[g], ModelQ_WLS$rhs[g], sep = "") %in% modeltest2$write.test.k) & ModelQ_WLS$est[g] == 0, 1, 0)
          }
          
          ##identify components of saturated model that were already estimated in user model and label as 2
          for (g in 1:nrow(ModelQ_WLS)){
            ModelQ_WLS$free[g]<-ifelse((paste(ModelQ_WLS$lhs[g], ModelQ_WLS$op[g], ModelQ_WLS$rhs[g], sep = "") %in% modeltest2$write.test.k) & ModelQ_WLS$est[g] != 0, 2, ModelQ_WLS$free[g])
          }
          
          tester<-vector(mode="list",length=nrow(ModelQ_WLS))
          
          ##replace components of saturdated model already estimated with residaul factor estimates
          for(g in 1:nrow(ModelQ_WLS)){
            if(ModelQ_WLS$free[g] == 2) { 
              t<-paste(ModelQ_WLS$lhs[g], ModelQ_WLS$op[g], ModelQ_WLS$rhs[g], sep = "")
              t2<-gsub("V", "VF", t)
              tester[[g]]<-t2}else{}
          }
          
          test2<-Filter(Negate(is.null), tester)
          
          for (g in 1:nrow(ModelQ_WLS)){
            ModelQ_WLS$free[g]<-ifelse((paste(ModelQ_WLS$lhs[g], ModelQ_WLS$op[g], ModelQ_WLS$rhs[g], sep = "") %in% test2), 1, ModelQ_WLS$free[g])
          }
          
          ModelQ_WLS$free<-ifelse(ModelQ_WLS$free != 1, 0, ModelQ_WLS$free)
          
          #want to freely estimate the residual factor variances and the residual covariances
          p<-length(ModelQ_WLS$free)-z
          
          ModelQ_WLS <- ModelQ_WLS[order(ModelQ_WLS$free),] 
          
          ModelQ_WLS$free <- c(rep(0, p),1:z)
          
          ModelQ_WLS$ustart <- ModelQ_WLS$est
          ModelQ_WLS$ustart<-ifelse(ModelQ_WLS$free > 0, .05, ModelQ_WLS$ustart)
          
          testQ<-tryCatch.W.E(ModelQ_Results_WLS <- sem(model = ModelQ_WLS, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs=2, start = ModelQ_WLS$ustart, optim.dx.tol = +Inf)) 
          testQ$warning$message[1]<-ifelse(is.null(testQ$warning$message), testQ$warning$message[1]<-"Safe", testQ$warning$message[1])
          testQ$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_WLS, "se")$theta[1,2]) == TRUE, testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ$warning$message[1])
          
          if(as.character(testQ$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
            
            ModelQ_WLS$ustart<-ifelse(ModelQ_WLS$free > 0, .01, ModelQ_WLS$ustart)
            
            testQ2<-tryCatch.W.E(ModelQ_Results_WLS <- sem(model = ModelQ_WLS, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2, start=ModelQ$ustart, optim.dx.tol = +Inf)) 
          }else{testQ2<-testQ}
          
          testQ2$warning$message[1]<-ifelse(is.null(testQ2$warning$message), testQ2$warning$message[1]<-"Safe", testQ2$warning$message[1])
          testQ2$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_WLS, "se")$theta[1,2]) == TRUE, testQ2$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ2$warning$message[1])
          
          if(as.character(testQ2$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
            
            ModelQ_WLS$ustart<-ifelse(ModelQ_WLS$free > 0, .1, ModelQ_WLS$ustart)
            
            testQ3<-tryCatch.W.E(ModelQ_Results_WLS <- sem(model = ModelQ_WLS, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, sample.nobs = 2, start=ModelQ$ustart, optim.dx.tol = +Inf)) 
          }else{testQ3<-testQ2}
          
          testQ3$warning$message[1]<-ifelse(is.null(testQ3$warning$message), testQ3$warning$message[1]<-"Safe", testQ3$warning$message[1])
          testQ3$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_WLS, "se")$theta[1,2]) == TRUE, testQ3$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ3$warning$message[1])
          
          if(as.character(testQ3$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
            
            #pull the delta matrix (this doesn't depend on N)
            S2.delt_Q <- lavInspect(ModelQ_Results_WLS, "delta")
            
            ##weight matrix from stage 2
            S2.W_Q <- lavInspect(ModelQ_Results_WLS, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
            bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q) 
            
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
            eta_test<-parTable(ModelQ_Results_WLS)
            eta_test<-subset(eta_test, eta_test$free != 0)
            eta<-cbind(eta_test[,14])
            
            #Ronald's magic combining all the pieces from above:
            Q_WLS<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta}else{Q_WLS<-NA}
        }
        
        ##remove parameter constraints, ghost parameters, and fixed effects from output to merge with SEs
        unstand<-subset(Model_WLS, Model_WLS$plabel != "" & Model_WLS$free > 0)[,c(2:4,8,11,14)]
        
        ##combine ghost parameters with rest of output
        if(exists("ghost2") == "TRUE"){
          unstand2<-rbind(cbind(unstand,SE),ghost2)
        }else{unstand2<-cbind(unstand,SE)}
        
        ##add in fixed effects and parameter constraints to output
        other<-subset(Model_WLS, (Model_WLS$plabel == "" & Model_WLS$op != ":=") | (Model_WLS$free == 0 & Model_WLS$plabel != ""))[,c(2:4,8,11,14)]
        other$SE<-rep(NA, nrow(other))
        
        ##combine fixed effects and parameter constraints with output if there are any
        if(nrow(other) > 0){
          final<-rbind(unstand2,other)
        }else{final<-unstand2}
        
        #reorder based on row numbers so it is in order the user provided
        final$index <- as.numeric(row.names(final))
        final<-final[order(final$index), ]
        final$index<-NULL
        
        final$Z_Estimate<-final$est/final$SE
        final$Pval_Estimate<-2*pnorm(abs(final$Z_Estimate),lower.tail=FALSE)
        
        if(modelchi == TRUE){
          ##replace V1-VX general form in output with user provided trait names
          for(g in 1:nrow(final)){
            for(p in 1:length(traits)){
              final$lhs[[g]]<-ifelse(final$lhs[[g]] %in% S_names[[p]], gsub(final$lhs[[g]], traits[[p]], final$lhs[[g]]), final$lhs[[g]])
              final$rhs[[g]]<-ifelse(final$rhs[[g]] %in% S_names[[p]], gsub(final$rhs[[g]], traits[[p]], final$rhs[[g]]), final$rhs[[g]])
            }}
          ##subest to only the pieces the user specified
          Modeltest <- parTable(ReorderModel1)
          for(t in 1:nrow(final)){
            final$test[t]<-ifelse(paste(final$lhs[t], final$op[t], final$rhs[t], sep = "") %in% paste(Modeltest$lhs, Modeltest$op, Modeltest$rhs, sep = ""), 1, 0)
            final$test[t]<-ifelse(paste(final$lhs[t], final$op[t], final$rhs[t], sep = "") %in% paste(Modeltest$rhs, Modeltest$op, Modeltest$lhs, sep = ""), 1, final$test[t])
          }
          final<-subset(final, final$test==1)
          final$test<-NULL
          
          ##add in model fit components to each row
          if(!(is.na(Q_WLS))){
            final$chisq<-rep(Q_WLS,nrow(final))
            final$chisq_df<-lavInspect(Model1_Results, "fit")["df"]
            final$chisq_pval<-pchisq(final$chisq,final$chisq_df,lower.tail=FALSE)
            final$AIC<-rep(Q_WLS + 2*lavInspect(Model1_Results, "fit")["npar"],nrow(final))}else{final$chisq<-rep(NA, nrow(final))
            final$chisq_df<-rep(NA,nrow(final))
            final$chisq_pval<-rep(NA,nrow(final))
            final$AIC<-rep(NA, nrow(final))}
        }
        
        ##add in error and warning messages 
        if(printwarn == TRUE){
          final$error<-ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
          final$warning<-ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]}
        
        ##combine with rs-id, BP, CHR, etc.
        final2<-cbind(i,n,SNPs[[n]][i,],final,row.names=NULL)
        
        if(!(sub[[1]]==FALSE)){
          final2<-subset(final2, paste0(final2$lhs, final2$op, final2$rhs, sep = "") %in% sub)
        }else{##pull results 
          final2$est<-ifelse(final2$op == "<" | final2$op == ">" | final2$op == ">=" | final2$op == "<=", final2$est == NA, final2$est)}
        
        ##results to be put into the output
        final2
        
      }
  }
  
  ##ML estimation
  if(estimation=="ML"){
    results<-foreach(n = icount(int), .combine = 'rbind') %:% 
      
      foreach (i=1:length(V_Full[[n]]), .combine='rbind', .packages = "lavaan") %dopar% { 
        
        #reorder sampling covariance matrix based on what lavaan expects given the specified model
        V_Full_Reorder <- V_Full[[n]][[i]][order,order]
        
        #import the S_Full matrix for appropriate run
        S_Fullrun<-S_Full[[n]][[i]]
        
        if(modelchi == TRUE){
          ##name the columns and rows of the S matrix in general format V1-VX
          rownames(S_Fullrun) <- S_names
          colnames(S_Fullrun) <- S_names
        }
        
        ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
        test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML", sample.nobs = 200, optim.dx.tol = +Inf))
        
        Model_ML <- parTable(Model1_Results)
        
        if(NA %in% Model_ML$se){
          SE<-rep("SE could not be computed", max(Model_ML$free))}else{
            #pull the delta matrix (this doesn't depend on N)
            S2.delt <- lavInspect(Model1_Results, "delta")
            
            ##weight matrix from stage 2
            S2.W <- lavInspect(Model1_Results, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
            bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt) 
            
            #create the "lettuce" part of the sandwich
            lettuce <- S2.W%*%S2.delt
            
            #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
            Ohtt <- bread %*% t(lettuce)%*%V_Full_Reorder%*%lettuce%*%bread  
            
            #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
            SE <- as.matrix(sqrt(diag(Ohtt)))}
        
        #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
        if(":=" %in% Model_ML$op){
          
          print("SEs of ghost parameters are not available for ML estimation")
          
          #pull the ghost parameter point estiamte
          ghost<-subset(Model_ML, Model_ML$op == ":=")[,c(2:4,8,11,14)]
          se.ghost<-rep("SE of ghost parameters not available for ML estimation", count(":=" %in% Model_ML$op)$freq)
          
          ##combine with delta method SE
          ghost2<-cbind(ghost,se.ghost)
          colnames(ghost2)[7]<-"SE"
        }else{
          if(":=" %in% Model_ML$op & (NA %in% Model_ML$se)){
            se.ghost<-rep("SE of ghost parameters not available for ML estimation", count(":=" %in% Model_ML$op)$freq)
            ghost<-subset(Model_ML, Model_ML$op == ":=")[,c(2:4,8,11,14)]
            ghost2<-cbind(ghost,se.ghost)
            colnames(ghost2)[7]<-"SE"}else{}}
        
        if(modelchi == TRUE){
          
          ModelQ_ML <- parTable(Model1_Results)
          
          ##remove any parameter constraint labels or ghost parameters
          ModelQ_ML<-subset(ModelQ_ML, ModelQ_ML$plabel != "")
          
          ##identify components of saturated model not already estimated and label as 1
          for (g in 1:nrow(ModelQ_ML)){
            ModelQ_ML$free[g]<-ifelse((paste(ModelQ_ML$lhs[g], ModelQ_ML$op[g], ModelQ_ML$rhs[g], sep = "") %in% modeltest2$write.test.k) & ModelQ_ML$est[g] == 0, 1, 0)
          }
          
          ##identify components of saturated model that were already estimated in user model and label as 2
          for (g in 1:nrow(ModelQ_ML)){
            ModelQ_ML$free[g]<-ifelse((paste(ModelQ_ML$lhs[g], ModelQ_ML$op[g], ModelQ_ML$rhs[g], sep = "") %in% modeltest2$write.test.k) & ModelQ_ML$est[g] != 0, 2, ModelQ_ML$free[g])
          }
          
          tester<-vector(mode="list",length=nrow(ModelQ_ML))
          
          ##replace components of saturdated model already estimated with residaul factor estimates
          for(g in 1:nrow(ModelQ_ML)){
            if(ModelQ_ML$free[g] == 2) { 
              t<-paste(ModelQ_ML$lhs[g], ModelQ_ML$op[g], ModelQ_ML$rhs[g], sep = "")
              t2<-gsub("V", "VF", t)
              tester[[g]]<-t2}else{}
          }
          
          test2<-Filter(Negate(is.null), tester)
          
          for (g in 1:nrow(ModelQ_ML)){
            ModelQ_ML$free[g]<-ifelse((paste(ModelQ_ML$lhs[g], ModelQ_ML$op[g], ModelQ_ML$rhs[g], sep = "") %in% test2), 1, ModelQ_ML$free[g])
          }
          
          ModelQ_ML$free<-ifelse(ModelQ_ML$free != 1, 0, ModelQ_ML$free)
          
          #want to freely estimate the residual factor variances and the residual covariances
          p<-length(ModelQ_ML$free)-z
          
          ModelQ_ML <- ModelQ_ML[order(ModelQ_ML$free),] 
          
          ModelQ_ML$free <- c(rep(0, p),1:z)
          
          ModelQ_ML$ustart <- ModelQ_ML$est
          ModelQ_ML$ustart<-ifelse(ModelQ_ML$free > 0, .05, ModelQ_ML$ustart)
          testQ<-tryCatch.W.E(ModelQ_Results_ML <- sem(model = ModelQ_ML, sample.cov = S_Fullrun, estimator = "ML", sample.nobs=200, start =  ModelQ_ML$ustart, optim.dx.tol = +Inf)) 
          testQ$warning$message[1]<-ifelse(is.null(testQ$warning$message), testQ$warning$message[1]<-"Safe", testQ$warning$message[1])
          
          testQ$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_ML, "se")$theta[1,2]) == TRUE, testQ$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ$warning$message[1])
          
          if(as.character(testQ$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
            
            ModelQ_ML$ustart<-ifelse(ModelQ_ML$free > 0, .01, ModelQ_ML$ustart)
            
            testQ2<-tryCatch.W.E(ModelQ_Results_ML <- sem(model = ModelQ_ML, sample.cov = S_Fullrun, estimator = "ML", sample.nobs=200, start =  ModelQ_ML$ustart, optim.dx.tol = +Inf)) 
          }else{testQ2<-testQ}
          
          
          testQ2$warning$message[1]<-ifelse(is.null(testQ2$warning$message), testQ2$warning$message[1]<-"Safe", testQ2$warning$message[1])
          testQ2$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_ML, "se")$theta[1,2]) == TRUE, testQ2$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ2$warning$message[1])
          
          if(as.character(testQ2$warning$message)[1] == "lavaan WARNING: model has NOT converged!"){
            
            ModelQ_ML$ustart<-ifelse(ModelQ_ML$free > 0, .1, ModelQ_ML$ustart)
            
            testQ3<-tryCatch.W.E(ModelQ_Results_ML <- sem(model = ModelQ_ML, sample.cov = S_Fullrun, estimator = "ML", sample.nobs=200, start =  ModelQ_ML$ustart, optim.dx.tol = +Inf)) 
          }else{testQ3<-testQ2}
          
          testQ3$warning$message[1]<-ifelse(is.null(testQ3$warning$message), testQ3$warning$message[1]<-"Safe", testQ3$warning$message[1])
          testQ3$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_ML, "se")$theta[1,2]) == TRUE, testQ3$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testQ3$warning$message[1])
          
          if(as.character(testQ3$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
            
            #pull the delta matrix (this doesn't depend on N)
            S2.delt_Q <- lavInspect(ModelQ_Results_ML, "delta")
            
            ##weight matrix from stage 2
            S2.W_Q <- lavInspect(ModelQ_Results_ML, "WLS.V") 
            
            #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
            bread_Q <- solve(t(S2.delt_Q)%*%S2.W_Q%*%S2.delt_Q) 
            
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
            eta_test<-parTable(ModelQ_Results_ML)
            eta_test<-subset(eta_test, eta_test$free != 0)
            eta<-cbind(eta_test[,14])
            
            #Combining all the pieces from above:
            Q_ML<-t(eta)%*%P1%*%solve(Eig)%*%t(P1)%*%eta}else{Q_ML<-NA}
        }
        
        ##remove parameter constraints, ghost parameters, and fixed effects from output to merge with SEs
        unstand<-subset(Model_ML, Model_ML$plabel != "" & Model_ML$free > 0)[,c(2:4,8,11,14)]
        
        ##combine ghost parameters with rest of output
        if(exists("ghost2") == "TRUE"){
          unstand2<-rbind(cbind(unstand,SE),ghost2)
        }else{unstand2<-cbind(unstand,SE)}
        
        ##add in fixed effects and parameter constraints to output
        other<-subset(Model_ML, (Model_ML$plabel == "" & Model_ML$op != ":=") | (Model_ML$free == 0 & Model_ML$plabel != ""))[,c(2:4,8,11,14)]
        other$SE<-rep(NA, nrow(other))
        
        ##combine fixed effects and parameter constraints with output if there are any
        if(nrow(other) > 0){
          final<-rbind(unstand2,other)
        }else{final<-unstand2}
        
        ##add in p-values
        final$Z_Estimate<-final$est/final$SE
        final$Pval_Estimate<-2*pnorm(abs(final$Z_Estimate),lower.tail=FALSE)
        
        #reorder based on row numbers so it is in order the user provided
        final$index <- as.numeric(row.names(final))
        final<-final[order(final$index), ]
        final$index<-NULL
        
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
          if(!(is.na(Q_ML))){
            final$chisq<-rep(Q_ML,nrow(final))
            final$chisq_df<-lavInspect(Model1_Results, "fit")["df"]
            final$chisq_pval<-pchisq(final$chisq,final$chisq_df,lower.tail=FALSE)
            final$AIC<-rep(Q_ML + 2*lavInspect(Model1_Results, "fit")["npar"],nrow(final))}else{final$chisq<-rep(NA, nrow(final))
            final$chisq_df<-rep(NA,nrow(final))
            final$chisq_pval<-rep(NA,nrow(final))
            final$AIC<-rep(NA, nrow(final))}
        }
        
        ##add in error and warning messages 
        if(printwarn == TRUE){
          final$error<-ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
          final$warning<-ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]}
        
        ##combine with rs-id, BP, CHR, etc.
        final2<-cbind(i,n,SNPs[[n]][i,],final,row.names=NULL)
        
        if(!(sub[[1]]==FALSE)){
          final2<-subset(final2, paste0(final2$lhs, final2$op, final2$rhs, sep = "") %in% sub)
        }else{##pull results 
          final2$est<-ifelse(final2$op == "<" | final2$op == ">" | final2$op == ">=" | final2$op == "<=", final2$est == NA, final2$est)}
        
        ##results to be put into the output
        final2
      }
  }
  
  
  ##sort results so it is in order of the output lists provided for the function
  results<- results[order(results$i, results$n),] 
  results$n<-NULL
  
  if(!sub[[1]] == FALSE){
    results$i<-NULL
    Results_List<-vector(mode="list",length=length(sub))
    for(y in 1:length(sub)){
      Results_List[[y]]<-as.data.frame(matrix(NA,ncol=ncol(results),nrow=nrow(results)/length(sub)))
      colnames(Results_List[[y]])<-colnames(results)
      Results_List[[y]]<-subset(results, paste0(results$lhs, results$op, results$rhs, sep = "") %in% sub[[y]])
    }
    rm(results)
  }
  
  
  if(sub[[1]]==FALSE){
    x<-count(results$chisq == results$chisq[[1]])[2,2]
    results$i<-rep(1:(nrow(results)/x),each=x)
    colnames(results)[1]<-"Model_Number"
    Results_List<-vector(mode="list",length=length(Output[[1]]))
    for(y in 1:max(factorPa$Model_Number)){
      Results_List[[y]]<-subset(results, results$Model_Number == y)
      Results_List[[y]]$Model_Number<-NULL
    }
    rm(results)
  }
  
  time_all<-proc.time()-time
  print(time_all[3])
  
  return(Results_List)
  
}
