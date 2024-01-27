usermodel <-function(covstruc,estimation="DWLS", model = "", CFIcalc=TRUE, std.lv=FALSE, imp_cov=FALSE,fix_resid=TRUE,toler=NULL){ 
  time<-proc.time()
  ##determine if the model is likely being listed in quotes and print warning if so
  test<-c(str_detect(model, "~"),str_detect(model, "="),str_detect(model, "\\+"))
  if(any(test) != TRUE){
    warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned an error about not locating the ReorderModel.")
  }
  
  ##read in the LD portion of the V (sampling covariance) matrix
  V_LD<-as.matrix(covstruc[[1]])
  
  ##read in the LD portion of the S (covariance) matrix
  S_LD<-as.matrix(covstruc[[2]])
  
  ##k = number of phenotypes in dataset (i.e., number of columns in LD portion of S matrix)
  k<-ncol(S_LD)
  
  ##size of V matrix used later in code to create diagonal V matrix
  z<-(k*(k+1))/2
  
  Model1<-model
  
  ##pull the column names specified in the munge function
  S_names<-colnames(S_LD)
  rownames(S_LD)<-colnames(S_LD)
  
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
  
  for(i in 1:length(S_names)){
    b<-grepl(S_names[i], model)
    if(b == FALSE){
      remove<-paste0("\\b", colnames(S_LD)[i],"\\b",sep="")
      remove2[w]<-i
      V_LD <- V_LD[-grep(pattern=remove[1],row.names(V_LD)),-grep(pattern=remove[1],colnames(V_LD))]
      w<-w+1
      if (!(is.matrix(V_LD))) {
        stop("None of the trait names in the LDSC output match names in the model")
      }
    }else{}
  }
  
  if(is.null(remove2) == FALSE){
    S_LD<-S_LD[-remove2,-remove2]
  }
  
  ##redefine k and z and model names after removing non-used variables
  k<-ncol(S_LD)
  z<-(k*(k+1))/2

  ##smooth to near positive definite if either V or S are non-positive definite
  S_LDb<-S_LD
  smooth1<-ifelse(eigen(S_LD)$values[nrow(S_LD)] <= 0, S_LD<-as.matrix((nearPD(S_LD, corr = FALSE))$mat), S_LD<-S_LD)
  LD_sdiff<-max(abs(S_LD-S_LDb))
  
  V_LDb<-V_LD
  smooth2<-ifelse(eigen(V_LD)$values[nrow(V_LD)] <= 0, V_LD<-as.matrix((nearPD(V_LD, corr = FALSE))$mat), V_LD<-V_LD)
  LD_sdiff2<-max(abs(V_LD-V_LDb))
  
  SE_pre<-matrix(0, k, k)
  SE_pre[lower.tri(SE_pre,diag=TRUE)] <-sqrt(diag(V_LDb))
  
  SE_post<-matrix(0, k, k)
  SE_post[lower.tri(SE_post,diag=TRUE)] <-sqrt(diag(V_LD))
  
  Z_pre<-S_LDb/SE_pre
  Z_post<-S_LD/SE_post
  Z_diff<-(Z_pre-Z_post)
  Z_diff[which(!is.finite(Z_diff))]<-0
  Z_diff<-max(Z_diff)
  rm(V_LDb,S_LDb,Z_pre,Z_post)
  
  ##run model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  #transform V_LD matrix into a weight matrix: 
  W <- solve(V_LD)
  
  if(CFIcalc==TRUE){
    
    ##code to write null model for calculation of CFI
    write.null<-function(k, label = "V", label2 = "VF") {
      Model3<-""
      for (p in 1:k) {
        linestart3 <- paste(colnames(S_LD)[p], " ~~ ", colnames(S_LD)[p], sep = "")
        Model3<-paste(Model3, linestart3, " \n ", sep = "")}
      
      Model2<-""
      for (p in 1:k) {
        linestart2 <- paste(label2, p, " =~ 1*", colnames(S_LD)[p], sep = "")
        Model2<-paste(Model2, linestart2, " \n ", sep = "")}
      
      Modelsat<-""
      for (i in 1:(k-1)) {
        linestartc <- paste(colnames(S_LD)[i], " ~~ 0*", colnames(S_LD)[i+1],  sep = "")
        if (k-i >= 2) { 
          linemidc <- ""
          for (j in (i+2):k) {
            linemidc <- paste(linemidc, " + 0*", colnames(S_LD)[j], sep = "")
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
    
    ##create inependence model for calculation of CFI
    modelCFI<-write.null(k)
    
    ##run CFI model so it knows the reordering for the independence model
    empty<-.tryCatch.W.E(fitCFI <- sem(modelCFI, sample.cov = S_LD, estimator = "DWLS", WLS.V = W,sample.nobs=2, optim.dx.tol = +Inf,optim.force.converged=TRUE,control=list(iter.max=1)))
   
    orderCFI <- .rearrange(k = k, fit =  fitCFI, names =  rownames(S_LD))
    
    ##reorder matrix for independence (i.e., null) model for CFI calculation
    V_Reorder2 <- V_LD[orderCFI,orderCFI]
    W_CFI<-diag(z)
    diag(W_CFI)<-diag(V_Reorder2)
    W_CFI<-solve(W_CFI)

  }
                         
  empty3<-.tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", WLS.V = W, sample.nobs = 2,warn=FALSE,std.lv=std.lv, optim.dx.tol = +Inf,optim.force.converged=TRUE,control=list(iter.max=1)))

  r<-nrow(lavInspect(ReorderModel, "cor.lv"))
  
  if(class(empty3$value) != "lavaan"){
    warning(paste("The function has stopped due to convergence issues for your primary model. Please contact us with your specific model and variables used or try specifying an alternative model"))
  }
  
  ##save the ordering
  order <- .rearrange(k = k, fit = ReorderModel, names = rownames(S_LD))
  
  ##reorder the weight (inverted V_LD) matrix
  V_Reorder<-V_LD[order,order]
  W_Reorder<-diag(z)
  diag(W_Reorder)<-diag(V_Reorder)
  W_Reorder<-solve(W_Reorder)

    print("Running primary model")
    
    if(estimation == "DWLS"){
    ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
    empty4<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "DWLS", std.lv=std.lv,WLS.V = W_Reorder, sample.nobs = 2,optim.dx.tol = +Inf))
    }
    
    if(estimation == "ML"){
        empty4<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "ML", sample.nobs = 200,std.lv=std.lv, optim.dx.tol = +Inf,sample.cov.rescale=FALSE))
    }
    
    empty4$warning$message[1]<-ifelse(is.null(empty4$warning$message), empty4$warning$message[1]<-0, empty4$warning$message[1])
    
    if(fix_resid == TRUE){
      if(class(empty4$value)[1] == "simpleError" | lavInspect(Model1_Results,"converged") == FALSE){
      
          #create unique combination of letters for residual variance parameter labels
          n<-combn(letters,4)[,sample(1:14000, k, replace=FALSE)]
          
      Model3<-""
      for (p in 1:k) {
        linestart3a <- paste(colnames(S_LD)[p], " ~~ ",  paste(n[,p],collapse=""), "*", colnames(S_LD)[p], sep = "")
        linestart3b <- paste(paste(n[,p],collapse=""), " > .0001", sep = "")
        Model3<-paste(Model3, " \n ", linestart3a, " \n ", linestart3b, " \n ", sep = "")}
      
      Model1<-paste(Model1,Model3)
        
        if(estimation == "DWLS"){
          empty4<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "DWLS",std.lv=std.lv, WLS.V = W_Reorder, sample.nobs = 2, optim.dx.tol = +Inf))
        }
        
        if(estimation == "ML"){
            empty4<-.tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_LD, estimator = "ML", std.lv=std.lv, sample.nobs = 200,optim.dx.tol = +Inf,sample.cov.rescale=FALSE))
        }
        
        #if adding in residuals fixed above 0 is duplicating user provided arguments then revert to original model
       if(class(empty4$value)[1] != "lavaan"){
        if(grepl("duplicate", as.character(empty4$value)[1]) == TRUE){
        Model1<-model
      }else{print("The model as initially specified failed to converge. A lower bound of 0 on residual variances was automatically added to try and troubleshoot this. This behavior can be toggled off by setting the fix_resid argument to FALSE.")
      }
      }else{print("The model as initially specified failed to converge. A lower bound of 0 on residual variances was automatically added to try and troubleshoot this. This behavior can be toggled off by setting the fix_resid argument to FALSE.")}
     
      }
    }
    
    if(class(empty4$value)[1] == "simpleError"){
      warning("The model failed to converge on a solution. Please try specifying an alternative model")}
    
    #pull the delta matrix (this doesn't depend on N)
    ##note that while the delta matrix is reordered based on the ordering in the model specification
    ##that the lavaan output is also reordered so that this actually ensures that the results match up 
    S2.delt <- lavInspect(Model1_Results, "delta")
    
    ##weight matrix from stage 2. S2.W is not reordered by including something like model constraints
    S2.W <- lavInspect(Model1_Results, "WLS.V") 
    
    #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
  bread2<-.tryCatch.W.E(bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt,tol=toler))
     
    if(!(is.null(empty4$warning))){
      if(lavInspect(Model1_Results,"converged") == FALSE){
        warning("The model failed to converge on a solution. Please try specifying an alternative model.")
      }}
    
    if(class(bread2$value)[1] != "matrix"){
      warning("Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
            The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
            that these results should not be interpreted.")
      check<-1
      unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
      unstand<-subset(unstand, unstand$free != 0)                    
      unstand$free<-NULL
      results<-unstand
      colnames(results)=c("lhs","op","rhs","Unstandardized_Estimate")
    
      print(results)  
    }
    
    if(class(bread2$value)[1] == "matrix"){
      #create the "lettuce" part of the sandwich
      lettuce <- S2.W%*%S2.delt
      
      #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
      Ohtt <- bread %*% t(lettuce)%*%V_Reorder%*%lettuce%*%bread  
      
      #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
      SE <- as.matrix(sqrt(diag(Ohtt)))
      
      Model_Output <- parTable(Model1_Results)
      
      constraints<-subset(Model_Output$label, Model_Output$label != "")
      constraints2<-duplicated(constraints)
      
      #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
      if(estimation == "DWLS"){
      if(":=" %in% Model_Output$op & !(is.na(SE[1]))){
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
        ghost<-subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
        
        ##combine with delta method SE
        ghost2<-cbind(ghost,se.ghost)
        colnames(ghost2)[7]<-"SE"
        
      }else{se.ghost<-NA
      if(":=" %in% Model_Output$op & is.na(se.ghost[1])){
        se.ghost<-rep("SE could not be computed", count(":=" %in% Model_Output$op)$freq)
        ghost<-subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
        ghost2<-cbind(ghost,se.ghost)
        colnames(ghost2)[7]<-"SE"}else{}}
      }
      
      if(estimation == "ML"){
        if(":=" %in% Model_Output$op){
          
          print("SEs of ghost parameters are not available for ML estimation")
          
          #pull the ghost parameter point estiamte
          ghost<-subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
          se.ghost<-rep(NA, sum(":=" %in% Model_Output$op))
          warning("SE for ghost parameter not available for ML")
          ##combine with delta method SE
          ghost2<-cbind(ghost,se.ghost)
          colnames(ghost2)[7]<-"SE"
        }
      }
      
      
      ##check whether correlations among latent variables is positive definite
      if(r > 1){
        empty<-.tryCatch.W.E(check<-lowerTriangle(lavInspect(Model1_Results,"cor.lv")[1:r,1:r]))
        t<-max(check)
        t2<-min(check)}else{
          t<-1
          t2<--1
        }
      
      if(t > 1 | t2 < -1  | t == "NaN"){
        print("Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. 
              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below 
              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested.")
        
        unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
        unstand<-subset(unstand, unstand$free != 0)                    
        unstand$free<-NULL
        results<-unstand
        colnames(results)=c("lhs","op","rhs","Unstandardized_Estimate")
        
        if(exists("ghost2") == "TRUE"){
          ghost2$free<-NULL
          ghost2$label<-NULL
          unstand2<-rbind(cbind(results,SE),ghost2)
        }else{unstand2<-cbind(results,SE)}
        
        print(unstand2)
        check<-1
      }else{
        check<-2
        
        #calculate model chi-square
        Eig<-as.matrix(eigen(V_LD)$values)
        Eig2<-diag(z)
        diag(Eig2)<-Eig
        
        #Pull P1 (the eigen vectors of V_eta)
        P1<-eigen(V_LD)$vectors
        
        implied<-as.matrix(fitted(Model1_Results))[1]
        implied_order<-colnames(S_LD)
        implied[[1]]<-implied[[1]][implied_order,implied_order]
        implied2<-S_LD-implied[[1]]
        eta<-as.vector(lowerTriangle(implied2,diag=TRUE))
        Q<-t(eta)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta

        if(CFIcalc == TRUE){
          print("Calculating CFI")
          ##now CFI
          ##run independence model
          if(estimation == "DWLS"){
          testCFI<-.tryCatch.W.E(fitCFI <- sem(modelCFI, sample.cov =  S_LD, estimator = "DWLS", WLS.V = W_CFI, sample.nobs=2, optim.dx.tol = +Inf))
          }
          
          if(estimation == "ML"){
            testCFI<-.tryCatch.W.E(fitCFI <- sem(modelCFI, sample.cov =  S_LD, estimator = "ML",sample.nobs=200, optim.dx.tol = +Inf,sample.cov.rescale=FALSE))
          }
          testCFI$warning$message[1]<-ifelse(is.null(testCFI$warning$message), testCFI$warning$message[1]<-"Safe", testCFI$warning$message[1])
          testCFI$warning$message[1]<-ifelse(is.na(inspect(fitCFI, "se")$theta[1,2]) == TRUE, testCFI$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testCFI$warning$message[1])
          
          if(as.character(testCFI$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
            
            ##code to estimate chi-square of independence model#
            #First pull the estimates from Step 2
            ModelQ_CFI <- parTable(fitCFI)
            p2<-length(ModelQ_CFI$free)-z
            
            ##fix variances and freely estimate covariances
            ModelQ_CFI$free <- c(rep(0, p2), 1:z)
            ModelQ_CFI$ustart <- ModelQ_CFI$est
            
            if(estimation == "DWLS"){
            testCFI2<-.tryCatch.W.E(ModelQ_Results_CFI <- sem(model = ModelQ_CFI, sample.cov = S_LD, estimator = "DWLS", WLS.V = W_CFI, sample.nobs=2, optim.dx.tol = +Inf))
            }
            
            if(estimation == "ML"){
              testCFI2<-.tryCatch.W.E(ModelQ_Results_CFI <- sem(model = ModelQ_CFI, sample.cov = S_LD, estimator = "ML", sample.nobs=200, optim.dx.tol = +Inf,sample.cov.rescale=FALSE))
            }
            
            testCFI2$warning$message[1]<-ifelse(is.null(testCFI2$warning$message), testCFI2$warning$message[1]<-"Safe", testCFI2$warning$message[1])
            testCFI2$warning$message[1]<-ifelse(is.na(inspect(ModelQ_Results_CFI , "se")$theta[1,2]) == TRUE, testCFI2$warning$message[1]<-"lavaan WARNING: model has NOT converged!", testCFI2$warning$message[1])
            
            if(as.character(testCFI2$warning$message)[1] != "lavaan WARNING: model has NOT converged!"){
              
              #pull the delta matrix (this doesn't depend on N)
              S2.delt_Q_CFI <- lavInspect(ModelQ_Results_CFI, "delta")
              
              ##weight matrix from stage 2
              S2.W_Q_CFI <- lavInspect(ModelQ_Results_CFI, "WLS.V") 
              
              #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
              bread_Q_CFI <- solve(t(S2.delt_Q_CFI)%*%S2.W_Q_CFI%*%S2.delt_Q_CFI) 
              
              #create the "lettuce" part of the sandwich
              lettuce_Q_CFI <- S2.W_Q_CFI%*%S2.delt_Q_CFI
              
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
              eta_test_CFI<-parTable(ModelQ_Results_CFI)
              eta_test_CFI<-subset(eta_test_CFI, eta_test_CFI$free != 0)
              eta_CFI<-cbind(eta_test_CFI[,14])
              
              #Combining all the pieces from above:
              Q_CFI<-t(eta_CFI)%*%P1_CFI%*%solve(Eig_CFI)%*%t(P1_CFI)%*%eta_CFI}else{Q_CFI<-"The null (i.e. independence) model did not converge"}}
          
          ##df of independence Model
          dfCFI <- (((k * (k + 1))/2) - k)
          
          ##df of user model
          df <- lavInspect(Model1_Results, "fit")["df"]
          
          if(!(is.character(Q_CFI)) & !(is.character(Q))){
            CFI<-as.numeric(((Q_CFI-dfCFI)-(Q-df))/(Q_CFI-dfCFI))
            CFI<-ifelse(CFI > 1, 1, CFI)
          }else{CFI<-"Either the chi-square or null (i.e. independence) model did not converge"}
          
        }
        
        print("Calculating Standardized Results")
        ##transform the S covariance matrix to S correlation matrix
        D=sqrt(diag(diag(S_LD)))
        S_Stand=solve(D)%*%S_LD%*%solve(D)
        rownames(S_Stand)<-rownames(S_LD)
        colnames(S_Stand)<-colnames(S_Stand)
        
        #obtain diagonals of the original V matrix and take their sqrt to get SE's
        Dvcov<-sqrt(diag(V_LD))
        
        #calculate the ratio of the rescaled and original S matrices
        scaleO=as.vector(lowerTriangle((S_Stand/S_LD),diag=T))
        
        ## MAke sure that if ratio in NaN (devision by zero) we put the zero back in: ### TEMP STUPID MICHEL FIX!
        scaleO[is.nan(scaleO)] <- 0
        
        #rescale the SEs by the same multiples that the S matrix was rescaled by
        Dvcovl<-as.vector(Dvcov*t(scaleO))
        
        #obtain the sampling correlation matrix by standardizing the original V matrix
        Vcor<-cov2cor(V_LD)
        
        #rescale the sampling correlation matrix by the appropriate diagonals
        V_stand<-diag(Dvcovl)%*%Vcor%*%diag(Dvcovl)
        V_stand2<-diag(z)
        diag(V_stand2)<-diag(V_stand)
        
        ### make sure no value on the diagonal of V is 0 
        diag(V_stand2)[diag(V_stand2) == 0] <- 2e-9
        
        W_stand<-solve(V_stand2[order,order])
        
        if(estimation == "DWLS"){
          emptystand<-.tryCatch.W.E(Fit_stand <- sem(Model1, sample.cov = S_Stand, estimator = "DWLS", WLS.V = W_stand, std.lv=std.lv,sample.nobs = 2, optim.dx.tol = +Inf))
          if(is.null(emptystand$warning$message[1])) {
             emptystand$warning$message[1] <- 0
             }
        }
        
        if(estimation == "ML"){
            emptystand<-.tryCatch.W.E(Fit_stand <- sem(Model1, sample.cov = S_Stand, estimator = "ML",  sample.nobs = 200, std.lv=std.lv, optim.dx.tol = +Inf,sample.cov.rescale=FALSE))
            if(is.null(emptystand$warning$message[1])) {
             emptystand$warning$message[1] <- 0
             }
        }
      
        ##perform same procedures for sandwich correction as in the unstandardized case
        delt_stand <- lavInspect(Fit_stand, "delta") 
        
        W_stand <- lavInspect(Fit_stand, "WLS.V") 

        bread_stand2<-.tryCatch.W.E(bread_stand <- solve(t(delt_stand)%*%W_stand %*%delt_stand,tol=toler))
   
        if(class(bread_stand2$value)[1] != "matrix" | lavInspect(Fit_stand,"converged") == FALSE | class(emptystand)[1] == "simpleError"){
          warning("The standardized model failed to converge. This likely indicates more general problems with the model solution. Unstandardized results are printed below but this should be interpreted with caution.")
          
          unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
          unstand<-subset(unstand, unstand$free != 0)                    
          unstand$free<-NULL
          results<-unstand
          colnames(results)=c("lhs","op","rhs","Unstandardized_Estimate")
          
          if(exists("ghost2") == "TRUE"){
            ghost2$free<-NULL
            ghost2$label<-NULL
            unstand2<-rbind(cbind(results,SE),ghost2)
          }else{unstand2<-cbind(results,SE)}
          
          print(unstand2)
          check<-1
        }else{
          
          lettuce_stand <- W_stand%*%delt_stand
          Vcov_stand<-as.matrix(V_stand[order,order])
          Ohtt_stand <- bread_stand %*% t(lettuce_stand)%*%Vcov_stand%*%lettuce_stand%*%bread_stand
          SE_stand <- as.matrix(sqrt(diag(Ohtt_stand)))
          
          Model_Stand <- parTable(Fit_stand)
          
          if(estimation == "DWLS"){
          #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
          if(":=" %in% Model_Stand$op & !(NA %in% Model_Stand$se)){
            #variance-covariance matrix of parameter estimates, q-by-q (this is the naive one)
            vcov <- lavInspect(Fit_stand, "vcov") 
            
            #internal lavaan representation of the model
            lavmodel <- Fit_stand@Model 
            
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
            ghost_stand<-subset(Model_Stand,  Model_Stand$op == ":=")[,c(2:4,8,11,14)]
            
            ##combine with delta method SE
            ghost2_stand<-cbind(ghost_stand,se.ghost_stand)
            colnames(ghost2_stand)[7]<-"SE_stand"
          }else{
            if(":=" %in% Model_Stand$op & (NA %in% Model_Stand$se)){
              se.ghost_stand<-rep("SE could not be computed", count(":=" %in% Model_Stand$op)$freq)
              ghost_stand<-subset(Model_Stand, Model_Stand$op == ":=")[,c(2:4,8,11,14)]
              ghost2_stand<-cbind(ghost_stand,se.ghost_stand)
              colnames(ghost2_stand)[7]<-"SE_stand"}else{}}
          }
          
          if(estimation == "ML"){
            if(":=" %in% Model_Stand$op){
              
              print("SEs of ghost parameters are not available for ML estimation")
              
              #pull the ghost parameter point estiamte
              ghost_stand<-subset(Model_Stand, Model_Stand$op == ":=")[,c(2:4,8,11,14)]
              se.ghost_stand<-rep(NA, sum(":=" %in% Model_Stand$op))
              warning("SEs for ghost parameters are not available for ML estimation")
              ##combine with delta method SE
              ghost2<-cbind(ghost_stand,se.ghost_stand)
              colnames(ghost2_stand)[7]<-"SE_stand"
            }
          }
          
          unstand<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
          unstand<-subset(unstand, unstand$free != 0)                    
          unstand$free<-NULL
          
          ##combine ghost parameters with rest of output
          if(exists("ghost2") == "TRUE"){
            ghost2$free<-NULL
            ghost2$label<-NULL
            unstand2<-rbind(cbind(unstand,SE),ghost2)
          }else{unstand2<-cbind(unstand,SE)}
          
          stand<-data.frame(inspect(Fit_stand,"list")[,c(8,14)])
          stand<-subset(stand, stand$free != 0)
          stand$free<-NULL
          
          ##combine ghost parameters with rest of output
          if(exists("ghost2_stand") == "TRUE"){
            ghost2_stand[,1:5]<-NULL
            stand2<-rbind(cbind(stand,SE_stand),ghost2_stand)
          }else{stand2<-cbind(stand,SE_stand)}
          
          colnames(stand2)<-c("est_stand","se_stand")
          
          ##df of user model
          df<-lavInspect(Model1_Results, "fit")["df"]
          
          if(!(is.character(Q))){
            chisq<-Q
            AIC<-(Q + 2*lavInspect(Model1_Results, "fit")["npar"])}else{chisq<-Q
            AIC<-NA}
          
          print("Calculating SRMR")
          
          SRMR<-lavInspect(Model1_Results, "fit")["srmr"]
          
          if(CFIcalc == TRUE){
            modelfit<-cbind(chisq,df,AIC,CFI,SRMR)}else{modelfit<-cbind(chisq,df,AIC,SRMR)}
          
          std_all<-standardizedSolution(Fit_stand)
          std_all<-subset(std_all, !(is.na(std_all$pvalue)))
          
          results<-cbind(unstand2, stand2)
          
          ##add in fixed effects
          base_model<-data.frame(inspect(ReorderModel, "list")[,c(2:4,8,14)])
          base_model<-subset(base_model,  !(paste0(base_model$lhs, base_model$op,base_model$rhs) %in% paste0(unstand2$lhs, unstand2$op, unstand2$rhs)))
          base_model<-subset(base_model, base_model$op == "=~" | base_model$op == "~~" | base_model$op == "~")
          if(nrow(base_model) > 0){
            base_model$free<-NULL
            base_model$SE<-""
            base_model[6]<-base_model$est
            base_model$SE_stand<-""
            colnames(base_model)<-colnames(results)
            results<-rbind(results,base_model)
          }
          std_all<-subset(std_all,  paste0(std_all$lhs, std_all$op, std_all$rhs) %in% paste0(results$lhs, results$op, results$rhs))
          std_all$order<-paste0(std_all$lhs, std_all$op, std_all$rhs)
          std_all<-data.frame(std_all$est.std,std_all$order)
          colnames(std_all)<-c("est.std","order")
          results$order<-paste0(results$lhs,results$op,results$rhs)
          results$order2<-1:nrow(results)
          results<-suppressWarnings(merge(results,std_all,by="order",all=T))
          results$est.std<-ifelse(is.na(results$est.std), results$est_stand, results$est.std)
          results<-results[order(results$order2),]
          results$order<-NULL
          results$order2<-NULL
        } 
      }
    }
  
  
  if(class(bread2$value)[1] == "matrix" & check == 2){  
    ##name the columns of the results file
    colnames(results)=c("lhs","op","rhs","Unstand_Est","Unstand_SE","STD_Genotype","STD_Genotype_SE", "STD_All")
   
    ##name model fit columns
    if(CFIcalc == TRUE){
      colnames(modelfit)=c("chisq","df","AIC","CFI","SRMR")}else{colnames(modelfit)=c("chisq","df","AIC","SRMR")}
    
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
            warning(paste("CFI estimates below 0 should not be trusted, and indicate that the other model fit estimates should be interpreted with caution. A negative CFI estimates typically appears due to negative residual variances."))
          }}
        modelfit$CFI<-ifelse(modelfit$df == 0, modelfit$CFI == NA, modelfit$CFI)
      }else{order<-c(1,2,5,3,4)
      modelfit<-modelfit[,order]
      }}
    
    time_all<-proc.time()-time
    print(time_all[3])
    
    if(modelfit$df == 0){
      print("Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)")
    }
    
    
    if(LD_sdiff > 0){
      print(paste("The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was ", LD_sdiff, "As a result of the smoothing, the largest Z-statistic change for the genetic covariances was ", Z_diff, ". We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS.", sep = " "))
    }
    
    if(LD_sdiff > .025){
      warning("A difference greater than .025 was observed pre- and post-smoothing in the genetic covariance matrix. This reflects a large difference and results should be interpreted with caution!! This can often result from including low powered traits, and you might consider removing those traits from the model. If you are going to run a multivariate GWAS we strongly recommend setting the smooth_check argument to true to check smoothing for each SNP.")
    }
    
    if(Z_diff > .025){
      warning("A difference greater than .025 was observed pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a large difference and results should be interpreted with caution!! This can often result from including low powered traits, and you might consider removing those traits from the model. If you are going to run a multivariate GWAS we strongly recommend setting the smooth_check argument to true to check smoothing for each SNP.")
    }
    
    if(LD_sdiff2 > 0){
      print(paste("The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was ", LD_sdiff2, "As a result of the smoothing, the largest Z-statistic change for the genetic covariances was ", Z_diff,  ". We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS.", sep = " "))
    }
    
    if(any(constraints2 == TRUE)){
      print("Please note that when equality constraints are used in the current version of Genomic SEM that the standardized output will also impose the same constraint.")
    }
    results$p_value<-2*pnorm(abs(as.numeric(results$Unstand_Est)/as.numeric(results$Unstand_SE)),lower.tail=FALSE)
    results$p_value<-ifelse(results$p_value == 0, "< 5e-300", results$p_value)

   if(empty4$warning$message[1] != 0 & !grepl("not recommended for continuous data", empty4$warning$message[1])){
      warning(paste0("The unstandardized model produced the following warning: ", empty4$warning$message[1],sep=""))
    }  
    
    if(emptystand$warning$message[1] != 0 & !grepl("not recommended for continuous data", emptystand$warning$message[1])){
      warning(paste0("The standardized model produced the following warning: ", emptystand$warning$message[1],sep=""))
    }  
    
    if(imp_cov == FALSE){
      return(list(modelfit=modelfit,results=results))
    }
    
    if(imp_cov == TRUE){
      resid_cov<-list()
      resid_cov[[1]]<-implied[[1]]  
      resid_cov[[2]]<-implied2
      names(resid_cov) <- c("Model Implied Covariance Matrix", "Residual Covariance Matrix: Calculated as Observed Cov - Model Implied Cov")
      return(list(modelfit=modelfit,results=results,resid_cov=resid_cov))
    }
  }
  
}
