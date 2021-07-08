userGWAS<-function(covstruc=NULL,SNPs=NULL,estimation="DWLS",model="",modelchi=TRUE,printwarn=TRUE,sub=FALSE,cores=NULL,toler=FALSE,SNPSE=FALSE,parallel=TRUE,GC="standard",MPI=FALSE,smooth_check=FALSE){ 
  time<-proc.time()
  
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
  
  Operating<-Sys.info()[['sysname']]
  
  #remove white spacing on subset argument so will exact match lavaan representation of parameter
  if(!(sub[[1]])==FALSE){
    sub<-str_replace_all(sub, fixed(" "), "")
  }
  
  ##make sure SNP and A1/A2 are character columns to avoid being shown as integers in ouput
  SNPs<-data.frame(SNPs)
  SNPs$A1<-as.character(SNPs$A1)
  SNPs$A2<-as.character(SNPs$A2)
  SNPs$SNP<-as.character(SNPs$SNP)
  
  #SNP variance
  varSNP=2*SNPs$MAF*(1-SNPs$MAF)  
  
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
  colnames(S_Full)<-c("SNP", colnames(S_LD))
  
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
    
    test2<-tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Full, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf,optim.force.converged=TRUE,control=list(iter.max=1)))
    
    order <- rearrange(k = k2, fit = ReorderModel, names = rownames(S_Full))
    
    suppressWarnings(df<-lavInspect(ReorderModel, "fit")["df"])
    suppressWarnings(npar<-lavInspect(ReorderModel, "fit")["npar"])
    
  }
  
  SNPs2<-SNPs[,1:6]
  rm(SNPs)
  
  if(parallel==FALSE){
    
    #f = number of SNPs in dataset
    f=nrow(beta_SNP) 
    
    #make empty list object for model results if not saving specific model parameter
    if(sub[[1]]==FALSE){
      Results_List<-vector(mode="list",length=f)}
    
    print("Starting GWAS Estimation")
    for (i in 1:f) { 
      
      #create empty shell of V_SNP matrix
      V_SNP<-diag(k)
      
      #loops to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
      
      #double GC correction using univariate LDSC intercepts
      if(GC == "conserv"){
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[i]^2)}
          if (x == y) {
            V_SNP[x,x]<-(SE_SNP[i,x]*I_LD[x,x]*varSNP[i])^2
          }
        }
      }
      
      #single GC correction using sqrt of univariate LDSC intercepts
      if(GC == "standard"){
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*varSNP[i]^2)}
          if (x == y) {
            V_SNP[x,x]<-(SE_SNP[i,x]*sqrt(I_LD[x,x])*varSNP[i])^2
          }
        }
      }
      
      #no GC correction
      if(GC == "none"){
        for (p in 1:nrow(coords)) { 
          x<-coords[p,1]
          y<-coords[p,2]
          if (x != y) { 
            V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*varSNP[i]^2)}
          if (x == y) {
            V_SNP[x,x]<-(SE_SNP[i,x]*varSNP[i])^2
          }
        }
      }
      
      if(smooth_check == TRUE){
        if(GC == "conserv"){
          Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*diag(I_LD))
        }
        if(GC=="standard"){
          Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*sqrt(diag(I_LD))) 
        }
        if(GC=="none"){
          Z_pre<-beta_SNP[i,]/SE_SNP[i,] 
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
      
      if(eigen(V_Full)$values[kv] <= 0){
        V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat)
        V_smooth<-1
      }
      
      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Full_Reorder <- V_Full[order,order]
      u<-nrow(V_Full_Reorder)
      W<-diag(u)
      diag(W)<-diag(V_Full_Reorder)
      
      ##invert the reordered sampling covariance matrix to create a weight matrix 
      if(toler==FALSE){
        W<- solve(W)
      }
      
      if(toler!=FALSE){
        W <- solve(W,tol=toler)
      }
      
      #create empty vector for S_SNP
      S_SNP<-vector(mode="numeric",length=k+1)
      
      #enter SNP variance from reference panel as first observation
      S_SNP[1]<-varSNP[i]
      
      #enter SNP covariances (standardized beta * SNP variance from refference panel)
      for (p in 1:k) {
        S_SNP[p+1]<-varSNP[i]*beta_SNP[i,p]
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
      
      if(eigen(S_Fullrun)$values[ks] <= 0){
        S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat)
        S_smooth<-1
      }
      
      if(smooth_check == TRUE){
        if(exists("S_smooth") | exists("V_smooth")){
          SE_smooth<-matrix(0, ks, ks)
          SE_smooth[lower.tri(SE_smooth,diag=TRUE)] <-sqrt(diag(V_Full))
          Z_smooth<-(S_Fullrun/SE_smooth)[2:ks,1]
          Z_smooth<-max(abs(Z_smooth-Z_pre))
        }else{Z_smooth<-0}
      }
      
      #name the columns
      colnames(S_Fullrun)<-c("SNP", colnames(S_LD))
      
      ##name rows like columns
      rownames(S_Fullrun)<-colnames(S_Fullrun)
      
      ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
      if(estimation == "DWLS"){
        test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf))
      }
      
      if(estimation == "ML"){
        test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML", sample.nobs = 200, optim.dx.tol = +Inf,sample.cov.rescale=FALSE))
      }
      
      test$warning$message[1]<-ifelse(is.null(test$warning$message), test$warning$message[1]<-0, test$warning$message[1])
      
      if(class(test$value)[1] == "lavaan" & grepl("solution has NOT",  as.character(test$warning)) != TRUE){
        Model_Output <- parTable(Model1_Results)
        
        resid_var1<-subset(Model_Output, Model_Output$op == "~~" & Model_Output$free != 0 & Model_Output$lhs == Model_Output$rhs)
        
        resid_var2<-min(resid_var1$est)}else{resid_var2<--9}
      
      if(resid_var2 > 0){
        
        #pull the delta matrix (this doesn't depend on N)
        S2.delt <- lavInspect(Model1_Results, "delta")
        
        ##weight matrix from stage 2
        S2.W <- lavInspect(Model1_Results, "WLS.V") 
        
        #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
        bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt, tol=toler) 
        
        #create the "lettuce" part of the sandwich
        lettuce <- S2.W%*%S2.delt
        
        #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
        Ohtt <- bread %*% t(lettuce)%*%V_Full_Reorder%*%lettuce%*%bread  
        
        #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
        SE <- as.matrix(sqrt(diag(Ohtt)))
        
        if(estimation == "DWLS"){
          #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
          if(":=" %in% Model_Output$op & !(NA %in% Model_Output$se)){
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
          }else{
            if(":=" %in% Model_Output$op & (NA %in% Model_Output$se)){
              se.ghost<-rep(NA, sum(":=" %in% Model_Output$op))
              warning("SE for ghost parameter could not be computed")
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
        
        #calculate model chi-square
        Eig<-as.matrix(eigen(V_Full)$values)
        Eig2<-diag((ncol(S_Fullrun)*(ncol(S_Fullrun)+1))/2)
        diag(Eig2)<-Eig
        
        #Pull P1 (the eigen vectors of V_eta)
        P1<-eigen(V_Full)$vectors
        
        implied<-as.matrix(fitted(Model1_Results))[1]
        implied_order<-colnames(S_Fullrun)
        implied[[1]]<-implied[[1]][implied_order,implied_order]
        implied2<-S_Fullrun-implied[[1]]
        eta<-as.vector(lowerTriangle(implied2,diag=TRUE))
        Q<-t(eta)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta
        
        ##remove parameter constraints, ghost parameters, and fixed effects from output to merge with SEs
        unstand<-subset(Model_Output, Model_Output$plabel != "" & Model_Output$free > 0)[,c(2:4,8,11,14)]
        
        ##combine ghost parameters with rest of output
        if(exists("ghost2") == "TRUE"){
          unstand2<-rbind(cbind(unstand,SE),ghost2)
        }else{unstand2<-cbind(unstand,SE)}
        
        ##add in fixed effects and parameter constraints to output
        other<-subset(Model_Output, (Model_Output$plabel == "" & Model_Output$op != ":=") | (Model_Output$free == 0 & Model_Output$plabel != ""))[,c(2:4,8,11,14)]
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
        
        ##add in model fit components to each row
        if(!(is.na(Q))){
          final$chisq<-rep(Q,nrow(final))
          final$chisq_df<-df
          final$chisq_pval<-pchisq(final$chisq,final$chisq_df,lower.tail=FALSE)
          final$AIC<-rep(Q + 2*npar,nrow(final))}else{final$chisq<-rep(NA, nrow(final))
          final$chisq_df<-rep(NA,nrow(final))
          final$chisq_pval<-rep(NA,nrow(final))
          final$AIC<-rep(NA, nrow(final))}
        
        ##add in error and warning messages 
        if(printwarn == TRUE){
          final$error<-ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
          final$warning<-ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]}
        
        ##combine results with SNP, CHR, BP, A1, A2 for particular model
        final2<-cbind(SNPs2[i,],final,row.names=NULL)
        
        if(smooth_check==TRUE){
          final2<-cbind(final2,Z_smooth)
        }
        
        if(!(sub[[1]])==FALSE){
          final2<-subset(final2, paste0(final2$lhs, final2$op, final2$rhs, sep = "") %in% sub)
          if(i == 1){
            Results_List<-vector(mode="list",length=nrow(final2))
            for(y in 1:nrow(final2)){
              Results_List[[y]]<-as.data.frame(matrix(NA,ncol=ncol(final2),nrow=f))
              colnames(Results_List[[y]])<-colnames(final2)
              Results_List[[y]][1,]<-final2[y,]
            }
          }else{
            for(y in 1:nrow(final2)){
              Results_List[[y]][i,]<-final2[y,]
            }
          }
          
        }else{##pull results and put into list object
          final2$est<-ifelse(final2$op == "<" | final2$op == ">" | final2$op == ">=" | final2$op == "<=", final2$est == NA, final2$est)
          Results_List[[i]]<-final2}
      }else{
        final<-data.frame(t(rep(NA, 13)))
        if(printwarn == TRUE){
          final$error<-ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
          if(resid_var2 != -9){
            final$error<-c("This particular run produced negative (residual) variances for either your latent or observed variables. You may discard the run for this SNP, re-run the model with constraints to keep variances above 0, or specify an alternative model.")
          }
          final$warning<-ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]}
        
        ##combine results with SNP, CHR, BP, A1, A2 for particular model
        final2<-cbind(SNPs2[i,],final,row.names=NULL)
        
        if(smooth_check==TRUE){
          final2<-cbind(final2,Z_smooth)
        }
        
        if(!(sub[[1]])==FALSE){
          final3<-as.data.frame(matrix(NA,ncol=ncol(final2),nrow=length(sub)))
          final3[1:length(sub),]<-final2[1,]
          if(smooth_check == TRUE){
            colnames(final3)<-c("SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC","error","warning","Z_smooth")
          }else{colnames(final3)<-c("SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC","error","warning")}
          
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
      
      if(i == 1){
        cat(paste0("Running Model: ", i, "\n"))
      }else{
        if(i %% 1000==0) {
          cat(paste0("Running Model: ", i, "\n"))
        }}
      
    }
    
    time_all<-proc.time()-time
    print(time_all[3])
    
    return(Results_List)
    
  }
  
  if(parallel == TRUE & Operating != "Windows"){
    
    if(is.null(cores)){
      ##if no default provided use 1 less than the total number of cores available so your computer will still function
      int <- detectCores() - 1
    }else{int<-cores}
    
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
    
    print("Starting GWAS Estimation")
    
    results<-foreach(n = icount(int), .combine = 'rbind') %:% 
      
      foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = "lavaan") %dopar% { 
        
        #create empty shell of V_SNP matrix
        V_SNP<-diag(k)
        
        #loop to add in the GWAS SEs, correct them for univariate and bivariate intercepts, and multiply by SNP variance from reference panel
        if(GC == "conserv"){
          for (p in 1:nrow(coords)) { 
            x<-coords[p,1]
            y<-coords[p,2]
            if (x != y) { 
              V_SNP[x,y]<-(SE_SNP[[n]][i,y]*SE_SNP[[n]][i,x]*I_LD[x,y]*I_LD[x,x]*I_LD[y,y]*varSNP[[n]][i]^2)}
            if (x == y) {
              V_SNP[x,x]<-(SE_SNP[[n]][i,x]*I_LD[x,x]*varSNP[[n]][i])^2
            }
          }
        }
        
        if(GC == "standard"){
          for (p in 1:nrow(coords)) { 
            x<-coords[p,1]
            y<-coords[p,2]
            if (x != y) { 
              V_SNP[x,y]<-(SE_SNP[[n]][i,y]*SE_SNP[[n]][i,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*varSNP[[n]][i]^2)}
            if (x == y) {
              V_SNP[x,x]<-(SE_SNP[[n]][i,x]*sqrt(I_LD[x,x])*varSNP[[n]][i])^2
            }
          }
        }
        
        if(GC == "none"){
          for (p in 1:nrow(coords)) { 
            x<-coords[p,1]
            y<-coords[p,2]
            if (x != y) { 
              V_SNP[x,y]<-(SE_SNP[[n]][i,y]*SE_SNP[[n]][i,x]*I_LD[x,y]*varSNP[[n]][i]^2)}
            if (x == y) {
              V_SNP[x,x]<-(SE_SNP[[n]][i,x]*varSNP[[n]][i])^2
            }
          }
        }
        
        if(smooth_check == TRUE){
          if(GC == "conserv"){
            Z_pre<-beta_SNP[[n]][i,]/(SE_SNP[[n]][i,]*diag(I_LD))
          }
          if(GC=="standard"){
            Z_pre<-beta_SNP[[n]][i,]/(SE_SNP[[n]][i,]*sqrt(diag(I_LD))) 
          }
          if(GC=="none"){
            Z_pre<-beta_SNP[[n]][i,]/SE_SNP[[n]][i,] 
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
        if(eigen(V_Full)$values[kv] <= 0){
          V_Full<-as.matrix((nearPD(V_Full, corr = FALSE))$mat)
          V_smooth<-1
        }
  
        #reorder sampling covariance matrix based on what lavaan expects given the specified model
        V_Full_Reorder <- V_Full[order,order]
        u<-nrow(V_Full_Reorder)
        W<-diag(u)
        diag(W)<-diag(V_Full_Reorder)
        
        ##invert the reordered sampling covariance matrix to create a weight matrix 
        if(toler==FALSE){
          W<- solve(W)
        }
        
        if(toler!=FALSE){
          W <- solve(W,tol=toler)
        }
        
        #create empty vector for S_SNP
        S_SNP<-vector(mode="numeric",length=k+1)
        
        #enter SNP variance from reference panel as first observation
        S_SNP[1]<-varSNP[[n]][i]
        
        #enter SNP covariances (standardized beta * SNP variance from refference panel)
        for (p in 1:k) {
          S_SNP[p+1]<-varSNP[[n]][i]*beta_SNP[[n]][i,p]
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
        if(eigen(S_Fullrun)$values[ks] <= 0){
          S_Fullrun<-as.matrix((nearPD(S_Fullrun, corr = FALSE))$mat)
          S_smooth<-1
        }
        
        if(smooth_check == TRUE){
          if(exists("S_smooth") | exists("V_smooth")){
            SE_smooth<-matrix(0, ks, ks)
            SE_smooth[lower.tri(SE_smooth,diag=TRUE)] <-sqrt(diag(V_Full))
            Z_smooth<-(S_Fullrun/SE_smooth)[2:ks,1]
            Z_smooth<-max(abs(Z_smooth-Z_pre))
          }else{Z_smooth<-0}
        }
        
        #name the columns
        colnames(S_Fullrun)<-c("SNP", colnames(S_LD))
        
        ##name rows like columns
        rownames(S_Fullrun)<-colnames(S_Fullrun)
        
        ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
        if(estimation == "DWLS"){
          test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf))
        }
        if(estimation == "ML"){
          test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML",sample.nobs = 200, optim.dx.tol = +Inf,sample.cov.rescale=FALSE))
        }
        
        test$warning$message[1]<-ifelse(is.null(test$warning$message), test$warning$message[1]<-0, test$warning$message[1])
        
        if(class(test$value)[1] == "lavaan" & grepl("solution has NOT",  as.character(test$warning)) != TRUE){
          Model_Output <- parTable(Model1_Results)
          
          resid_var1<-subset(Model_Output, Model_Output$op == "~~" & Model_Output$free != 0 & Model_Output$lhs == Model_Output$rhs)
          
          resid_var2<-min(resid_var1$est)}else{resid_var2<--9}
        
        if(resid_var2 > 0){
          
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
          if(estimation == "DWLS"){
            if(":=" %in% Model_Output$op & !(NA %in% Model_Output$se)){
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
            }else{
              if(":=" %in% Model_Output$op & (NA %in% Model_Output$se)){
                se.ghost<-rep(NA, sum(":=" %in% Model_Output$op))
                warning("SE for ghost parameter not available for ML")
                ghost<-subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
                ghost2<-cbind(ghost,se.ghost)
                colnames(ghost2)[7]<-"SE"}else{}}
          }
          
          #code for computing SE of ghost parameter (e.g., indirect effect in mediation model)
          if(estimation == "ML"){
            if(":=" %in% Model_Output$op){
              
              #pull the ghost parameter point estiamte
              ghost<-subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
              se.ghost<-rep(NA, sum(":=" %in% Model_Output$op))
              warning("SE for ghost parameter not available for ML")
              ##combine with delta method SE
              ghost2<-cbind(ghost,se.ghost)
              colnames(ghost2)[7]<-"SE"
            }else{
              if(":=" %in% Model_Output$op & (NA %in% Model_Output$se)){
                se.ghost<-rep(NA, sum(":=" %in% Model_Output$op))
                warning("SE for ghost parameter not available for ML")
                ghost<-subset(Model_Output, Model_Output$op == ":=")[,c(2:4,8,11,14)]
                ghost2<-cbind(ghost,se.ghost)
                colnames(ghost2)[7]<-"SE"}else{}}
          }
          
          #calculate model chi-square
          Eig<-as.matrix(eigen(V_Full)$values)
          Eig2<-diag((ncol(S_Fullrun)*(ncol(S_Fullrun)+1))/2)
          diag(Eig2)<-Eig
          
          #Pull P1 (the eigen vectors of V_eta)
          P1<-eigen(V_Full)$vectors
          
          implied<-as.matrix(fitted(Model1_Results))[1]
          implied_order<-colnames(S_Fullrun)
          implied[[1]]<-implied[[1]][implied_order,implied_order]
          implied2<-S_Fullrun-implied[[1]]
          eta<-as.vector(lowerTriangle(implied2,diag=TRUE))
          Q<-t(eta)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta
          
          ##remove parameter constraints, ghost parameters, and fixed effects from output to merge with SEs
          unstand<-subset(Model_Output, Model_Output$plabel != "" & Model_Output$free > 0)[,c(2:4,8,11,14)]
          
          ##combine ghost parameters with rest of output
          if(exists("ghost2") == "TRUE"){
            unstand2<-rbind(cbind(unstand,SE),ghost2)
          }else{unstand2<-cbind(unstand,SE)}
          
          ##add in fixed effects and parameter constraints to output
          other<-subset(Model_Output, (Model_Output$plabel == "" & Model_Output$op != ":=") | (Model_Output$free == 0 & Model_Output$plabel != ""))[,c(2:4,8,11,14)]
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
          
          ##add in model fit components to each row
          if(!(is.na(Q))){
            final$chisq<-rep(Q,nrow(final))
            final$chisq_df<-df
            final$chisq_pval<-pchisq(final$chisq,final$chisq_df,lower.tail=FALSE)
            final$AIC<-rep(Q + 2*npar,nrow(final))}else{final$chisq<-rep(NA, nrow(final))
            final$chisq_df<-rep(NA,nrow(final))
            final$chisq_pval<-rep(NA,nrow(final))
            final$AIC<-rep(NA, nrow(final))}
          
          ##add in error and warning messages 
          if(printwarn == TRUE){
            final$error<-ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
            final$warning<-ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]}
          
          ##combine with rs-id, BP, CHR, etc.
          final2<-cbind(i,n,SNPs2[[n]][i,],final,row.names=NULL)
          
          if(smooth_check==TRUE){
            final2<-cbind(final2,Z_smooth)
          }
          
          if(!(sub[[1]])==FALSE){
            final2<-subset(final2, paste0(final2$lhs, final2$op, final2$rhs, sep = "") %in% sub)
          }else{##pull results 
            final2$est<-ifelse(final2$op == "<" | final2$op == ">" | final2$op == ">=" | final2$op == "<=", final2$est == NA, final2$est)}
          
          ##results to be put into the output
          final2
          
        }else{ 
          
          final<-data.frame(t(rep(NA, 13)))
          if(printwarn == TRUE){
            final$error<-ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
            if(resid_var2 != -9){
              final$error<-c("This particular run produced negative (residual) variances for either your latent or observed variables. You may discard the run for this SNP, re-run the model with constraints to keep variances above 0, or specify an alternative model.")
            }
            final$warning<-ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]}
          
          ##combine results with SNP, CHR, BP, A1, A2 for particular model
          final2<-cbind(i,n,SNPs2[[n]][i,],final,row.names=NULL)
          if(smooth_check==TRUE){
            final2<-cbind(final2,Z_smooth)
            colnames(final2)<-c("i", "n", "SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC","error","warning","Z_smooth")
          }else{colnames(final2)<-c("i", "n", "SNP", "CHR", "BP", "MAF", "A1", "A2", "lhs", "op", "rhs", "free", "label", "est", "SE", "Z_Estimate", "Pval_Estimate","chisq","chisq_df","chisq_pval", "AIC","error","warning")}
        
          final2
        }
        
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
  
  if(parallel == TRUE & Operating == "Windows"){
    stop("Parallel processing is not currently available for Windows operating systems. Please set the parallel argument to FALSE, or switch to a Linux or Mac operating system.")
  }
  
}
