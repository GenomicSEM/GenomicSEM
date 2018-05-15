

userGWAS<-function(Output,estimation="DWLS",model=""){ 
  time<-proc.time()
  
  ##split the V and S matrices into as many (cores - 1) as are aviailable on the local computer
  V_Full<-(Output[[1]])
  S_Full<-(Output[[2]])
  
  #enter in k for number of phenotypes 
  k<-ncol(S_Full[[1]])-1
  
  f<-length(Output[[1]])
  
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
  
  ##run one model that specifies the factor structure so that lavaan knows how to rearrange the V (i.e., sampling covariance) matrix
  for (i in 1) {
    
    #transform sampling covariance matrix into a weight matrix: 
    W <- solve(V_Full[[i]])
    
    S_Fullrun<-S_Full[[i]]
    
    test2<-tryCatch.W.E(ReorderModel <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2))
    
    order <- rearrange(k = k+1, fit = ReorderModel, names = rownames(S_Full[[1]]))
  }
  
  #make empty list object for model results
  Results_List<-vector(mode="list",length=f)
  
  ##estimation for WLS
  if(estimation=="DWLS"){

    for (i in 1:f) { 
      
      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Full_Reorder <- V_Full[[i]][order,order]
      u<-nrow(V_Full_Reorder)
      V_Full_Reorderb<-diag(u)
      diag(V_Full_Reorderb)<-diag(V_Full_Reorder)
      
      ##invert the reordered sampling covariance matrix to create a weight matrix 
      W <- solve(V_Full_Reorderb) 
      
      #import the S_Full matrix for appropriate run
      S_Fullrun<-S_Full[[i]]
      
      ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
      test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "DWLS", WLS.V = W, sample.nobs = 2))
      
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
     
      ifelse(":=" %in% Model_WLS$op & !(NA %in% Model_WLS$se), 1, 0)
      
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
        ghost<-subset(Model_WLS, Model_WLS$op == ":=")[,c(2:4,8,,11,14)]
        ghost2<-cbind(ghost,se.ghost)
        colnames(ghost2)[7]<-"SE"}else{}} 
        
      ##remove parameter constraints, ghost parameters, and fixed effects from output to merge with SEs
      unstand<-subset(Model_WLS, Model_WLS$plabel != "" & Model_WLS$free > 0)[,c(2:4,8,11,14)]
      
      ##combine ghost parameters with rest of output
      if(exists("ghost2") == "TRUE"){
        unstand2<-rbind(cbind(unstand,SE),ghost2)
      }else{unstand2<-cbind(unstand,SE)}
      
      ##add in fixed effects and parameter constraints to output
      other<-subset(Model_WLS, (Model_WLS$plabel == "" & Model_WLS$op != ":=") | (Model_WLS$free == 0 & Model_WLS$plabel != ""))[,c(2:4,8,14)]
      other$SE<-rep(NA, nrow(other))
     
      ##combine fixed effects and parameter constraints with output if there are any
      if(nrow(other) > 0){
      final<-rbind(unstand2,other)
      }else{final<-unstand2}
      
      #reorder based on row numbers so it is in order the user provided
      final$index <- as.numeric(row.names(final))
      final<-final[order(final$index), ]
      final$index<-NULL
      
      ##add in error and warning messages 
      final$error<-ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
      final$warning<-ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]
      
      ##combine results with SNP, CHR, BP, A1, A2 for particular model
      final2<-cbind(Output[[3]][i,],final,row.names=NULL)
       
      ##pull results and put into list object
      Results_List[[i]]<-final2
      print(i)

    }
  }
  
  ##ML estimation
  if(estimation=="ML"){
    
    for (i in 1:f) { 
      
      #reorder sampling covariance matrix based on what lavaan expects given the specified model
      V_Full_Reorder <- V_Full[[i]][order,order]
      
      #import the S_Full matrix for appropriate run
      S_Fullrun<-S_Full[[i]]
      
      ##run the model. save failed runs and run model. warning and error functions prevent loop from breaking if there is an error. 
      test<-tryCatch.W.E(Model1_Results <- sem(Model1, sample.cov = S_Fullrun, estimator = "ML", sample.nobs = 200))
      
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
      if(":=" %in% Model_ML$op & !(NA %in% Model_ML$se)){
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
        ghost<-subset(Model_ML, Model_ML$op == ":=")[,c(2:4,8,11,14)]
        
        ##combine with delta method SE
        ghost2<-cbind(ghost,se.ghost)
        colnames(ghost2)[7]<-"SE"
      }else{
        if(":=" %in% Model_ML$op & (NA %in% Model_ML$se)){
          se.ghost<-rep("SE could not be computed", count(":=" %in% Model_ML$op)$freq)
          ghost<-subset(Model_ML, Model_ML$op == ":=")[,c(2:4,8,11,14)]
          ghost2<-cbind(ghost,se.ghost)
          colnames(ghost2)[7]<-"SE"}else{}}
      
      ##remove parameter constraints, ghost parameters, and fixed effects from output to merge with SEs
      unstand<-subset(Model_ML, Model_ML$plabel != "" & Model_ML$free > 0)[,c(2:4,8,11,14)]
      
      ##combine ghost parameters with rest of output
      if(exists("ghost2") == "TRUE"){
        unstand2<-rbind(cbind(unstand,SE),ghost2)
      }else{unstand2<-cbind(unstand,SE)}
      
      ##add in fixed effects and parameter constraints to output
      other<-subset(Model_ML, (Model_ML$plabel == "" & Model_ML$op != ":=") | (Model_ML$free == 0 & Model_ML$plabel != ""))[,c(2:4,8,14)]
      other$SE<-rep(NA, nrow(other))
      
      ##combine fixed effects and parameter constraints with output if there are any
      if(nrow(other) > 0){
        final<-rbind(unstand2,other)
      }else{final<-unstand2}
      
      #reorder based on row numbers so it is in order the user provided
      final$index <- as.numeric(row.names(final))
      final<-final[order(final$index), ]
      final$index<-NULL
      
      ##add in error and warning messages 
      final$error<-ifelse(class(test$value) == "lavaan", 0, as.character(test$value$message))[1]
      final$warning<-ifelse(class(test$warning) == 'NULL', 0, as.character(test$warning$message))[1]
      
      ##combine results with SNP, CHR, BP, A1, A2 for particular model
      final2<-cbind(Output[[3]][i,],final,row.names=NULL)
      
      ##pull results and put into list object
      Results_List[[i]]<-final2
      print(i)
    }
  }
  
  time_all<-proc.time()-time
  print(time_all[3])
  
  return(Results_List)
  
}
