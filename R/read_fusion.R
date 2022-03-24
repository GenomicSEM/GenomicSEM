
read_fusion <- function(files,trait.names=NULL,binary=NULL,N=NULL,perm=FALSE){
  
  length <- length(files)
  
  if(is.null(trait.names)){
    
    names.beta <- paste0("beta.",1:length)
    names.se <- paste0("se.",1:length)
    
  }else{
    
    names.beta <- paste0("beta.",trait.names)
    names.se <- paste0("se.",trait.names)
    
  }
  
  files = lapply(files, read.table, header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))

  print("Done reading in FUSION files")
  
  if(is.null(binary)){
    print("Running read.fusion assuming all traits being analyzed are binary.
          If this is not true please use the binary argument to specify which traits and binary or continuous.")
    binary<-rep(TRUE,length)
  }
  
  for(i in 1:length){
    
    if(perm == TRUE){
    #set p-values of 1 to (1-permutations)/permutations = upper bound of possible empirical p-value estimates.
    #dont leave at 1 so backing out SE is possible below
    files[[i]]$PERM.PV<-ifelse(files[[i]]$PERM.PV == 1, (files[[i]]$PERM.N-1)/files[[i]]$PERM.N,files[[i]]$PERM.PV)
    
    #set to lowest possible empirical p-value given number of iteratations if empirical p = 0
    files[[i]]$PERM.PV<-ifelse(files[[i]]$PERM.PV ==0, 1/files[[i]]$PERM.N,files[[i]]$PERM.PV)
    
    files[[i]]$Z.perm <- sign(files[[i]]$TWAS.Z) * sqrt(qchisq(files[[i]]$PERM.PV,1,lower=F))
    }
    
    if(binary[i] == TRUE){
      
      if(perm == FALSE){
      files[[i]]$effect <- files[[i]]$TWAS.Z /sqrt((N[i]/4) *files[[i]]$HSQ)
      files[[i]]$SE<-1/sqrt((N[i]/4)*files[[i]]$HSQ)
      }
      
      if(perm == TRUE){
        files[[i]]$effect <- files[[i]]$Z.perm/sqrt((N[i]/4) *files[[i]]$HSQ)
        files[[i]]$SE<-abs(files[[i]]$effect/files[[i]]$Z.perm)
        
      }
      
      output <- cbind.data.frame(files[[i]]$FILE,files[[i]]$ID,
                                 (files[[i]]$effect)/((files[[i]]$effect^2) * files[[i]]$HSQ + (pi^2)/3)^.5,
                                 (files[[i]]$SE)/(((files[[i]]$effect)^2) * files[[i]]$HSQ + (pi^2)/3)^.5)  

    }
    
    if(binary[i] == FALSE){
      
      if(perm == FALSE){
        files[[i]]$effect <- files[[i]]$TWAS.Z /sqrt(N[i] * files[[i]]$HSQ)
        
        output <- cbind.data.frame(files[[i]]$FILE,files[[i]]$ID,
                                   files[[i]]$effect,
                                   abs(files[[i]]$effect/files[[i]]$TWAS.Z))}
      if(perm == TRUE){
        files[[i]]$effect <- files[[i]]$Z.perm /sqrt(N[i] * files[[i]]$HSQ)
        output <- cbind.data.frame(files[[i]]$FILE,files[[i]]$ID,
                                   files[[i]]$effect,
                                   abs(files[[i]]$effect/files[[i]]$Z.perm))
      }
      
    }
    
    colnames(output) <- c("Panel","Gene",names.beta[i],names.se[i])
    
    #remove weights folder naming of Panel automatically added by FUSION
    output$Panel<-gsub(".*//","", output$Panel)
 
    if(i ==1){
      data.frame.out <- cbind.data.frame(files[[i]]$HSQ,output)
      
      colnames(data.frame.out ) <- c("HSQ","Panel","Gene",names.beta[i],names.se[i])
      data.frame.out <- na.omit(data.frame.out)
      
    }else{
      data.frame.out <- merge(data.frame.out,output,by=c("Gene","Panel"),all.x=F,all.y=F) 
      data.frame.out <- na.omit(data.frame.out)
    }
    
    
  }
  data.frame.out <- unique(data.frame.out)
  data.frame.out
  
}
