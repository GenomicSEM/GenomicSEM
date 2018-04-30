read.fusion <- function(files,trait.names=NULL,OLS=FALSE,linprob=FALSE,prop=FALSE,N=NA){
  
  length <- length(files)
  
  if(is.null(trait.names)){
    
    names.beta <- paste0("beta.",1:length)
    names.se <- paste0("se.",1:length)
    
  }else{
    
    names.beta <- paste0("beta.",trait.names)
    names.se <- paste0("se.",trait.names)
    
  }
  
  
  files = lapply(files, read.table, header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
  
  

  
  for(i in 1:length){


    if(linprob[i] == T){
      
       files[[i]]$effect <- files[[i]]$TWAS.Z /sqrt((prop[i]*(1-prop[i])*N[i] *files[[i]]$HSQ))
      files[[i]]$SE<-1/sqrt((prop[i]*(1-prop[i])*(N[i]*files[[i]]$HSQ)))}
    
    output <- cbind.data.frame(files[[i]]$ID,
                               files[[i]]$effect,
                               files[[i]]$SE )
    
    colnames(output) <- c("Gene",names.beta[i],names.se[i])
  
  }

    
    if(OLS[i] == T){
      
      
      files[[i]]$effect <- files[[i]]$TWAS.Z / sqrt(N[i] * files[[i]]$HSQ )
      files[[i]]$SE
      
      output <- cbind.data.frame(files[[i]]$ID,
                                 files[[i]]$effect,
                                 abs(files[[i]]$effect/files[[i]]$TWAS.Z ))
                    
      colnames(output) <- c("Gene",names.beta[i],names.se[i])                           
      
      
    }
    
    
    if(linprob[i] == F){
      if(OLS[i] == F){                                     
       
        # placeholder, effective N from Willier et al. 2010.
        N[i] <-4/(1/(prop[i]*N[i]) + 1/((1-prop[i])*N[i]))
        
        files[[i]]$effect <- files[[i]]$TWAS.Z /sqrt(N[i] *files[[i]]$HSQ)
        files[[i]]$SE <-1/sqrt(N[i]*files[[i]]$HSQ)
         
        
          output <- cbind.data.frame(files[[i]]$ID,
                                 files[[i]]$effect,
                                 files[[i]]$SE )
          
          colnames(output) <- c("Gene",names.beta[i],names.se[i])}}
    
  
    
    if(i ==1){
      data.frame.out <- cbind.data.frane(output,files[[i]]$HSQ)
      colnames(data.frame.out ) <- c("Gene",names.beta[i],names.se[i],"HSQ")
      data.frame.out <- na.omit(data.frame.out)
      
    }else{
      data.frame.out <- merge(data.frame.out,output,by=1,all.x=F,all.y=F) 
      data.frame.out <- na.omit(data.frame.out)
    }
    
  }
  
  data.frame.out
  
}
