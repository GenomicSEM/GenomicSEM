write.model<-function(Loadings,S_LD,cutoff,fix_resid=TRUE,bifactor=FALSE,mustload=FALSE,common=FALSE){
  Model<-""
  if(common == TRUE){
    for(f in 1){
      u<-1
      Model1<-""
      for(i in 1:nrow(S_LD)){
        if(u == 1){
          linestart<-paste("F", f, "=~",  colnames(S_LD)[i], sep = "")
          u<-u+1
          linemid<-""
        }else{
          linemid<-paste(linemid, " + ", colnames(S_LD)[i], sep = "")
        }
      }
    }
    Model<-paste(Model,linestart, linemid, " \n ", sep="")
  }else{
    if(mustload == TRUE){
      Mins<-apply(abs(Loadings),1,max)
      for(i in 1:nrow(Loadings)){
        for(f in 1:ncol(Loadings)){
          if(Mins[i] == abs(Loadings[i,f])){
            Loadings[i,f]<-cutoff+.01
          }
        }
      }
    }
    
    
    for(f in 1:ncol(Loadings)){
      u<-1
      Model1<-""
      for(i in 1:nrow(Loadings)){
        if(abs(Loadings[i,f]) > cutoff){
          if(u == 1){
            linestart<-paste("F", f, "=~",  colnames(S_LD)[i], sep = "")
            u<-u+1
            linemid<-""
          }else{
            linemid<-paste(linemid, " + ", colnames(S_LD)[i], sep = "")
          }
        }
      }
      Model<-paste(Model,linestart, linemid, " \n ", sep="")
      linestart<-""
      linemid<-""
    }
    
    if(bifactor==TRUE){
      Model_bi<-""
      u<-1
      for(i in 1:ncol(S_LD)){
        b<-grepl(colnames(S_LD)[i], Model)
        if(b == TRUE){
          if(u == 1){
            linestart_bi<-paste("Common_F", "=~",  colnames(S_LD)[i], sep = "")
            u<-u+1
            linemid_bi<-""
          }else{
            linemid_bi<-paste(linemid_bi, " + ", colnames(S_LD)[i], sep = "")
          }
        }
      }
      Model_bi<-paste(linestart_bi,linemid_bi," \n ", sep="")
      r<-1
      Factor_bi<-""
      for(i in 1:ncol(Loadings)){
        Factor_bi_start<-paste("Common_F~~0*", "F", i, " \n ", sep = "")
        Factor_bi<-paste(Factor_bi,Factor_bi_start,sep="")
      }
      
      Modelsat<-""
      for (i in 1:(ncol(Loadings)-1)) {
        linestartc <- paste("", "F", i, "~~0*F", i+1, sep = "")
        if (ncol(Loadings)-i >= 2) {
          linemidc <- ""
          for (j in (i+2):ncol(Loadings)) {
            linemidc <- paste("", linemidc, "F", i, "~~0*F", j, " \n ", sep="")
            
          }
        } else {linemidc <- ""}
        Modelsat <- paste(Modelsat, linestartc, " \n ", linemidc, sep = "")
      } 
      
      Model<-paste(Model,Model_bi,Factor_bi,Modelsat)
    }
    
  
  }
  
  if(fix_resid==TRUE){
    Model3<-""
    #create unique combination of letters for residual variance parameter labels
    n<-combn(letters,4)[,sample(1:14000, ncol(S_LD), replace=FALSE)]
    for(i in 1:ncol(S_LD)){
      if(grepl(colnames(S_LD)[i],Model) == TRUE){
        linestart3a <- paste(colnames(S_LD)[i], " ~~ ",  paste(n[,i],collapse=""), "*", colnames(S_LD)[i], sep = "")
        linestart3b <- paste(paste(n[,i],collapse=""), " > .0001", sep = "")
        Model3<-paste(Model3, linestart3a, " \n ", linestart3b, " \n ", sep = "")
      }
    }
    Model<-paste(Model,Model3,sep="")
  }
  
  
  return(Model)
}
