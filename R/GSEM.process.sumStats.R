GSEM.process.sumStats <-function(sumstats,type){
  
  if(type="OLS"){
  
cat("Setting unified headers...") 
x<-lapply(sumstats,setNames,nm=c("SNP","CHR","BP","A1","A2","EAF","Beta","se","Z","N","PVAL"))
cat("done\n")

cat("Removing duplicate entries...")
x<-lapply(x,function(x){return(x[duplicated(x[,2])==F,])})
cat("done\n")

cat("Changing alleles to upper-case...")
x<-lapply(x,function(x){
  dum<-cbind(x[,1:4],apply(x[,5:6],2,toupper),x[,7:ncol(x)])
  dum$A1<-as.character(dum$A1)
  dum$A2<-as.character(dum$A2)
  return(dum)
}
)
cat("done\n")

cat("Extracting the intersection of SNPs in x...")
extract_SNPs<-function(x,column){
  SNPlist.all<-unique(unlist(lapply(x,"[[",column)))
  SNPlist.mat<-as.data.frame(matrix(0,nrow=length(SNPlist.all),ncol=length(x),dimnames=list(SNPlist.all,names(x))),stringsAsFactors=F)
  for(i in 1:ncol(SNPlist.mat)){
    SNPlist.mat[rownames(SNPlist.mat)%in%x[[i]][,column],i]<-1
  }
  SNPlist.mat$TOTAL<-rowSums(SNPlist.mat)
  return(SNPlist.mat)
}
SNPlist.mat<-extract_SNPs(x,2)
SNPs<-data.frame(MarkerName=rownames(SNPlist.mat[SNPlist.mat$TOTAL==length(x),]),stringsAsFactors=F)
cat("done\n",format(nrow(SNPs),big.mark=",")," SNPs were found in all files",sep="")
#willen we dit ook opslaan (i.e. later komt dit nog een keer voor de SNPs die echt in het model gaan)?
if(save.matrix){
  cat("\nAn overview is saved at ",out,".N_weighted_GWAMA.SNP_list.matrix",sep="")
  SNPlist.mat.out<-as.data.frame(cbind(rownames(SNPlist.mat),as.data.frame(SNPlist.mat)))
  colnames(SNPlist.mat.out)<-c("SNP",colnames(SNPlist.mat))
  if(output.gz){
    cat(".gz\n")
    gz1<-gzfile(paste0(out,".N_weighted_GWAMA.SNP_list.matrix.gz"),"w")
  }else{
    cat("\n")
    gz1<-paste0(out,".N_weighted_GWAMA.SNP_list.matrix")
  }
  write.table(x=SNPlist.mat.out,file=gz1,quote=F,sep="\t",row.names=F)
  if(is.character(gz1)==F){
    close(gz1)
  }
  rm(SNPlist.mat.out,gz1)
  gc()
}
cat("\nKeeping only common SNPs between the different inputs...")
x_common<-lapply(x,function(x,y){return(x[match(y$MarkerName,x$MarkerName),])},y=SNPs)
cat("done\n")
rm(x,SNPs,size.range,SNPlist.mat)
gc()

#SANITY CHECKS
cat("Checking whether alleles are aligned between the inputs...")
A1A2_file1<-subset(x_common[[1]],select=c(A1,A2))
colnames(A1A2_file1)[1]<-"A1.ref"
colnames(A1A2_file1)[2]<-"A2.ref"
cat("done (NOTE: The first input in x was used as reference file)\n")

#Check whether the alleles of all data frames in your list are aligned properly
identicalYesNo<-lapply(x_common,function(x,y){identical(x$A1,y$A1.ref) & identical(x$A2,y$A2.ref)},y=A1A2_file1)
print(identicalYesNo)

if(sum(unlist(identicalYesNo))!=length(x_common)){
  cat("Unaligned SNPs were found\n")
  #cbind with the "ref file"-> first data frame in list will be considered the ref allele
  x_reffile<-lapply(x_common,function(x){cbind(x,A1A2_file1)})
  rm(x_common,A1A2_file1)
  gc()
  
  #Check which alleles are aligned between files and which alleles are not aligned
  aligned<-lapply(x_reffile,function(x){subset(x,A1==A1.ref & A2==A2.ref)})
  
  flipped<-lapply(x_reffile,function(x){subset(x,A1==A2.ref & A2==A1.ref)})
  N_flipped<-sapply(flipped,nrow)
  cat("How many markers were flipped?\n")
  print(noquote(format(N_flipped,big.mark=",")))
  
  if(sum(unlist(N_flipped))!=0){
    cat("Flipping Z-scores for SNPs that are not aligned...")
    Zlist<-lapply(flipped,"[","Z")
    Zlist<-lapply(Zlist,function(x){x$Z * -1})
    #make data frame again
    Zlist<-lapply(Zlist,function(x){as.data.frame(x)})
    #cbind newbetas with not aligned file
    flipped_new_Z<-Map(cbind,flipped,Zlist) 
    #format right column names
    flipped_new_Z<-lapply(flipped_new_Z,setNames,nm=c("SNP","CHR","BP","A1","A2","EAF","Beta","se","Z_NOTUSE","N","PVAL","A1.ref","A2.ref","Z"))
    cat("done\n")
    
    cat("Flipping Beta's for SNPs that are not aligned...")
    Blist<-lapply(flipped,"[","Beta")
    Blist<-lapply(Blist,function(x){x$Beta * -1})
    #make data frame again
    Blist<-lapply(Blist,function(x){as.data.frame(x)})
    #cbind newbetas with not aligned file
    flipped_new_Z<-Map(cbind,flipped_new_Z,Blist) 
    #format right column names
    flipped_new_Z<-lapply(flipped_new_Z,setNames,nm=c("SNP","CHR","BP","A1","A2","EAF","Beta_NOTUSE","se","Z_NOTUSR","N","PVAL","A1.ref","A2.ref","Z","Beta"))
    cat("done\n")
    
    cat("Flipping EAF for SNPs that are not aligned...")
    EAFlist<-lapply(flipped,"[","EAF")
    EAFlist<-lapply(EAFlist,function(x){1-x$EAF})
    EAFlist<-lapply(EAFlist,function(x){as.data.frame(x)})
    flipped_new_Z_EAF<-Map(cbind,flipped_new_Z,EAFlist) 
    flipped_new_Z_EAF<-lapply(flipped_new_Z_EAF,setNames,nm=c("SNP","CHR","BP","A1","A2","EAF_NOTUSE","Beta_NOTUSE","se","Z_NOTUSE","N","PVAL","A1.ref","A2.ref","Z","BEta","EAF"))
    cat("done\n")
    
    rm(flipped_new_Z)
    gc()
    
    #select right columns for rbind with data frames in list that are not complete due to misalignment
    flipped_reorder<-lapply(flipped_new_Z_EAF,function(x){subset(x,select=(c(SNP,MarkerName,CHR,BP,A2,A1,EAF,N,Z,PVAL,A1.ref,A2.ref)))})
    flipped_reorder<-lapply(flipped_reorder,setNames,nm=c("SNP","MarkerName","CHR","BP","A1","A2","EAF","N","Z","PVAL","A1.ref","A2.ref"))
    
    rm(flipped_new_Z_EAF)
    gc()
    
    cat("Combining the newly aligned SNPs with SNPs that were already combined...")
    x_common_aligned<-Map(rbind,aligned,flipped_reorder)
    cat("done\n")
    
    rm(aligned,flipped_reorder,Zlist,EAFlist)
    gc()
  }
  rm(flipped,N_flipped)
  gc()
  
  if(exists("x_common_aligned")==F){
    x_common_aligned<-aligned
    rm(aligned)
    gc()
  }
  
  all_aligned<-lapply(x_common_aligned,function(x){subset(x,A1==A1.ref & A2==A2.ref)})
  rm(x_common_aligned)
  gc()
  
  unalignable<-lapply(x_reffile,function(x){subset(x,(A1!=A1.ref & A1!=A2.ref)|(A2!=A1.ref & A2!=A2.ref))})
  N_unalignable<-sapply(unalignable,nrow)
  if(sum(N_unalignable)!=0){
    cat("On rare occasions, SNPs cannot be aligned (e.g. A-C in the reference file vs A-T in the input). This might lead to summary statistics of different length\n")
    cat("How many markers could not be aligned?\n")
    print(noquote(format(N_unalignable,big.mark=",")))
    cat("These were removed from the data frame(s)\nA list was saved at ",out,".not_aligned",sep="")
    unalignable_SNPs<-unique(unlist(lapply(unalignable,"[[","MarkerName")))
    if(output.gz){
      cat(".gz\n")
      gz1<-gzfile(paste0(out,".not_aligned.gz"),"w")
    }else{
      cat("\n")
      gz1<-paste0(out,".not_aligned")
    }
    write.table(x=data.frame(MarkerName=unalignable_SNPs,stringsAsFactors=F),file=gz1,quote=F,sep="\t",row.names=F)
    if(is.character(gz1)==F){
      close(gz1)
    }
    rm(unalignable_SNPs,gz1)
    gc()
  }
  rm(unalignable,N_unalignable,x_reffile)
  gc()
  
  cat("Selecting only SNPs that are common again and extracting them from the data frames...")
  SNPlist.mat<-extract_SNPs(all_aligned,2)
  SNPs<-data.frame(MarkerName=rownames(SNPlist.mat[SNPlist.mat$TOTAL==length(all_aligned),]),stringsAsFactors=F)
  all_aligned<-lapply(all_aligned,function(x,y){x[match(y$MarkerName,x$MarkerName),]},y=SNPs)
  cat("done\n")
  rm(SNPlist.mat,SNPs)
  gc()
}else{
  cat("All SNPs were aligned\n")
  all_aligned<-x_common
  rm(x_common,A1A2_file1)
  gc()
}
rm(identicalYesNo)
gc()

cat(format(nrow(all_aligned[[1]]),big.mark=",")," SNPs will be used in the N-weighted GWAMA\n",sep="")
if(save.matrix){
  cat("A list was saved at ",out,".N_weighted_GWAMA.input",sep="")
  if(output.gz){
    cat(".gz\n")
    gz1<-gzfile(paste0(out,".N_weighted_GWAMA.input.gz"),"w")
  }else{
    cat("\n")
    gz1<-paste0(out,".N_weighted_GWAMA.input")
  }		
  write.table(x=all_aligned[[1]][,1],file=gz1,quote=F,sep="\t",row.names=F,col.names=c("MarkerName"))
  if(is.character(gz1)==F){
    close(gz1)
  }
  rm(gz1)
}


cat("Extracting Z-scores from aligned summary statistics...")
Zlist<-lapply(all_aligned,"[","Z")
cat("done\n")

  }
  
  if(type="Logistic"){
    cat("Setting unified headers...") 
    x<-lapply(sumstats,setNames,nm=c("SNP","CHR","BP","A1","A2","EAF","Beta","se","Z","N","PVAL"))
    cat("done\n")
    
    cat("Removing duplicate entries...")
    x<-lapply(x,function(x){return(x[duplicated(x[,2])==F,])})
    cat("done\n")
    
    cat("Changing alleles to upper-case...")
    x<-lapply(x,function(x){
      dum<-cbind(x[,1:4],apply(x[,5:6],2,toupper),x[,7:ncol(x)])
      dum$A1<-as.character(dum$A1)
      dum$A2<-as.character(dum$A2)
      return(dum)
    }
    )
    cat("done\n")
    
    cat("Extracting the intersection of SNPs in x...")
    extract_SNPs<-function(x,column){
      SNPlist.all<-unique(unlist(lapply(x,"[[",column)))
      SNPlist.mat<-as.data.frame(matrix(0,nrow=length(SNPlist.all),ncol=length(x),dimnames=list(SNPlist.all,names(x))),stringsAsFactors=F)
      for(i in 1:ncol(SNPlist.mat)){
        SNPlist.mat[rownames(SNPlist.mat)%in%x[[i]][,column],i]<-1
      }
      SNPlist.mat$TOTAL<-rowSums(SNPlist.mat)
      return(SNPlist.mat)
    }
    SNPlist.mat<-extract_SNPs(x,2)
    SNPs<-data.frame(MarkerName=rownames(SNPlist.mat[SNPlist.mat$TOTAL==length(x),]),stringsAsFactors=F)
    cat("done\n",format(nrow(SNPs),big.mark=",")," SNPs were found in all files",sep="")
    #willen we dit ook opslaan (i.e. later komt dit nog een keer voor de SNPs die echt in het model gaan)?
    if(save.matrix){
      cat("\nAn overview is saved at ",out,".N_weighted_GWAMA.SNP_list.matrix",sep="")
      SNPlist.mat.out<-as.data.frame(cbind(rownames(SNPlist.mat),as.data.frame(SNPlist.mat)))
      colnames(SNPlist.mat.out)<-c("SNP",colnames(SNPlist.mat))
      if(output.gz){
        cat(".gz\n")
        gz1<-gzfile(paste0(out,".N_weighted_GWAMA.SNP_list.matrix.gz"),"w")
      }else{
        cat("\n")
        gz1<-paste0(out,".N_weighted_GWAMA.SNP_list.matrix")
      }
      write.table(x=SNPlist.mat.out,file=gz1,quote=F,sep="\t",row.names=F)
      if(is.character(gz1)==F){
        close(gz1)
      }
      rm(SNPlist.mat.out,gz1)
      gc()
    }
    cat("\nKeeping only common SNPs between the different inputs...")
    x_common<-lapply(x,function(x,y){return(x[match(y$MarkerName,x$MarkerName),])},y=SNPs)
    cat("done\n")
    rm(x,SNPs,size.range,SNPlist.mat)
    gc()
    
    #SANITY CHECKS
    cat("Checking whether alleles are aligned between the inputs...")
    A1A2_file1<-subset(x_common[[1]],select=c(A1,A2))
    colnames(A1A2_file1)[1]<-"A1.ref"
    colnames(A1A2_file1)[2]<-"A2.ref"
    cat("done (NOTE: The first input in x was used as reference file)\n")
    
    #Check whether the alleles of all data frames in your list are aligned properly
    identicalYesNo<-lapply(x_common,function(x,y){identical(x$A1,y$A1.ref) & identical(x$A2,y$A2.ref)},y=A1A2_file1)
    print(identicalYesNo)
    
    if(sum(unlist(identicalYesNo))!=length(x_common)){
      cat("Unaligned SNPs were found\n")
      #cbind with the "ref file"-> first data frame in list will be considered the ref allele
      x_reffile<-lapply(x_common,function(x){cbind(x,A1A2_file1)})
      rm(x_common,A1A2_file1)
      gc()
      
      #Check which alleles are aligned between files and which alleles are not aligned
      aligned<-lapply(x_reffile,function(x){subset(x,A1==A1.ref & A2==A2.ref)})
      
      flipped<-lapply(x_reffile,function(x){subset(x,A1==A2.ref & A2==A1.ref)})
      N_flipped<-sapply(flipped,nrow)
      cat("How many markers were flipped?\n")
      print(noquote(format(N_flipped,big.mark=",")))
      
      if(sum(unlist(N_flipped))!=0){
        cat("Flipping Z-scores for SNPs that are not aligned...")
        Zlist<-lapply(flipped,"[","Z")
        Zlist<-lapply(Zlist,function(x){x$Z * -1})
        #make data frame again
        Zlist<-lapply(Zlist,function(x){as.data.frame(x)})
        #cbind newbetas with not aligned file
        flipped_new_Z<-Map(cbind,flipped,Zlist) 
        #format right column names
        flipped_new_Z<-lapply(flipped_new_Z,setNames,nm=c("SNP","CHR","BP","A1","A2","EAF","Beta","se","Z_NOTUSE","N","PVAL","A1.ref","A2.ref","Z"))
        cat("done\n")
        
        cat("Flipping Beta's for SNPs that are not aligned...")
        Blist<-lapply(flipped,"[","Beta")
        Blist<-lapply(Blist,function(x){x$Beta * -1})
        #make data frame again
        Blist<-lapply(Blist,function(x){as.data.frame(x)})
        #cbind newbetas with not aligned file
        flipped_new_Z<-Map(cbind,flipped_new_Z,Blist) 
        #format right column names
        flipped_new_Z<-lapply(flipped_new_Z,setNames,nm=c("SNP","CHR","BP","A1","A2","EAF","Beta_NOTUSE","se","Z_NOTUSR","N","PVAL","A1.ref","A2.ref","Z","Beta"))
        cat("done\n")
        
        cat("Flipping EAF for SNPs that are not aligned...")
        EAFlist<-lapply(flipped,"[","EAF")
        EAFlist<-lapply(EAFlist,function(x){1-x$EAF})
        EAFlist<-lapply(EAFlist,function(x){as.data.frame(x)})
        flipped_new_Z_EAF<-Map(cbind,flipped_new_Z,EAFlist) 
        flipped_new_Z_EAF<-lapply(flipped_new_Z_EAF,setNames,nm=c("SNP","CHR","BP","A1","A2","EAF_NOTUSE","Beta_NOTUSE","se","Z_NOTUSE","N","PVAL","A1.ref","A2.ref","Z","BEta","EAF"))
        cat("done\n")
        
        rm(flipped_new_Z)
        gc()
        
        #select right columns for rbind with data frames in list that are not complete due to misalignment
        flipped_reorder<-lapply(flipped_new_Z_EAF,function(x){subset(x,select=(c(SNP,MarkerName,CHR,BP,A2,A1,EAF,N,Z,PVAL,A1.ref,A2.ref)))})
        flipped_reorder<-lapply(flipped_reorder,setNames,nm=c("SNP","MarkerName","CHR","BP","A1","A2","EAF","N","Z","PVAL","A1.ref","A2.ref"))
        
        rm(flipped_new_Z_EAF)
        gc()
        
        cat("Combining the newly aligned SNPs with SNPs that were already combined...")
        x_common_aligned<-Map(rbind,aligned,flipped_reorder)
        cat("done\n")
        
        rm(aligned,flipped_reorder,Zlist,EAFlist)
        gc()
      }
      rm(flipped,N_flipped)
      gc()
      
      if(exists("x_common_aligned")==F){
        x_common_aligned<-aligned
        rm(aligned)
        gc()
      }
      
      all_aligned<-lapply(x_common_aligned,function(x){subset(x,A1==A1.ref & A2==A2.ref)})
      rm(x_common_aligned)
      gc()
      
      unalignable<-lapply(x_reffile,function(x){subset(x,(A1!=A1.ref & A1!=A2.ref)|(A2!=A1.ref & A2!=A2.ref))})
      N_unalignable<-sapply(unalignable,nrow)
      if(sum(N_unalignable)!=0){
        cat("On rare occasions, SNPs cannot be aligned (e.g. A-C in the reference file vs A-T in the input). This might lead to summary statistics of different length\n")
        cat("How many markers could not be aligned?\n")
        print(noquote(format(N_unalignable,big.mark=",")))
        cat("These were removed from the data frame(s)\nA list was saved at ",out,".not_aligned",sep="")
        unalignable_SNPs<-unique(unlist(lapply(unalignable,"[[","MarkerName")))
        if(output.gz){
          cat(".gz\n")
          gz1<-gzfile(paste0(out,".not_aligned.gz"),"w")
        }else{
          cat("\n")
          gz1<-paste0(out,".not_aligned")
        }
        write.table(x=data.frame(MarkerName=unalignable_SNPs,stringsAsFactors=F),file=gz1,quote=F,sep="\t",row.names=F)
        if(is.character(gz1)==F){
          close(gz1)
        }
        rm(unalignable_SNPs,gz1)
        gc()
      }
      rm(unalignable,N_unalignable,x_reffile)
      gc()
      
      cat("Selecting only SNPs that are common again and extracting them from the data frames...")
      SNPlist.mat<-extract_SNPs(all_aligned,2)
      SNPs<-data.frame(MarkerName=rownames(SNPlist.mat[SNPlist.mat$TOTAL==length(all_aligned),]),stringsAsFactors=F)
      all_aligned<-lapply(all_aligned,function(x,y){x[match(y$MarkerName,x$MarkerName),]},y=SNPs)
      cat("done\n")
      rm(SNPlist.mat,SNPs)
      gc()
    }else{
      cat("All SNPs were aligned\n")
      all_aligned<-x_common
      rm(x_common,A1A2_file1)
      gc()
    }
    rm(identicalYesNo)
    gc()
    
    cat(format(nrow(all_aligned[[1]]),big.mark=",")," SNPs will be used in the N-weighted GWAMA\n",sep="")
    if(save.matrix){
      cat("A list was saved at ",out,".N_weighted_GWAMA.input",sep="")
      if(output.gz){
        cat(".gz\n")
        gz1<-gzfile(paste0(out,".N_weighted_GWAMA.input.gz"),"w")
      }else{
        cat("\n")
        gz1<-paste0(out,".N_weighted_GWAMA.input")
      }		
      write.table(x=all_aligned[[1]][,1],file=gz1,quote=F,sep="\t",row.names=F,col.names=c("MarkerName"))
      if(is.character(gz1)==F){
        close(gz1)
      }
      rm(gz1)
    }
    
    
    cat("Extracting Z-scores from aligned summary statistics...")
    Zlist<-lapply(all_aligned,"[","Z")
    cat("done\n")
    
    cat("Extracting Beta's ...")
    Blist<-lapply(all_aligned,"[","Beta")
    
    Blist <- Blist 
    
    
    
  }
    
  }
  
  
}
