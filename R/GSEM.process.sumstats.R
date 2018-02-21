
GSEM.process.sumstats <- function(files,ref,trait.names=NULL,se.logit){
  
  length <- length(files)
  
  if(is.null(trait.names)){
    
    names.beta <- paste0("beta.",1:length)
    names.se <- paste0("se.",1:length)
    
  }else{
    
    names.beta <- paste0("beta.",trait.names)
    names.se <- paste0("se.",trait.names)
    
  }
  
  
  files = lapply(files, read.table, header=T, quote="\"")
  
  ref <- fread(ref,header=T,data.table=F)
  
  files_1 <- files
  files <- files_1
  
  data.frame.out <- ref
  
  
  for(i in 1:length){
    
    hold_names <- names(files[[i]])
    
    hold_names[hold_names %in%c("snp","SNP","snpid","rsid","RSID","RS_NUMBER","rs_number","RS_NUMBERS","rs_numbers")] <- "SNP"
    hold_names[hold_names %in%c("a1","A1","allele1","ALLELE1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA")] <- "A1"
    hold_names[hold_names %in%c("a2","A2","allele2","ALLELE2","OTHER_ALLELE","NON_EFFECT_ALLELE","DEC_ALLELE","NEA")]  <- "A2"
    hold_names[hold_names %in%c("OR","or","B","beta","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT")] <- "effect"
    hold_names[hold_names %in%c("se","StdErr","SE")] <- "SE"
    
    
    names(files[[i]]) <- hold_names
    
    
    files[[i]]$A1 <- factor(toupper(files[[i]]$A1), c("A", "C", "G", "T"))
    files[[i]]$A2 <- factor(toupper(files[[i]]$A2), c("A", "C", "G", "T"))
    
    files[[i]]$effect <- ifelse(rep(round(mean(files[[i]]$effect)) == 1,nrow(files[[i]])), log(files[[i]]$effect),files[[i]]$effect)  
    
    
    files[[i]] <- merge(ref,files[[i]],by="SNP",all.x=F,all.y=F)
    
    # Flip all alleles to the ref
    files[[i]]$effect <-  ifelse(files[[i]]$A1.x != (files[[i]]$A1.y) & files[[i]]$A1.x == (files[[i]]$A2.y)  ,files[[i]]$effect*-1,files[[i]]$effect)
    files[[i]]$A2 <- ifelse(files[[i]]$A2.x == (files[[i]]$A1.y),files[[i]]$A1.y,files[[i]]$A2.y)
    files[[i]]$A1 <-ifelse(files[[i]]$A1.x == (files[[i]]$A2.y),files[[i]]$A2.y,files[[i]]$A1.y)
    
    
    varSNP<-2*files[[i]]$MAF*(1-files[[i]]$MAF)  
    
    if(se.logit[i] == F){
      output <- cbind.data.frame(files[[i]]$SNP,
                      (files[[i]]$effect)/((files[[i]]$effect^2) * varSNP + (pi^2)/3)^.5,
                      (files[[i]]$SE/exp(files[[i]]$effect))/(((files[[i]]$SE/exp(files[[i]]$effect))^2) * varSNP + (pi^2)/3)^.5    )
    
      colnames(output) <- c("SNP",names.beta[i],names.se[i])
      
      } else{
      output <- cbind.data.frame(files[[i]]$SNP,
      (files[[i]]$effect)/((files[[i]]$effect^2) * varSNP + (pi^2)/3)^.5,
      (files[[i]]$SE)/(((files[[i]]$SE)^2) * varSNP + (pi^2)/3)^.5    )  
 
      colnames(output) <- c("SNP",names.beta[i],names.se[i])
    }
    
    
    
    
    if(i ==1){
      data.frame.out <- merge(data.frame.out,output,by=1,all.x=F,all.y=F)
    }else{
      data.frame.out <- merge(data.frame.out,output,by=1,all.x=F,all.y=F) 
    }
    
  }
  
  
  
  data.frame.out
  
}
