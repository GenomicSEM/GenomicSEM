



GSEM.process.munge <- function(files,hm3,trait.names=NULL,N,info.filter = .9,maf.filter=0.01){
  
  length <- length(files)
  
  
  
  files = lapply(files, read.table, header=T, quote="\"",fill=T)
  
  ref <- fread(hm3,header=T,data.table=F)
  
  
  
  data.frame.out <- ref
  
  for(i in 1:length){
    
    hold_names <- names(files[[i]])
    
    hold_names[hold_names %in%c("snp","SNP","snpid","SNPID","rsid","RSID","RS_NUMBER","rs_number","RS_NUMBERS","rs_numbers","MarkerName", "markername", "MARKERNAME")] <- "SNP"
    hold_names[hold_names %in%c("a1","A1","allele1","Allele1", "ALLELE1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA")] <- "A1"
    hold_names[hold_names %in%c("a2","A2","allele2","Allele2","ALLELE2","OTHER_ALLELE","NON_EFFECT_ALLELE","DEC_ALLELE","NEA")]  <- "A2"
    hold_names[hold_names %in%c("OR","or","B","beta","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT", "Effect")] <- "effect"
    hold_names[hold_names %in%c("se","StdErr","SE")] <- "SE"
    hold_names[hold_names %in%c("INFO","info")] <- "INFO"
    hold_names[hold_names %in%c("P","p","PVALUE","pvalue","P_VALUE","p_value","PVAL","pval","P_VAL","p_val","GC_PVALUE","gc_pvalue" )] <- "P"
    
    ##rename common MAF labels
    hold_names[hold_names %in%c("MAF","maf", "CEUaf", "Freq1")] <- "MAF"
    
    names(files[[i]]) <- hold_names
    
    ##make sure all alleles are upper case for matching
    files[[i]]$A1 <- factor(toupper(files[[i]]$A1), c("A", "C", "G", "T"))
    files[[i]]$A2 <- factor(toupper(files[[i]]$A2), c("A", "C", "G", "T"))
    
    ##merge with ref file
    files[[i]] <- merge(ref,files[[i]],by="SNP",all.x=F,all.y=F)
    
    ##determine whether it is OR or logistic/continuous effect based on median effect size 
    files[[i]]$effect<-ifelse(rep(round(median(files[[i]]$effect)) == 1,nrow(files[[i]])), log(files[[i]]$effect),files[[i]]$effect)
    
    # Flip effect to match ordering in ref file
    files[[i]]$effect <-  ifelse(files[[i]]$A1.x != (files[[i]]$A1.y) & files[[i]]$A1.x == (files[[i]]$A2.y),files[[i]]$effect*-1,files[[i]]$effect)
    
    ##remove SNPs that don't match A1 OR A2 in ref. confirm this works
    files[[i]]<-subset(files[[i]], !(files[[i]]$A2.x != (files[[i]]$A2.y)  & files[[i]]$A2.x !=  (files[[i]]$A1.y)))
    files[[i]]<-subset(files[[i]], !(files[[i]]$A1.x != (files[[i]]$A1.y)  & files[[i]]$A1.x != (files[[i]]$A2.y)))
    
    # Validity checks:
    
    if((sum(files[[i]]$P > 1) + sum(files[[i]]$P < 0)) > 100){
      
      print("in excess of 100 SNPs have P val above 1 or below 0, column may be mislabled!")
      
    }
    
    #Compute Z score
    files[[i]]$Z <- sign(files[[i]]$effect) * sqrt(qchisq(files[[i]]$P,1,lower=F))
    
    
    
    
    if("INFO" %in% colnames(files[[i]])) {
      files[[i]] <- files[[i]][files[[i]]$INFO >= info.filter,]
      
    }else{print("No INFO column, cant filter on INFO, may influence results")}
    
    if("MAF" %in% colnames(files[[i]])) {
      files[[i]] <- files[[i]][files[[i]]$MAF >= maf.filter,]
      
    }else{print("No MAF column, cant filter on MAF, may influence results")}
    
    
    
    
    output <- cbind.data.frame(files[[i]]$SNP,N[i],files[[i]]$Z,files[[i]]$A1.x,files[[i]]$A2.x)
    
    colnames(output) <- c("SNP","N","Z","A1","A2")
    
    
    write.table(x = output,file = paste0(trait.names[i],".sumstats"), quotes = FALSE)
    gzip(paste0(trait.names[i],".sumstats"))
    
  }
  
  
  
  
}
