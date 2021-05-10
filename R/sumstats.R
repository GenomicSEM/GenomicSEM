sumstats <- function(files,ref,trait.names=NULL,se.logit,OLS=NULL,linprob=NULL,prop=NULL,N=NULL,info.filter = .6,maf.filter=0.01,keep.indel=FALSE,parallel=FALSE,cores=NULL){
  
  begin.time <- Sys.time()
  
  length <- length(files)
  
  filenames <- as.vector(files)
  
  ref2<-ref
  
  if(is.null(OLS)){
    OLS<-rep(FALSE,length)
  }
  
  if(is.null(linprob)){
    linprob<-rep(FALSE,length)
  }
  
  if(is.null(trait.names)){
    names.beta <- paste0("beta.",1:length)
    names.se <- paste0("se.",1:length)
  }else{
    names.beta <- paste0("beta.",trait.names)
    names.se <- paste0("se.",trait.names)
  }
  
  if(parallel == FALSE){
    
    log2<-paste(trait.names,collapse="_")
    
    log.file <- file(paste0(log2, "_sumstats.log"),open="wt")
    
    cat(print(paste0("The preparation of ", length(trait.names), " summary statistics for use in Genomic SEM began at: ",begin.time), sep = ""),file=log.file,sep="\n",append=TRUE)
    cat(print(paste0("Please note that the files should be in the same order that they were listed for the ldsc function"), sep = ""),file=log.file,sep="\n",append=TRUE)
  
    cat(print("Reading in reference file"),file=log.file,sep="\n",append=TRUE)
    ref <- fread(ref,header=T,data.table=F)
    
    ##filter ref file on user provided maf.filter
    cat(print(paste("Applying MAF filer of", maf.filter, "to the refernece file.")),file=log.file,sep="\n",append=TRUE)
    ref<-subset(ref, ref$MAF >= maf.filter)
    
    data.frame.out <- ref
    
    cat(print(paste("Reading summary statistics for", paste(files,collapse=" "), ". Please note that this step usually takes a few minutes due to the size of summary statistic files.")),file=log.file,sep="\n",append=TRUE)
    
    ##note that fread is not used here as we have observed different formatting for column headers causing mismatched columns
    files = lapply(files, read.table, header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
    
    cat(print("All files loaded into R!"),file=log.file,sep="\n",append=TRUE)
    
    for(i in 1:length){
      
      cat(paste("     "),file=log.file,sep="\n",append=TRUE)
      cat(paste("     "),file=log.file,sep="\n",append=TRUE)
      
      cat(print(paste("Preparing summary statistics for file:", filenames[i])),file=log.file,sep="\n",append=TRUE)
      
      hold_names <- toupper(names(files[[i]]))
      names1<-hold_names
      
      if("SNP" %in% hold_names) cat(print(paste("Interpreting the SNP column as the SNP column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in% c("SNP","SNPID","RSID","RS_NUMBER","RS_NUMBERS", "MARKERNAME", "ID", "SNP_ID")] <- "SNP"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the SNP column.")),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("A1" %in% hold_names) cat(print(paste("Interpreting the A1 column as the A1 column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("A1", "ALLELE1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA","REF")] <- "A1"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the A1 column.")),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("A2" %in% hold_names) cat(print(paste("Interpreting the A2 column as the A2 column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("A2","ALLELE2","ALLELE0","OTHER_ALLELE","REF","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA", "ALT")]  <- "A2"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the A2 column.")),file=log.file,sep="\n",append=TRUE)
      
      if(linprob[i] == F){
        if(OLS[i] == F){ 
          if("Z" %in% names1 | "ZSCORE" %in% names1 | "Z-SCORE" %in% names1 | "ZSTATISTIC" %in% names1 | "Z-STATISTIC" %in% names1){
            warning(paste0("There appears to be a Z-statistic column in the summary statistic file for ", trait.names[i], ". Transformations for case/control traits require either an OR or logistic beta column. Please remove/replace the Z-statistic column"))
            cat(print(paste("WARNING: There appears to be a Z-statistic column in the summary statistic file for ", trait.names[i], ". Transformations for case/control traits require either an OR or logistic beta column. Please remove/replace the Z-statistic column"),file=log.file,sep="\n",append=TRUE))
          }
        }}
      
      
      names1<-hold_names
      if("effect" %in% hold_names) cat(print(paste("Interpreting the effect column as the effect column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("OR","B","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT", "Z","ZSCORE","EST","ZSTAT","ZSTATISTIC", "BETA1")] <- "effect"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the effect column.")),file=log.file,sep="\n",append=TRUE)
    
     
      names1<-hold_names
      if("INFO" %in% hold_names) cat(print(paste("Interpreting the INFO column as the INFO column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("INFO")] <- "INFO"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the INFO column.")),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("SE" %in% hold_names) cat(print(paste("Interpreting the SE column as the SE column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("STDERR","SE")] <- "SE"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the SE column.")),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("P" %in% hold_names) cat(print(paste("Interpreting the P column as the P column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE")] <- "P"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the", setdiff(names1, hold_names), "column as the P column.")),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("N" %in% hold_names) cat(print(paste("Interpreting the N column as the N (sample size) column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("N","WEIGHT","NCOMPLETESAMPLES", "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N","N_COMPLETE_SAMPLES", "SAMPLESIZE")] <- "N"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the ", setdiff(names1, hold_names), " column as the N (sample size) column.")),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("N_CAS" %in% hold_names) cat(print(paste("Interpreting the N_CAS column as the N_CAS (sample size for cases) column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("NCASE","N_CASE","N_CASES","N_CAS", "NCAS", "NCA")] <- "N_CAS"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the ", setdiff(names1, hold_names), " column as the N_CAS (sample size for cases) column.")),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("N_CON" %in% hold_names) cat(print(paste("Interpreting the N_CON column as the N_CON (sample size for controls) column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("NCONTROL","N_CONTROL","N_CONTROLS","N_CON","CONTROLS_N", "NCON", "NCO")] <- "N_CON"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(print(paste("Interpreting the ", setdiff(names1, hold_names), " column as the N_CON (sample size for controls) column.")),file=log.file,sep="\n",append=TRUE)
      
       
      
      # Print a message for misisng P value, rs, effect or allele column
      if(sum(hold_names %in% "P") == 0) cat(print(paste0('Cannot find P-value column, try renaming it P in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "A1") == 0) cat(print(paste0('Cannot find effect allele column, try renaming it A1 in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "A2") == 0) cat(print(paste0('Cannot find other allele column, try renaming it A2 in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "effect") == 0) cat(print(paste0('Cannot find beta or effect column, try renaming it effect in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "SNP") == 0) cat(print(paste0('Cannot find rs-id column, try renaming it SNP in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
      
      # Print a warning message when multiple columns interpreted as P-values, rsID, effect or allele columns
      if(sum(hold_names %in% "P") > 1) cat(print(paste0('Multiple columns are being interpreted as the P-value column. Try renaming the column you dont want interpreted as P to P2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "A1") > 1) cat(print(paste0('Multiple columns are being interpreted as the effect allele column. Try renaming the column you dont want interpreted as effect allele column to A1_2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "A2") > 1) cat(print(paste0('Multiple columns are being interpreted as the other allele column. Try renaming the column you dont want interpreted as the other allele column to A2_2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "effect") > 1) cat(print(paste0('Multiple columns are being interpreted as the beta or effect column. Try renaming the column you dont want interpreted as the beta or effect column to effect2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "SNP") > 1) cat(print(paste0('Multiple columns are being interpreted as the rs-id column. Try renaming the column you dont want interpreted as rs-id to SNP2 for:',filenames[i])),file=log.file,sep="\n",append=TRUE)
      
      # Throw warnings for misisng P valuue, rs, effect or allele columns
      if(sum(hold_names %in% "P") == 0) warning(paste0('Cannot find P-value column, try renaming it P in the summary statistics file for:',trait.names[i]))
      if(sum(hold_names %in% "A1") == 0) warning(paste0('Cannot find effect allele column, try renaming it A1 in the summary statistics file for:',trait.names[i]))
      if(sum(hold_names %in% "A2") == 0) warning(paste0('Cannot find other allele column, try renaming it A2 in the summary statistics file for:',trait.names[i]))
      if(sum(hold_names %in% "effect") == 0) warning(paste0('Cannot find beta or effect column, try renaming it effect in the summary statistics file for:',trait.names[i]))
      if(sum(hold_names %in% "SNP") == 0) warning(paste0('Cannot find rs-id column, try renaming it SNP in the summary statistics file for:',trait.names[i]))
      
      # Print a warning message when multiple columns interpreted as P-values, rsID, effect or allele columns
      if(sum(hold_names %in% "P") > 1) warning(paste0('Multiple columns are being interpreted as the P-value column. Try renaming the column you dont want interpreted as P to P2 for:',filenames[i]))
      if(sum(hold_names %in% "A1") > 1) warning(paste0('Multiple columns are being interpreted as the effect allele column. Try renaming the column you dont want interpreted as effect allele column to A1_2 for:',filenames[i]))
      if(sum(hold_names %in% "A2") > 1) warning(paste0('Multiple columns are being interpreted as the other allele column. Try renaming the column you dont want interpreted as the other allele column to A2_2 for:',filenames[i]))
      if(sum(hold_names %in% "effect") > 1) warning(paste0('Multiple columns are being interpreted as the beta or effect column. Try renaming the column you dont want interpreted as the beta or effect column to effect2 for:',filenames[i]))
      if(sum(hold_names %in% "SNP") > 1) warning(paste0('Multiple columns are being interpreted as the rs-id column. Try renaming the column you dont want interpreted as rs-id to SNP2 for:',filenames[i]))
      
      
      ##rename common MAF labels to MAF_Other so MAF from ref file is used across traits for conversions
      hold_names[hold_names %in%c("MAF", "CEUAF", "FREQ1", "EAF", "FREQ1.HAPMAP", "FREQALLELE1HAPMAPCEU", "FREQ.ALLELE1.HAPMAPCEU", "EFFECT_ALLELE_FREQ", "FREQ.A1")] <- "MAF_Other"
      
      names(files[[i]]) <- hold_names
      
      # Compute N as N cases and N control if reported:
      if("N_CAS" %in% colnames(files[[i]]) & "N_CON" %in% colnames(files[[i]])){
        files[[i]]$N <- files[[i]]$N_CAS + files[[i]]$N_CON
        cat(print(paste("As the file includes both N_CAS and N_CON columns, the summation of these two columns will be used as the total sample size")),file=log.file,sep="\n",append=TRUE)
      }
      
      if("N" %in% colnames(files[[i]]) & !(is.null(N))){
        if(!(is.na(N[i]))){
          cat(print(paste("As the summary statistics file includes a sample size column, this is being used in place of the user provided sample size. If the user wishes to still use the provided sample size, as opposed to the sample size listed in the summary statistics file, please change the sample size column header to N2.However, we note that sample size is only used for LPM and OLS conversions.")),file=log.file,sep="\n",append=TRUE)
        }}
      
      if(!(is.null(N)) & !("N" %in% colnames(files[[i]]))){
        if(!(is.na(N[i]))){
          files[[i]]$N<-N[i]
        }}
      
      if(keep.indel == TRUE){
        files[[i]]$A1 <- factor(toupper(files[[i]]$A1))
        files[[i]]$A2 <- factor(toupper(files[[i]]$A2))
        cat(print(paste0("Keeping variants other than SNPs, this may cause problems when alligning allle's across traits and the reference file")),file=log.file,sep="\n",append=TRUE)
      }
      
      if(keep.indel == FALSE){
        ##make sure all alleles are upper case for matching
        files[[i]]$A1 <- factor(toupper(files[[i]]$A1), c("A", "C", "G", "T"))
        files[[i]]$A2 <- factor(toupper(files[[i]]$A2), c("A", "C", "G", "T"))
      }
      
      ##merge with ref file
      cat(print(paste0("Merging file: ", filenames[i], "with the reference file: ", ref2)),file=log.file,sep="\n",append=TRUE)
      b<-nrow(files[[i]])
      cat(print(paste0(b, " rows present in the full ", filenames[i], " summary statistics file.")),file=log.file,sep="\n",append=TRUE)
      files[[i]] <- suppressWarnings(inner_join(ref,files[[i]],by="SNP",all.x=F,all.y=F))
      cat(print(paste((b-nrow(files[[i]])), "rows were removed from the", filenames[i], "summary statistics file as the rsIDs for these SNPs were not present in the reference file.")),file=log.file,sep="\n",append=TRUE)
      
      ##remove any rows with missing p-values
      b<-nrow(files[[i]])
      if("P" %in% colnames(files[[i]])) {
        files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$P)))
      }
      if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to missing values in the P-value column")),file=log.file,sep="\n",append=TRUE)
      
      ##remove any rows with missing effects
      b<-nrow(files[[i]])
      if("effect" %in% colnames(files[[i]])) {
        files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$effect)))
      }
      if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to missing values in the effect column")),file=log.file,sep="\n",append=TRUE)
      
      ##determine whether it is OR or logistic/continuous effect based on median effect size 
      a1<-files[[i]]$effect[[1]]
      files[[i]]$effect<-ifelse(rep(round(median(files[[i]]$effect,na.rm=T)) == 1,nrow(files[[i]])), log(files[[i]]$effect),files[[i]]$effect)
      a2<-files[[i]]$effect[[1]]
      if(a1 != a2) cat(print(paste("The effect column was determined to be coded as an odds ratio (OR) for the", filenames[i], "summary statistics file based on the median of the effect column being close to 1. Please ensure the interpretation of this column as an OR is correct.")),file=log.file,sep="\n",append=TRUE)
      if(a1 == a2) cat(print(paste("The effect column was determined NOT to be coded as an odds ratio (OR) for the", filenames[i], "summary statistics file based on the median of the effect column being close to 0.")),file=log.file,sep="\n",append=TRUE)
     
      ##remove any rows printed as exactly 0
      b<-nrow(files[[i]])
      if("effect" %in% colnames(files[[i]])) {
        files[[i]]<-subset(files[[i]], files[[i]]$effect != 0)
      }
      if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to effect values estimated at exactly 0 as this causes problems for matrix inversion necessary for later Genomic SEM analyses.")),file=log.file,sep="\n",append=TRUE)
      
      
      if(OLS[i] == T){
        cat(print(paste("An OLS transformation is being used for file:", filenames[i])),file=log.file,sep="\n",append=TRUE)
        
        files[[i]]$Z <- sign(files[[i]]$effect) * sqrt(qchisq(files[[i]]$P,1,lower=F))
        
        if("N" %in% colnames(files[[i]])){
          files[[i]]$effect <- files[[i]]$Z/ sqrt(files[[i]]$N * 2 * (files[[i]]$MAF *(1-files[[i]]$MAF)))}else{cat(print("ERROR: A Sample Size (N) is needed for OLS Standardization. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N."),file=log.file,sep="\n",append=TRUE)}}
      
      
      if(linprob[i] == T){
        cat(print(paste("An LPM transformation is being used for file:", filenames[i])),file=log.file,sep="\n",append=TRUE)
        
        files[[i]]$Z <- sign(files[[i]]$effect) * sqrt(qchisq(files[[i]]$P,1,lower=F))
        
        if("N" %in% colnames(files[[i]])){
          files[[i]]$effect <- files[[i]]$Z/sqrt((prop[i]*(1-prop[i])*(2*files[[i]]$N*files[[i]]$MAF*(1-files[[i]]$MAF))))
          files[[i]]$SE<-1/sqrt((prop[i]*(1-prop[i])*(2*files[[i]]$N*files[[i]]$MAF*(1-files[[i]]$MAF))))}else{cat(print("ERROR: A Sample Size (N) is needed for LPM Standardization. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N."),file=log.file,sep="\n",append=TRUE)}}
      
      # Flip effect to match ordering in ref file
      files[[i]]$effect <-  ifelse(files[[i]]$A1.x != (files[[i]]$A1.y) & files[[i]]$A1.x == (files[[i]]$A2.y),files[[i]]$effect*-1,files[[i]]$effect)
      
      ##remove SNPs that don't match A1 OR A2 in ref. 
      b<-nrow(files[[i]])
      files[[i]]<-subset(files[[i]], !(files[[i]]$A1.x != (files[[i]]$A1.y)  & files[[i]]$A1.x != (files[[i]]$A2.y)))
      if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "row(s) were removed from the", filenames[i], "summary statistics file due to the effect allele (A1) column not matching A1 or A2 in the reference file.")),file=log.file,sep="\n",append=TRUE)
      
      b<-nrow(files[[i]])
      files[[i]]<-subset(files[[i]], !(files[[i]]$A2.x != (files[[i]]$A2.y)  & files[[i]]$A2.x !=  (files[[i]]$A1.y)))
      if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "row(s) were removed from the", filenames[i], "summary statistics file due to the other allele (A2) column not matching A1 or A2 in the reference file.")),file=log.file,sep="\n",append=TRUE)
      
      #Check that p-value column does not contain an excess of 1s/0s
      if((sum(files[[i]]$P > 1) + sum(files[[i]]$P < 0)) > 100){
        cat(print("In excess of 100 SNPs have P val above 1 or below 0. The P column may be mislabled!"),file=log.file,sep="\n",append=TRUE)
      }
      
      if("INFO" %in% colnames(files[[i]])) {
        b<-nrow(files[[i]])
        files[[i]] <- files[[i]][files[[i]]$INFO >= info.filter,]
        cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to INFO values below the designated threshold of", info.filter)),file=log.file,sep="\n",append=TRUE)
      }else{cat(print("No INFO column, cannot filter on INFO, which may influence results"),file=log.file,sep="\n",append=TRUE)}
      
      varSNP<-2*files[[i]]$MAF*(1-files[[i]]$MAF)  
      
      if(OLS[i] == T){
        output <- cbind.data.frame(files[[i]]$SNP,
                                   files[[i]]$effect,
                                   abs(files[[i]]$effect/files[[i]]$Z))
        output<-na.omit(output)                           
        colnames(output) <- c("SNP",names.beta[i],names.se[i])                           
      }
      
      if(linprob[i] == T){
        output<-cbind.data.frame(files[[i]]$SNP,
                                 (files[[i]]$effect)/((files[[i]]$effect^2) * varSNP + (pi^2)/3)^.5,
                                 (files[[i]]$SE)/(((files[[i]]$effect)^2) * varSNP + (pi^2)/3)^.5)  
        output<-na.omit(output)
        output<-output[apply(output!=0, 1, all),]
        colnames(output) <- c("SNP",names.beta[i],names.se[i])                                         
      }
      
      if(linprob[i] == F){
        if(OLS[i] == F){                                     
          if(se.logit[i] == F){
            cat(print(paste("Performing transformation under the assumption that the effect column is either an odds ratio or logistic beta (please see output above to determine whether it was interpreted as an odds ratio) and the SE column is the SE of the odds ratio (i.e., NOT on the logistic scale) for:", filenames[i])),file=log.file,sep="\n",append=TRUE)
            
            if(sum(hold_names %in% "SE") == 0) cat(print(paste0('Cannot find SE column, try renaming it SE in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
            if(sum(hold_names %in% "SE") == 0) warning(paste0('Cannot find SE column, try renaming it SE in the summary statistics file for:',trait.names[i]))
            
            output <- cbind.data.frame(files[[i]]$SNP,
                                       (files[[i]]$effect)/((files[[i]]$effect^2) * varSNP + (pi^2)/3)^.5,
                                       (files[[i]]$SE/exp(files[[i]]$effect))/(((files[[i]]$effect)^2 * varSNP + (pi^2)/3)^.5))
            output<-na.omit(output)  
            colnames(output) <- c("SNP",names.beta[i],names.se[i])}}}
      
      if(se.logit[i]== T){
        cat(print(paste("Performing transformation under the assumption that the effect column is either an odds ratio or logistic beta (please see output above to determine whether it was interpreted as an odds ratio) and the SE column is a logistic SE (i.e., NOT the SE of the odds ratio) for:", filenames[i])),file=log.file,sep="\n",append=TRUE)
        
        if(sum(hold_names %in% "SE") == 0) cat(print(paste0('Cannot find SE column, try renaming it SE in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
        if(sum(hold_names %in% "SE") == 0) warning(paste0('Cannot find SE column, try renaming it SE in the summary statistics file for:',trait.names[i]))
        
        output <- cbind.data.frame(files[[i]]$SNP,
                                   (files[[i]]$effect)/((files[[i]]$effect^2) * varSNP + (pi^2)/3)^.5,
                                   (files[[i]]$SE)/(((files[[i]]$effect)^2) * varSNP + (pi^2)/3)^.5)  
        output<-na.omit(output)  
        colnames(output) <- c("SNP",names.beta[i],names.se[i])}
      
      cat(print(paste(nrow(output), "SNPs are left in the summary statistics file", filenames[i], "after QC and merging with the reference file.")),file=log.file,sep="\n",append=TRUE)
      
      if(mean(abs(output[,2]/output[,3])) > 5){
        cat(print(paste0('WARNING: The average value of estimate over standard error (i.e., Z) is > 5 for ',trait.names[i], ". This suggests a column was misinterpreted or arguments were misspecified. Please post on the google group if you are unable to figure out the issue.")),file=log.file,sep="\n",append=TRUE)
        warning(paste0('The average value of estimate over standard error (i.e., Z) is > 5 for ',trait.names[i], ". This suggests a column was misinterpreted or arguments were misspecified. Please post on the google group if you are unable to figure out the issue."))
      }
      
      if(i ==1){
        data.frame.out <- suppressWarnings(inner_join(data.frame.out,output,by="SNP",all.x=F,all.y=F))
      }else{
        b<-nrow(output)
        data.frame.out <- suppressWarnings(inner_join(data.frame.out,output,by="SNP",all.x=F,all.y=F)) 
        cat(print(paste((b-nrow(data.frame.out)), "rows were removed from the", filenames[i], "summary statistics file as the rsIDs for these SNPs were not present for the other summary statistics.")),file=log.file,sep="\n",append=TRUE)
      }
    }
  }
  
  if(parallel == TRUE){
    
    print("Please note that the files should be in the same order that they were listed for the ldsc function.")
    print("Reading in reference file")
    ref <- fread(ref,header=T,data.table=F)
    
    ##filter ref file on user provided maf.filter
    print(paste("Applying MAF filer of", maf.filter, "to the reference file."))
    ref<-subset(ref, ref$MAF >= maf.filter)
    
    data.frame.out <- ref
    
   if(is.null(cores)){
  ##if no default provided use 1 less than the total number of cores available so your computer will still function
  int <- detectCores() - 1
}else{int<-cores}

#if more cores than traits then set to number of traits
if(int > length){
  int<-length
}
    
    print("Performing conversions of individual summary statistics using parallel processing. Please note this step typically takes 10-20 minutes due to the size of the files.")
    Output<-mclapply(X=1:length,FUN=function(X){
      i<-X
      
      log.file <- file(paste0(trait.names[i], "_sumstats.log"),open="wt")  
      files2<-read.table(files[[i]], header = T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
      
      cat(paste("     "),file=log.file,sep="\n",append=TRUE)
      cat(paste("     "),file=log.file,sep="\n",append=TRUE)
      
      cat(print(paste0("Please note that the files should be in the same order that they were listed for the ldsc function"), sep = ""),file=log.file,sep="\n",append=TRUE)
      cat(print(paste("Preparing summary statistics for file:", filenames[i])),file=log.file,sep="\n",append=TRUE)
      
      hold_names <- toupper(names(files2))
      names1<-hold_names
      
      if("SNP" %in% hold_names) cat(print(paste("Interpreting the SNP column as the SNP column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in% c("SNP","SNPID","RSID","RS_NUMBER","RS_NUMBERS", "MARKERNAME", "ID")] <- "SNP"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the", setdiff(names1, hold_names), "column as the SNP column."),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("A1" %in% hold_names) cat(print(paste("Interpreting the A1 column as the A1 column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("A1", "ALLELE1","EFFECT_ALLELE","INC_ALLELE","REFERENCE_ALLELE","EA","REF")] <- "A1"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the", setdiff(names1, hold_names), "column as the A1 column."),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("A2" %in% hold_names) cat(print(paste("Interpreting the A2 column as the A2 column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("A2","ALLELE2","ALLELE0","OTHER_ALLELE","REF","NON_EFFECT_ALLELE","DEC_ALLELE","OA","NEA", "ALT")]  <- "A2"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the", setdiff(names1, hold_names), "column as the A2 column."),file=log.file,sep="\n",append=TRUE)
      
      if(linprob[i] == F){
        if(OLS[i] == F){ 
          if("Z" %in% names1 | "ZSCORE" %in% names1 | "Z-SCORE" %in% names1 | "ZSTATISTIC" %in% names1 | "Z-STATISTIC" %in% names1){
            cat(paste("WARNING: There appears to be a Z-statistic column in the summary statistic file for ", trait.names[i], ". Transformations for case/control traits require either an OR or logistic beta column. Please remove/replace the Z-statistic column"),file=log.file,sep="\n",append=TRUE)
          }
        }}
      
      
      names1<-hold_names
      if("effect" %in% hold_names) cat(print(paste("Interpreting the effect column as the effect column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("OR","B","BETA","LOG_ODDS","EFFECTS","EFFECT","SIGNED_SUMSTAT", "Z","ZSCORE","EST","ZSTAT","ZSTATISTIC")] <- "effect"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the", setdiff(names1, hold_names), "column as the effect column."),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("INFO" %in% hold_names) cat(print(paste("Interpreting the INFO column as the INFO column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("INFO")] <- "INFO"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the", setdiff(names1, hold_names), "column as the INFO column."),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("SE" %in% hold_names) cat(print(paste("Interpreting the SE column as the SE column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("STDERR","SE")] <- "SE"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the", setdiff(names1, hold_names), "column as the SE column."),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("P" %in% hold_names) cat(print(paste("Interpreting the P column as the P column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("P","PVALUE","PVAL","P_VALUE","P-VALUE","P.VALUE","P_VAL","GC_PVALUE")] <- "P"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the", setdiff(names1, hold_names), "column as the P column."),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("N" %in% hold_names) cat(print(paste("Interpreting the N column as the N (sample size) column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("N","WEIGHT","NCOMPLETESAMPLES", "TOTALSAMPLESIZE", "TOTALN", "TOTAL_N","N_COMPLETE_SAMPLES")] <- "N"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the ", setdiff(names1, hold_names), " column as the N (sample size) column."),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("N_CAS" %in% hold_names) cat(print(paste("Interpreting the N_CAS column as the N_CAS (sample size for cases) column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("NCASE","N_CASE","N_CASES","N_CAS", "NCAS", "NCA")] <- "N_CAS"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the ", setdiff(names1, hold_names), " column as the N_CAS (sample size for cases) column."),file=log.file,sep="\n",append=TRUE)
      
      names1<-hold_names
      if("N_CON" %in% hold_names) cat(print(paste("Interpreting the N_CON column as the N_CON (sample size for controls) column.")),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in%c("NCONTROL","N_CONTROL","N_CONTROLS","N_CON","CONTROLS_N", "NCON", "NCO")] <- "N_CON"
      if(length(base::setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the ", setdiff(names1, hold_names), " column as the N_CON (sample size for controls) column."),file=log.file,sep="\n",append=TRUE)
      
      # Print a message for misisng P value, rs, effect or allele column
      if(sum(hold_names %in% "P") == 0) cat(paste0('Cannot find P-value column, try renaming it P in the summary statistics file for:',trait.names[i]),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "A1") == 0) cat(paste0('Cannot find effect allele column, try renaming it A1 in the summary statistics file for:',trait.names[i]),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "A2") == 0) cat(paste0('Cannot find other allele column, try renaming it A2 in the summary statistics file for:',trait.names[i]),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "effect") == 0) cat(paste0('Cannot find beta or effect column, try renaming it effect in the summary statistics file for:',trait.names[i]),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "SNP") == 0) cat(paste0('Cannot find rs-id column, try renaming it SNP in the summary statistics file for:',trait.names[i]),file=log.file,sep="\n",append=TRUE)
      
      # Print a warning message when multiple columns interpreted as P-values, rsID, effect or allele columns
      if(sum(hold_names %in% "P") > 1) cat(paste0('Multiple columns are being interpreted as the P-value column. Try renaming the column you dont want interpreted as P to P2 for:',filenames[i]),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "A1") > 1) cat(paste0('Multiple columns are being interpreted as the effect allele column. Try renaming the column you dont want interpreted as effect allele column to A1_2 for:',filenames[i]),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "A2") > 1) cat(paste0('Multiple columns are being interpreted as the other allele column. Try renaming the column you dont want interpreted as the other allele column to A2_2 for:',filenames[i]),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "effect") > 1) cat(paste0('Multiple columns are being interpreted as the beta or effect column. Try renaming the column you dont want interpreted as the beta or effect column to effect2 for:',filenames[i]),file=log.file,sep="\n",append=TRUE)
      if(sum(hold_names %in% "SNP") > 1) cat(paste0('Multiple columns are being interpreted as the rs-id column. Try renaming the column you dont want interpreted as rs-id to SNP2 for:',filenames[i]),file=log.file,sep="\n",append=TRUE)
      
      # Throw warnings for misisng P valuue, rs, effect or allele columns
      if(sum(hold_names %in% "P") == 0) warning(paste0('Cannot find P-value column, try renaming it P in the summary statistics file for:',trait.names[i]))
      if(sum(hold_names %in% "A1") == 0) warning(paste0('Cannot find effect allele column, try renaming it A1 in the summary statistics file for:',trait.names[i]))
      if(sum(hold_names %in% "A2") == 0) warning(paste0('Cannot find other allele column, try renaming it A2 in the summary statistics file for:',trait.names[i]))
      if(sum(hold_names %in% "effect") == 0) warning(paste0('Cannot find beta or effect column, try renaming it effect in the summary statistics file for:',trait.names[i]))
      if(sum(hold_names %in% "SNP") == 0) warning(paste0('Cannot find rs-id column, try renaming it SNP in the summary statistics file for:',trait.names[i]))
      
      # Print a warning message when multiple columns interpreted as P-values, rsID, effect or allele columns
      if(sum(hold_names %in% "P") > 1) warning(paste0('Multiple columns are being interpreted as the P-value column. Try renaming the column you dont want interpreted as P to P2 for:',filenames[i]))
      if(sum(hold_names %in% "A1") > 1) warning(paste0('Multiple columns are being interpreted as the effect allele column. Try renaming the column you dont want interpreted as effect allele column to A1_2 for:',filenames[i]))
      if(sum(hold_names %in% "A2") > 1) warning(paste0('Multiple columns are being interpreted as the other allele column. Try renaming the column you dont want interpreted as the other allele column to A2_2 for:',filenames[i]))
      if(sum(hold_names %in% "effect") > 1) warning(paste0('Multiple columns are being interpreted as the beta or effect column. Try renaming the column you dont want interpreted as the beta or effect column to effect2 for:',filenames[i]))
      if(sum(hold_names %in% "SNP") > 1) warning(paste0('Multiple columns are being interpreted as the rs-id column. Try renaming the column you dont want interpreted as rs-id to SNP2 for:',filenames[i]))
      
      
      ##rename common MAF labels
      names1<-hold_names
      if("MAF" %in% hold_names) cat(paste("Interpreting the MAF column as the MAF (minor allele frequency) column."),file=log.file,sep="\n",append=TRUE)
      hold_names[hold_names %in% c("MAF","maf", "CEUaf", "Freq1", "EAF", "Freq1.Hapmap", "FreqAllele1HapMapCEU", "Freq.Allele1.HapMapCEU", "EFFECT_ALLELE_FREQ", "Freq.A1")] <- "MAF"
      if(length(setdiff(names1,hold_names)) > 0) cat(paste("Interpreting the ", setdiff(names1, hold_names), " column as the MAF (minor allele frequency) column."),file=log.file,sep="\n",append=TRUE)
      
      ##rename common MAF labels to MAF_Other so MAF from ref file is used across traits for conversions
      hold_names[hold_names %in%c("MAF","maf", "CEUaf", "Freq1", "EAF", "Freq1.Hapmap", "FreqAllele1HapMapCEU", "Freq.Allele1.HapMapCEU", "EFFECT_ALLELE_FREQ", "Freq.A1")] <- "MAF_Other"
      
      names(files2) <- hold_names
      
      # Compute N as N cases and N control if reported:
      if("N_CAS" %in% colnames(files2) & "N_CON" %in% colnames(files2)){
        files2$N <- files2$N_CAS + files2$N_CON
        cat(print(paste("As the file includes both N_CAS and N_CON columns, the summation of these two columns will be used as the total sample size")),file=log.file,sep="\n",append=TRUE)
      }
      
      if("N" %in% colnames(files2) & !(is.null(N))){
        if(!(is.na(N[i]))){
          cat(print(paste("As the summary statistics file includes a sample size column, this is being used in place of the user provided sample size. If the user wishes to still use the provided sample size, as opposed to the sample size listed in the summary statistics file, please change the sample size column header to N2. However, we note that sample size is only used for LPM and OLS conversions.")),file=log.file,sep="\n",append=TRUE)
        }}
      
      if(!(is.null(N)) & !("N" %in% colnames(files2))){
        if(!(is.na(N[i]))){
          files2$N<-N[i]
        }}
      
      ##make sure all alleles are upper case for matching
      files2$A1 <- factor(toupper(files2$A1), c("A", "C", "G", "T"))
      files2$A2 <- factor(toupper(files2$A2), c("A", "C", "G", "T"))
      
      ##merge with ref file
      cat(print(paste("Merging file:", filenames[i], "with the reference file:", ref2)),file=log.file,sep="\n",append=TRUE)
      b<-nrow(files2)
      cat(print(paste(b, "rows present in the full", filenames[i], "summary statistics file.")),file=log.file,sep="\n",append=TRUE)
      files2 <- suppressWarnings(inner_join(ref,files2,by="SNP",all.x=F,all.y=F))
      cat(print(paste((b-nrow(files2)), "rows were removed from the", filenames[i], "summary statistics file as the rsIDs for these SNPs were not present in the reference file.")),file=log.file,sep="\n",append=TRUE)
      
      ##remove any rows with missing p-values
      b<-nrow(files2)
      if("P" %in% colnames(files2)) {
        files2<-subset(files2, !(is.na(files2$P)))
      }
      if(b-nrow(files2) > 0) cat(print(paste(b-nrow(files2), "rows were removed from the", filenames[i], "summary statistics file due to missing values in the P-value column")),file=log.file,sep="\n",append=TRUE)
      
      ##remove any rows with missing effects
      b<-nrow(files2)
      if("effect" %in% colnames(files2)) {
        files2<-subset(files2, !(is.na(files2$effect)))
      }
      if(b-nrow(files2) > 0) cat(print(paste(b-nrow(files2), "rows were removed from the", filenames[i], "summary statistics file due to missing values in the effect column")),file=log.file,sep="\n",append=TRUE)
      
      ##determine whether it is OR or logistic/continuous effect based on median effect size 
      a1<-files2$effect[[1]]
      files2$effect<-ifelse(rep(round(median(files2$effect,na.rm=T)) == 1,nrow(files2)), log(files2$effect),files2$effect)
      a2<-files2$effect[[1]]
      if(a1 != a2) cat(print(paste("The effect column was determined to be coded as an odds ratio (OR) for the", filenames[i], "summary statistics file based on the median of the effect column being close to 1. Please ensure the interpretation of this column as an OR is correct.")),file=log.file,sep="\n",append=TRUE)
      if(a1 == a2) cat(print(paste("The effect column was determined NOT to be coded as an odds ratio (OR) for the", filenames[i], "summary statistics file based on the median of the effect column being close to 0.")),file=log.file,sep="\n",append=TRUE)
      
      
      ##remove any rows printed as exactly 0
      b<-nrow(files2)
      if("effect" %in% colnames(files2)){
        files2<-subset(files2, files2$effect != 0)
      }
      if(b-nrow(files2) > 0) cat(print(paste(b-nrow(files2), "rows were removed from the", filenames[i], "summary statistics file due to effect values estimated at exactly 0 as this causes problems for matrix inversion necessary for later Genomic SEM analyses.")),file=log.file,sep="\n",append=TRUE)
      
      if(OLS[i] == T){
        cat(print(paste("An OLS transformation is being used for file:", filenames[i])),file=log.file,sep="\n",append=TRUE)
        
        files2$Z <- sign(files2$effect) * sqrt(qchisq(files2$P,1,lower=F))
        
        if("N" %in% colnames(files2)){
          files2$effect <- files2$Z/ sqrt(files2$N * 2 * (files2$MAF *(1-files2$MAF)))}else{cat(print("ERROR: A Sample Size (N) is needed for OLS Standardization. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N."),file=log.file,sep="\n",append=TRUE)}}
      
      
      if(linprob[i] == T){
        cat(print(paste("An LPM transformation is being used for file:", filenames[i])),file=log.file,sep="\n",append=TRUE)
        
        files2$Z <- sign(files2$effect) * sqrt(qchisq(files2$P,1,lower=F))
        
        if("N" %in% colnames(files2)){
          files2$effect <- files2$Z/sqrt((prop[i]*(1-prop[i])*(2*files2$N*files2$MAF*(1-files2$MAF))))
          files2$SE<-1/sqrt((prop[i]*(1-prop[i])*(2*files2$N*files2$MAF*(1-files2$MAF))))}else{cat(print("ERROR: A Sample Size (N) is needed for LPM Standardization. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N."),file=log.file,sep="\n",append=TRUE)}}
      
      # Flip effect to match ordering in ref file
      files2$effect <-  ifelse(files2$A1.x != (files2$A1.y) & files2$A1.x == (files2$A2.y),files2$effect*-1,files2$effect)
      
      ##remove SNPs that don't match A1 or A2 in reference file.
      b<-nrow(files2)
      files2<-subset(files2, !(files2$A1.x != (files2$A1.y)  & files2$A1.x != (files2$A2.y)))
      if(b-nrow(files2) > 0) cat(print(paste(b-nrow(files2), "row(s) were removed from the", filenames[i], "summary statistics file due to the effect allele (A1) column not matching A1 or A2 in the reference file.")),file=log.file,sep="\n",append=TRUE)
      
      b<-nrow(files2)
      files2<-subset(files2, !(files2$A2.x != (files2$A2.y)  & files2$A2.x !=  (files2$A1.y)))
      if(b-nrow(files2) > 0) cat(print(paste(b-nrow(files2), "row(s) were removed from the", filenames[i], "summary statistics file due to the other allele (A2) column not matching A1 or A2 in the reference file.")),file=log.file,sep="\n",append=TRUE)
      
      #Check that p-value column does not contain an excess of 1s/0s
      if((sum(files2$P > 1) + sum(files2$P < 0)) > 100){
        cat(print("In excess of 100 SNPs have P val above 1 or below 0. The P column may be mislabled!"),file=log.file,sep="\n",append=TRUE)
      }
      
      if("INFO" %in% colnames(files2)) {
        b<-nrow(files2)
        files2 <- files2[files2$INFO >= info.filter,]
        cat(print(paste(b-nrow(files2), "rows were removed from the", filenames[i], "summary statistics file due to INFO values below the designated threshold of", info.filter)),file=log.file,sep="\n",append=TRUE)
      }else{cat(print("No INFO column, cannot filter on INFO, which may influence results"),file=log.file,sep="\n",append=TRUE)}
      
      varSNP<-2*files2$MAF*(1-files2$MAF)  
      
      if(OLS[i] == T){
        output <- cbind.data.frame(files2$SNP,
                                   files2$effect,
                                   abs(files2$effect/files2$Z))
        output<-na.omit(output)                           
        colnames(output) <- c("SNP",names.beta[i],names.se[i])                           
      }
      
      if(linprob[i] == T){
        output<-cbind.data.frame(files2$SNP,
                                 (files2$effect)/((files2$effect^2) * varSNP + (pi^2)/3)^.5,
                                 (files2$SE)/(((files2$effect)^2) * varSNP + (pi^2)/3)^.5)  
        output<-na.omit(output)
        output<-output[apply(output!=0, 1, all),]
        colnames(output) <- c("SNP",names.beta[i],names.se[i])                                         
      }
      
      if(linprob[i] == F){
        if(OLS[i] == F){                                     
          if(se.logit[i] == F){
            cat(print(paste("Performing transformation under the assumption that the effect column is either an odds ratio or logistic beta (please see output above to determine whether it was interpreted as an odds ratio) and the SE column is the SE of the odds ratio (i.e., NOT on the logistic scale) for:", filenames[i])),file=log.file,sep="\n",append=TRUE)
            
            if(sum(hold_names %in% "SE") == 0) cat(print(paste0('Cannot find SE column, try renaming it SE in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
            if(sum(hold_names %in% "SE") == 0) warning(paste0('Cannot find SE column, try renaming it SE in the summary statistics file for:',trait.names[i]))
            
            output <- cbind.data.frame(files2$SNP,
                                       (files2$effect)/((files2$effect^2) * varSNP + (pi^2)/3)^.5,
                                       (files2$SE/exp(files2$effect))/(((files2$effect)^2 * varSNP + (pi^2)/3)^.5))
            output<-na.omit(output)  
            colnames(output) <- c("SNP",names.beta[i],names.se[i])}}}
      
      if(se.logit[i]== T){
        cat(print(paste("Performing transformation under the assumption that the effect column is either an odds ratio or logistic beta (please see output above to determine whether it was interpreted as an odds ratio) and the SE column is a logistic SE (i.e., NOT the SE of the odds ratio) for:", filenames[i])),file=log.file,sep="\n",append=TRUE)
        
        if(sum(hold_names %in% "SE") == 0) cat(print(paste0('Cannot find SE column, try renaming it SE in the summary statistics file for:',trait.names[i])),file=log.file,sep="\n",append=TRUE)
        if(sum(hold_names %in% "SE") == 0) warning(paste0('Cannot find SE column, try renaming it SE in the summary statistics file for:',trait.names[i]))
        
        output <- cbind.data.frame(files2$SNP,
                                   (files2$effect)/((files2$effect^2) * varSNP + (pi^2)/3)^.5,
                                   (files2$SE)/(((files2$effect)^2) * varSNP + (pi^2)/3)^.5)  
        output<-na.omit(output)  
        colnames(output) <- c("SNP",names.beta[i],names.se[i])}
      
      cat(print(paste(nrow(output), "SNPs are left in the summary statistics file", filenames[i], "after QC and merging with the reference file.")),file=log.file,sep="\n",append=TRUE)
      
      if(mean(abs(output[,2]/output[,3])) > 5){
        cat(print(paste0('WARNING: The average value of estimate over standard error (i.e., Z) is > 5 for ',trait.names[i], ". This suggests a column was misinterpreted or arguments were misspecified. Please post on the google group if you are unable to figure out the issue.")),file=log.file,sep="\n",append=TRUE)
      }
      
      output
      
    },mc.cores=int)
    
    for(i in 1:length){
      if(i == 1){
        data.frame.out <- suppressWarnings(inner_join(data.frame.out,Output[[i]],by="SNP",all.x=F,all.y=F))
      }else{
        data.frame.out <- suppressWarnings(inner_join(data.frame.out,Output[[i]],by="SNP",all.x=F,all.y=F)) 
      }
    }
  }
  
  b<-nrow(data.frame.out)
  #data.frame.out<-data.frame.out[!duplicated(data.frame.out$BP),]
  
  end.time <- Sys.time()
  
  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time-mins*60
  
  if(parallel == FALSE){
    cat(paste("     "),file=log.file,sep="\n",append=TRUE)
    #cat(print(paste(b-nrow(data.frame.out), "rows were removed from the final summary statistics file due to duplicated base pair (BP) values")),file=log.file,sep="\n",append=TRUE)
    cat(print(paste0("After merging across all summary statistics using listwise deletion, performing QC, and merging with the reference file, there are ",nrow(data.frame.out), " SNPs left in the final multivariate summary statistics file"), sep = ""),file=log.file,sep="\n",append=TRUE)
    cat(print(paste0("Sumstats finished running at ",end.time), sep = ""),file=log.file,sep="\n",append=TRUE)
    cat(print(paste0("Running sumstats for all files took ",mins," minutes and ",secs," seconds"), sep = ""),file=log.file,sep="\n",append=TRUE)
    cat(print(paste("Please check the log file", paste0(log2, "_sumstats.log"), "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files.")),file=log.file,sep="\n",append=TRUE)
    flush(log.file)
    close(log.file)
  }
  
  if(parallel == TRUE){
    #print(paste(b-nrow(data.frame.out), "rows were removed from the final summary statistics file due to duplicated base pair (BP) values"))
    print(paste0("After merging across all summary statistics using listwise deletion, performing QC, and merging with the reference file, there are ",nrow(data.frame.out), " SNPs left in the final multivariate summary statistics file"), sep = "") 
    print(paste0("Sumstats finished running at ",end.time), sep = "")
    print(paste0("Running sumstats for all files took ",mins," minutes and ",secs," seconds"), sep = "")
    print(paste0("Please check the log files to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files."))
  }


  
  data.frame.out
  
}
