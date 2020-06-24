munge <- function(files,reference,trait.names=NULL,N,info.filter=.9,maf.filter=0.01,munge.log='munge.log') {
  
  sink(munge.log, append=FALSE, split=TRUE)
  
  length <- length(files)
  filenames <- as.vector(files)

  log.file <- file(munge.log,open="wt")
  
  begin.time <- Sys.time()
  
  cat( paste0( "The munging of ", length(trait.names), " summary statistics started at ", begin.time ) )
  
  cat( paste0( "Reading summary statistics for ", paste(files,collapse=" "), ".\n",
    "Please note that this step usually takes a few minutes due to the size of summary statistic files."
  ) )
  
  ##note that fread is not used here due to formatting differences across summary statistic files
  files = lapply(files, read.table,header=T, quote="\"",fill=T,na.string=c(".",NA,"NA",""))
  cat("Reading in reference file\n")
  ref <- fread(reference,header=T,data.table=F)
  cat("All files loaded into R!\n")
 
  for(i in 1:length) {
    
    cat(paste("\n\nMunging file:", filenames[i]))

    names(files[[i]]) <- GenomicSEM:::.sumstats.parse.header( files[[i]] )

    # Print a warning message when multiple columns are interpreted as p-values, effects, rsIDs or alleles columns
    if(sum(names(files[[i]]) %in% "P") > 1) warning(paste0(
      'Multiple columns are being interpreted as the P-value column. Try renaming the column you dont want interpreted"
      "as P, e.g., to P2 for:', filenames[i]))
    if(sum(names(files[[i]]) %in% "effect") > 1) warning(paste0(
      'Multiple columns are being interpreted as the effect column. Try renaming the column you dont want interpreted"
      "as effect, e.g., to effect2 for:', filenames[i]))
    if(sum(names(files[[i]]) %in% "SNP") > 1) warning(paste0(
      'Multiple columns are being interpreted as the variant ID column. Try renaming the column you dont want"
      "interpreted as variant ID, e.g., to SNP2 for:', filenames[i]))
    if(sum(names(files[[i]]) %in% "A1") > 1) warning(paste0(
      'Multiple columns are being interpreted as the effect allele column. Try renaming the column you dont want"
      "interpreted as effect allele, e.g., to A1_2 for:', filenames[i]))
    if(sum(names(files[[i]]) %in% "A2") > 1) warning(paste0(
      'Multiple columns are being interpreted as the non-effect allele column. Try renaming the column you dont want"
      "interpreted as non-effect allele, e.g., to A2_2 for:', filenames[i]))

    if("MAF" %in% colnames(files[[i]])) {
      ##make sure MAF is actually MAF (i.e., max value is .5 or less)
      files[[i]]$MAF<-ifelse(files[[i]]$MAF <= .5, files[[i]]$MAF, (1-files[[i]]$MAF))
    }
    
    # Compute N is N cases and N control is reported:
    if("N_CAS" %in% colnames(files[[i]])) {
      files[[i]]$N <- files[[i]]$N_CAS + files[[i]]$N_CON
      writeLines( strwrap( paste0(
        "As the file includes both N_CAS and N_CON columns, the summation of these two columns will be used as the",
        "total sample size"
      ) ) )
    }
    
    ##make sure all alleles are upper case for matching to reference file
    files[[i]]$A1 <- factor(toupper(files[[i]]$A1), c("A", "C", "G", "T"))
    files[[i]]$A2 <- factor(toupper(files[[i]]$A2), c("A", "C", "G", "T"))
    
    ##merge with ref file
    cat( "Merging file:", filenames[i], "with the reference file:", reference, "\n" )
    b<-nrow(files[[i]])
    cat( b, "rows present in the full", filenames[i], "summary statistics file.\n" )
    files[[i]] <- merge(ref,files[[i]],by="SNP",all.x=F,all.y=F)
    writeLines( strwrap( paste0(
      b-nrow(files[[i]]), " rows were removed from the ", filenames[i], " summary statistics file as the rs-ids for",
      "these rows were not present in the reference file."
    ) ) )
    
    ##remove any rows with missing p-values
    b<-nrow(files[[i]])
    if("P" %in% colnames(files[[i]])) {
      files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$P)))
    }
    if(b-nrow(files[[i]]) > 0) writeLines( strwrap( paste0(
      b-nrow(files[[i]]), " rows were removed from the ", filenames[i], " summary statistics file due to missing",
      "values in the P-value column."
    ) ) )
    
    ##remove any rows with missing effects
    b<-nrow(files[[i]])
    if("EFFECT" %in% colnames(files[[i]])) {
      files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$EFFECT)))
    }
    if(b-nrow(files[[i]]) > 0) writeLines( strwrap( paste0(
      b-nrow(files[[i]]), " rows were removed from the ", filenames[i], " summary statistics file due to missing",
      "values in the effect column."
    ) ) )
    
    ##determine whether it is OR or logistic/continuous effect based on median effect size
    EFFECT.untransf<-files[[i]]$EFFECT[[1]]
    files[[i]]$EFFECT<-ifelse(
      rep(round(median(files[[i]]$EFFECT,na.rm=T)) == 1,nrow(files[[i]])),
      log(files[[i]]$EFFECT),
      files[[i]]$EFFECT
    )
    EFFECT.transf<-files[[i]]$EFFECT[[1]]
    if(EFFECT.untransf != EFFECT.transf) writeLines( strwrap( paste0(
      "Values in the effect column identified as odds ratios (ORs) for the ", filenames[i], " summary statistics",
      "file. Please ensure that this is correct.",
    ) ) )
    
    # Flip effect to match ordering in ref file
    files[[i]]$EFFECT<-ifelse(
      files[[i]]$A1.x != (files[[i]]$A1.y) & files[[i]]$A1.x == (files[[i]]$A2.y),
      -files[[i]]$EFFECT,
      files[[i]]$EFFECT
    )
    
    ##remove SNPs that don't match A1 OR A2 in reference file.
    b<-nrow(files[[i]])
    files[[i]]<-subset(files[[i]], !(files[[i]]$A1.x != (files[[i]]$A1.y) & files[[i]]$A1.x != (files[[i]]$A2.y)))
    if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "row(s) were removed from the", filenames[i], "summary statistics file due to the effect allele (A1) column not matching A1 or A2 in the reference file.")),file=log.file,sep="\n",append=TRUE)
  
    b<-nrow(files[[i]])
    files[[i]]<-subset(files[[i]], !(files[[i]]$A2.x != (files[[i]]$A2.y) & files[[i]]$A2.x !=  (files[[i]]$A1.y)))
    if(b-nrow(files[[i]]) > 0) cat(print(paste(b-nrow(files[[i]]), "row(s) were removed from the", filenames[i], "summary statistics file due to the other allele (A2) column not matching A1 or A2 in the reference file.")),file=log.file,sep="\n",append=TRUE)
    
    ####VALIDITY CHECKS#####
    
    #Check that p-value column does not contain an excess of 1s/0s
    if((sum(files[[i]]$P > 1) + sum(files[[i]]$P < 0)) > 100){
      cat(print("In excess of 100 SNPs have P val above 1 or below 0. The P column may be mislabled!"),file=log.file,sep="\n",append=TRUE)
    }
   
    #Compute Z score
    files[[i]]$Z <- sign(files[[i]]$EFFECT) * sqrt(qchisq(files[[i]]$P,1,lower=F))
    
    ##filter on INFO column at designated threshold provided for the info.filter argument (default = 0.9)
    if("INFO" %in% colnames(files[[i]])) {
      b<-nrow(files[[i]])
      files[[i]] <- files[[i]][files[[i]]$INFO >= info.filter,]
      cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to INFO values below the designated threshold of", info.filter)),file=log.file,sep="\n",append=TRUE)
    }else{cat(print("No INFO column, cannot filter on INFO, which may influence results"),file=log.file,sep="\n",append=TRUE)}
    
    ##filter on MAF filter at designated threshold provided for the maf.filter argument (default = 0.01)
    if("MAF" %in% colnames(files[[i]])) {
      files[[i]]$MAF <- as.numeric(as.character(files[[i]]$MAF))
      b<-nrow(files[[i]])
      files[[i]] <- files[[i]][files[[i]]$MAF >= maf.filter,]
      files[[i]] <- subset(files[[i]], !(is.na(files[[i]]$MAF)))
      cat(print(paste(b-nrow(files[[i]]), "rows were removed from the", filenames[i], "summary statistics file due to missing MAF information or MAFs below the designated threshold of", maf.filter)),file=log.file,sep="\n",append=TRUE)
    }else{cat(print("No MAF column, cannot filter on MAF, which may influence results"),file=log.file,sep="\n",append=TRUE)}

    if("N" %in% colnames(files[[i]])) {
      output <- cbind.data.frame(files[[i]]$SNP,files[[i]]$N,files[[i]]$Z,files[[i]]$A1.x,files[[i]]$A2.x)
    }else{output <- cbind.data.frame(files[[i]]$SNP,N[i],files[[i]]$Z,files[[i]]$A1.x,files[[i]]$A2.x) }
    
    if(!("N" %in% names(files[[i]])) & (exists("N") == FALSE)) cat(warning(paste0('Cannot find sample size column for',filenames[i], " and a sample size was not provided for the N argument. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N.")),file=log.file,sep="\n",append=TRUE)

    colnames(output) <- c("SNP","N","Z","A1","A2")
    cat( paste(nrow(output), "SNPs are left from file", filenames[i], "after QC."), sep="\n" )
    
    write.table(x = output, file = paste0(trait.names[i],".sumstats"), sep = "\t", quote = FALSE, row.names = F)
    gzip(paste0(trait.names[i],".sumstats"))
    cat( paste("I am done munging file:", filenames[i]), sep="\n" )
    cat(
      paste("The munged summary statistics are saved as", paste0(trait.names[i],".sumstats.gz")),
      "in the current working directory.",
      sep="\n"
    )
  }
  
  end.time <- Sys.time()
  
  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time-mins*60
  
  cat( "     ", paste( "Munging completed at", end.time), sep="\n" )
  cat( paste( "Munging of all files took", mins, "minutes and", secs, "seconds"), sep="\n" )
  cat(
    paste( "Please check the log file", munge.log, "to ensure that all columns were interpreted" ),
    "correctly and no warnings were issued for any of the summary statistics files.",
    sep="\n"
  )
  
  sink()

}

