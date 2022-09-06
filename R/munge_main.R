.munge_main <- function(i, utilfuncs, file, filename, trait.name, N, ref, hm3, info.filter, maf.filter, column.names, overwrite, log.file=NULL) {
  if (is.null(log.file)) {
    log.file <- file(paste0(trait.name, "_munge.log"),open="wt")
    on.exit(flush(log.file))
    on.exit(close(log.file))
  } else {
    .LOG("\n\n",file=log.file, print=FALSE)
  }
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
        assign(j, utilfuncs[[j]], envir=environment())
    }
  }
  if (is.null(file)) {
    file <- read.table(filename, header=T, quote="\"", fill=T, na.string=c(".", NA, "NA", ""))
  }
  .LOG("Munging file: ", filename,file=log.file, print=TRUE)

  N_provided <- (!is.na(N))
  
  colnames(file)<-toupper(colnames(file))
if("NEFFDIV2" %in% colnames(file)) {
  .LOG("Found an NEFFDIV2 column for sample size. \n
Please note that this is likely effective sample size cut in half. The function is automatically doubling this value. This should only be used for liability h^2 conversion for binary traits and that it should reflect the sum of effective sample sizes across cohorts.", file=log.file)
  file$NEFF<-file$NEFFDIV2*2
  file$NEFFDIV2<-NULL
}

if("NEFF_HALF" %in% colnames(file)) {
  .LOG("Found an NEFF_HALF column for sample size. \n
Please note that this is likely effective sample size cut in half. The function is automatically doubling this value. This should only be used for liability h^2 conversion for binary traits and that it should reflect the sum of effective sample sizes across cohorts.", file=log.file)
  file$NEFF<-file$NEFF_HALF*2
  file$NEFF_HALF<-NULL
}
  hold_names <- .get_renamed_colnames(toupper(names(file)),
                                      column.names, c("P", "A1", "A2", "effect", "SNP"), filename,
                                      N_provided=N_provided, log.file=log.file,
                                      warn_for_missing=c("P", "A1", "A2", "effect", "SNP", "N"),
                                      utilfuncs=utilfuncs)
  colnames(file) <- hold_names
  if (N_provided) {
    file$N <- N
    .LOG("Using provided N (",N,") for file:",filename, file=log.file)
  }

  if("MAF" %in% colnames(file)) {
    ##make sure MAF is actually MAF (i.e., max value is .5 or less)
    file$MAF<-ifelse(file$MAF <= .5, file$MAF, (1-file$MAF))
  }

  ##make sure all alleles are upper case for matching to reference file
  file$A1 <- factor(toupper(file$A1), c("A", "C", "G", "T"))
  file$A2 <- factor(toupper(file$A2), c("A", "C", "G", "T"))
  
  ##merge with ref file
  .LOG("Merging file:", filename, " with the reference file:", hm3,file=log.file)
  b <- nrow(file)
  .LOG(b, " rows present in the full ", filename, " summary statistics file.",file=log.file)
  file <- merge(ref,file,by="SNP",all.x=F,all.y=F)
  .LOG((b-nrow(file)), " rows were removed from the ", filename, " summary statistics file as the rs-ids for these rows were not present in the reference file.",file=log.file)
  
  ##remove any rows with missing p-values
  b<-nrow(file)
  if("P" %in% colnames(file)) {
    file<-subset(file, !(is.na(file$P)))
  }
  if(b-nrow(file) > 0) .LOG(b-nrow(file), " rows were removed from the ", filename, " summary statistics file due to missing values in the P-value column",file=log.file)
  
  ##remove any rows with missing effects
  b<-nrow(file)
  if("effect" %in% colnames(file)) {
    file<-subset(file, !(is.na(file$effect)))
  }
  if(b-nrow(file) > 0) .LOG(b-nrow(file), " rows were removed from the ", filename, " summary statistics file due to missing values in the effect column",file=log.file)
  
  ##determine whether it is OR or logistic/continuous effect based on median effect size 
  a1<-file$effect[[1]]
  file$effect<-ifelse(rep(round(median(file$effect,na.rm=T)) == 1,nrow(file)), log(file$effect),file$effect)
  a2<-file$effect[[1]]
  if(a1 != a2) .LOG("The effect column was determined to be coded as an odds ratio (OR) for the ", filename, " summary statistics file. Please ensure this is correct.",file=log.file)
  
  # Flip effect to match ordering in ref file
  file$effect<-ifelse(file$A1.x != (file$A1.y) & file$A1.x == (file$A2.y),file$effect*-1,file$effect)
  
  ##remove SNPs that don't match A1 OR A2 in reference file.
  b<-nrow(file)
  file<-subset(file, !(file$A1.x != (file$A1.y)  & file$A1.x != (file$A2.y)))
  if(b-nrow(file) > 0) .LOG(b-nrow(file), " row(s) were removed from the ", filename, " summary statistics file due to the effect allele (A1) column not matching A1 or A2 in the reference file.",file=log.file)

  b<-nrow(file)
  file<-subset(file, !(file$A2.x != (file$A2.y)  & file$A2.x !=  (file$A1.y)))
  if(b-nrow(file) > 0) .LOG(b-nrow(file), " row(s) were removed from the ", filename, " summary statistics file due to the other allele (A2) column not matching A1 or A2 in the reference file.",file=log.file)
  
  ####VALIDITY CHECKS#####
  #Check that p-value column does not contain an excess of 1s/0s
  if((sum(file$P > 1) + sum(file$P < 0)) > 100){
    .LOG("In excess of 100 SNPs have P val above 1 or below 0. The P column may be mislabled!",file=log.file)
  }
 
  #Compute Z score
  file$Z <- sign(file$effect) * sqrt(qchisq(file$P,1,lower=F))
  
  ##filter on INFO column at designated threshold provided for the info.filter argument (default = 0.9)
  if("INFO" %in% colnames(file)) {
    b<-nrow(file)
    file <- file[file$INFO >= info.filter,]
    .LOG(b-nrow(file), " rows were removed from the ", filename, " summary statistics file due to INFO values below the designated threshold of", info.filter,file=log.file)
  }else{.LOG("No INFO column, cannot filter on INFO, which may influence results",file=log.file)}
  
  ##filter on MAF filter at designated threshold provided for the maf.filter argument (default = 0.01)
  if("MAF" %in% colnames(file)) {
    file$MAF<-as.numeric(as.character(file$MAF))
    b<-nrow(file)
    file <- file[file$MAF >= maf.filter,]
    file<-subset(file, !(is.na(file$MAF)))
    .LOG(b-nrow(file), " rows were removed from the ", filename, " summary statistics file due to missing MAF information or MAFs below the designated threshold of", maf.filter,file=log.file)
  }else{
    .LOG("No MAF column, cannot filter on MAF, which may influence results",file=log.file)
  }

  if("N" %in% colnames(file)) {
    output <- cbind.data.frame(file$SNP,file$N,file$Z,file$A1.x,file$A2.x)
  }else{output <- cbind.data.frame(file$SNP,N,file$Z,file$A1.x,file$A2.x) }
  
  if(!("N" %in% names(file)) & (exists("N") == FALSE)) .LOG('Cannot find sample size column for',filename, " and a sample size was not provided for the N argument. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N.",file=log.file)

  colnames(output) <- c("SNP","N","Z","A1","A2")
  .LOG(nrow(output), "SNPs are left in the summary statistics file ", filename, " after QC.",file=log.file)
  
  #remove spaces in trait.names file to avoid errors with fread functionality used for s_ldsc
  trait.name<-str_replace_all(trait.name, fixed(" "), "") 
  
  write.table(x = output,file = paste0(trait.name,".sumstats"),sep="\t", quote = FALSE, row.names = F)
  gzip(paste0(trait.name,".sumstats"), overwrite=overwrite)
  .LOG("I am done munging file: ", filename,file=log.file)
  .LOG("The file is saved as ", paste0(trait.name,".sumstats.gz"), " in the current working directory.",file=log.file)
  return()
}
