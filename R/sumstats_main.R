.sumstats_main <- function(X, utilfuncs, filename, trait.name, N, keep.indel, OLS, beta, info.filter, linprob,
                           se.logit, name.beta, name.se, ref, ref2, file=NULL, log.file=NULL,direct.filter) {
  if (!is.null(file)) {
    file <- data.frame(file)
  } else{
    file <- data.frame(read.table(filename, header = T, quote="\"",fill=T,na.string=c(".",NA,"NA","")))
  }
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
      assign(j, utilfuncs[[j]], envir=environment())
    }
  }
  if (is.null(log.file)) {
    log.file <- base::file(paste0(trait.name, "_sumstats.log"),open="wt")
    on.exit(flush(log.file))
    on.exit(close(log.file))
  } else {
    .LOG("\n\n",file=log.file, print=FALSE)
  }
  
  .LOG("Preparing summary statistics for file: ", filename,file=log.file)
  N_provided <- (!is.na(N))
  if (linprob){
    stop_on_missing <- c("effect", "SNP")
  }  else if ((se.logit) | (all(!(c(linprob, OLS, se.logit))))) {
    # if se.logit == T, or linprob,OLS,se.logit are all FALSE, SE is required.
    stop_on_missing <- c("effect", "SNP", "SE")
  } else {
    stop_on_missing <- c("effect", "SNP")
  }
  hold_names <- .get_renamed_colnames(toupper(names(file)),
                                      userprovided=list(), checkforsingle=c("P", "A1", "A2", "effect", "SNP"),
                                      N_provided=FALSE, filename, log.file, warnz=((!linprob) & (!OLS)),
                                      warn_for_missing = c("P", "A1", "A2", "N"),
                                      stop_on_missing = stop_on_missing,
                                      utilfuncs)
  colnames(file) <- hold_names
  if (N_provided) {
    file$N <- N
    if(OLS){
      .LOG("Using user provided N of ", N, " for ", filename, " . Please note that this should reflect the total sample size.",file=log.file)
    } else {
      .LOG("Using user provided N of ", N, " for ", filename, " . Please note that this should reflect the sum of effective sample sizes if the linprob argument is being used to back out logistic betas.",file=log.file)
    }
  }
  
  b<-nrow(file)
  file<-file[!(duplicated(file$SNP) | duplicated(file$SNP, fromLast = TRUE)), ]
  .LOG((b-nrow(file)), " rows were removed from the ", filename, " summary statistics file due to entries that were duplicated for rsID. These are removed as they likely reflect multiallelic variants.",file=log.file)
  
  if(keep.indel){
    file$A1 <- factor(toupper(file$A1))
    file$A2 <- factor(toupper(file$A2))
    .LOG("Keeping variants other than SNPs, this may cause problems when alligning alleles across traits and the reference file",file=log.file)
  } else {
    ##make sure all alleles are upper case for matching
    file$A1 <- factor(toupper(file$A1), c("A", "C", "G", "T"))
    file$A2 <- factor(toupper(file$A2), c("A", "C", "G", "T"))
  }
  
  if(direct.filter){
    if("DIRECTION" %in% colnames(file)) {
      b<-nrow(file)
      file$DIRECTION<-as.character(file$DIRECTION)
      file<-subset(file, (str_count(file$DIRECTION, "\\?")/nchar(file$DIRECTION)) < .5)
      .LOG((b-nrow(file)), "rows were removed from the", filename, "summary statistics file due to missingness across 50% or more of contributing cohorts as determined by missing info in the direction column.",file=log.file)
    }else{.LOG("No DIRECTION column, cannot filter on missingness by conributing cohorts, which may influence results",file=log.file)}
  }
  
  ##merge with ref file
  .LOG("Merging file: ", filename, " with the reference file: ", ref2,file=log.file)
  b<-nrow(file)
  .LOG(b, " rows present in the full ", filename, " summary statistics file.",file=log.file)
  file <- suppressWarnings(inner_join(ref,file,by="SNP"))
  .LOG((b-nrow(file)), " rows were removed from the ", filename, " summary statistics file as the rsIDs for these SNPs were not present in the reference file.",file=log.file)
  
  ##remove any rows with missing p-values
  b<-nrow(file)
  if("P" %in% colnames(file)) {
    file<-subset(file, !(is.na(file$P)))
  }
  if(b-nrow(file) > 0) .LOG(b-nrow(file), "rows were removed from the ", filename, " summary statistics file due to missing values in the P-value column",file=log.file)
  
  ##remove any rows with missing effects
  b<-nrow(file)
  if("effect" %in% colnames(file)) {
    file<-subset(file, !(is.na(file$effect)))
  }
  if(b-nrow(file) > 0) .LOG(b-nrow(file), "rows were removed from the ", filename, " summary statistics file due to missing values in the effect column",file=log.file)
  
  #use sample specific MAF for later conversions when possible; otherwise use ref MAF
  if("MAF.y" %in% colnames(file)){
    file$MAF.y<-ifelse(file$MAF.y > .5, 1-file$MAF.y, file$MAF.y)
    b<-nrow(file)
    file<-subset(file, file$MAF.y != 0 & file$MAF.y != 1)
    if(b-nrow(file) > 0) .LOG(b-nrow(file), " rows were removed from the ", filename, " summary statistics file due to allele frequencies printed as exactly 1 or 0", file=log.file)
    file$varSNP<-2*file$MAF.y*(1-file$MAF.y)
  }else{
    file$varSNP<-2*file$MAF*(1-file$MAF)
  }
  
  ##determine whether it is OR or logistic/continuous effect based on median effect size
  a1<-file$effect[[1]]
  file$effect <- ifelse(rep(round(median(file$effect,na.rm=T)) == 1,nrow(file)), log(file$effect),file$effect)
  a2<-file$effect[[1]]
  if(a1 != a2) .LOG("The effect column was determined to be coded as an odds ratio (OR) for the ", filename, " summary statistics file based on the median of the effect column being close to 1. Please ensure the interpretation of this column as an OR is correct.",file=log.file)
  if(a1 == a2) .LOG("The effect column was determined NOT to be coded as an odds ratio (OR) for the ", filename, " summary statistics file based on the median of the effect column being close to 0.",file=log.file)
  
  ##remove any rows printed as exactly 0
  b<-nrow(file)
  if("effect" %in% colnames(file)) {
    file<-subset(file, file$effect != 0)
  }
  if(b-nrow(file) > 0) .LOG(b-nrow(file), "rows were removed from the", filename, "summary statistics file due to effect values estimated at exactly 0 as this causes problems for matrix inversion necessary for later Genomic SEM analyses.",file=log.file)
  
  file$Z <- sign(file$effect) * sqrt(qchisq(file$P,1,lower=F))
  
  if(OLS & beta){
    .LOG("User provided arguments indicate that a GWAS of a continuous trait with already standardized betas is being provided for: ", filename,file=log.file)
  }
  
  if(OLS & !beta){
    if("N" %in% colnames(file)){
      file$effect <- file$Z/sqrt(file$N * file$varSNP)
    }else{
      .LOG("ERROR: A Sample Size (N) is needed for OLS Standardization. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N.",file=log.file, print=FALSE)
    }
  }
  
  if(linprob){
    .LOG("An transformation used to back out logistic betas for binary traits is being applied for: ", filename,file=log.file)
    if("N" %in% colnames(file)){
      file$effect <- file$Z/sqrt((file$N/4)*file$varSNP)
      file$SE<-1/sqrt((file$N/4)*file$varSNP)
    }else{
      .LOG("ERROR: An effective sample Size (N) is needed for backing out betas for binary traits. Please provide the sum of effective sample sizes to the N argument.",file=log.file)
    }
  }
  
  # Flip effect to match ordering in ref file
  file$effect <-  ifelse(file$A1.x != (file$A1.y) & file$A1.x == (file$A2.y),file$effect*-1,file$effect)
  
  ##remove SNPs that don't match A1 OR A2 in ref.
  b<-nrow(file)
  file<-subset(file, !(file$A1.x != (file$A1.y)  & file$A1.x != (file$A2.y)))
  if(b-nrow(file) > 0) .LOG(b-nrow(file), " row(s) were removed from the" , filename, " summary statistics file due to the effect allele (A1) column not matching A1 or A2 in the reference file.",file=log.file)
  
  b<-nrow(file)
  file<-subset(file, !(file$A2.x != (file$A2.y)  & file$A2.x !=  (file$A1.y)))
  if(b-nrow(file) > 0) .LOG(b-nrow(file), " row(s) were removed from the ", filename, " summary statistics file due to the other allele (A2) column not matching A1 or A2 in the reference file.",file=log.file)
  
  #Check that p-value column does not contain an excess of 1s/0s
  if((sum(file$P > 1) + sum(file$P < 0)) > 100){
    .LOG("In excess of 100 SNPs have P val above 1 or below 0. The P column may be mislabled!",file=log.file)
  }
  
  if("INFO" %in% colnames(file)) {
    b<-nrow(file)
    file <- file[file$INFO >= info.filter,]
    .LOG(b-nrow(file), "rows were removed from the ", filename, " summary statistics file due to INFO values below the designated threshold of ", info.filter,file=log.file)
  }else{.LOG("No INFO column, cannot filter on INFO, which may influence results",file=log.file)}
  
  if(OLS){
    output <- cbind.data.frame(file$SNP,
                               file$effect,
                               abs(file$effect/file$Z))
    output<-na.omit(output)
    colnames(output) <- c("SNP",name.beta,name.se)
  }
  
  if(linprob){
    output<-cbind.data.frame(file$SNP,
                             (file$effect)/((file$effect^2) * file$varSNP + (pi^2)/3)^.5,
                             (file$SE)/(((file$effect)^2) * file$varSNP + (pi^2)/3)^.5)
    output<-na.omit(output)
    output<-output[apply(output!=0, 1, all),]
    colnames(output) <- c("SNP",name.beta,name.se)
  }
  
  if(!linprob & !OLS & !se.logit){
    .LOG("Performing transformation under the assumption that the effect column is either an odds ratio or logistic beta (please see output above to determine whether it was interpreted as an odds ratio) and the SE column is the SE of the odds ratio (i.e., NOT on the logistic scale) for:", filename,file=log.file)
    
    output <- cbind.data.frame(file$SNP,
                               (file$effect)/((file$effect^2) * file$varSNP + (pi^2)/3)^.5,
                               (file$SE/exp(file$effect))/(((file$effect)^2 * file$varSNP + (pi^2)/3)^.5))
    output<-na.omit(output)
    colnames(output) <- c("SNP",name.beta,name.se)
  }
  if(se.logit){
    .LOG("Performing transformation under the assumption that the effect column is either an odds ratio or logistic beta (please see output above to determine whether it was interpreted as an odds ratio) and the SE column is a logistic SE (i.e., NOT the SE of the odds ratio) for:", filename,file=log.file)
    
    output <- cbind.data.frame(file$SNP,
                               (file$effect)/((file$effect^2) * file$varSNP + (pi^2)/3)^.5,
                               (file$SE)/(((file$effect)^2) * file$varSNP + (pi^2)/3)^.5)
    output<-na.omit(output)
    colnames(output) <- c("SNP",name.beta,name.se)
  }
  
  .LOG(nrow(output), " SNPs are left in the summary statistics file ", filename, " after QC and merging with the reference file.",file=log.file)
  
  if(mean(abs(output[,2]/output[,3])) > 5){
    .LOG('WARNING: The average value of estimate over standard error (i.e., Z) is > 5 for ',trait.name, ". This suggests a column was misinterpreted or arguments were misspecified. Please post on the google group if you are unable to figure out the issue.",file=log.file, print=FALSE)
    warning(paste0('The average value of estimate over standard error (i.e., Z) is > 5 for ',trait.name, ". This suggests a column was misinterpreted or arguments were misspecified. Please post on the google group if you are unable to figure out the issue."))
  }
  return(output)
}
