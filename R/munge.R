munge <- function(files,hm3,trait.names=NULL,N=NULL,info.filter = .9,maf.filter=0.01,log.name=NULL, column.names=list(), overwrite=TRUE){
  # Sanity checks
  .check_file_exists(files)
  .check_file_exists(hm3)
  .check_equal_length(files, trait.names)
  .check_equal_length(files, N)
  if (!is.null(N)) {.check_range(N, min=0, max=Inf, allowNA=TRUE)}
  .check_range(info.filter)
  .check_range(maf.filter)
  if (any(!(names(column.names) %in% c("SNP", "A1", "A2", "effect", "INFO", "P", "N", "N_CAS", "N_CON", "MAF", "Z")))) {
    stop(paste0("Names in column.names not recognized. Please use the following keys:\n        ",
                paste(c("SNP", "A1", "A2", "effect", "INFO", "P", "N", "N_CAS", "N_CON", "MAF", "Z"), collapse=", ")))
  }
  .check_boolean(overwrite)
  # Sanity checks finished

  filenames <- as.vector(files)
  if (is.null(N))  {
    N <- rep(NA, length(files))
  }
  if(is.null(log.name)){
    log2<-paste(trait.names,collapse="_")
      if(nchar(log2) > 200){
      log2<-substr(log2,1,100)}
  log.file <- file(paste0(log2, "_munge.log"),open="wt")
  }
  
  if(!is.null(log.name)){ 
   log.file <- file(paste0(log.name, "_munge.log"),open="wt") 
  }
  
  begin.time <- Sys.time()
  
  .LOG("The munging of ", length(trait.names), " summary statistics started at ", begin.time, file=log.file)
  if (overwrite) {
    existing_files <- c()
    for (trait.name in trait.names) {
      if (file.exists(paste0(trait.name, ".sumstats")))
        existing_files <- c(existing_files, paste0(trait.name, ".sumstats"))
    }
    if (length(existing_files) > 0)
      .LOG("File(s) ", paste0(existing_files, collapse = ", "), " already exist and will be overwritten", file=log.file)
  }
  .LOG("Reading summary statistics for ", paste(files,collapse=" "), ". Please note that this step usually takes a few minutes due to the size of summary statistic files.", file=log.file)
  
  ##note that fread is not used here due to formatting differences across summary statistic files
  files <- lapply(files, read.table, header=T, quote="\"", fill=T, na.string=c(".", NA, "NA", ""))
  .LOG("Reading in reference file",file=log.file)
  ref <- fread(hm3,header=T,data.table=F)
  .LOG("All files loaded into R!",file=log.file)
 
  for(i in 1:length(files)){
    
    .LOG("\n\n",file=log.file, print=FALSE)
    
    .LOG("Munging file: ", filenames[i],file=log.file, print=TRUE)
    if (!is.na(N[i])) {
      N_provided <- TRUE
    }
    hold_names <- .get_renamed_colnames(toupper(names(files[[i]])),
                                        column.names, c("P", "A1", "A2", "effect", "SNP"), filenames[i], N_provided=N_provided, log.file)
    colnames(files[[i]]) <- hold_names
    if (N_provided) {
      files[[i]]$N <- N[i]
      .LOG("Using provided N (",N,") for file:",filenames[i], file=log.file)
    }

    if("MAF" %in% colnames(files[[i]])) {
      ##make sure MAF is actually MAF (i.e., max value is .5 or less)
      files[[i]]$MAF<-ifelse(files[[i]]$MAF <= .5, files[[i]]$MAF, (1-files[[i]]$MAF))
    }

    ##make sure all alleles are upper case for matching to reference file
    files[[i]]$A1 <- factor(toupper(files[[i]]$A1), c("A", "C", "G", "T"))
    files[[i]]$A2 <- factor(toupper(files[[i]]$A2), c("A", "C", "G", "T"))
    
    ##merge with ref file
    .LOG("Merging file:", filenames[i], " with the reference file:", hm3,file=log.file)
    b<-nrow(files[[i]])
    .LOG(b, " rows present in the full ", filenames[i], " summary statistics file.",file=log.file)
    files[[i]] <- merge(ref,files[[i]],by="SNP",all.x=F,all.y=F)
    .LOG((b-nrow(files[[i]])), " rows were removed from the ", filenames[i], " summary statistics file as the rs-ids for these rows were not present in the reference file.",file=log.file)
    
    ##remove any rows with missing p-values
    b<-nrow(files[[i]])
    if("P" %in% colnames(files[[i]])) {
      files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$P)))
    }
    if(b-nrow(files[[i]]) > 0) .LOG(b-nrow(files[[i]]), " rows were removed from the ", filenames[i], " summary statistics file due to missing values in the P-value column",file=log.file)
    
    ##remove any rows with missing effects
    b<-nrow(files[[i]])
    if("effect" %in% colnames(files[[i]])) {
      files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$effect)))
    }
    if(b-nrow(files[[i]]) > 0) .LOG(b-nrow(files[[i]]), " rows were removed from the ", filenames[i], " summary statistics file due to missing values in the effect column",file=log.file)
    
    ##determine whether it is OR or logistic/continuous effect based on median effect size 
    a1<-files[[i]]$effect[[1]]
    files[[i]]$effect<-ifelse(rep(round(median(files[[i]]$effect,na.rm=T)) == 1,nrow(files[[i]])), log(files[[i]]$effect),files[[i]]$effect)
    a2<-files[[i]]$effect[[1]]
    if(a1 != a2) .LOG("The effect column was determined to be coded as an odds ratio (OR) for the ", filenames[i], " summary statistics file. Please ensure this is correct.",file=log.file)
    
    # Flip effect to match ordering in ref file
    files[[i]]$effect<-ifelse(files[[i]]$A1.x != (files[[i]]$A1.y) & files[[i]]$A1.x == (files[[i]]$A2.y),files[[i]]$effect*-1,files[[i]]$effect)
    
    ##remove SNPs that don't match A1 OR A2 in reference file.
    b<-nrow(files[[i]])
    files[[i]]<-subset(files[[i]], !(files[[i]]$A1.x != (files[[i]]$A1.y)  & files[[i]]$A1.x != (files[[i]]$A2.y)))
    if(b-nrow(files[[i]]) > 0) .LOG(b-nrow(files[[i]]), " row(s) were removed from the ", filenames[i], " summary statistics file due to the effect allele (A1) column not matching A1 or A2 in the reference file.",file=log.file)
  
    b<-nrow(files[[i]])
    files[[i]]<-subset(files[[i]], !(files[[i]]$A2.x != (files[[i]]$A2.y)  & files[[i]]$A2.x !=  (files[[i]]$A1.y)))
    if(b-nrow(files[[i]]) > 0) .LOG(b-nrow(files[[i]]), " row(s) were removed from the ", filenames[i], " summary statistics file due to the other allele (A2) column not matching A1 or A2 in the reference file.",file=log.file)
    
    ####VALIDITY CHECKS#####
    
    #Check that p-value column does not contain an excess of 1s/0s
    if((sum(files[[i]]$P > 1) + sum(files[[i]]$P < 0)) > 100){
      .LOG("In excess of 100 SNPs have P val above 1 or below 0. The P column may be mislabled!",file=log.file)
    }
   
    #Compute Z score
    files[[i]]$Z <- sign(files[[i]]$effect) * sqrt(qchisq(files[[i]]$P,1,lower=F))
    
    ##filter on INFO column at designated threshold provided for the info.filter argument (default = 0.9)
    if("INFO" %in% colnames(files[[i]])) {
      b<-nrow(files[[i]])
      files[[i]] <- files[[i]][files[[i]]$INFO >= info.filter,]
      .LOG(b-nrow(files[[i]]), " rows were removed from the ", filenames[i], " summary statistics file due to INFO values below the designated threshold of", info.filter,file=log.file)
    }else{.LOG("No INFO column, cannot filter on INFO, which may influence results",file=log.file)}
    
    ##filter on MAF filter at designated threshold provided for the maf.filter argument (default = 0.01)
    if("MAF" %in% colnames(files[[i]])) {
      files[[i]]$MAF<-as.numeric(as.character(files[[i]]$MAF))
      b<-nrow(files[[i]])
      files[[i]] <- files[[i]][files[[i]]$MAF >= maf.filter,]
      files[[i]]<-subset(files[[i]], !(is.na(files[[i]]$MAF)))
      .LOG(b-nrow(files[[i]]), " rows were removed from the ", filenames[i], " summary statistics file due to missing MAF information or MAFs below the designated threshold of", maf.filter,file=log.file)
    }else{
      .LOG("No MAF column, cannot filter on MAF, which may influence results",file=log.file)
    }

    if("N" %in% colnames(files[[i]])) {
      output <- cbind.data.frame(files[[i]]$SNP,files[[i]]$N,files[[i]]$Z,files[[i]]$A1.x,files[[i]]$A2.x)
    }else{output <- cbind.data.frame(files[[i]]$SNP,N[i],files[[i]]$Z,files[[i]]$A1.x,files[[i]]$A2.x) }
    
    if(!("N" %in% names(files[[i]])) & (exists("N") == FALSE)) .LOG('Cannot find sample size column for',filenames[i], " and a sample size was not provided for the N argument. Please either provide a total sample size to the N argument or try changing the name of the sample size column to N.",file=log.file)

    colnames(output) <- c("SNP","N","Z","A1","A2")
    .LOG(nrow(output), "SNPs are left in the summary statistics file", filenames[i], "after QC.",file=log.file)
    
    #remove spaces in trait.names file to avoid errors with fread functionality used for s_ldsc
    trait.names[i]<-str_replace_all(trait.names[i], fixed(" "), "") 
    
    write.table(x = output,file = paste0(trait.names[i],".sumstats"),sep="\t", quote = FALSE, row.names = F)
    gzip(paste0(trait.names[i],".sumstats"), overwrite=overwrite)
    .LOG("I am done munging file:", filenames[i],file=log.file)
    .LOG("The file is saved as", paste0(trait.names[i],".sumstats.gz"), "in the current working directory.",file=log.file)
  }
  
  end.time <- Sys.time()
  
  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time-mins*60
  
  .LOG("     ",file=log.file)
  .LOG("Munging was completed at ",end.time,file=log.file)
  
  .LOG("The munging of all files took ",mins," minutes and ",secs," seconds",file=log.file)
  .LOG("Please check the .log file to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files",file=log.file)
    
  flush(log.file)
  close(log.file)
    
}
