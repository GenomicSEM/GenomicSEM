sumstats <- function(files,ref,trait.names=NULL,se.logit,OLS=NULL,linprob=NULL,N=NULL,betas=NULL,
                     info.filter = .6,maf.filter=0.01,keep.indel=FALSE,parallel=FALSE,cores=NULL){
  len <- length(files)
  if (is.null(N)) N <- rep(NA, len)
  if(is.null(OLS)){
    OLS<-rep(FALSE,len)
  }
  if(is.null(linprob)){
    linprob<-rep(FALSE,len)
  }
  if(is.null(betas)){
    betas<-rep(FALSE,len)
  }
  # Sanity checks
  .check_file_exists(files)
  .check_file_exists(ref)
  .check_equal_length(files, trait.names)
  .check_equal_length(files, OLS)
  .check_equal_length(files, linprob)
  .check_equal_length(files, N)
  .check_equal_length(files, betas)
  .check_range(info.filter)
  .check_range(maf.filter)
  if (!is.null(N)) {.check_range(N, min=0, max=Inf, allowNA=TRUE)}
  .check_boolean(keep.indel)
  .check_boolean(parallel)
  .check_range(cores, min=0, max=Inf)
  # Sanity checks finished


  begin.time <- Sys.time()

  filenames <- as.vector(files)

  ref2<-ref

  if(is.null(trait.names)){
    names.beta <- paste0("beta.",1:len)
    names.se <- paste0("se.",1:len)
  }else{
    names.beta <- paste0("beta.",trait.names)
    names.se <- paste0("se.",trait.names)
  }

  log2<-paste(trait.names,collapse="_")
  log2<-str_remove_all(log2, "/")

  #subset log name to first 200 characters
  if(object.size(log2) > 200){
    log2<-substr(log2,1,100)
  }

  log.file <- file(paste0(log2, "_sumstats.log"),open="wt")

  .LOG("The preparation of ", length(trait.names), " summary statistics for use in Genomic SEM began at: ",begin.time,file=log.file)
  .LOG("Please note that the files should be in the same order that they were listed for the ldsc function",file=log.file)

  .LOG("Reading in reference file",file=log.file,append=TRUE)
  ref <- fread(ref,header=T,data.table=F)

  ##filter ref file on user provided maf.filter
  .LOG("Applying MAF filer of", maf.filter, "to the reference file.",file=log.file,append=TRUE)
  ref<-subset(ref, ref$MAF >= maf.filter)

  data.frame.out <- ref


  if(parallel == FALSE){
    ##note that fread is not used here as we have observed different formatting for column headers causing mismatched columns
    files <- lapply(files, read.table, header=T, quote="\"", fill=T, na.string=c(".", NA, "NA", ""))

    .LOG("All files loaded into R!",file=log.file,append=TRUE)
    Output <- list()
    for(i in 1:len){
      Output[[i]] <- .sumstats_main(i, filenames, trait.names, N, keep.indel, OLS, betas, info.filter, linprob, se.logit, names.beta, names.se, ref, ref2, files, log.file)
    }
  }

  if(parallel == TRUE){
    if(is.null(cores)){
      ##if no default provided use 1 less than the total number of cores available so your computer will still function
      int <- detectCores() - 1
    }else{int<-cores}
    #if more cores than traits then set to number of traits
    if(int > len){
      int<-len
    }

    print("Performing conversions of individual summary statistics using parallel processing. Please note this step typically takes 10-20 minutes due to the size of the files.")
    Output<-mclapply(X=1:len,FUN=.sumstats_main,filenames, trait.names, N, keep.indel, OLS, betas, info.filter, linprob, se.logit, names.beta, names.se, ref, ref2,
                     mc.cores=int)
  }
  for(i in 1:len){
    data.frame.out <- suppressWarnings(inner_join(data.frame.out,Output[[i]],by="SNP",all.x=F,all.y=F))
  }

  b<-nrow(data.frame.out)
  #data.frame.out<-data.frame.out[!duplicated(data.frame.out$BP),]

  end.time <- Sys.time()

  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time-mins*60

  .LOG("     ",file=log.file,append=TRUE, print=FALSE)
  #cat(print(paste(b-nrow(data.frame.out), "rows were removed from the final summary statistics file due to duplicated base pair (BP) values")),file=log.file,sep="\n",append=TRUE)
  .LOG("After merging across all summary statistics using listwise deletion, performing QC, and merging with the reference file, there are ",nrow(data.frame.out), " SNPs left in the final multivariate summary statistics file",file=log.file,append=TRUE)
  .LOG("Sumstats finished running at ",end.time,file=log.file,append=TRUE)
  .LOG("Running sumstats for all files took ",mins," minutes and ",secs," seconds",file=log.file,append=TRUE)
  .LOG("Please check the log file", paste0(log2, "_sumstats.log"), "to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files.",file=log.file,append=TRUE)
  flush(log.file)
  close(log.file)

  data.frame.out

}