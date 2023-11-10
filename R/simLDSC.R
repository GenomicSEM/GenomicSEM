
simLDSC <- function(covmat,N,seed,ld,rPheno=NULL,int=NULL,N_overlap=.99,r=1,gzip_output=TRUE,parallel=FALSE,cores=NULL) {
  if (!(is.character(covmat) | is.matrix(covmat))) stop("argument covmat should be character or matrix")
  if (!(is.numeric(N) | is.matrix(N))) stop("argument N should be integer or matrix")
  if (is.matrix(N)){
    if(TRUE %in% (diag(N) <= 0)) stop("argument N should be greater than 0")
  }
  if (!is.matrix(N)){
    if (is.numeric(N)) {
      if (N <= 0) stop("argument N should be greater than 0")
      if ((N_overlap < 0) | (N_overlap > 1)) stop("argument N_overlap should be between 0 and 1 (inclusive)")
    }
  }
  if (!dir.exists(ld)) stop("directory for argument ld (", ld, ") not found")
  if (!is.numeric(r)) stop("argument r should be numeric")
  if (is.numeric(r)) {if (r < 1) stop("argument r should be equal to or greater than 1")}
  if (!is.numeric(seed)) stop("argument seed should be numeric")
  if (!(gzip_output %in% c(TRUE, FALSE))) stop("argument gzip_output should be TRUE or FALSE")
  list.of.packages <- c("data.table","readr","GenomicSEM","R.utils","MASS","tictoc","lavaan","doParallel","tictoc")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) stop("Missing package(s) ", paste0(new.packages, collapse=" and "))
  lapply(list.of.packages, library,character.only = TRUE)
  ld_path<-ld
  chr <- 22
  cat("--------------------","\n"," Reading LD scores","\n","--------------------","\n", sep ="")
  x <- do.call("rbind", lapply(1:chr, function(i) {
    suppressMessages(read_delim(
      file.path(ld, paste0(i, ".l2.ldscore.gz")),
      delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
  }))
  cat("--------------------","\n","Filtering MHC region","\n","--------------------","\n", sep ="")
  x<-data.frame(x)
  if(class(covmat)[1]=="character"){
    model <- covmat
    A <- simulateData(model = model, model.type = "cfa", std.lv = F, sample.nobs = 5, empirical = T, return.fit = TRUE)
    names <- names(A)
    covMatrix <- as.data.frame(simulateData(model = model, model.type = "cfa",std.lv = F,
                                            seed = seed,return.type = "cov", empirical = T))
    correctorder <- paste0("Pheno", 1:ncol(covMatrix))
    fixOrder <- order(names, correctorder)
    covMatrix <- covMatrix[fixOrder, fixOrder]
  }else{
    covMatrix <- as.data.frame(covmat)
  }
  if (!is.matrix(N)) {
    N <- diag(N,ncol(covMatrix),ncol(covMatrix))
    N[lower.tri(N,diag = F)] <- N[1, 1]*N_overlap
    N[upper.tri(N,diag = F)] <- N[1, 1]*N_overlap
  }
  phenos <- nrow(covMatrix)
  rownames(covMatrix) <- paste("Pheno_", 1:phenos, sep = "")
  colnames(covMatrix) <- paste("Pheno_", 1:phenos, sep = "")
  rG <- cov2cor(as.matrix(covMatrix))
  if (is.null(int)){
    int <- rep(1,phenos)
  }
  if (is.null(rPheno)){
    rPheno <- rG
  }
  if (!is.matrix(rPheno)) {
    rPheno <- matrix(rPheno, ncol(covMatrix), ncol(covMatrix))
    diag(rPheno) <- 1
  }
  colnames(rPheno) <- colnames(rG)
  rownames(rPheno) <- colnames(rG)
  ld <- x$L2 
  M <- do.call("rbind", lapply(1:chr, function(i) {
    suppressMessages(read_csv(file.path(ld_path, paste0(i, ".l2.M_5_50")), col_names = FALSE))
  }))
  M <-  sum(M)
  colnames(N) <- paste("Pheno_",1:phenos,sep = "")
  rownames(N) <- colnames(N)
  mu <- rep(0,phenos)
  varZ <- list()
  for(i in 1:phenos){
    varZ[[i]] <- (diag(N)[i]*covMatrix[i,i]/M)*ld+int[i]
    names(varZ)[i] <- paste(colnames(covMatrix)[i],",",colnames(covMatrix)[i],sep = "")
  }
  covZ <- list()
  covs <- t(unique(combn(colnames(rG), 2)))
  for(i in 1:nrow(covs)){
    P1<-covs[i,1]
    P2<-covs[i,2]
    covZ[[i]] <- (sqrt(diag(N)[P1]*diag(N)[P2])*covMatrix[P1,P2]/M)*ld+rPheno[P1,P2]*N[P1,P2]/sqrt(diag(N)[P1]*diag(N)[P2])
    names(covZ)[i] <- paste0(P1,",",P2)
  }
  Z <- matrix(NA,nrow=M,ncol=phenos)
  SigmaNames <- matrix(NA,phenos,phenos)
  for (k in 1:phenos){
    for(c in 1:phenos){
      SigmaNames[k,c] <- paste(rownames(covMatrix)[k],",",colnames(covMatrix)[c],sep = "")
    }
  }
  makeSymm <- function(m) {
    m[upper.tri(m)] <- t(m)[upper.tri(m)]
    return(m)
  }
  SigmaNames <- makeSymm(t(SigmaNames))
  varcovarZ <- c(varZ,covZ)
  Sigma <- list()
  for (p in 1:length(c(SigmaNames))){
    Sigma[[p]] <-  varcovarZ[[c(SigmaNames)[[p]]]]
    names(Sigma)[p] <- c(SigmaNames)[[p]]
  }
  Sigma <- do.call(cbind,Sigma)
  Sigmalist <- lapply(asplit(Sigma,1),matrix,ncol = phenos)
  cat("--------------------","\n","Starting simulation","\n","--------------------","\n", sep ="")
  set.seed(seed)
  if(parallel == FALSE){
    for (r in 1:r){ 
      tic("Iteration")
      cat("Iteration number",r,"\n")
      cat("--------------------","\n",sep ="")
      cat("\n","Generating Z statistics","\n","\n", sep ="")
      Z <- t(mapply(function(x, y) mvrnorm(n = 1, mu = x, Sigma=y), list(mu), Sigmalist))
      sumstats <- list()
      for (i in 1:phenos){
        GWAS <- x
        GWAS$Z=Z[,i]
        GWAS$N=diag(N)[i]
        GWAS$A1="A"
        sumstats[[i]] <- GWAS
        names(sumstats)[[i]] <- rownames(rG)[i]
        cat("Writing sumstats for Phenotype",i,"\n","\n", sep ="")
        write.table(x = sumstats[[i]],file = paste0("iter",r,"GWAS",i,".sumstats"),
                    sep="\t", quote = FALSE, row.names = F)
        if (gzip_output) gzip(paste0("iter",r,"GWAS",i,".sumstats"), overwrite=TRUE)
      }
      toc()
      cat("-------------------------------------","\n",sep ="")
    }
  }
  if(parallel == TRUE){
    tic("Iteration")
    if(is.null(cores)){
      int <- detectCores() - 1
    }else{int<-cores}
    
    registerDoParallel(int)
    makeCluster(int, type="FORK")
      foreach(r=1:r) %dopar% {
      Z <- t(mapply(function(x, y) mvrnorm(n = 1, mu = x, Sigma=y), list(mu), Sigmalist))
      sumstats <- list()
      for (i in 1:phenos){
        GWAS <- x
        GWAS$Z=Z[,i]
        GWAS$N=diag(N)[i]
        GWAS$A1="A"
        sumstats[[i]] <- GWAS
        names(sumstats)[[i]] <- rownames(rG)[i]
        cat("Writing sumstats for Phenotype",i,"\n","\n", sep ="")
        write.table(x = sumstats[[i]],file = paste0("iter",r,"GWAS",i,".sumstats"),
                    sep="\t", quote = FALSE, row.names = F)
        if (gzip_output) gzip(paste0("iter",r,"GWAS",i,".sumstats"), overwrite=TRUE)
      }
    }
  }
  toc()
  cat("-------------------------------------","\n",sep ="")
  save(covMatrix,file = "PopCovMat.RData")
}
