#### GenomicSEM multivariable HDL function, based on the amazing work by Ning, Pawitan and Shen, Nature Genetics (2020)

hdl <- function(traits,sample.prev=NA,population.prev=NA,trait.names=NULL,LD.path,Nref = 335265,method="piecewise"){

  ### Do some data wrangling for the LD files:
  
  cat("GenomicSEM multivariable HDL function, based on the original implmentation of HDL please cite: Ning, Pawitan & Shen, Nature Genetics (2020)")
  
  if(is.null(trait.names)){
    traits2 <- paste0("V",1:length(traits))
    trat.names<-(traits2)
  }
  
  
  LD.files <- list.files(LD.path)
  
  if(any(grepl(x = LD.files, pattern = "UKB_snp_counter.*"))){
    snp_counter_file <- LD.files[grep(x = LD.files, pattern = "UKB_snp_counter.*")]
    snp_list_file <- LD.files[grep(x = LD.files, pattern = "UKB_snp_list.*")]
    load(file=paste(LD.path, snp_counter_file, sep = "/"))
    load(file=paste(LD.path, snp_list_file, sep = "/"))
    if("nsnps.list.imputed" %in% ls()){
      snps.name.list <- snps.list.imputed.vector
      nsnps.list <- nsnps.list.imputed
    }
  } else{
    error.message <- "It seems this directory does not contain all files needed for HDL. Please check your LD.path again. The version of HDL implementerd in GenomicSEM only support pre-computed LD reference panels."
    stop(error.message)
  }
  
  num.pieces <- length(unlist(nsnps.list))
  
  
  
  # Dimensions
  n.traits <- length(traits)
  n.V <- (n.traits^2 / 2) + .5*n.traits
  
  if(!(is.null(trait.names))){
    check_names<-str_detect(trait.names, "-")
    if(any(check_names==TRUE)){warning("Your trait names specified include mathematical arguments (e.g., + or -) that will be misread by lavaan. Please rename the traits using the trait.names argument.")}
  }
  
  if(length(traits)==1){warning("Our version of hdl requires 2 or more traits. Please include an additional trait.")}
  
  
  # Storage:
  S <- cov <- matrix(NA,nrow=n.traits,ncol=n.traits)
  V.hold <- matrix(NA,nrow=num.pieces,ncol=n.V)
  N.vec <- matrix(NA,nrow=1,ncol=n.V)
  Liab.S <- matrix(1,nrow=1,ncol=n.traits)
  I <- matrix(NA,nrow=n.traits,ncol=n.traits)
  complete <- matrix(1,nrow=1,ncol=n.traits)
# Start with general utility functions needed for hdl and standard errors


#### Define the liklihood funcrion to be optimized for h2:


llfun <-  function(param, N, M,Nref=1000, lam, bstar, lim=exp(-10)){
    h2 <- param[1]
    int <- param[2]
    lamh2 <- h2/M*lam^2 - h2*lam/Nref + int*lam/N
    lamh2 <- ifelse(lamh2<lim, lim, lamh2)
    ll <- sum(log(lamh2)) + sum(bstar^2/(lamh2))
    return(ll)
  }

#### Define the liklihood funcrion to be optimized for genetic covariance:

llfun.gcov.part.2 <- function(param, h11, h22, rho12, M, N1, N2, N0, Nref, lam0, lam1, lam2, bstar1, bstar2, lim=exp(-10)){
  h12 <- param[1]
  int <- param[2]
  ## sample fractions
  p1 <- N0/N1; p2 <- N0/N2
  ## must follow the formula for lamh2 used in llfun4
  lam11 <- h11[1]/M*lam1^2 - h11[1]*lam1/Nref + h11[2]*lam1/N1
  lam11 <- ifelse(lam11<lim, lim, lam11)
  lam22 <- h22[1]/M*lam2^2 - h22[1]*lam2/Nref + h22[2]*lam2/N2
  lam22 <- ifelse(lam22<lim, lim, lam22)
  #lam12 = h12/M*lam1*lam2 - p1*p2*h12*lam1/Nref + p1*p2*int*lam1/N0
  if (N0>0) lam12 <- h12/M*lam1*lam2 + p1*p2*int*lam1/N0  ## key change here
  if (N0==0) lam12 <- h12/M*lam1*lam2
  ##  resid of bstar2 ~bstar1
  ustar <- bstar2 - lam12/lam11*bstar1  ## see note
  lam22.1 <- lam22 - lam12^2/lam11
  lam22.1 <- ifelse(lam22.1<lim, lim, lam22.1)
  ll <- sum(log(lam22.1)) + sum(ustar^2/(lam22.1))
  return(ll)
}


################# Begin the loop to compute S and V cell by cell using HDL

s<-1 # coubnt for elements in V
for(j in 1:n.traits){
  for(d in j:n.traits){
cat("\n")
  if(j == d){
    chi1 <- traits[j]
    gwas.df <- as.data.frame(suppressMessages(read_delim(chi1, "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))
    
    gwas.df <- gwas.df %>% filter(SNP %in% snps.name.list)
    gwas.df$A1 <- as.character(gwas.df$A1)
    gwas.df$A2 <- as.character(gwas.df$A2)
    if (!("Z" %in% colnames(gwas.df))) {
      if (("b" %in% colnames(gwas.df)) && ("se" %in% colnames(gwas.df))) {
        gwas.df$Z <- gwas.df$b/gwas.df$se
      }
      else {
        error.message <- "Z is not available, meanwhile either b or se is missing. Please check."

        stop(error.message)
      }
    }
    gwas.df <- gwas.df %>% filter(!is.na(Z))
    k1 <- nrow(gwas.df)
    k1.percent <- paste("(", round(100 * k1/length(snps.name.list), 
                                   2), "%)", sep = "")
    cat(k1, "out of", length(snps.name.list), k1.percent, "SNPs in reference panel are available in the GWAS of",trait.names[j] , 
        " \n")
  
    if (k1 < length(snps.name.list) * 0.99) {
      error.message <- "Warning: More than 1% SNPs in reference panel are missed in the GWAS. This may generate bias in estimation. Please make sure that you are using correct reference panel.  \n"
      cat(error.message)
    }
    complete[j] <- k1.percent
    
    N1 <- median(gwas.df[, "N"])
    N <- N1
    bstar1.v <- lam.v <- list()
    HDL11.df <- names.row <- NULL
    logL.df <- NULL
    counter <- 0
    message <- ""
    num.pieces <- length(unlist(nsnps.list))
    for (chr in 1:22) {
      k <- length(nsnps.list[[chr]])
      for (piece in 1:k) {
        LD_rda_file <- LD.files[grep(x = LD.files, pattern = paste0("chr", chr, ".", piece, ".*rda"))]
        LD_bim_file <- LD.files[grep(x = LD.files, pattern = paste0("chr", chr, ".", piece, ".*bim"))]
        load(file = paste(LD.path, LD_rda_file, sep = "/"))
        snps.ref.df <- read.table(paste(LD.path, LD_bim_file, 
                                        sep = "/"))
        colnames(snps.ref.df) <- c("chr", "id", "non", "pos", "A1", "A2")
        snps.ref <- snps.ref.df$id
        A2.ref <- snps.ref.df$A2
        names(A2.ref) <- snps.ref
        gwas.df.subset <- gwas.df %>% filter(SNP %in% snps.ref)
        bhat1.raw <- gwas.df.subset[, "Z"]/sqrt(gwas.df.subset[, "N"])
        A2.gwas1 <- gwas.df.subset[, "A2"]
        names(bhat1.raw) <- names(A2.gwas1) <- gwas.df.subset$SNP
        idx.sign1 <- A2.gwas1 == A2.ref[names(A2.gwas1)]
        bhat1.raw <- bhat1.raw * (2 * as.numeric(idx.sign1) -  1)
        M <- length(LDsc)
        bhat1 <- numeric(M)
        names(bhat1) <- snps.ref
        bhat1[names(bhat1.raw)] <- bhat1.raw
        a11 <- bhat1^2
        reg <- lm(a11 ~ LDsc)
        h11.ols <- c(summary(reg)$coef[1:2, 1:2] * c(N1,  M))
        h11v <- (h11.ols[2] * LDsc/M + 1/N1)^2
        reg <- lm(a11 ~ LDsc, weight = 1/h11v)
        h11.wls <- c(summary(reg)$coef[1:2, 1:2] * c(N1,  M))
        bstar1 <- crossprod(V, bhat1)
        opt <- optim(c(h11.wls[2], 1), llfun, N = N1, Nref = Nref,
                     lam = lam, bstar = bstar1, M = M, lim = exp(-18),
                     method = "L-BFGS-B", lower = c(0, 0), upper = c(1,10))
        h11.hdl <- opt$par
        logL <- opt$value
        
        HDL11.df <- rbind(HDL11.df, h11.hdl)
        bstar1.v <- c(bstar1.v, list(bstar1))
        lam.v <- c(lam.v, list(lam))
        counter <- counter + 1
        value <- round(counter/num.pieces * 100)
        backspaces <- paste(rep("\b", nchar(message)), collapse = "")
        message <- paste("Estimation for cell: ",s," out of: ",n.V," cells is ongoing ... ", value,  "%", sep = "", collapse = "")
        cat(backspaces, message, sep = "")
      }
    }
   
    
     # This runs the optimisation at once for the entire genome, and computes the V via jackknive. Its HDL defaiult behaviour, but not GenomicSEM defaiult behaviour.
    
    if(method=="jackknife"){ 
      
      M.ref <- sum(unlist(nsnps.list))
      
      opt <- optim(c(sum(HDL11.df[, 1]), 1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v), bstar=unlist(bstar1.v), M=M.ref,
                   lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
      h11.hdl <- opt$par
     
      cat("Continuing computing standard error with jackknife \n")

      counter <- 0
      message <- ""
       h11.jackknife <- numeric(length(lam.v))
      for(i in 1:length(lam.v)){
        opt <- optim(h11.hdl, llfun, N=N1, Nref=Nref, lam=unlist(lam.v[-i]), bstar=unlist(bstar1.v[-i]), M=M.ref,
                     lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
        h11.hdl.jackknife <- opt$par
        

      
          V.hold[i,s] <- h11.hdl.jackknife[1]
        
        ## Report progress ##
        
        counter <- counter + 1
        value <- round(counter/length(lam.v)*100)
        backspaces <- paste(rep("\b", nchar(message)), collapse = "")
        message <- paste("Progress... ", value, "%", sep = "", 
                         collapse = "")
        cat(backspaces, message, sep = "")
      }
      
      
      I[j,d] <- I[d,j] <- h11.hdl[2]
      S[j,d] <- S[d,j] <- h11.hdl[1]
      
    }

    if(method=="piecewise"){
    
    rownames(HDL11.df) <- names.row
    h1_2 <- sum(HDL11.df[, 1])
  
    
    S[j,j] <- h1_2
    I[j,j] <- mean(HDL11.df[,2])

    for(i in 1:num.pieces){
    V.hold[i,s] <- (num.pieces/(num.pieces-1))*sum(HDL11.df[-i,1])
    }
    }
    
    s <- s+1
  }  
    
  if (j != d){  
  chi1 <- traits[j]
  chi2 <- traits[d]
  gwas1.df <- as.data.frame(suppressMessages(read_delim(chi1, "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))
  gwas2.df <- as.data.frame(suppressMessages(read_delim(chi2, "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))

  gwas1.df <- gwas1.df %>% filter(SNP %in% snps.name.list)
  gwas2.df <- gwas2.df %>% filter(SNP %in% snps.name.list)
  gwas1.df$A1 <- as.character(gwas1.df$A1)
  gwas1.df$A2 <- as.character(gwas1.df$A2)
  gwas2.df$A1 <- as.character(gwas2.df$A1)
  gwas2.df$A2 <- as.character(gwas2.df$A2)
  
  
  N0 <- min(gwas1.df$N, gwas2.df$N)
  
  if (!("Z" %in% colnames(gwas1.df))) {
    if (("b" %in% colnames(gwas1.df)) && ("se" %in% colnames(gwas1.df))) {
      gwas1.df$Z <- gwas1.df$b/gwas1.df$se
    }
    else {
      error.message <- "Z is not available, meanwhile either b or se is missing. Please check."
      stop(error.message)
    }
  }
  if (!("Z" %in% colnames(gwas2.df))) {
    if (("b" %in% colnames(gwas2.df)) && ("se" %in% colnames(gwas2.df))) {
      gwas2.df$Z <- gwas2.df$b/gwas2.df$se
    }
    else {
      error.message <- "Z is not available, meanwhile either b or se is missing. Please check."
      stop(error.message)
    }
  }
  gwas1.df <- gwas1.df %>% filter(!is.na(Z))
  gwas2.df <- gwas2.df %>% filter(!is.na(Z))

  N1 <- median(gwas1.df[, "N"])
  N2 <- median(gwas2.df[, "N"])
  N <- sqrt(N1) * sqrt(N2)
  p1 <- N0/N1
  p2 <- N0/N2
  rho12 <- suppressWarnings(inner_join(gwas1.df , gwas2.df , by = "SNP") %>% summarise(x = cor(Z.x, Z.y, use = "complete.obs")) %>% unlist)
  bstar1.v <- bstar2.v <- lam.v <- list()
  HDL11.df <- HDL12.df <- HDL22.df <- names.row <- NULL
  counter <- 0
  message <- ""
  num.pieces <- length(unlist(nsnps.list))
  
  for (chr in 1:22) {
    k <- length(nsnps.list[[chr]])
    for (piece in 1:k) {
      LD_rda_file <- LD.files[grep(x = LD.files, pattern = paste0("chr", chr, ".", piece, ".*rda"))]
      LD_bim_file <- LD.files[grep(x = LD.files, pattern = paste0("chr", chr, ".", piece, ".*bim"))]
      load(file = paste(LD.path, LD_rda_file, sep = "/"))
      snps.ref.df <- read.table(paste(LD.path, LD_bim_file, sep = "/"))
      colnames(snps.ref.df) <- c("chr", "id", "non", "pos", "A1", "A2")
      snps.ref <- snps.ref.df$id
      A2.ref <- snps.ref.df$A2
      names(A2.ref) <- snps.ref
      gwas1.df.subset <- gwas1.df %>% filter(SNP %in%  snps.ref)
      bhat1.raw <- gwas1.df.subset[, "Z"]/sqrt(gwas1.df.subset[, "N"])
      A2.gwas1 <- gwas1.df.subset[, "A2"]
      names(bhat1.raw) <- names(A2.gwas1) <- gwas1.df.subset$SNP
      idx.sign1 <- A2.gwas1 == A2.ref[names(A2.gwas1)]
      bhat1.raw <- bhat1.raw * (2 * as.numeric(idx.sign1) -  1)
      gwas2.df.subset <- gwas2.df %>% filter(SNP %in% 
                                               snps.ref)
      bhat2.raw <- gwas2.df.subset[, "Z"]/sqrt(gwas2.df.subset[, "N"])
      A2.gwas2 <- gwas2.df.subset[, "A2"]
      names(bhat2.raw) <- names(A2.gwas2) <- gwas2.df.subset$SNP
      idx.sign2 <- A2.gwas2 == A2.ref[names(A2.gwas2)]
      bhat2.raw <- bhat2.raw * (2 * as.numeric(idx.sign2) -  1)
      M <- length(LDsc)
      bhat1 <- bhat2 <- numeric(M)
      names(bhat1) <- names(bhat2) <- snps.ref
      bhat1[names(bhat1.raw)] <- bhat1.raw
      bhat2[names(bhat2.raw)] <- bhat2.raw
      a11 <- bhat1^2
      a22 <- bhat2^2
      a12 <- bhat1 * bhat2
      reg <- lm(a11 ~ LDsc)
      h11.ols <- c(summary(reg)$coef[1:2, 1:2] * c(N1,M))
      reg <- lm(a22 ~ LDsc)
      h22.ols <- c(summary(reg)$coef[1:2, 1:2] * c(N2,   M))
      reg <- lm(a12 ~ LDsc)
      if (N0 > 0) 
        h12.ols <- c(summary(reg)$coef[1:2, 1:2] * c((N0/p1/p2), M))
      if (N0 == 0) 
        h12.ols <- c(summary(reg)$coef[1:2, 1:2] * c(N, M))
      h11v <- (h11.ols[2] * LDsc/M + 1/N1)^2
      h22v <- (h22.ols[2] * LDsc/M + 1/N2)^2
      reg <- lm(a11 ~ LDsc, weight = 1/h11v)
      h11.wls <- c(summary(reg)$coef[1:2, 1:2] * c(N1, M))
      reg <- lm(a22 ~ LDsc, weight = 1/h22v)
      h22.wls <- c(summary(reg)$coef[1:2, 1:2] * c(N2,M))
      if (N0 > 0) 
        h12v <- sqrt(h11v * h22v) + (h12.ols[2] * LDsc/M + p1 * p2 * rho12/N0)^2
      if (N0 == 0) 
        h12v <- sqrt(h11v * h22v) + (h12.ols[2] * LDsc/M)^2
      reg <- lm(a12 ~ LDsc, weight = 1/h12v)
      if (N0 > 0) 
        h12.wls <- c(summary(reg)$coef[1:2, 1:2] * c((N0/p1/p2), M))
      if (N0 == 0) 
        h12.wls <- c(summary(reg)$coef[1:2, 1:2] * c(N, M))
      bstar1 <- crossprod(V, bhat1)
      bstar2 <- crossprod(V, bhat2)
      opt <- optim(c(h11.wls[2], 1), llfun, N = N1, Nref = Nref,
                   lam = lam, bstar = bstar1, M = M, lim = exp(-18),
                   method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 10))
      h11.hdl <- opt$par
      opt <- optim(c(h22.wls[2], 1), llfun, N = N2, Nref = Nref,
                   lam = lam, bstar = bstar2, M = M, lim = exp(-18),
                   method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 10))
      h22.hdl <- opt$par
      opt <- optim(c(h12.wls[2], rho12), llfun.gcov.part.2,
                   h11 = h11.hdl, h22 = h22.hdl, rho12 = rho12,
                   M = M, N1 = N1, N2 = N2, N0 = N0, Nref = Nref,
                   lam0 = lam, lam1 = lam, lam2 = lam, bstar1 = bstar1,
                   bstar2 = bstar2, lim = exp(-18), method = "L-BFGS-B",
                   lower = c(-1, -10), upper = c(1, 10))
      h12.hdl <- opt$par

      HDL11.df <- rbind(HDL11.df, h11.hdl)
      HDL22.df <- rbind(HDL22.df, h22.hdl)
      HDL12.df <- rbind(HDL12.df, h12.hdl)
      bstar1.v <- c(bstar1.v, list(bstar1))
      bstar2.v <- c(bstar2.v, list(bstar2))
      lam.v <- c(lam.v, list(lam))
      counter <- counter + 1
      value <- round(counter/num.pieces * 100)
      backspaces <- paste(rep("\b", nchar(message)), collapse = "")
      message <- paste("Estimation for cell: ",s," out of: ",n.V," cells is ongoing ... ", value, 
                       "%", sep = "", collapse = "")
      cat(backspaces, message, sep = "")
    }
    
  }
  
  if(method=="jackknife"){ # This runs the optimisation at once for the entire genome, and computes the V via jackknive. Its HDL defaiult behaviour, but not GenomicSEM defaiult behaviour.
    
    
  M.ref <- sum(unlist(nsnps.list))
    
  opt <- optim(c(sum(HDL11.df[, 1]), 1), llfun, N=N1, Nref=Nref, lam=unlist(lam.v), bstar=unlist(bstar1.v), M=M.ref,
               lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
  h11.hdl <- opt$par
  opt <- optim(c(sum(HDL22.df[, 1]), 1), llfun, N=N2, Nref=Nref, lam=unlist(lam.v), bstar=unlist(bstar2.v), M=M.ref,
               lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
  h22.hdl <- opt$par
  
  
 
  opt <- optim(c(sum(HDL12.df[, 1]), rho12), llfun.gcov.part.2, h11=h11.hdl, h22=h22.hdl,
               rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref,
               lam0=unlist(lam.v), lam1=unlist(lam.v), lam2=unlist(lam.v),
               bstar1=unlist(bstar1.v), bstar2=unlist(bstar2.v),
               lim=exp(-18), method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
  if(opt$convergence != 0){
    starting.value.v <- c(0,-sqrt(h11*h22)*0.5, sqrt(h11*h22)*0.5)
    k <- 1
    while(opt$convergence != 0){
      starting.value <- starting.value.v[k]
      opt <- optim(c(starting.value, rho12), llfun.gcov.part.2, h11=h11.hdl, h22=h22.hdl,
                   rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref,
                   lam0=unlist(lam.v), lam1=unlist(lam.v), lam2=unlist(lam.v),
                   bstar1=unlist(bstar1.v), bstar2=unlist(bstar2.v),
                   lim=exp(-18), method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
      k <- k + 1
      if(k > length(starting.value.v)){
        error.message <- "Algorithm failed to converge after trying different initial values. \n"
        stop(error.message)
      }
    }}
  h12.hdl <- opt$par
  
  
  cat("Continuing computing standard error with jackknife \n")
  counter <- 0
  message <- ""
  rg.jackknife <- h11.jackknife <- h12.jackknife <- h22.jackknife <- numeric(length(lam.v))
  for(i in 1:length(lam.v)){
    opt <- optim(h11.hdl, llfun, N=N1, Nref=Nref, lam=unlist(lam.v[-i]), bstar=unlist(bstar1.v[-i]), M=M.ref,
                 lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
    h11.hdl.jackknife <- opt$par
    
    
    opt <- optim(h22.hdl, llfun, N=N2, Nref=Nref, lam=unlist(lam.v[-i]), bstar=unlist(bstar2.v[-i]), M=M.ref,
                 lim=exp(-18), method ='L-BFGS-B', lower=c(0,0), upper=c(1,10))
    h22.hdl.jackknife <- opt$par
    
    opt <- optim(h12.hdl, llfun.gcov.part.2, h11=h11.hdl, h22=h22.hdl,
                 rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref,
                 lam0=unlist(lam.v[-i]), lam1=unlist(lam.v[-i]), lam2=unlist(lam.v[-i]),
                 bstar1=unlist(bstar1.v[-i]), bstar2=unlist(bstar2.v[-i]),
                 lim=exp(-18), method ='L-BFGS-B', lower=c(-1,-10), upper=c(1,10))
    h12.hdl.jackknife <- opt$par
    
    
    
    V.hold[i,s] <- h12.hdl.jackknife[1]
    
    
    ## Report progress ##
    
    counter <- counter + 1
    value <- round(counter/length(lam.v)*100)
    backspaces <- paste(rep("\b", nchar(message)), collapse = "")
    message <- paste("Progress... ", value, "%", sep = "", 
                     collapse = "")
    cat(backspaces, message, sep = "")
  }
  
 
  I[j,d] <- I[d,j] <- h12.hdl[2]
  S[j,d] <- S[d,j] <- h12.hdl[1]
  
  }
  
  if(method=="piecewise"){
  
  rownames(HDL11.df) <- rownames(HDL22.df) <- rownames(HDL12.df) <- names.row
  h1_2 <- sum(HDL11.df[, 1])
  h2_2 <- sum(HDL22.df[, 1])
  gen.cov <- sum(HDL12.df[, 1])

    
  I[j,d] <- I[d,j] <- mean(HDL12.df[,2])
  S[j,d] <- S[d,j] <- gen.cov
  for(i in 1:num.pieces){
    V.hold[i,s] <- (num.pieces/(num.pieces-1))*sum(HDL12.df[-i,1])
  }
  
  }
  
  s <- s+ 1
  }
    
  }
}

Liab.S <- matrix(1,nrow=1,ncol=n.traits)

for(z in 1:n.traits){
  pop.prev <- population.prev[z]
  samp.prev <- sample.prev[z]
  
if(is.na(pop.prev)==F & is.na(samp.prev)==F){
  conversion.factor <- (pop.prev^2*(1-pop.prev)^2)/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
  Liab.S[,z] <- conversion.factor
}}
  



V <- cov(V.hold)*(num.pieces-1)  


S2 <- S


### Scale S and V to liability:
S <- diag(as.vector(sqrt(Liab.S))) %*% S %*% diag(as.vector(sqrt(Liab.S)))

#calculate the ratio of the rescaled and original S matrices
scaleO <- as.vector(lowerTriangle((S/S2), diag=T))

#obtain diagonals of the original V matrix and take their sqrt to get SE's
Dvcov<-sqrt(diag(V))

#rescale the SEs by the same multiples that the S matrix was rescaled by
Dvcovl<-as.vector(Dvcov*t(scaleO))

#obtain the sampling correlation matrix by standardizing the original V matrix
vcor<-cov2cor(V)

#rescale the sampling correlation matrix by the appropriate diagonals
V<-diag(Dvcovl)%*%vcor%*%diag(Dvcovl)

colnames(S) <- trait.names  

return(list(V = V,S = S,I = I,complete=complete))
}

