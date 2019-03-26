par_ldsc <- function(ld,traits,wld,frq,n.blocks=200,population.prev=NULL,sample.prev=NULL,out,trait.names=NULL){
 
  
 # Functions doesnt need these includes right?
 # require(plyr)
 # require(e1071)
 # require(data.table)
 # require(readr)
 # require(gdata)
  
  ##create log file
  log2<-paste(trait.names,collapse="_")
  log.file <- file(paste0(log2, "_Partitioned.log"),open="wt")
  sink(log.file,split=T,type="output")
  error.file <- file(paste0(out,".error"),open="wt")
  sink(error.file,type="message")
  begin.time <- Sys.time()
  
  ##print start time
  cat(paste("Analysis started at",begin.time),"\n")

  if(!is.null(traits)){
    cat("The following traits were analyzed:",trait.names,"\n")
  }
  
  cat("The following annotations were added to the model: ")
  cat(ld,sep=",")
  
  ld2 <- NULL
  if(length(grep(pattern="baseline$",x=ld,perl=T))==1){
    ld1 <- ld[grep(pattern="baseline$",x=ld,perl=T)]
    if(length(ld)>1){
      ld2 <- ld[grep(pattern="baseline$",x=ld,perl=T,invert=T)]
    }
  }else{
    ld1 <- ld[1]
    if(length(ld)>1){
      ld2 <- ld[2:length(ld)]
    }
  }
  
  cat("\n","Reading in LD scores from ",paste0(ld1,".[1-22]"),"\n",sep="")
  
  ##determine name of LD score files
  x.files <- sort(Sys.glob(paste0(ld1,"*l2.ldscore*")))
  
  ##function to read in the files. used with ldply below
  readLdFunc <- function(LD.in){
    if(substr(x=LD.in,start=nchar(LD.in)-1,stop=nchar(LD.in))=="gz"){
      dum=fread(input=paste("gzcat",LD.in),header=T,showProgress=F,data.table=F)
    }else{
      dum=fread(input=LD.in,header=T,showProgress=F,data.table=F)
    }
  }
  ##^changed this to gzcat for non-Linux users [AG]
  
  x <- ldply(.data=x.files,.fun=readLdFunc)
  x$CM <- NULL
  x$MAF <- NULL
  
  ##read in the M_5_50 files (number of SNPs in annotation by chromosome)
  m.files <- sort(Sys.glob(paste0(ld1,"*M_5_50")))
  readMFunc <- function(x){dum=read.table(file=x,header=F)}
  m <- ldply(.data=m.files,.fun=readMFunc)
  
  ##read in additional annotations on top of baseline if relevant
  if(!is.null(ld2)){
    for(i in 1:length(ld2)){
      cat("Reading in LD scores from ",paste0(ld2[i],".[1-22]"),"\n",sep="")
      extra.x.files <- sort(Sys.glob(paste0(ld2[i],"*l2.ldscore*")))
      extra.ldscore <- ldply(.data=extra.x.files,.fun=readLdFunc)
      extra.ldscore$CHR <- NULL
      extra.ldscore$BP <- NULL
      if(ncol(extra.ldscore)==2){colnames(extra.ldscore)[2] <- c(ld2[i])}
      extra.m.files <- sort(Sys.glob(paste0(ld2[i],"*M_5_50")))
      extra.m <- ldply(.data=extra.m.files,.fun=readMFunc)
      if((("primary" %in% ld2 && "bloodInBrain" %in% ld2)||("brain_all" %in% ld2 && "brainInBlood" %in% ld2))&&(ld2[i]=="brainInBlood" || ld2[i]=="bloodInBrain")){
        extra.ldscore <- extra.ldscore[,1:2]
        extra.m <- extra.m[,1]
      }
      if(identical(as.character(x$SNP),as.character(extra.ldscore$SNP))==T){
        colnames.x <- colnames(x)
        colnames.extra.ldscore <- colnames(extra.ldscore)[2:ncol(extra.ldscore)]
        x <- cbind(x,extra.ldscore[,2:ncol(extra.ldscore)])
        colnames(x) <- c(colnames.x,colnames.extra.ldscore)
      }else{
        x <- merge(x,extra.ldscore,by="SNP")
      }
      m <- cbind(m,extra.m)
      rm(list=ls(pattern="extra."))
      gc()
    }
  }
  
  if(ncol(m)!=(ncol(x)-3)){
    cat("ERROR:number of annotations does not match between LD scores and .M_5_50 files")
    sink()
    stop()
  }
  
  n.annot <- ncol(m)
  if(n.annot < 1){
    cat("ERROR:The files do not contain any annotations")
    sink()
    stop()
  }else if(n.annot==1){
    cat("ERROR:The function cannot handle single annotation LDSR yet")
    sink()
    stop()
  }
  
  colnames(m) <- tail(colnames(x),n=n.annot)
  m <- as.matrix(colSums(m))
  M.tot <- sum(m)
  
  cat("LD scores contain",nrow(x),"SNPs and",n.annot,"annotations","\n")
  
  cat("Reading in weighted LD scores from",paste0(wld,".[1-22]"),"\n")
  
  w.files <- sort(Sys.glob(paste0(wld,"*l2.ldscor*")))
  w <- ldply(w.files,readLdFunc)
  w$CM <- NULL
  w$MAF <- NULL
  colnames(w)[ncol(w)] <- "wLD"
  
  cat("Weighted LD scores contain ",nrow(w)," SNPs","\n","Reading in summary statistics from ",traits,"\n",sep=" ")
  
  ##ADDED FOR COVARIANCE [AG]
  n.traits <- length(traits)
  n.V <- (n.traits^2 / 2) + .5*n.traits
  
  ##create lists of S and V with length = number of annotations [AG]
  ##in baseline model, this includes annotation with all SNPs for creating overall S and V 
  S_List<-vector(mode="list",length=n.annot)
  Tau_List<-vector(mode="list",length=n.annot)
  
  for (i in 1:n.annot) { 
    S_List[[i]] = matrix(NA,nrow=n.traits,ncol=n.traits)
    Tau_List[[i]] = matrix(NA,nrow=n.traits,ncol=n.traits)
  } 

  total_pseudo<-matrix(NA,nrow=200,ncol=n.V*n.annot)
  total_pseudo_tau<-matrix(NA,nrow=200,ncol=n.V*n.annot)
  
  # Storage for N and Intercept matrix [added]:
  N.vec <- matrix(NA,nrow=1,ncol=n.V)
  I <- matrix(NA,nrow=n.traits,ncol=n.traits)
  
  ##matrix of 1s so that if it is non-liability it will just matrix multiply by 
  ##1s during liability conversion [AG]
  Liab.S <- matrix(1,nrow=1,ncol=n.traits)
  
  #moved up so only happens once; my understanding is this is general for all traits [AG]
  annot.files <- sort(Sys.glob(paste0(ld1,"*annot.gz")))
  frq.files <- sort(Sys.glob(paste0(frq,"*.frq")))
  
  M.tot.annot <- 0
  first.annot.matrix <- matrix(data=NA,nrow=(n.annot*length(annot.files)),ncol=n.annot)
  replace.annot.from <- seq(from=1,to=nrow(first.annot.matrix),by=n.annot)
  replace.annot.to <- seq(from=n.annot,to=nrow(first.annot.matrix),by=n.annot)
  
  for(i in 1:length(annot.files)){
    header <- read.table(annot.files[i], header = TRUE, nrow = 1)
    annot <- fread(input=paste("gzcat",annot.files[i]),skip=1,header=FALSE,showProgress=F,data.table=F)
    setnames(annot, colnames(header))
    annot$CM <- NULL
    if(is.null(ld2)==F){
      for(j in 1:length(ld2)){
        extra.annot.files <- sort(Sys.glob(paste0(ld2[j],"*annot.gz")))
        extra.annot <- fread(input=paste("gzcat",extra.annot.files[i]),header=T,showProgress=F,data.table=F)
        extra.annot$CHR <- NULL
        extra.annot$BP <- NULL
        extra.annot$CM <- NULL
        if(ncol(extra.annot)==2){colnames(extra.annot)[2] <- ld2[j]}
        if((("primary" %in% ld && "bloodInBrain" %in% ld)||("brain_all" %in% ld && "brainInBlood" %in% ld))&&(ld2[j]=="brainInBlood" || ld2[j]=="bloodInBrain")){
          extra.annot <- extra.annot[,1:2]
        }
        if(identical(as.character(annot$SNP),as.character(extra.annot$SNP))==T){
          colnames.annot <- colnames(annot)
          colnames.extra.annot <- colnames(extra.annot)[2:ncol(extra.annot)]
          annot <- cbind(annot,extra.annot[,2:ncol(extra.annot)])
          colnames(annot) <- c(colnames.annot,colnames.extra.annot)
        }else{
          annot <- merge(x=annot,y=extra.annot,by="SNP")
        }
        rm(extra.annot)
        gc()
      }
    }
    
    frq <- fread(input=frq.files[i],header=T,showProgress=F,data.table=F)
    
    frq <- frq[frq$MAF > 0.05 & frq$MAF < 0.95,]
    selected.annot <- merge(x=frq[,c("SNP","MAF")],y=annot,by="SNP")
    
    rm(frq,annot)
    gc()
    M.tot.annot <- M.tot.annot+nrow(selected.annot)
    first.annot.matrix[replace.annot.from[i]:replace.annot.to[i],] <- t(as.matrix(selected.annot[,5:ncol(selected.annot)]))%*% as.matrix(selected.annot[,5:ncol(selected.annot)])
  }
  
  annot.matrix <- matrix(data=NA,nrow=n.annot,ncol=n.annot)
  for(i in 1:n.annot){annot.matrix[i,] <- t(colSums(first.annot.matrix[seq(from=i,to=nrow(first.annot.matrix),by=n.annot),]))}
  rm(first.annot.matrix)
  gc()
  
  overlap.matrix <- matrix(data=NA,nrow=n.annot,ncol=n.annot)
  colnames(overlap.matrix) <- rownames(overlap.matrix) <- rownames(m)
  for(i in 1:n.annot){overlap.matrix[i,] <- annot.matrix[i,]/m}
  #########
  
  ##added##
  j<-1
  k<-1
  # count the total nummer of runs, both loops
  s <- 1
  
  for(j in 1:n.traits){
    
    chi1 <- traits[j]
    
    if(substr(x=chi1,start=nchar(chi1)-1,stop=nchar(chi1))=="gz"){
      y1 <- fread(paste("gzcat",chi1),header=T,showProgress=F,data.table=F)
    }else{
      y1 <- fread(chi1,header=T,showProgress=F,data.table=F)
    }
    y1 <- na.omit(y1)
    y1$chi <- y1$Z^2
    
    cat("Read in summary statistics for",nrow(y1),"SNPs","\n")
    
    for(k in j:n.traits){
      
      ##HERITABILITY
      if(j == k){
        samp.prev <- sample.prev[j]
        pop.prev <- population.prev[j]
        
        cat("Merging the files","\n")
        
        merged <- merge(x=y1[,c("SNP","chi","N")],y=w[,c("SNP","wLD")],by="SNP")
        merged <- merge(x=merged,y=x,by="SNP")
        merged <- merged[with(merged,order(CHR,BP)),]
        remaining.snps <- nrow(merged)
        
        cat(remaining.snps,"SNPs remaining","after merging the files","\n")
        
        chisq.max <- max(0.001*max(merged$N),80)
        merged <- merged[merged$chi < chisq.max,]
        n.snps <- nrow(merged)
        removed.snps <- remaining.snps-n.snps
        gc()
        
        cat("Removed",removed.snps,"SNPs with Chi^2 >",chisq.max,paste0("(",n.snps," SNPs remain)"),"\n")
        
        merged$intercept <- 1
        merged$x.tot <- rowSums(merged[,c(rownames(m))])
        merged$x.tot.intercept <- 1
        
        tot.agg <- (M.tot*(mean(merged$chi)-1))/mean(merged$x.tot*merged$N)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- as.numeric(lapply(X=merged$x.tot,function(x){max(x,1)}))
        merged$w.ld <- as.numeric(lapply(X=merged$wLD,function(x){max(x,1)}))
        merged$c <- tot.agg*merged$N/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)
        merged$weights <- merged$initial.w/sum(merged$initial.w)
        
        N.bar <- mean(merged$N)
        merged[,(ncol(merged)-10-n.annot):(ncol(merged)-11)] <- (merged[,(ncol(merged)-10-n.annot):(ncol(merged)-11)]*merged$N)/N.bar
        
        LD.scores <- as.matrix(merged[,(ncol(merged)-10-n.annot):(ncol(merged)-10)])
        weighted.LD <- as.matrix(LD.scores*merged$weights)
        weighted.chi <- as.matrix(merged$chi*merged$weights)
        
        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:200]-1,n.snps)
        
        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow =((n.annot+1)* n.blocks),ncol =(n.annot+1))
        colnames(xty.block.values) <- colnames(xtx.block.values) <- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by =(n.annot+1))
        replace.to <- seq(from =(n.annot+1),to=nrow(xtx.block.values),by =(n.annot+1))
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
        colnames(xtx) <- colnames(weighted.LD)
        for(i in 1:nrow(xtx)){xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))}
        
        reg <- solve(xtx)%*% xty
        intercept <- reg[length(reg)]
        coefs <- reg[1:(length(reg)-1)]/N.bar
        cats <- coefs*m
        reg.tot <- sum(cats)
        est <- as.matrix(cats/reg.tot)

        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values) <- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete)%*% xty.delete
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        
        pseudo.values_tau <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values_tau) <- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values_tau[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])} 
       
        jackknife.cov <- cov(pseudo.values_tau)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        
        Coefficient.std.error=data.frame(sqrt(diag(coef.cov)))
        Coefficient.Z<-data.frame(coefs/sqrt(diag(coef.cov)))
  
        cat.cov <- coef.cov*(m %*% t(m))
        #Tau.std.error=data.frame(sqrt(diag(cat.cov)))
        ##^^Added 3.24.19; think this is Tau SE, but wouldn't provide off-diagonal elements in current form
        
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
       
        ##for overall heritability
        if(is.null(c(pop.prev,samp.prev))){
          cat("h2:",round(reg.tot,4),"(",round(tot.se,4),")","\n")
        }else if(!is.null(c(pop.prev,samp.prev))){
          conversion.factor <- pop.prev^2*(1-pop.prev)^2/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
          h2.liab <- reg.tot*conversion.factor
          h2.liab.se <- tot.se*conversion.factor
          cat("Liability h2:",round(h2.liab,4),"(",round(h2.liab.se,4),")","\n")
        }else{
          cat("ERROR:one of the prevalence values was not supplied")
          sink()
          stop()
        }
        
        lambda.gc <- median(merged$chi)/0.4549
        mean.Chi <- mean(merged$chi)
        ratio <- (intercept-1)/(mean(merged$chi)-1)
        ratio.se <- intercept.se/(mean(merged$chi)-1)
        
        cat("Lambda GC:",round(lambda.gc,4),"\n")
        cat("Mean Chi^2:",round(mean.Chi,4),"\n")
        cat("Intercept: ",round(intercept,4),"(",round(intercept.se,4),")","\n",sep="")
        cat("Ratio: ",round(ratio,4),"(",round(ratio.se,4),")","\n",sep="")
        
        ####NEW HILL CODE######
        ########## begin part with 2 surroundings blocks removed ##########
        #delete.values.K3 <- matrix(data=NA,nrow=n.blocks,ncol=ncol(weighted.LD))
        #colnames(delete.values.K3)<- colnames(weighted.LD)
        #first block
        #xty.delete.K3 <- xty-xty.block.values[1,]-xty.block.values[2,]
        #xtx.delete.K3 <- xtx-xtx.block.values[delete.from[1]:delete.to[1],]-xtx.block.values[delete.from[2]:delete.to[2],]
        #delete.values.K3[1,] <- solve(xtx.delete.K3)%*% xty.delete.K3
        #last block
        #xty.delete.K3 <- xty-xty.block.values[(n.blocks-1),]-xty.block.values[n.blocks,]
        #xtx.delete.K3 <- xtx-xtx.block.values[delete.from[(n.blocks-1)]:delete.to[(n.blocks-1)],]-xtx.block.values[delete.from[n.blocks]:delete.to[n.blocks],]
        #delete.values.K3[n.blocks,] <- solve(xtx.delete.K3)%*% xty.delete.K3
        #remaining blocks
        #for(i in 2:(n.blocks-1)){
        #  xty.delete.K3 <- xty-xty.block.values[(i-1),]-xty.block.values[i,]-xty.block.values[(i+1),]
        #  xtx.delete.K3 <- xtx-xtx.block.values[delete.from[(i-1)]:delete.to[(i-1)],]-xtx.block.values[delete.from[i]:delete.to[i],]-xtx.block.values[delete.from[(i+1)]:delete.to[(i+1)],]
        #  delete.values.K3[i,] <- solve(xtx.delete.K3)%*% xty.delete.K3
        #}
        ########## end part with 2 surrounding blocks removed ##########
        
        ########## begin part about cross-validation ##########
        #fitted.values <- merged[,c("CHR","BP","SNP")]
        #fitted.values$Z <- y[y$SNP %in% fitted.values$SNP,]$Z
        #fitted.values$observed.chi <- merged$chi
        #fitted.values$weighted.chi <- weighted.chi
        #fitted.values$predicted.weighted.chi <- NA
        #fitted.values$predicted.weighted.chi.K <- NA
        #fitted.values$predicted.weighted.chi.K3 <- NA
        #fitted.values$predicted.unweighted.chi <- NA
    
        ##using unweighted LD scores excluding the intercept
        #for(i in 1:n.blocks){
        #  fitted.values[select.from[i]:select.to[i],]$predicted.weighted.chi <- LD.scores[select.from[i]:select.to[i],-length(reg)] %*% reg[-length(reg)]
        #  fitted.values[select.from[i]:select.to[i],]$predicted.weighted.chi.K <-  LD.scores[select.from[i]:select.to[i],-length(reg)] %*% delete.values[i,-length(reg)]
        #  fitted.values[select.from[i]:select.to[i],]$predicted.weighted.chi.K3 <- LD.scores[select.from[i]:select.to[i],-length(reg)] %*% delete.values.K3[i,-length(reg)]
       # }

        #fitted.values$predicted.weighted.chi.K_tau<-(fitted.values$predicted.weighted.chi)*merged$weights
        
        cat("Partitioning the heritability over the annotations","\n")
        
        ##Michels added code for standard heritability
        HSQ.MICHEL <- overlap.matrix %*% cats
        
        ###new delete values with N.bar added to denominator for standard heritability [AG]

        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values) <- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- ((solve(xtx.delete)%*% xty.delete)/N.bar)
        }
        
        #updated to multiply by m in delete.values to keep scaling consistent [AG]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(cats))
        colnames(pseudo.values) <- colnames(weighted.LD)[1:(ncol(weighted.LD)-1)]
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*cats)-((n.blocks-1)* delete.values[i,1:(ncol(delete.values)-1)]*m)}
        
        if(is.null(population.prev)==F & is.null(sample.prev)==F){
          conversion.factor <- (pop.prev^2*(1-pop.prev)^2)/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
          Liab.S[,j] <- conversion.factor
        }
        
        N.vec[1,s] <- N.bar 
        I[j,j] <- intercept
        
        ##loop saving relevant pieces to S and V
        f<-1
        for(f in 1:n.annot){
          S_List[[f]][j,j] <- HSQ.MICHEL[f]
          Tau_List[[f]][j,j]<- cats[f]
        }
        
        if(s == 1){
          t<-1
          d<-n.annot
        }else{t<-t+n.annot
        d<-(s*n.annot)}
        
        total_pseudo[,t:d]<-pseudo.values
      
        ### Total count
        s <- s+1
      }
      
      ##GENETIC COVARIANCE CODE
      if(j != k){
        
        chi2 <- traits[k]
        
        if(substr(x=chi2,start=nchar(chi2)-1,stop=nchar(chi2))=="gz"){
          y2 <- fread(paste("gzcat",chi2),header=T,showProgress=F,data.table=F)
        }else{
          y2 <- fread(chi2,header=T,showProgress=F,data.table=F)
        }
        
        y2 <- na.omit(y2)
        y2$chi2 <- y2$Z^2
        cat("Read in summary statistics from",chi2,"\n")
        cat("Read in summary statistics for",nrow(y2),"SNPs","\n")
        
        y <- merge(y1,y2,by="SNP")
        
        y[y$A1.y == y$A1.x,]$Z.x <- y[y$A1.y == y$A1.x,]$Z.x
        y[y$A1.y != y$A1.x,]$Z.x <-  -1* y[y$A1.y != y$A1.x,]$Z.x
        
        y$ZZ <- y$Z.y * y$Z.x
        y <- na.omit(y)
        
        cat("After merging",chi1,"and",chi2,"summary statistics for",nrow(y),"SNPs","\n")
        
        merged <- merge(x=y[,c("SNP","chi","chi2","N.x","N.y","ZZ")],y=w[,c("SNP","wLD")],by="SNP")
        
        merged <- merge(x=merged,y=x,by="SNP")
        
        merged <- merged[with(merged,order(CHR,BP)),]
        remaining.snps <- nrow(merged)
        cat(remaining.snps,"SNPs remaining","after merging the files","\n")
        
        ## REMOVE SNPS with excess chi-square:
        chisq.max1 <- max(0.001*max(merged$N.x),80)
        merged <- merged[merged$chi < chisq.max1,]
        
        chisq.max2 <- max(0.001*max(merged$N.y),80)
        merged <- merged[merged$chi2 < chisq.max2,]
        
        n.snps <- nrow(merged)
        removed.snps <- remaining.snps-n.snps
        
        cat("Removed",removed.snps,"SNPs with Chi^2 >",chisq.max2,paste0("(",n.snps," SNPs remain)"),"\n")
        
        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- rowSums(merged[,c(rownames(m))])
        merged$x.tot.intercept <- 1
        
        #### MAKE WEIGHTS for genetic covariance (Phenotype 1): 
        tot.agg <- (M.tot*(mean(merged$chi)-1))/mean(merged$x.tot*merged$N.x)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- as.numeric(lapply(X=merged$x.tot,function(x){max(x,1)}))
        merged$w.ld <- as.numeric(lapply(X=merged$wLD,function(x){max(x,1)}))
        merged$c <- tot.agg*merged$N.x/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)
        
        #### MAKE WEIGHTS for genetic covariance (Phenotype 2): 
        tot.agg2 <- (M.tot*(mean(merged$chi2)-1))/mean(merged$x.tot*merged$N.y)
        tot.agg2 <- max(tot.agg2,0)
        tot.agg2 <- min(tot.agg2,1)
        merged$ld2 <- as.numeric(lapply(X=merged$x.tot,function(x){max(x,1)}))
        merged$w.ld2 <- as.numeric(lapply(X=merged$wLD,function(x){max(x,1)}))
        merged$c2 <- tot.agg2*merged$N.y/M.tot
        merged$het.w2 <- 1/(2*(1+(merged$c2*merged$ld2))^2)
        merged$oc.w2 <- 1/merged$w.ld2
        merged$w2 <- merged$het.w2*merged$oc.w2
        merged$initial.w2 <- sqrt(merged$w2)
        
        merged$weights_cov <- (merged$initial.w + merged$initial.w2)/sum(merged$initial.w + merged$initial.w2 )
        
        N.bar <- sqrt(mean(merged$N.x)*mean(merged$N.y))
        
        merged$N <- sqrt((merged$N.x)*(merged$N.y))
        
        a<-which(colnames(merged)=="baseL2")
        b<-which(colnames(merged)=="intercept")
        LD.scores <- as.matrix(merged[,a:b])
        weighted.LD <- as.matrix(LD.scores*merged$weights_cov)
        weighted.chi <- as.matrix(merged$ZZ *merged$weights_cov)
        
        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:200]-1,n.snps)
        
        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow =((n.annot+1)* n.blocks),ncol =(n.annot+1))
        colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by =(n.annot+1))
        replace.to <- seq(from =(n.annot+1),to=nrow(xtx.block.values),by =(n.annot+1))
        
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
        colnames(xtx)<- colnames(weighted.LD)
        for(i in 1:nrow(xtx)){xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))}
        
        reg <- solve(xtx)%*% xty
        intercept <- reg[length(reg)]
        coefs <- reg[1:(length(reg)-1)]/N.bar
        cats <- coefs*m
        reg.tot <- sum(cats)
        
        est <- as.matrix(cats/reg.tot)
        
       delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
       delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
       delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        colnames(delete.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete)%*% xty.delete
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values_tau <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values_tau)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values_tau[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}
        
        jackknife.cov <- cov(pseudo.values_tau)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        cat.cov <- coef.cov*(m %*% t(m))
        #Tau.std.error=data.frame(sqrt(diag(cat.cov)))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        ##Added
        Coefficient.std.error=data.frame(sqrt(diag(coef.cov)))
        Coefficient.Z<-data.frame(coefs/sqrt(diag(coef.cov)))

        
        mean.ZZ <- mean(merged$ZZ)
        cat("Results for covariance between:",chi1,"and",chi2,"\n")
        cat("Mean Z*Z:",round(mean.ZZ,4),"\n")
        cat("Cross trait Intercept: ",round(intercept,4),"(",round(intercept.se,4),")","\n",sep="")
        cat("cov_g:",round(reg.tot,4),"(",round(tot.se,4),")","\n")
        
        cat("Partitioning the genetic covariance over the annotations","\n")
        
        COV.MICHEL <- overlap.matrix %*% cats
        
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values) <- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- ((solve(xtx.delete)%*% xty.delete)/N.bar)
        }
        
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(cats))
        colnames(pseudo.values) <- colnames(weighted.LD)[1:(ncol(weighted.LD)-1)]
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*cats)-((n.blocks-1)* delete.values[i,1:(ncol(delete.values)-1)]*m)}

        N.vec[1,s] <- N.bar 
        
        ##add in bivariate intercept to both sides of the intercept [I] matrix
        I[k,j] <- intercept
        I[j,k] <- intercept
        
        ##loop saving relevant pieces to S and V
        f<-1
        for(f in 1:n.annot){
          ##add in partitioned genetic covariances to both sides of the matrix
          S_List[[f]][k,j] <- COV.MICHEL[f]
          S_List[[f]][j,k] <- COV.MICHEL[f]
          Tau_List[[f]][k,j]<- cats[f]
          Tau_List[[f]][j,k]<- cats[f]
        }
        
        if(s == 1){
          t<-1
          d<-n.annot
        }else{t<-t+n.annot
        d<-(s*n.annot)}
        
        total_pseudo[,t:d]<-pseudo.values
        
        ### Total count
        s <- s+1
      }
    }
  }
  
  S<-vector(mode="list",length=n.annot)
  S_Tau<-vector(mode="list",length=n.annot)
  V<-vector(mode="list",length=n.annot)
  V_Tau<-vector(mode="list",length=n.annot)
  v.out<-vector(mode="list",length=n.annot)
  v.out_Tau<-vector(mode="list",length=n.annot)
  
  for (i in 1:n.annot) { 
    S[[i]] = matrix(NA,nrow=n.traits,ncol=n.traits)
    S_Tau[[i]] = matrix(NA,nrow=n.traits,ncol=n.traits)
    V[[i]] = matrix(NA, nrow=n.V,ncol=n.V)
    V_Tau[[i]] = matrix(NA, nrow=n.V,ncol=n.V)
    v.out[[i]] = matrix(NA, nrow=n.V,ncol=n.V)
    v.out_Tau[[i]] = matrix(NA, nrow=n.V,ncol=n.V)
  } 
  
  total_pseudo2<-(cov(total_pseudo)/n.blocks)
  Small_V<-vector(mode="list",length=n.V*n.V)
  i<-1
  for (i in 1:length(Small_V)) { 
    Small_V[[i]] = matrix(NA,nrow=n.annot,ncol=n.annot)
  } 
  
  #loop pulling chunks of total_pseudo2 
  r<-1
  u<-1
  p<-1
  for(u in 1:n.V){
    for(p in 1:n.V){
      if(u + p == 2){
        Small_V[[r]]<-total_pseudo2[1:n.annot,1:n.annot]
        }
      if(u == 1 & p != 1){
        n<-n.annot*(p-1)+1  
        v<-n+n.annot-1
        Small_V[[r]]<-total_pseudo2[1:n.annot,n:v]
      }
      if(u != 1 & p == 1){
        n<-n.annot*(u-1)+1  
        v<-n+n.annot-1
        Small_V[[r]]<-total_pseudo2[n:v,1:n.annot]
      }
      if(u != 1 & p != 1){
        n1<-n.annot*(u-1)+1  
        v1<-n1+n.annot-1
        n2<-n.annot*(p-1)+1  
        v2<-n2+n.annot-1
        Small_V[[r]]<-total_pseudo2[n1:v1,n2:v2]
      }
      ##total count
      r<-r+1
    }}
  
  ##loop calculating the sampling variances and covariances across partitions
  i<-1  
  SampleVar<-vector(mode="list",length=n.V*n.V)
  SampleVar_Tau<-vector(mode="list",length=n.V*n.V)
  for (i in 1:length(SampleVar)){ 
    SampleVar[[i]] = matrix(NA,nrow=n.annot,ncol=1)
    SampleVar_Tau[[i]] = matrix(NA,nrow=n.annot,ncol=1)
  } 
  
  i<-1
  for(i in 1:length(Small_V)){
    SampleVar[[i]]<-diag((overlap.matrix %*% Small_V[[i]])%*% t(overlap.matrix))
    ##the non-diagonal elements are the cross-partition sampling covariances
    ##^^Note for future cross-partition analyses
    
    SampleVar_Tau[[i]]<-diag(Small_V[[i]])
  }
  
  #loop inputting SampleVar into separate Vs for each partition
  r<-1
  u<-1
  i<-1
  s<-1
  for(i in 1:n.annot){
    for(u in 1:n.V){
      for(r in 1:n.V){
        v.out[[i]][u,r]<-SampleVar[[s]][i]
        v.out_Tau[[i]][u,r]<-SampleVar_Tau[[s]][i]
        s<-s+1
      }}
    s<-1}
  
  f<-1
  length<-ncol(S[[1]])
  for(f in 1:n.annot){
 
    ### Scale S and V to liability (Liab.S = matrix of 1s when no pop or samp prev provided)
    ##cov is what changes across partitions
    S[[f]] <- diag(as.vector(sqrt(Liab.S))) %*% S_List[[f]] %*% diag(as.vector(sqrt(Liab.S)))
    S_Tau[[f]]<-diag(as.vector(sqrt(Liab.S))) %*% Tau_List[[f]] %*% diag(as.vector(sqrt(Liab.S)))
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO=as.vector(lowerTriangle((S[[f]]/S_List[[f]]),diag=T))
    
    #calculate the ratio of the rescaled and original S_Tau matrices
    scaleO_Tau=as.vector(lowerTriangle((S_Tau[[f]]/Tau_List[[f]]),diag=T))
    
    #obtain diagonals of the original V matrix and take their sqrt to get SE's
    Dvcov<-sqrt(diag(v.out[[f]]))
    
    #obtain diagonals of the original V matrix and take their sqrt to get SE's
    Dvcov_Tau<-sqrt(diag(v.out_Tau[[f]]))
    
    #rescale the SEs by the same multiples that the S matrix was rescaled by
    Dvcovl<-as.vector(Dvcov*t(scaleO))
    
    #rescale the SEs by the same multiples that the S_Tau matrix was rescaled by
    Dvcovl_Tau<-as.vector(Dvcov_Tau*t(scaleO_Tau))
  
    #obtain the sampling correlation matrix by standardizing the original V matrix
    vcor<-cov2cor(v.out[[f]])
    
    #obtain the sampling correlation matrix by standardizing the original V_Tau matrix
    vcor_Tau<-cov2cor(v.out_Tau[[f]])
    
    #rescale the sampling correlation matrix by the appropriate diagonals
    V[[f]]<-diag(Dvcovl)%*%vcor%*%diag(Dvcovl)
   
    #rescale the sampling correlation matrix by the appropriate diagonals
    V_Tau[[f]]<-diag(Dvcovl_Tau)%*%vcor_Tau%*%diag(Dvcovl_Tau)

    ##name columns of S based on trait names. Provide general form V1:VX if none provided
    if(is.null(trait.names)){
      traits <- paste0("V",1:length)
      colnames(S[[f]])<-(traits)
      colnames(S_Tau[[f]])<-(traits)
      }else{colnames(S[[f]])<-(trait.names)
      colnames(S_Tau[[f]])<-(trait.names)
      }
    
    ##name list object by partition
    names(V)[[f]]<-names(SampleVar[[1]][f])
    names(S)[[f]]<-names(SampleVar[[1]][f])
    names(S_Tau)[[f]]<-names(SampleVar[[1]][f])
    names(V_Tau)[[f]]<-names(SampleVar[[1]][f])
  }
   
  Tau_Flag = matrix(NA,nrow=length(S_Tau),ncol=1)
  for(i in 1:length(S_Tau)){
    if(any(diag(S_Tau[[i]]) < 0)==TRUE){
      Tau_Flag[[i,1]]<-1
    }else{Tau_Flag[[i,1]]<-0}
  }
  rownames(Tau_Flag)<-names(S_Tau)

  Results<-(list(V=V,S=S,S_Tau=S_Tau,V_Tau=V_Tau,Tau_Flag=Tau_Flag,I=I,N=N.vec,m=m))
  save(Results, file = paste0(log2,".RData",sep=""))
  
  end.time <- Sys.time()
  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time-mins*60
  cat(paste("Analysis ended at",end.time),"\n")
  cat("Analysis took ",mins,"minutes and",secs,"seconds","\n")
  rm(list=setdiff(ls(),lsf.str()))
  gc()
  sink()
}
