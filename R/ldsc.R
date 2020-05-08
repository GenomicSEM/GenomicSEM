


ldsc <- function(traits,sample.prev,population.prev,ld,wld,trait.names=NULL,sep_weights = FALSE,chr=22,n.blocks=200,ldsc.log=NULL,stand=FALSE){
  time <- proc.time()
  
  begin.time <- Sys.time()
  
  if(is.null(ldsc.log)){
    log2<-paste(traits,collapse="_")
    if(object.size(log2) > 200){
      log2<-substr(log2,1,100)
    }
    log.file <- file(paste0(log2, "_ldsc.log"),open="wt")
  }else{log.file<-file(paste0(ldsc.log, "_ldsc.log"),open="wt")}
  
  log3<-paste(traits,collapse=" ")
  cat(print(paste0("Multivariate ld-score regression of ", length(traits), " traits ", "(", log3, ")", " began at: ",begin.time), sep = ""),file=log.file,sep="\n",append=TRUE)
  
  
  # Dimensions
  n.traits <- length(traits)
  n.V <- (n.traits^2 / 2) + .5*n.traits
  
  if(!(is.null(trait.names))){
    check_names<-str_detect(trait.names, "-")
    if(any(check_names==TRUE)){warning("Your trait names specified include mathematical arguments (e.g., + or -) that will be misread by lavaan. Please rename the traits using the trait.names argument.")}
  }
  
  if(length(traits)==1){warning("Our version of ldsc requires 2 or more traits. Please include an additional trait.")}
  
  
  # Storage:
  cov <- matrix(NA,nrow=n.traits,ncol=n.traits)
  V.hold <- matrix(NA,nrow=n.blocks,ncol=n.V)
  N.vec <- matrix(NA,nrow=1,ncol=n.V)
  Liab.S <- matrix(1,nrow=1,ncol=n.traits)
  I <- matrix(NA,nrow=n.traits,ncol=n.traits)
  
  # Current working directory
  curwd = getwd()
  
  
  
  #########  READ LD SCORES:
  cat(print("Reading in LD scores"),file=log.file,sep="\n",append=TRUE)
  setwd(ld)
  
  x <- suppressMessages(read_delim("1.l2.ldscore.gz", "\t", escape_double = FALSE, trim_ws = TRUE,progress = F))
  for(i in 2:chr){
    x <- rbind(x,suppressMessages(read_delim(paste0(i,".l2.ldscore.gz"), "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))
  }
  
  x$CM <- NULL
  x$MAF <- NULL
  
  setwd(curwd)
  
  ######### READ weights:
  
  
  setwd(wld)
  
  if(sep_weights==T){
    w <- suppressMessages(read_delim("1.l2.ldscore.gz", "\t", escape_double = FALSE, trim_ws = TRUE,progress = F))
    for(i in 2:chr){
      w <- rbind(w,suppressMessages(read_delim(paste0(i,".l2.ldscore.gz"), "\t", escape_double = FALSE, trim_ws = TRUE,progress = F)))
    }   }
  
  w <- x
  
  w$CM <- NULL
  w$MAF <- NULL
  
  setwd(curwd)
  
  colnames(w)[ncol(w)] <- "wLD"
  
  ### READ M
  setwd(ld)
  m  <- suppressMessages(read_csv("1.l2.M_5_50",  col_names = FALSE))
  for(i in 2:chr){
    m <- rbind(m,suppressMessages(read_csv(paste0(i,".l2.M_5_50"),  col_names = FALSE)))
  }
  setwd(curwd)
  M.tot <- sum(m)
  m <- M.tot
  
  # count the total nummer of runs, both loops
  s <- 1
  
  for(j in 1:n.traits){
    
    chi1 <- traits[j]
    ######### READ chi2
    
    cat(paste("     "),file=log.file,sep="\n",append=TRUE)
    cat(paste("     "),file=log.file,sep="\n",append=TRUE)
    
    cat(print(paste("Estimating heritability for:", traits[j])),file=log.file,sep="\n",append=TRUE)
    
    y1 <- suppressMessages(read_delim(chi1, "\t", escape_double = FALSE, trim_ws = TRUE,progress = F))
    
    y1 <- na.omit(y1)
    y1$chi1 <- y1$Z^2
    
    cat(print(paste("Read in summary statistics from:", chi1)),file=log.file,sep="\n",append=TRUE)
    cat(print(paste("Read in summary statistics for", nrow(y1), "SNPs")),file=log.file,sep="\n",append=TRUE)
    
    for(k in j:length(traits)){
      
      ##### HERITABILITY code
      
      if(j == k){
        
        samp.prev <- sample.prev[j]
        pop.prev <- population.prev[j]
        
        ######## Merge files
        
        merged <- merge(x=y1[,c("SNP","chi1","N")],y=w[,c("SNP","wLD")],by="SNP")
        merged <- merge(x=merged,y=x,by="SNP")
        merged <- merged[with(merged,order(CHR,BP)),]
        remaining.snps <- nrow(merged)
        
        cat(print(paste(remaining.snps, "SNPs remaining after merging file with LD-score files")),file=log.file,sep="\n",append=TRUE)
        
        
        ## REMOVE SNPS with excess chi-square:
        chisq.max <- max(0.001*max(merged$N),80)
        merged <- merged[merged$chi1 < chisq.max,]
        n.snps <- nrow(merged)
        removed.snps <- remaining.snps-n.snps
        
        cat(print(paste("Removed",removed.snps,"SNPs with Chi^2 >",chisq.max)),file=log.file,sep="\n",append=TRUE)
        cat(print(paste(n.snps, "SNPs remain")),file=log.file,sep="\n",append=TRUE)
        
        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1
        
        
        #### MAKE WEIGHTS:
        
        tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- as.numeric(lapply(X=merged$L2,function(x){max(x,1)}))
        merged$w.ld <- as.numeric(lapply(X=merged$wLD,function(x){max(x,1)}))
        merged$c <- tot.agg*merged$N/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)
        merged$weights <- merged$initial.w/sum(merged$initial.w)
        
        N.bar <- mean(merged$N)
        
        
        ## preweight LD and chi:
        
        weighted.LD <- as.matrix(cbind(merged$L2,merged$intercept)*merged$weights)
        weighted.chi <- as.matrix(merged$chi1*merged$weights)
        
        
        ## Perfrom analysis:
        
        n.annot <- 1
        
        
        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)
        
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
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m
        
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete)%*% xty.delete
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        V.hold[,s] <- pseudo.values[,1]
        N.vec[1,s] <- N.bar
        
        if(is.na(pop.prev)==F & is.na(samp.prev)==F){
          conversion.factor <- (pop.prev^2*(1-pop.prev)^2)/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
          Liab.S[,j] <- conversion.factor
          cat(paste("     "),file=log.file,sep="\n",append=TRUE)
          
          cat(print(paste("Please note that the results initially printed to the screen and log file reflect the NON-liability h2 and cov_g. However, a liability conversion is being used for trait", chi1, "when creating the genetic covariance matrix used as input for Genomic SEM and liability scale results are printed at the end of the log file.")),file=log.file,sep="\n",append=TRUE)
          cat(paste("     "),file=log.file,sep="\n",append=TRUE)
          t<-1
        }
        
        cov[j,j] <- reg.tot
        I[j,j] <- intercept
        
        lambda.gc <- median(merged$chi1)/0.4549
        mean.Chi <- mean(merged$chi1)
        ratio <- (intercept-1)/(mean(merged$chi1)-1)
        ratio.se <- intercept.se/(mean(merged$chi1)-1)
        
        cat(print(paste("Heritability Results for trait:",chi1)),file=log.file,sep="\n",append=TRUE)
        cat(print(paste("Lambda GC:",round(lambda.gc,4))),file=log.file,sep="\n",append=TRUE)
        cat(print(paste("Mean Chi^2 across remaining SNPs:",round(mean.Chi,4))),file=log.file,sep="\n",append=TRUE)
        cat(print(paste0("Intercept: ",round(intercept,4)," (",round(intercept.se,4),")")),file=log.file,sep="\n",append=TRUE)
        cat(print(paste("Lambda GC:",round(lambda.gc,4))),file=log.file,sep="\n",append=TRUE)
        cat(print(paste0("Ratio: ",round(ratio,4)," (",round(ratio.se,4),")")),file=log.file,sep="\n",append=TRUE)
        cat(print(paste0("Total Observed Scale h2: ",round(reg.tot,4)," (",round(tot.se,4),")")),file=log.file,sep="\n",append=TRUE)
        cat(print(paste0("h2 Z: ", format(reg.tot/tot.se),digits=3)),file=log.file,sep="\n",append=TRUE)
        
        ### Total count
        s <- s+1
        
      }
      
      
      ##### GENETIC COVARIANCE code
      
      if(j != k)
      {cat(paste("     "),file=log.file,sep="\n",append=TRUE)
        
        
        chi2 <- traits[k]
        ######### READ chi2
        cat(print(paste("Calculating genetic covariance for traits:",chi1, "and", chi2)),file=log.file,sep="\n",append=TRUE)
        
        # Reuse the data read in for heritability
        y2 <- suppressMessages(read_delim(chi2, "\t", escape_double = FALSE, trim_ws = TRUE,progress = F))
        
        y2 <- na.omit(y2)
        y2$chi2 <- y2$Z^2
        
        cat(print(paste("Read in summary statistics from",chi2)),file=log.file,sep="\n",append=TRUE)
        cat(print(paste("Read in summary statistics for",nrow(y2),"SNPs")),file=log.file,sep="\n",append=TRUE)
        
        
        y <- merge(y1,y2,by="SNP")
        
        y[y$A1.y == y$A1.x,]$Z.x <- y[y$A1.y == y$A1.x,]$Z.x
        y[y$A1.y != y$A1.x,]$Z.x <-  -1* y[y$A1.y != y$A1.x,]$Z.x
        
        
        y$ZZ <- y$Z.y * y$Z.x
        y <- na.omit(y)
        
        cat(print(paste("After merging",chi1,"and",chi2,"summary statistics for",nrow(y),"SNPs remain")),file=log.file,sep="\n",append=TRUE)
        
        ######## Merge files
        
        merged <- merge(x=y[,c("SNP","chi1","chi2","N.x","N.y","ZZ")],y=w[,c("SNP","wLD")],by="SNP")
        merged <- merge(x=merged,y=x,by="SNP")
        merged <- merged[with(merged,order(CHR,BP)),]
        remaining.snps <- nrow(merged)
        
        cat(print(paste(remaining.snps,"SNPs remaining after merging the files with LD-score files")),file=log.file,sep="\n",append=TRUE)
        
        ## REMOVE SNPS with excess chi-square:
        chisq.max1 <- max(0.001*max(merged$N.x),80)
        merged <- merged[merged$chi1 < chisq.max1,]
        
        chisq.max2 <- max(0.001*max(merged$N.y),80)
        merged <- merged[merged$chi2 < chisq.max2,]
        
        n.snps <- nrow(merged)
        removed.snps <- remaining.snps-n.snps
        
        cat(print(paste("Removed",removed.snps,"SNPs with Chi^2 >",chisq.max1, "for", chi1, "or SNPs with Chi^2 >",chisq.max2, "for", chi2, paste0("(",n.snps," SNPs remain)"))),file=log.file,sep="\n",append=TRUE)
        
        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1
        
        
        #### MAKE WEIGHTS:
        
        tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N.x)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- as.numeric(lapply(X=merged$L2,function(x){max(x,1)}))
        merged$w.ld <- as.numeric(lapply(X=merged$wLD,function(x){max(x,1)}))
        merged$c <- tot.agg*merged$N.x/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)
        
        tot.agg2 <- (M.tot*(mean(merged$chi2)-1))/mean(merged$L2*merged$N.y)
        tot.agg2 <- max(tot.agg2,0)
        tot.agg2 <- min(tot.agg2,1)
        merged$ld2 <- as.numeric(lapply(X=merged$L2,function(x){max(x,1)}))
        merged$w.ld2 <- as.numeric(lapply(X=merged$wLD,function(x){max(x,1)}))
        merged$c2 <- tot.agg2*merged$N.y/M.tot
        merged$het.w2 <- 1/(2*(1+(merged$c2*merged$ld))^2)
        merged$oc.w2 <- 1/merged$w.ld2
        merged$w2 <- merged$het.w2*merged$oc.w2
        merged$initial.w2 <- sqrt(merged$w2)
        
        
        merged$weights_cov <- (merged$initial.w + merged$initial.w2)/sum(merged$initial.w + merged$initial.w2 )
        
        N.bar <- sqrt(mean(merged$N.x)*mean(merged$N.y))
        
        ## preweight LD and chi:
        
        weighted.LD <- as.matrix(cbind(merged$L2,merged$intercept)*merged$weights)
        weighted.chi <- as.matrix(merged$ZZ *merged$weights_cov)
        
        ## Perfrom analysis:
        
        
        n.annot <- 1
        
        
        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)
        
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
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m
        
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete)%*% xty.delete
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        V.hold[,s] <- (pseudo.values[,1])
        N.vec[1,s] <- N.bar
        
        cov[k,j] <- reg.tot
        cov[j,k] <- reg.tot
        I[k,j] <- intercept
        I[j,k] <- intercept
        
        mean.ZZ <- mean(merged$ZZ)
        
        
        cat(print(paste("Results for genetic covariance between:",chi1,"and",chi2)),file=log.file,sep="\n",append=TRUE)
        cat(print(paste("Mean Z*Z:",round(mean.ZZ,4))),file=log.file,sep="\n",append=TRUE)
        cat(print(paste0("Cross trait Intercept: ",round(intercept,4)," (",round(intercept.se,4),")")),file=log.file,sep="\n",append=TRUE)
        cat(print(paste0("Total Observed Scale Genetic Covariance (g_cov): ",round(reg.tot,4)," (",round(tot.se,4),")")),file=log.file,sep="\n",append=TRUE)
        cat(print(paste0("g_cov Z: ", format(reg.tot/tot.se),digits=3)),file=log.file,sep="\n",append=TRUE)
        cat(print(paste0("g_cov P-value: ", format(2*pnorm(abs(reg.tot/tot.se),lower.tail=FALSE),digits=5))),file=log.file,sep="\n",append=TRUE)
        
        
        ### Total count
        s <- s+1
        
      }
      
    }
  }
  
  
  ## Scale V to N per study (assume m constant)
  v.out <- ((cov(V.hold)/n.blocks) /t(N.vec) %*% N.vec) * m^2
  
  ### Scale S and V to liability:
  S <- diag(as.vector(sqrt(Liab.S))) %*% cov %*% diag(as.vector(sqrt(Liab.S)))
  
  #calculate the ratio of the rescaled and original S matrices
  scaleO=as.vector(lowerTriangle((S/cov),diag=T))
  
  #obtain diagonals of the original V matrix and take their sqrt to get SE's
  Dvcov<-sqrt(diag(v.out))
  
  #rescale the SEs by the same multiples that the S matrix was rescaled by
  Dvcovl<-as.vector(Dvcov*t(scaleO))
  
  #obtain the sampling correlation matrix by standardizing the original V matrix
  vcor<-cov2cor(v.out)
  
  #rescale the sampling correlation matrix by the appropriate diagonals
  V<-diag(Dvcovl)%*%vcor%*%diag(Dvcovl)
  
  
  #name traits according to trait.names argument
  #use general format of V1-VX if no names provided
  if(is.null(trait.names)){
    traits2 <- paste0("V",1:ncol(S))
    colnames(S)<-(traits2)
  }else{
    colnames(S)<-(trait.names)
  }
  
  

  
  if(mean(Liab.S)!=1){
    r<-nrow(S)
    SE<-matrix(0, r, r)
    SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(V))
    
    cat(paste("     "),file=log.file,sep="\n",append=TRUE)
    cat(paste("     "),file=log.file,sep="\n",append=TRUE)
    cat(print(paste("Liability Scale Results")),file=log.file,sep="\n",append=TRUE)
    
    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j == k){
          cat(paste("     "),file=log.file,sep="\n",append=TRUE)
          cat(print(paste("Liability scale results for:", chi1)),file=log.file,sep="\n",append=TRUE)
          cat(print(paste0("Total Liability Scale h2: ",round(S[j,j],4)," (",round(SE[j,j],4),")")),file=log.file,sep="\n",append=TRUE)
        }
        
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          cat(print(paste0("Total Liability Scale Genetic Covariance between ", chi1, " and ",chi2, ": ", round(S[k,j],4)," (",round(SE[k,j],4),")")),file=log.file,sep="",append=TRUE)
          cat(paste("     "),file=log.file,sep="\n",append=TRUE)
        }
      }
    }
  }
  
  
  if(all(diag(S) > 0)){
    ##calculate standardized results to print genetic correlations to log and screen
    D<-sqrt(diag(diag(S)))
    S_Stand=solve(D)%*%S%*%solve(D)
    
    #obtain diagonals of the original V matrix and take their sqrt to get SE's
    Dvcov<-sqrt(diag(V))
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO=as.vector(lowerTriangle((S_Stand/S),diag=T))
    
    ## MAke sure that if ratio in NaN (devision by zero) we put the zero back in
    scaleO[is.nan(scaleO)] <- 0
    
    #rescale the SEs by the same multiples that the S matrix was rescaled by
    Dvcovl<-as.vector(Dvcov*t(scaleO))
    
    #obtain the sampling correlation matrix by standardizing the original V matrix
    Vcor<-cov2cor(V)
    
    #rescale the sampling correlation matrix by the appropriate diagonals
    V_Stand<-diag(Dvcovl)%*%Vcor%*%diag(Dvcovl)
    
    #enter SEs from diagonal of standardized V
    r<-nrow(S)
    SE_Stand<-matrix(0, r, r)
    SE_Stand[lower.tri(SE_Stand,diag=TRUE)] <-sqrt(diag(V_Stand))
 
    
    cat(paste("     "),file=log.file,sep="\n",append=TRUE)
    cat(paste("     "),file=log.file,sep="\n",append=TRUE)
    cat(print(paste("Genetic Correlation Results")),file=log.file,sep="\n",append=TRUE)
    
    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          cat(print(paste0("Genetic Correlation between ", chi1, " and ",chi2, ": ", round(S_Stand[k,j],4)," (",round(SE_Stand[k,j],4),")")),file=log.file,sep="",append=TRUE)
          cat(paste("     "),file=log.file,sep="\n",append=TRUE)
        }
      }
    }
  }else{
    warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
    cat(paste0("Your genetic covariance matrix includes traits estimated to have a negative heritability."),file=log.file,sep="",append=TRUE)
    cat(print(paste0("Genetic correlation results could not be computed due to negative heritability estimates.")),file=log.file,sep="",append=TRUE)
  }
  
  end.time <- Sys.time()
  
  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time-mins*60
  
  cat(paste("     "),file=log.file,sep="\n",append=TRUE)
  cat(print(paste0("LDSC finished running at ",end.time), sep = ""),file=log.file,sep="\n",append=TRUE)
  cat(print(paste0("Running LDSC for all files took ",mins," minutes and ",secs," seconds"), sep = ""),file=log.file,sep="\n",append=TRUE)
  cat(paste("     "),file=log.file,sep="\n",append=TRUE)
  
  if(stand == FALSE){
  return(list(V=V,S=S,I=I,N=N.vec,m=m))
  }
  
  if(stand == TRUE){
    return(list(V=V,S=S,I=I,N=N.vec,m=m,V_Stand=V_Stand,S_Stand=S_Stand))
  }
  
  
}
