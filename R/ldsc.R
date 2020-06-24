ldsc <- function(traits,sample.prev,population.prev,ld,wld,trait.names=NULL,sep_weights=FALSE,chr=22,n.blocks=200,log.file='ldsc.log') {

  begin.time <- Sys.time()
  
  sink(log.file, append=FALSE, split=TRUE)

  writeLines( strwrap( paste0(
    "Multivariate LD-score regression of ", length(traits), " traits ", "(", paste(trait.names,collapse=", "), ")",
    " began at: ", begin.time, "\n"
  ) ) )

  # Dimensions
  n.traits <- length(traits)
  n.V <- (n.traits^2 / 2) + .5*n.traits
  
  if(!(is.null(trait.names))) {
    check_names<-str_detect(trait.names, "-")
    if(any(check_names==TRUE)) {
      warning( paste(
        "The trait names you specified include mathematical arguments (e.g., + or -) that will be misread by lavaan.",
        "Please rename the traits using the trait.names argument.", sep='\n'
      ) )
    }
  }
  
  if(length(traits)==1) {
    warning( "Our version of ldsc requires 2 or more traits to give meaningful results." )
  }
  
  # Storage:
  Gcov.mat <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I.mat <- matrix(NA,nrow=n.traits,ncol=n.traits)
  Gcov_p.mat <- matrix(NA,nrow=n.traits,ncol=n.traits)
  I_p.mat <- matrix(NA,nrow=n.traits,ncol=n.traits)
#  Gcov_pasympt.mat <- matrix(NA,nrow=n.traits,ncol=n.traits)
#  I_pasympt.mat <- matrix(NA,nrow=n.traits,ncol=n.traits)
  N.vec <- matrix(NA,nrow=1,ncol=n.V)
  V.hold <- matrix(NA,nrow=n.blocks,ncol=n.V)
  liab.scale.conv.fact <- matrix(1,nrow=1,ncol=n.traits)

  # Current working directory
  curwd = getwd()
  
  
  
  #########  READ LD SCORES:

  cat("Reading in LD scores.\n")
  
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
  m <- suppressMessages(read_csv("1.l2.M_5_50", col_names = FALSE))
  for(i in 2:chr){
    m <- rbind(m, suppressMessages(read_csv(paste0(i,".l2.M_5_50"), col_names = FALSE)))
  }
  setwd(curwd)
  M.tot <- sum(m)
  m <- M.tot  #is this right? in a couple of occasions below, 'm' is treated as an array
  
  # count the total nummer of runs, both loops
  s <- 1
  
  for(j in 1:n.traits) {
    
    chi1 <- traits[j]

    ######### READ chi2
    
    cat('\n\n')
    
    cat("Estimating heritability for:", traits[j], "\n")
    
    y1 <- suppressMessages(read_delim(chi1, "\t", escape_double = FALSE, trim_ws = TRUE, progress = F))
    
    y1 <- na.omit(y1)
    y1$chi1 <- y1$Z^2

    cat("Read in summary statistics from",chi1,"\n")
    cat("Read in summary statistics for",nrow(y1),"SNPs","\n")
    
    for(k in j:length(traits)) {
      
      ##### HERITABILITY code
      
      if(j == k) {
        
        samp.prev <- sample.prev[j]
        pop.prev <- population.prev[j]
        
        ######## Merge files
        
        merged <- merge(x=y1[,c("SNP","chi1","N")],y=w[,c("SNP","wLD")],by="SNP")
        merged <- merge(x=merged,y=x,by="SNP")
        merged <- merged[with(merged,order(CHR,BP)),]
        remaining.snps <- nrow(merged)

        cat(remaining.snps,"SNPs remaining","after merging the files","\n")
        
        ## REMOVE SNPS with excess chi-square:
        chisq.max <- max(0.001*max(merged$N),80)
        merged <- merged[merged$chi1 < chisq.max,]
        n.snps <- nrow(merged)
        removed.snps <- remaining.snps-n.snps
        
        cat("Removed",removed.snps,"SNPs with Chi^2 >",chisq.max,paste0("(",n.snps," SNPs remain)"),"\n")
        cat(n.snps, "SNPs remain\n")
        
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
        
        ## Perform analysis:
        
        n.annot <- 1
        
        select.from <- floor(seq(from=1,to=n.snps,length.out=(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)
        
        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol=(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow=((n.annot+1)*n.blocks),ncol=(n.annot+1))
        colnames(xty.block.values) <- colnames(xtx.block.values) <- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by=(n.annot+1))
        replace.to <- seq(from=(n.annot+1),to=nrow(xtx.block.values),by=(n.annot+1))
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],]) %*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],]) %*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow=(n.annot+1),ncol=(n.annot+1))
        colnames(xtx) <- colnames(weighted.LD)
        for(i in 1:nrow(xtx)) {
          xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))
        }
        
        reg <- solve(xtx) %*% xty
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m
        
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol=(n.annot+1))
        colnames(delete.values) <- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete) %*% xty.delete
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values) <- colnames(weighted.LD)
        for(i in 1:n.blocks) pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1) * delete.values[i,])
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        V.hold[,s] <- pseudo.values[,1]
        N.vec[1,s] <- N.bar
        
        if ( !is.na(pop.prev) && !is.na(samp.prev) ) {
          liab.scale.conv.fact[j] <- (pop.prev^2*(1-pop.prev)^2) / (samp.prev*(1-samp.prev)*dnorm(qnorm(1-pop.prev))^2)
          cat("\n")
          writeLines( strwrap( paste(
            "Please note that the results initially printed to screen and log file reflect the observed dichotomous",
            "phenotype's h2 and cov_g. However, a conversion to the liability scale is being used for dichotomous",
            "traits when creating the genetic covariance matrix used as input for Genomic SEM, and liability scale",
            "results are printed at the end of the log file."
          ) ) )
          cat("\n")
#           t<-1  # what is this?
        }
        
        Gcov.mat[j,j] <- reg.tot
        I.mat[j,j] <- intercept
        Gcov_p.mat[j,j] <- 2*pt( abs(reg.tot) / tot.se, df=n.blocks-2, lower.tail=F )
        I_p.mat[j,j] <- 2*pt( abs(intercept - 1) / intercept.se, df=n.blocks-2, lower.tail=F )
#         Gcov_pasympt.mat[j,j] <- pnorm( reg.tot / tot.se, lower.tail=F )
#         I_pasympt.mat[j,j] <- pnorm( intercept / intercept.se, lower.tail=F )

        lambda.gc <- median(merged$chi1)/0.4549
        mean.Chi <- mean(merged$chi1)
        ratio <- (intercept-1)/(mean(merged$chi1)-1)
        ratio.se <- intercept.se/(mean(merged$chi1)-1)
         
        cat("Heritability results for trait",chi1,"\n")
        cat("Mean Chi^2:",round(mean.Chi,4),"\n")
        cat("Intercept:",round(intercept,4),"(",round(intercept.se,4),")","\n")
        cat("Lambda GC:",round(lambda.gc,4),"\n")
        cat("Ratio:",round(ratio,4),"(",round(ratio.se,4),")","\n")
        cat("Total observed scale h2:",round(reg.tot,4),"(",round(tot.se,4),")","\n")
        cat("Total standardized h2:", format(reg.tot/tot.se,digits=3),"\n")
        cat("Total h2 P-value:", format(2*pnorm(abs(reg.tot/tot.se),lower.tail=FALSE),digits=5),"\n")
        
        ### Total count
        s <- s+1
        
      }
      
      
      ##### GENETIC COVARIANCE code
      
      if(j != k) {
        
        cat('\n')

        chi2 <- traits[k]

        ######### READ chi2

        cat("Calculating genetic covariance for traits:", chi1, "and", chi2, "\n")
        
        # Reuse the data read in for heritability
        
        y2 <- suppressMessages(read_delim(chi2, "\t", escape_double = FALSE, trim_ws = TRUE,progress = F))
        
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
        
        ######## Merge files
        
        merged <- merge(x=y[,c("SNP","chi1","chi2","N.x","N.y","ZZ")],y=w[,c("SNP","wLD")],by="SNP")
        merged <- merge(x=merged,y=x,by="SNP")
        merged <- merged[with(merged,order(CHR,BP)),]
        remaining.snps <- nrow(merged)
        
        writeLines( strwrap( paste( remaining.snps, "SNPs remaining after merging the files with LD-score files." ) ) )
        
        ## REMOVE SNPS with excess chi-square:
        chisq.max1 <- max(0.001*max(merged$N.x),80)
        merged <- merged[merged$chi1 < chisq.max1,]
        
        chisq.max2 <- max(0.001*max(merged$N.y),80)
        merged <- merged[merged$chi2 < chisq.max2,]
        
        n.snps <- nrow(merged)
        removed.snps <- remaining.snps-n.snps
        
        writeLines( strwrap( paste(
          "Removed", removed.snps, "SNPs with Chi^2 >", chisq.max1, "for", chi1, "or Chi^2 >", chisq.max2, "for", chi2,
          "(", n.snps, "SNPs remaining)"
        ) ) )
        
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
        merged$weights_cov <- (merged$initial.w + merged$initial.w2)/sum(merged$initial.w + merged$initial.w2)
        
        N.bar <- sqrt(mean(merged$N.x)*mean(merged$N.y))
        
        ## preweight LD and chi:
        
        weighted.LD <- as.matrix(cbind(merged$L2,merged$intercept)*merged$weights)
        weighted.chi <- as.matrix(merged$ZZ*merged$weights_cov)
        
        ## Perform analysis:
        
        n.annot <- 1
        
        select.from <- floor(seq(from=1,to=n.snps,length.out=(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)
        
        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol=(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow=((n.annot+1)*n.blocks),ncol=(n.annot+1))
        colnames(xty.block.values) <- colnames(xtx.block.values) <- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by=(n.annot+1))
        replace.to <- seq(from=(n.annot+1),to=nrow(xtx.block.values),by=(n.annot+1))
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],]) %*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],]) %*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow=(n.annot+1),ncol=(n.annot+1))
        colnames(xtx) <- colnames(weighted.LD)
        for(i in 1:nrow(xtx)) {
          xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))
        }
        
        reg <- solve(xtx) %*% xty
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m
        
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol=(n.annot+1))
        colnames(delete.values) <- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete) %*% xty.delete
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values) <- colnames(weighted.LD)
        for(i in 1:n.blocks) pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1) * delete.values[i,])
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        V.hold[,s] <- pseudo.values[,1]
        N.vec[1,s] <- N.bar
        
        Gcov.mat[k,j] <- reg.tot
        Gcov.mat[j,k] <- reg.tot
        I.mat[k,j] <- intercept
        I.mat[j,k] <- intercept
        Gcov_p.mat[k,j] <- 2*pt( abs(reg.tot) / tot.se, df=n.blocks-2, lower.tail=F )
        Gcov_p.mat[j,k] <- Gcov_p.mat[k,j]
        I_p.mat[k,j] <- 2*pt( abs(intercept) / intercept.se, df=n.blocks-2, lower.tail=F )
        I_p.mat[j,k] <- I_p.mat[k,j]
#         Gcov_pasympt.mat[k,j] <- pnorm( reg.tot / tot.se, lower.tail=F )
#         Gcov_pasympt.mat[j,k] <- Gcov_pasympt.mat[k,j]
#         I_pasympt.mat[k,j] <- pnorm( intercept / intercept.se, lower.tail=F )
#         I_pasympt.mat[j,k] <- I_pasympt.mat[k,j]
        
        mean.ZZ <- mean(merged$ZZ)
        cat("Results for genetic covariance between",chi1,"and",chi2,"\n")
        cat("Mean Z*Z:",round(mean.ZZ,4),"\n")
        cat("Cross trait Intercept:",round(intercept,4),"(",round(intercept.se,4),")","\n")
        cat("Total observed scale genetic covariance (cov_g):",round(reg.tot,4),"(",round(tot.se,4),")","\n")
        cat("Total standardized genetic covariance:", format(reg.tot/tot.se,digits=3),"\n")
        cat("Total genetic covariance P-value:", format(2*pnorm(abs(reg.tot/tot.se),lower.tail=FALSE),digits=5),"\n")
      
        ### Total count
        s <- s+1
        
      }
      
    }

  }
  
  
  ## Scale V to N per study (assume m constant)
  v.out <- ((cov(V.hold)/n.blocks) / t(N.vec) %*% N.vec) * m^2
  
  ### Scale S and V to liability:
  S.mat <- diag(as.vector(sqrt(liab.scale.conv.fact))) %*% Gcov.mat %*% diag(as.vector(sqrt(liab.scale.conv.fact)))
  
  #calculate the ratio of the rescaled and original S matrices
  scaleO=as.vector(lowerTriangle((S.mat/Gcov.mat),diag=T))
  
  #obtain diagonals of the original V matrix and take their sqrt to get SE's
  Dvcov<-sqrt(diag(v.out))
  
  #rescale the SEs by the same multiples that the S matrix was rescaled by
  Dvcovl<-as.vector(Dvcov*t(scaleO))
  
  #obtain the sampling correlation matrix by standardizing the original V matrix
  vcor<-cov2cor(v.out)
  
  #rescale the sampling correlation matrix by the appropriate diagonals
  V.mat <- diag(Dvcovl)%*%vcor%*%diag(Dvcovl)
 
  n.traits<-ncol(S.mat)
  utlog<-upper.tri(S.mat, diag=T)
  if(is.null(trait.names)) {
    tmp.names <- paste0("V",1:n.traits)
    colnames(S.mat)<-tmp.names
    rownames(S.mat)<-tmp.names
    colnames(I.mat)<-tmp.names
    rownames(I.mat)<-tmp.names
    colnames(Gcov_p.mat)<-tmp.names
    rownames(Gcov_p.mat)<-tmp.names
    colnames(I_p.mat)<-tmp.names
    rownames(I_p.mat)<-tmp.names
    colnames(V.mat)<-c( matrix(
      unlist( lapply( tmp.names, paste, tmp.names, sep='.' ) ),
      length( tmp.names )
    )[utlog] )
    rownames(V.mat)<-colnames(V.mat)
    colnames(N.vec)<-colnames(V.mat)
  } else {
    colnames(S.mat)<-(trait.names)
    rownames(S.mat)<-(trait.names)
    colnames(I.mat)<-(trait.names)
    rownames(I.mat)<-(trait.names)
    colnames(Gcov_p.mat)<-(trait.names)
    rownames(Gcov_p.mat)<-(trait.names)
    colnames(I_p.mat)<-(trait.names)
    rownames(I_p.mat)<-(trait.names)
    colnames(V.mat)<-c( matrix(
      unlist( lapply(
        make.names(substr(trait.names,1,3), unique=T), paste,
        make.names(substr(trait.names,1,3), unique=T), sep='.'
      ) ),
      length( trait.names )
    )[utlog] )
    rownames(V.mat)<-colnames(V.mat)
    colnames(N.vec)<-colnames(V.mat)
  }
  
  if (mean(liab.scale.conv.fact)!=1) {
    r<-nrow(S.mat)
    SE<-matrix(0, r, r)
    SE[lower.tri(SE,diag=TRUE)]<-sqrt(diag(V.mat))
    
    cat("\n\nLiability scale results\n")
    cat( rep( '-', nchar("Liability scale results") ), sep='' )
    
    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j == k){
          cat("\n")
          cat("Liability scale results for", chi1, "\n")
          cat("Liability scale h2", ":", round(S.mat[j,j],4), "(", round(SE[j,j],4), ")\n")
        }
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          cat("Liability scale genetic covariance between", chi1, "and", chi2, ":", round(S.mat[k,j],4), "(", round(SE[k,j],4), ")\n")
        }
      }
    }
  }
  
  if(all(diag(S.mat) > 0)) {
    ##calculate standardized results to print genetic correlations to log and screen
    D<-sqrt(diag(diag(S.mat)))
    S_Stand=solve(D)%*%S.mat%*%solve(D)
    
    #obtain diagonals of the original V matrix and take their sqrt to get SE's
    Dvcov<-sqrt(diag(V.mat))
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO=as.vector(lowerTriangle((S_Stand/S.mat),diag=T))
    
    ## MAke sure that if ratio in NaN (devision by zero) we put the zero back in
    scaleO[is.nan(scaleO)] <- 0
    
    #rescale the SEs by the same multiples that the S matrix was rescaled by
    Dvcovl<-as.vector(Dvcov*t(scaleO))
    
    #obtain the sampling correlation matrix by standardizing the original V matrix
    Vcor<-cov2cor(V.mat)
    
    #rescale the sampling correlation matrix by the appropriate diagonals
    V_Stand<-diag(Dvcovl)%*%Vcor%*%diag(Dvcovl)
    
    #enter SEs from diagonal of standardized V
    r<-nrow(S.mat)
    SE_Stand<-matrix(0, r, r)
    SE_Stand[lower.tri(SE_Stand,diag=TRUE)]<-sqrt(diag(V_Stand))
    
    cat( "\n\n\nGenetic correlation results\n" )
    
    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          writeLines( strwrap( paste0(
            "Genetic correlation between ", chi1, " and ", chi2, ": ", round(S_Stand[k,j],4),
            " (", round(SE_Stand[k,j],4), ")\n"
          ) ) )
        }
      }
    }
  } else {
    warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
    cat("Your genetic covariance matrix includes traits estimated to have a negative heritability.\n")
    cat("Genetic correlations could not be computed due to negative heritability estimates.\n")
  }
  
  end.time <- Sys.time()
  
  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- total.time-mins*60
  
  writeLines( strwrap( paste0( 
    "\nFinished running at ", end.time, "."
  ) ) )
  writeLines( strwrap( paste0( 
    "Running ldsc for all files took ", mins," minutes and ", round(secs)," seconds."
  ) ) )

  sink()

  return( list(
    V=V.mat,
    S=S.mat,
    I=I.mat,
    N=N.vec,
    m=m,
    V_stand=V_Stand,
    S_stand=S_Stand,
    pS=Gcov_p.mat,
    pI=I_p.mat,
  #     pSalt=Gcov_pasympt.mat,
  #     pIalt=I_pasympt.mat,
    n=n.blocks
  ) )
  
}

