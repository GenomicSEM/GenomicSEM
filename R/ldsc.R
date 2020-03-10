
ldsc <- function(traits,sample.prev,population.prev,ld,wld,trait.names=NULL,sep_weights=FALSE,chr=22,n.blocks=200){
  time <- proc.time()
  
  # Dimensions
  n.traits <- length(traits)
  n.V <- (n.traits^2 / 2) + .5*n.traits
  
  # Storage:
  cov <- matrix(NA,nrow=n.traits,ncol=n.traits)
  V.hold <- matrix(NA,nrow=n.blocks,ncol=n.V)
  N.vec <- matrix(NA,nrow=1,ncol=n.V)
  Liab.S <- matrix(1,nrow=1,ncol=n.traits)
  I <- matrix(NA,nrow=n.traits,ncol=n.traits)

  # Current working directory
  curwd = getwd()
  
  
  
  #########  READ LD SCORES:
  
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
    
    y1 <- suppressMessages(read_delim(chi1, "\t", escape_double = FALSE, trim_ws = TRUE,progress = F))
    
    y1 <- na.omit(y1)
    y1$chi1 <- y1$Z^2
    cat("Read in summary statistics from",chi1,"\n")
    cat("Read in summary statistics for",nrow(y1),"SNPs","\n")
    
    
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
        cat(remaining.snps,"SNPs remaining","after merging the files","\n")
        
        ## REMOVE SNPS with excess chi-square:
        chisq.max <- max(0.001*max(merged$N),80)
        merged <- merged[merged$chi1 < chisq.max,]
        n.snps <- nrow(merged)
        removed.snps <- remaining.snps-n.snps
        
        cat("Removed",removed.snps,"SNPs with Chi^2 >",chisq.max,paste0("(",n.snps," SNPs remain)"),"\n")
        
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
        #end.delete.values <- (tot.delete.values %*% m)/N.bar
        #write.table(x=end.delete.values,file=paste0(out,".end.del.val.txt"),quote=F,sep="\t",row.names=F)
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        #write.table(x=coef.cov,file=paste0(out,".coef.cov.txt"),quote=F,sep="\t")
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        V.hold[,s] <- pseudo.values[,1]
        N.vec[1,s] <- N.bar
        
        
        
        if(is.na(pop.prev)==F & is.na(samp.prev)==F){
          conversion.factor <- (pop.prev^2*(1-pop.prev)^2)/(samp.prev*(1-samp.prev)* dnorm(qnorm(1-pop.prev))^2)
          Liab.S[,j] <- conversion.factor
          
        }
        
        
        cov[j,j] <- reg.tot
        I[j,j] <- intercept
        
        lambda.gc <- median(merged$chi1)/0.4549
        mean.Chi <- mean(merged$chi1)
        ratio <- (intercept-1)/(mean(merged$chi1)-1)
        ratio.se <- intercept.se/(mean(merged$chi1)-1)
        
        cat("Results for trait",chi1,"\n")
        cat("Lambda GC:",round(lambda.gc,4),"\n")
        cat("Mean Chi^2:",round(mean.Chi,4),"\n")
        cat("Intercept: ",round(intercept,4),"(",round(intercept.se,4),")","\n")
        cat("Ratio: ",round(ratio,4),"(",round(ratio.se,4),")","\n")
        cat("h2:",round(reg.tot,4),"(",round(tot.se,4),")","\n")
        
        ### Total count
        s <- s+1
        
      }
      
      
      ##### GENETIC COVARIANCE code
      
      if(j != k)
      {
        
        
        chi2 <- traits[k]
        ######### READ chi2
        
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
        cat(remaining.snps,"SNPs remaining","after merging the files","\n")
        
        ## REMOVE SNPS with excess chi-square:
        chisq.max1 <- max(0.001*max(merged$N.x),80)
        merged <- merged[merged$chi1 < chisq.max1,]
        
        chisq.max2 <- max(0.001*max(merged$N.y),80)
        merged <- merged[merged$chi2 < chisq.max2,]
        
        n.snps <- nrow(merged)
        removed.snps <- remaining.snps-n.snps
        
        cat("Removed",removed.snps,"SNPs with Chi^2 >",chisq.max1, "for", chi1, "or SNPs with Chi^2 >",chisq.max2, "for", chi2, paste0("(",n.snps," SNPs remain)"),"\n")
       
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
        #end.delete.values <- (tot.delete.values %*% m)/N.bar
        #write.table(x=end.delete.values,file=paste0(out,".end.del.val.txt"),quote=F,sep="\t",row.names=F)
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        #write.table(x=coef.cov,file=paste0(out,".coef.cov.txt"),quote=F,sep="\t")
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
        cat("Results for covariance between:",chi1,"and",chi2,"\n")
        cat("Mean Z*Z:",round(mean.ZZ,4),"\n")
        cat("Cross trait Intercept: ",round(intercept,4),"(",round(intercept.se,4),")","\n")
        cat("cov_g:",round(reg.tot,4),"(",round(tot.se,4),")","\n")
        
        ### Total count
        s <- s+1
        
      }
      
    }
  }
  
  time_all <- proc.time()-time
  print(time_all[3])
  
  ## Scale V to N per study (assume m constant)
  v.out <- ((cov(V.hold)/n.blocks) /t(N.vec) %*% N.vec) * m^2
  
  ### Scale S and V to liability:
  S <- diag(as.vector(sqrt(Liab.S))) %*% cov %*% diag(as.vector(sqrt(Liab.S)))
  
  ##Elliot Method
  
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
 
  length<-ncol(S)
  if(is.null(trait.names)){

    traits <- paste0("V",1:length)
    colnames(S)<-(traits)

  }else{
    colnames(S)<-(trait.names)
  }
   
  return(list(V=V,S=S,I=I,N=N.vec,m=m))
  
  
}
 
