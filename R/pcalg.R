makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

GSEM.suffStat <- function(suffStat,seed=NULL,rep=250){
  set.seed(seed=seed)
  if(is.null(suffStat$S )){
  baz <- mvrnorm(rep,suffStat$S_Full[lower.tri(suffStat$S_Full,diag = T)],Sigma = suffStat$V_Full) 
  suffStat$baz <- baz
  }
  if(!is.null(suffStat$S)){
    baz <- mvrnorm(rep,suffStat$S[lower.tri(suffStat$S,diag = T)],Sigma = suffStat$V) 
    suffStat$baz <- baz
  }

  suffStat

}

GSEM.indepTest <- function(x,y,S,suffStat){

  if(is.null(suffStat$S )){
  sig <- matrix(NA,nrow(suffStat$S_Full),ncol(suffStat$S_Full))
  }else{
  sig <- matrix(NA,nrow(suffStat$S),ncol(suffStat$S))
  }
  
pcor <- matrix(NA,250,1)

for(i in 1:250){

sig[lower.tri(sig,diag = T)] <- suffStat$baz[i,]
sig<- makeSymm(sig)
pcor[i,] <- -cor2pcor(sig[c(x,y,S),c(x,y,S)])[1,2]

}



p_val <- pchisq((mean(pcor)/sd(pcor))^2,1,lower=F)
p_val
}


