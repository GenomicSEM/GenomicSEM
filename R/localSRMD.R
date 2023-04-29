localSRMD <- function(unconstrained, constrained, lhsvar, rhsvar){
  
  print(paste0("calculating localSRMD across h = ",length(unconstrained)/length(lhsvar[[1]])," parameters in g = ",length(lhsvar[[1]])," groups (",length(unconstrained)," parameters total)"))
  
  #convert lists to data frames
  lhsframe <- t(data.frame(lhsvar))
  rhsframe <- t(data.frame(rhsvar))
  
  #pool standard deviation for each variable
  lhspooledsd <- sqrt(rowMeans(lhsframe))
  rhspooledsd <- sqrt(rowMeans(rhsframe))
  
  #concatenate all information 
  df <- as.data.frame(cbind(unconstrained, constrained, lhspooledsd, rhspooledsd))
  
  #calculate each parameter's contribution to localsrmd
  df$localdiff <- ((df$unconstrained - df$constrained)/(df$lhspooledsd*df$rhspooledsd))^2
  
  value <- sqrt(mean(df$localdiff, na.rm = TRUE))
  
  return(value)
}
