\name{userGWAS2}
\alias{userGWAS2}
\title{Estimate SNP effects for a user specified multivariate GWAS}
\description{
Function to obtain model estimates for a user-specified model across SNPs and create the SNP specific covariance matrices in one step. 
}
\usage{
userGWAS2(covstruc,SNPs,estimation="DWLS",model="",modelchi=FALSE,printwarn=TRUE,sub=FALSE,toler=FALSE,SNPSE=FALSE, \dots)

}
\arguments{
   \item{covstruc}{Output from Genomic SEM `ldsc` function}
   \item{SNPs}{Summary statistics file created using the 'sumstats' function}
   \item{estimation}{The estimation method to be used when running the factor model. The options are Diagonally Weighted Least Squares ("DWLS", this is the default) or Maximum Likelihood ("ML")}
   \item{model}{The user-specified model to use in model estimation using lavaan syntax. The SNP is referred to as 'SNP' in the model.}
   \item{modelchi}{Whether you want model fit indices for each SNP. This can be very time consuming when done across all SNPs, so the default is FALSE.}
   \item{printwarn}{Whether you want warnings and errors printed for each run. This can take up significant space across all SNPs, but the default is set to TRUE as these warnings may not be safe to ignore.} 
   \item{sub}{Whether you want to only output a piece of the model results (e.g., F1 ~ SNP). The argument takes a vector, as multiple pieces of the model result can be output.}  
   \item{toler}{The tolerance to use for matrix inversion.} 
   \item{SNPSE}{Whether the user wants to provide a different standard error (SE) of the SNP variance than the package default. The default is to use 1E-6 to reflect the fact that the SNP SE is assumed to be population fixed.}
}

\value{
  The function outputs results from the multivariate GWAS. If the sub argument is used, it will output as many list objects as there are sub objects requested. 
  If the sub argument is FALSE (as is the package default), the function will ouptut as many list objects as there are SNPs.
}



\examples{

}
