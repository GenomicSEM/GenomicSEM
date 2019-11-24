\name{userGWASpar}
\alias{userGWASpar}
\title{Estimate SNP effects for a user specified multivariate GWAS using parallel processing}
\description{
Function to obtain model estimates for a user-specified model across SNPs using parallel processing. 
}
\usage{
userGWASpar(Output,estimation="DWLS",model="",modelchi=FALSE,printwarn=TRUE,sub=FALSE,cores=NULL,toler=FALSE, \dots)

}
\arguments{
   \item{Output}{The Output obtained from using the addSNPs function. This includes as many S and V matrices as there are SNPs.}
   \item{estimation}{The estimation method to be used when running the factor model. The options are Diagonally Weighted Least Squares ("DWLS", this is the default) or Maximum Likelihood ("ML")}
   \item{model}{The user-specified model to use in model estimation using lavaan syntax. The SNP is referred to as 'SNP' in the model.}
   \item{modelchi}{Whether you want model fit indices for each SNP. This can be very time consuming when done across all SNPs, so the default is FALSE.}
   \item{printwarn}{Whether you want warnings and errors printed for each run. This can take up significant space across all SNPs, but the default is set to TRUE as these warnings may not be safe to ignore.} 
   \item{sub}{Whether you want to only output a piece of the model results (e.g., F1 ~ SNP). The argument takes a vector, as multiple pieces of the model result can be output.} 
   \item{cores}{The number of cores to use for parallel processing. If no number is provided, the function will automatically use 1 less than the total number of cores detected in the local environment.} 
   \item{toler}{The tolerance to use for matrix inversion.} 
}

\value{
  The function outputs a list of model outputs from multivariate GWAS using parallel processing. If the sub argument is used, it will output as many list objects as there are sub objects requested. 
  If the sub argument is FALSE (as is the package default), the function will ouptut as many list objects as there are SNPs.
}



\examples{

}
