\name{userGWAS}
\alias{userGWAS}
\title{Estimate SNP effects for a user specified multivariate GWAS}
\description{
Function to obtain model estimates for a user-specified model across SNPs. 
}
\usage{
userGWAS(Output,estimation="DWLS",model="",modelchi=FALSE, \dots)

}
\arguments{
   \item{Output}{The Output obtained from using the addSNPs function. This includes as many S and V matrices as there are SNPs.}
   \item{estimation}{The estimation method to be used when running the factor model. The options are Diagonally Weighted Least Squares ("DWLS", this is the default) or Maximum Likelihood ("ML")}
   \item{model}{The user-specified model to use in model estimation using lavaan syntax. The SNP is referred to as 'SNP' in the model.}
   \item{modelchi}{Whether you want model fit indices for each SNP. This can be very time consuming when done across all SNPs, so the default is FALSE.}
    
}

\value{
  The function outputs a list of model outputs from multivariate GWAS. As it is unclear what piece of the output the user is most interested in for the user-specified case, the function returns all elements of model results. 

}



\examples{

}
