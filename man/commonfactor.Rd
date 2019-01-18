\name{commonfactor}
\alias{commonfactor}
\title{Run common factor model on genetic covariance and sampling covariance matrix}
\description{
Function to run a common factor model based on output from multivariable LDSC
}
\usage{
commonfactor(covstruc,estimation="DWLS", \dots)

}
\arguments{
   \item{covstruc}{Output from the multivariable LDSC function of Genomic SEM}
   \item{estimation}{Options are either Diagonally Weighted Least Squares ("DWLS"; the default) or Maximum Likelihood ("ML")} 
    
}

\value{
  The function estimates a common factor model, along with model fit indices, using output from GenomicSEM LDSC. 
  The function returns a list with 2 named entries
   \item{modelfit}{The model fit results (e.g., model chi-square, AIC, CFI) from the specified model.}
   \item{results}{Parameter estimates and sandwich corrected standard errors from the specified model.}
}


\examples{

}
