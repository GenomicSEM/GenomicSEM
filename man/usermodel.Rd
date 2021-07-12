\name{usermodel}
\alias{usermodel}
\title{Run user specified model on LDSC output}
\description{
Function to run a user specified model based on output from multivariable LDSC
}
\usage{
usermodel(covstruc,estimation="DWLS", model = "", CFIcalc=TRUE,std.lv=FALSE,imp_cov=FALSE,fix_resid=TRUE,toler=FALSE, \dots)

}
\arguments{
   \item{covstruc}{Output from the multivariable LDSC function of Genomic SEM}
   \item{estimation}{Options are either Diagonally Weighted Least Squares ("DWLS"; the default) or Maximum Likelihood ("ML")} 
   \item{model}{Model to be specified using lavaan notation} 
   \item{CFIcalc}{Optional argument to denote whether CFI is being requested (default = TRUE). In some cases the estimation of the independent (i.e., Null) model for calculation of CFI can be time consuming. If the funciton seems to be stuck on this step, we would suggest re-running with this option set to FALSE} 
   \item{std.lv}{Optional argument to denote whether all latent variables are standardized using unit variance identification (default = FALSE)}
   \item{imp_cov}{Optional argument to denote whether the user wants the model implied and residual covariance matrix included in the usermodel output (default = FALSE)}
   \item{fix_resid}{Optional argument to denote whether the user wants Genomic SEM to try troubleshooting a model that does not converge by fixing residual variances to be above 0 (default = TRUE)}
   \item{toler}{Optional argument to set lower tolerance for matrix inversion used to produce sadnwich corrected standard errors. (default = FALSE)}
}

\value{
  The function estimates a user-specified model, along with model fit indices, using output from GenomicSEM LDSC. 

}



\examples{

}
