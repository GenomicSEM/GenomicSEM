\name{enrich}
\alias{enrich}
\title{Estimate enrichment of model parameter for a user specified model}
\description{
Function to take output from multivariable S-LDSC and estimate enrichment of model parameter for user specified model
}
\usage{
enrich(s_covstruc, model = "",params,fix= "regressions",std.lv=FALSE,rm_flank=TRUE,tau=FALSE,base=TRUE,toler=NULL,fixparam=NULL, \dots)

}
\arguments{
   \item{s_covstruc}{Output from the multivariable S-LDSC function of Genomic SEM (s_ldsc)}
   \item{model}{Model to be specified using lavaan notation} 
   \item{params}{Parameters of interest to be examined for enrichment (e.g., factor variances).} 
   \item{fix}{What components of the model should be fixed for follow-up enrichment models. Default = "regressions", which will fix all regression parameters.}
   \item{std.lv}{Optional argument to denote whether all latent variables are standardized using unit variance identification (default = FALSE)}
   \item{rm_flank}{Optional argument to denote whether flanking window annotations should automatically be removed from output (default = TRUE)}
   \item{tau}{Optional argument to denote whether the user wants to use the tau genetic covariance matrices, as opposed to the default zero-order matrices, for estimation of enrichment (default = FALSE)}
   \item{base}{Optional argument to denote whether the user wants to include the full model output from the genome-wide (i.e., baseline) matrix (default = TRUE)}
   \item{toler}{Optional argument to manually set tolerance for matrix inverison.}
   \item{fixparam}{Optional argument to manually fix paramters when estimating the model within annotations}
}

\value{
 Function to take output from multivariable S-LDSC and estimate enrichment of model parameter for user specified model 

}



\examples{

}
