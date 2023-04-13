\name{read_fusion}
\alias{read_fusion}
\title{Format univariate FUSION TWAS output across multiple traits for subsequent use in a multivariate TWAS [T-SEM]}
\description{
Function to format TWAS summary statistics for T-SEM
}
\usage{
read_fusion(files,trait.names=NULL,binary=NULL,N=NULL,perm=FALSE, \dots)

}
\arguments{
   \item{files}{List of FUSION files.}
   \item{trait.names}{What names to use when naming the beta and SE columns for each trait}
   \item{binary}{Vector specifying whether each trait is binary (TRUE) or continuous (FALSE)}
   \item{N}{The sample size to use when backing out the betas and SEs from FUSION Z-statistics. This should reflect the sum of effective sample size for binary traits}
   \item{perm}{Whether you want to use the permutation test p-values from FUSION. Default = FALSE}
   
   
}

\value{
   The function combines and formats FUSION TWAS summary statistics for subsequent T-SEM analyses.
   The function returns a single data.frame. 
}



\examples{

}
