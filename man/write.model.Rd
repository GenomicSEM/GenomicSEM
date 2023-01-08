\name{write.model}
\alias{write.model}
\title{Automate writing model syntax using EFA output}
\description{
Function to automate writing model syntax based on EFA loadings. This is most likely to be useful when examining larger numbers of traits (e.g., > 10).
}
\usage{
write.model(Loadings,S_LD,cutoff,fix_resid=TRUE,bifactor=FALSE,mustload=FALSE,common=FALSE, \dots)
}
\arguments{
   \item{Loadings}{The matrix of EFA loadings. Note that the number of columns in this matrix determines how many factors are specifeid in the model.}
   \item{S_LD}{The LDSC genetic covariance matrix}
   \item{cutoff}{The EFA standardized loadings cutoff to determine which traits should load on a factor}
   \item{fix_resid}{Whether to apply constraint on all variables to keep residual variances above .001. Default is TRUE.}
   \item{bifactor}{Whether to specify a bifactor model in which a general factor predicts all included traits and the remaining factors are specifided to be orthogonal of one another.} 
   \item{mustload}{Whether all variables should load on at least one factor, even if they dont meet the threshold specified using the cutoff argument.}
   \item{common}{Whether to specify a common factor model.}
}

\value{
  The function outputs model syntax that can be used to run the model using the usermodel function in Genomic SEM. 
}



\examples{

}
