\name{localSRMD}
\alias{localSRMD}
\title{localSRMD for Genomic measurement invariance models}
\description{
Given a set of parameter values from Genomic SEM models, calculate the extent to which groups differ on these parameters. 
}

\usage{
 localSRMD(unconstrained, constrained, lhsvar, rhsvar, \dots)
}

\arguments{  
   \item{unconstrained}{A vector of parameter values from an unconstrained structural equation model, where focal parameters are estimated freely in each group.}
   \item{constrained}{A vector of parameter values from a constrained structural equation model, where focal parameters are constrained to be equal across groups.}
   \item{lhsvar}{A list containing the variances, in each group, of the variables in the usermodel() results column lhs}
   \item{rhsvar}{A list containing the variances, in each group, of the variables in the usermodel() results column rhs}

}

\value{
 The function returns the average standardized extent to which estimates from a constrained set of structural equation model parameters differ from those obtained when the same set of parameters are freely estimated.
}



\examples{

}
