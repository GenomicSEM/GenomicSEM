\name{commonfactorGWAS}
\alias{commonfactorGWAS}
\title{Estimate SNP effects on a single common factor}
\description{
Function to obtain SNP effects on common factor along with index of SNP heterogeneity
}
\usage{
commonfactorGWAS(Output,estimation="DWLS", \dots)

}
\arguments{
   \item{Output}{The Output obtained from using the addSNPs function. This includes as many S and V matrices as there are SNPs.}
   \item{estimation}{The estimation method to be used when running the factor model. The options are Diagonally Weighted Least Squares ("DWLS", this is the default) or Maximum Likelihood ("ML")}   
}

\value{
  The function outputs a series of SNP effects with their SEs and estimate of QSNP (the heterogeneity index). The output is a single object.

}



\examples{

}
