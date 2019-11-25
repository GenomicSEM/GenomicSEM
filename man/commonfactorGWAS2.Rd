\name{commonfactorGWAS2}
\alias{commonfactorGWAS2}
\title{Estimate SNP effects and run addSNPs for a single common factor in one}
\description{
Function to obtain SNP effects on common factor along with index of SNP heterogeneity
}
\usage{
commonfactorGWAS(covstruc,SNPs,estimation="DWLS",toler=FALSE,SNPSE=FALSE, \dots)

}
\arguments{
   \item{covstruc}{Output from Genomic SEM `ldsc` function}
   \item{SNPs}{Summary statistics file created using the 'sumstats' function}
   \item{estimation}{The estimation method to be used when running the factor model. The options are Diagonally Weighted Least Squares ("DWLS", this is the default) or Maximum Likelihood ("ML")}
   \item{toler}{The tolerance to use for matrix inversion.} 
   \item{SNPSE}{Whether the user wants to provide a different standard error (SE) of the SNP variance than the package default. The default is to use 1E-6 to reflect the fact that the SNP SE is assumed to be population fixed.}
}

\value{
  The function outputs a series of SNP effects with their SEs and estimate of QSNP (the heterogeneity index). The output is a single object.

}



\examples{

}
