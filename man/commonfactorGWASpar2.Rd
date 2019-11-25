\name{commonfactorGWASpar2}
\alias{commonfactorGWASpar2}
\title{Run addSNPs and common Factor GWAS using parallel processing in one step}
\description{
Function to estimate SNP effects, SEs, and estimates of SNP level heterogeneity on common factor using parallel processing. 
}
\usage{
commonfactorGWASpar2(covstruc,SNPs,estimation="DWLS",cores=NULL,SNPSE=FALSE, \dots)

}
  \arguments{
   \item{covstruc}{Output from Genomic SEM `ldsc` function}
   \item{SNPs}{Summary statistics file created using the 'sumstats' function}
   \item{estimation}{The estimation method to be used when running the factor model. The options are Diagonally Weighted Least Squares ("DWLS", this is the default) or Maximum Likelihood ("ML")}
   \item{cores}{The number of cores to use on your computer for parallel processing. If no number is provided, the default is to use one less core then is available on your computer}
   \item{SNPSE}{Whether the user wants to provide a different standard error (SE) of the SNP variance than the package default. The default is to use 1E-6 to reflect the fact that the SNP SE is assumed to be population fixed.}
}

\value{
      The function outputs a series of SNP effects with their SEs and estimate of QSNP (the heterogeneity index) using parallel processing. The output is a single object.
}



\examples{

}
