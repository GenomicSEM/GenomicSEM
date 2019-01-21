\name{commonfactorGWASpar}
\alias{commonfactorGWASpar}
\title{Run Common Factor GWAS using parallel processing}
\description{
Function to estimate SNP effects, SEs, and estimates of SNP level heterogeneity on common factor using parallel processing. 
}
\usage{
commonfactorGWASpar(Output,estimation="DWLS",cores=NULL, \dots)

}
  \arguments{
   \item{Output}{The Output obtained from using the addSNPs function. This includes as many S and V matrices as there are SNPs.}
   \item{estimation}{The estimation method to be used when running the factor model. The options are Diagonally Weighted Least Squares ("DWLS", this is the default) or Maximum Likelihood ("ML")}
   \item{cores}{The number of cores to use on your computer for parallel processing. If no number is provided, the default is to use one less core then is available on your computer}
}

\value{
      The function outputs a series of SNP effects with their SEs and estimate of QSNP (the heterogeneity index) using parallel processing. The output is a single object.
}



\examples{

}
