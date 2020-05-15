\name{addSNPs}
\alias{addSNPs}
\title{Combine LDSC and summary statistic output for multivariate GWAS using GenomicSEM}
\description{
Function to expand the S and V matrices to include SNP effects for multivariate GWAS in GenomicSEM
}
\usage{
addSNPs(covstruc, SNPs,SNPSE=FALSE,parallel=TRUE,cores=NULL,GC="standard", \dots)

}
\arguments{
   \item{covstruc}{Output from Genomic SEM `ldsc` function}
   \item{SNPs}{Summary statistics file created using the 'sumstats' function}
   \item{SNPSE}{Whether the user wants to provide a different standard error (SE) of the SNP variance than the package default. The default is to use .0005 to reflect the fact that the SNP SE is assumed to be population fixed.}
   \item{parallel}{addSNPs automatically uses mclapply to create the S and V matrices in parallel. Sometimes running in parallel can cause memory issues within the computing cores. If this is the case, the parallel argument can be set to FALSE, and addSNPs will create the S and V matrices serially.}
   \item{cores}{addSNPs automatically uses mclapply to create the S and V matrices in parallel. If the user does not provide an argument to the cores option, then addSNPs will automatically use one less than the total number of cores available.}
   \item{GC}{Level of Genomic Control (GC) you want the function to use. The default is 'standard' which adjusts the univariate GWAS standard errors by multiplying them by the square root of the univariate LDSC intercept. Additional options include 'conserv' which corrects standard errors using the univariate LDSC intercept, and 'none' which does not correct the standard errors.}
   
}

\value{
   The function expands the S and V matrices to include SNP effects. As many S and V matrices will be created as there are rows in the summary statistics file (i.e., one S and V matrix per SNP).
   The function returns a list with 3 named entries:
   \item{V_Full}{variance covariance matrix of the parameter estimates in S that includes an individual SNP effect}
   \item{S_Full}{genetic covaraiance matrix including individual SNP effect}
   \item{RS}{A list containing relevant genetic information (e.g., rsID, basepair, A1/A2) to be appended to the output from other functions (e.g., userGWAS)}
}



\examples{

}
