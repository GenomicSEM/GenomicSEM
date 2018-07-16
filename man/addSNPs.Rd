\name{addSNPs}
\alias{addSNPs}
\title{Combined LDSC and summary statistic output for multivariate GWAS using GenomicSEM}
\description{
Function to expand the S and V matrices to include SNP effects for multivariate GWAS in GenomicSEM
}
\usage{
addSNPs(covstruc, SNPs, \dots)

}
\arguments{
   \item{covstruc}{Output from Genomic SEM multivariable LDSC}
   \item{SNPs}{Summary statistics file created using the sumstats function}
  
    
}

\value{
   The function expands the S and V matrices to include SNP effects. As many S and V matrices will be created as there are rows in the summary statistics file (i.e., one S and V matrix per SNP).
}



\examples{

}
