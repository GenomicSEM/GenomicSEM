\name{multiSNP}
\alias{multiSNP}
\title{Combine LDSC, summary statistic output, and LD information for models including multiple SNPs}
\description{
Function to expand the S and V matrices to include multiple SNP effects in a single matrix, along with LD information across these SNPs. 
}
\usage{
multiSNP (covstruc, SNPs, LD, SNPSE = FALSE, SNPlist = NA, \dots)

}
\arguments{
   \item{covstruc}{Output from Genomic SEM multivariable LDSC}
   \item{SNPs}{Summary statistics file created using the sumstats function}
   \item{LD}{Matrix of LD information across the SNPs. If only independent SNPs are being provided a matrix of 0s can be entered. Note that the function requires that A1 and A2 be included in the LD matrix column names (e.g., rs12345_A_T)}
   \item{SNPSE}{User provided SE of the SNP variance for entry in the V matrix. If no number is provided the package defaults to using .0005 to reflect a practically fixed population value taken from a reference panel}
   \item{SNPlist}{List of rsIDs if the user wishes to subset out a set of SNPs from a full set of summary statistics}
     
}

\value{
   The function expands the S and V matrices to include multiple SNP effects. These matrices include LD information across the SNPs.
}



\examples{

}
