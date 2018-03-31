\name{sumstats}
\alias{sumstats}
\title{Allign summary statistics from univariate GWAS for a GWAS in GenomicSEM}
\description{
Function to process GWAS summary statistis files and prepair them for a GWAS in genomicSEM
}
\usage{
sumstats(files,ref,trait.names=NULL,se.logit,info.filter = .6,maf.filter=0.01, \dots)

}
\arguments{
   \item{files}{a list of file names, fuiles must be located in the working directory, or a path must be provided.}
   \item{ref}{A reference file of SNPs to keep in your GWAS, one based on 1000 genomes phase 3 is provided. }
   \item{trait.names}{ a vector of trait names which will be used as names for the munged files}
     \item{se.logt}{ a logical vector indivating whether the standard errors in each set of summary statistics is on the logit scale}
   \item{info.filter}{ Numeric value which is used aas a lower bound for inputation quality (INFO) }
     \item{maf.filter}{ Num,eric value used as a lower bound for minor allel frequency  }
    
}

\value{
   The function ensures the SNPs in each file are aligned to the same reference allele, it attempts to filter strand issues, it retains SNPs present in the reference file. 
}



\examples{

}
