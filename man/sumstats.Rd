\name{sumstats}
\alias{sumstats}
\title{Allign summary statistics from univariate GWAS for a GWAS in GenomicSEM}
\description{
Function to process GWAS summary statistis files and prepair them for a GWAS in genomicSEM
}
\usage{
sumstats(files,ref,trait.names=NULL,se.logit,OLS=NULL,linprob=NULL,prop=NULL,N=NULL,info.filter = .6,maf.filter=0.01,parallel=TRUE,cores=NULL, \dots)
}
\arguments{
   \item{files}{a list of file names, fuiles must be located in the working directory, or a path must be provided.}
   \item{ref}{A reference file of SNPs to keep in your GWAS, one based on 1000 genomes phase 3 is provided. }
   \item{trait.names}{ a vector of trait names which will be used as names for the munged files}
   \item{se.logit}{ a logical vector indicating whether the standard errors in each set of summary statistics is on the logit scale}
   \item{OLS}{ a logical vector indicating whether the GWAS was for a continuous trait and used OLS (or a LMM)}
   \item{linprob}{a logical vector indicating whether the GWAS used a linear probability model i.e. a dichotomous trait using OLS (or a LMM)}
   \item{info.filter}{ Numeric value which is used as a lower bound for inputation quality (INFO) }
   \item{maf.filter}{ Numeric value used as a lower bound for minor allel frequency}
   \item{parallel}{Indicates whether sumstats should process the summary statistics files in parallel or serial fashion. Default is TRUE, indicating that it will run in parallel.}
   \item{cores}{Indicates how many cores to use when running in parallel. The default is NULL, in which case sumstats will use 1 less than the total number of cores available in the local environment.}
   
}

\value{
   The function ensures the SNPs in each file are aligned to the same reference allele, it attempts to filter strand issues, it retains SNPs present in the reference file. The function can deal with GWAS of continous traits, dichotomous traits using logistic regression and even dichotomous traits using (misspecified) OLS regression or a mixed model. 
   The function returns .log files that should be inspected to ensure that all column names were appropriately interpreted. 
      }



\examples{

}
