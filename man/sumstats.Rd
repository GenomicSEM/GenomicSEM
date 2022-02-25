\name{sumstats}
\alias{sumstats}
\title{Allign summary statistics from univariate GWAS for a GWAS in GenomicSEM}
\description{
Function to process GWAS summary statistics files and prepare them for a GWAS in genomicSEM
}
\usage{
sumstats(files,ref,trait.names=NULL,se.logit,OLS=NULL,linprob=NULL,N=NULL,betas=NULL,info.filter = .6,maf.filter=0.01,
         keep.indel=FALSE,parallel=FALSE,cores=NULL,ambig=FALSE,direct.filter=FALSE, \dots)
}
\arguments{
   \item{files}{a vector of file names, files must be located in the working directory, or a path must be provided.}
   \item{ref}{A reference file of SNPs to keep in your GWAS, one based on 1000 genomes phase 3 is provided. }
   \item{trait.names}{ a vector of trait names which will be used as names for the munged files}
   \item{se.logit}{ a logical vector indicating whether the standard errors in each set of summary statistics is on the logit scale}
   \item{OLS}{ a logical vector indicating whether the GWAS was for a continuous trait and used OLS (or a LMM)}
   \item{linprob}{a logical vector indicating whether the GWAS is a binary outcome with only Z-statistics or was analyzed using a linear probability model i.e. a dichotomous trait using OLS (or a LMM)}
   \item{N}{A vector of total sample sizes for continuous traits and the sum of effective sample sizes for binary traits}
   \item{betas}{A vector of column names of betas for continuous traits that are known to have been standardized prior to running the GWAS}
   \item{N}{A vector of sample size}
   \item{info.filter}{ Numeric value which is used as a lower bound for imputation quality (INFO) }
   \item{maf.filter}{ Numeric value used as a lower bound for minor allele frequency}
   \item{keep.indel}{Indicates whether insertion-deletion mutations (indels) should be included in your summary statistics. The default is FALSE.}
   \item{parallel}{Indicates whether sumstats should process the summary statistics files in parallel or serial fashion. Default is TRUE, indicating that it will run in parallel.}
   \item{cores}{Indicates how many cores to use when running in parallel. The default is NULL, in which case sumstats will use 1 less than the total number of cores available in the local environment.}
   \item{ambig}{Indicates whether strand ambiguous SNPs should be removed from output. }
   \item{direct.filter}{Indicates whether SNPs that have missing information for more than half of contributing cohorts, as indicated by missing information in the direction column, should be removed.}        
   
}

\value{
   The function ensures the SNPs in each file are aligned to the same reference allele, it attempts to filter strand issues, it retains SNPs present in the reference file. The function can deal with GWAS of continous traits, dichotomous traits using logistic regression and even dichotomous traits using (misspecified) OLS regression or a mixed model. 
   The function returns .log files that should be inspected to ensure that all column names were appropriately interpreted. 
      }



\examples{

}
