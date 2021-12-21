\name{munge}
\alias{munge}
\title{Clean and munge files to enable LD score regression}
\description{
Function to process GWAS summary statistis files and prepair them for LD score regression 
}
\usage{
munge(files,hm3,trait.names=NULL,N,info.filter = .9,maf.filter=0.01, column.names=list(),parallel=FALSE,cores=NULL,overwrite=TRUE \dots)

}
\arguments{
   \item{files}{A vector of file names, files must be located in the working directory, or a path must be provided.}
   \item{hm3}{A file of SNPs with A1, A2 and rsID used to allign alleles across traits. We suggest using an (UNZIPPED) file of HAPMAP3 SNPs with some basic cleaning applied (e.g., MHC region removed) that is supplied and created by the original LD score regression developers and available here: https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2: }
   \item{trait.names}{A vector of trait names which will be used as names for the munged files}
   \item{N}{A vector of sample size}
   \item{info.filter}{Numeric value which is used as a lower bound for inputation quality (INFO)}
   \item{maf.filter}{Numeric value used as a lower bound for minor allel frequency}
   \item{column.names}{Optional list detailing which columns represent, SNP, MAF, etc. e.g. list(SNP=my_snp_column)}
   \item{parallel}{Indicates whether munge should process the summary statistics files in parallel or serial fashion. Default is TRUE, indicating that it will run in parallel.}
   \item{cores}{Indicates how many cores to use when running in parallel. The default is NULL, in which case munge will use 1 less than the total number of cores available in the local environment.}
   \item{overwrite}{Indicates whether existing .sumstats.gz files should be overwritten}
}

\value{
   The function writes files of the ".sumstats" format, which can be used to estimate SNP heritability and genetic covariance using the ldsc() function.  
   The function will also output a .log file that should be examined to ensure that column names are being interpret correctly. 
}



\examples{

}
