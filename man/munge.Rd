\name{munge}
\alias{munge}
\title{Clean and munge files to enable LD score regression}
\description{
Function to process GWAS summary statistis files and prepare them for LD score regression
}
\usage{
munge(files,hm3,trait.names=NULL,N,info.filter=.9,maf.filter=0.01, \dots)

}
\arguments{
   \item{files}{A list of file names. the files are assumed to be located in the working directory unless an absolute path is provided.}
   \item{hm3}{A file (UNZIPPED) of HAPMAP3 SNPs with some basic cleaning applied. one was created and supplied by the original LD score regression developers and is available here: https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2}
   \item{trait.names}{A vector of trait names which will be used as names for the munged files}
   \item{info.filter}{Numeric value which is used as a lower bound for inputation quality (INFO)}
   \item{maf.filter}{Numeric value used as a lower bound for minor allele frequency}
    
}

\value{
   The function writes files of the ".sumstats" format, which can be used to estimate SNP heritability and genetic covariance using the ldsc() function.  
   The function will also output a .log file that should be examined to ensure that column names are being interpret correctly. 
}



\examples{

}

