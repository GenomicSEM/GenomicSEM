\name{munge}
\alias{munge}
\title{Clean and munge files to enable LD score regression}
\description{
Function to process GWAS summary statistis files and prepair them for LD score regression 
}
\usage{
munge(files,hm3,trait.names=NULL,N,info.filter = .6,maf.filter=0.05, \dots)

}
\arguments{
   \item{files}{A list of file names, fuiles must be located in the working directory, or a path must be provided.}
   \item{hm3}{A file (UNZIPPED) of HAPMAP3 SNPs with some basic cleaning applied it is supplied and created by the original LD score regression developers and available here: https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2: }
   \item{trait.names}{A vector of trait names which will be used as names for the munged files}
   \item{info.filter}{Numeric value which is used aas a lower bound for inputation quality (INFO)}
     \item{maf.filter}{Numeric value used as a lower bound for minor allel frequency}
    
}

\value{
   The function writes files of the ".sumstats" format, which can be used to estimate SNP heritability and genetic covariance using the ldsc() function.  

}



\examples{

}
