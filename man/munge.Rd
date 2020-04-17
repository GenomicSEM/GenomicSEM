\name{munge}
\alias{munge}
\title{Clean and munge files to enable LD score regression}
\description{
Function to process GWAS summary statistics files and prepare them for LD score regression
}
\usage{
munge(files,reference,trait.names=NULL,N,info.filter=.9,maf.filter=0.01,munge.log='munge.log' \dots)

}
\arguments{
  \item{files}{
    A list of file names. the files are assumed to be located in the working directory unless an
    absolute path is provided.
  }
  \item{reference}{
    A variant list file with A1, A2 and rsID used to align alleles across traits. We suggest using
    the (UNZIPPED) list of HAPMAP3 SNPs with some basic cleaning applied (e.g., MHC region removed)
    created and supplied by the original LD score regression developers
    (\url{https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2}).
  }
  \item{trait.names}{
    A vector of trait names which will be used as names for the munged files.
  }
  \item{info.filter}{
    Numeric value which is used as a lower bound for imputation quality (INFO).
  }
  \item{maf.filter}{
    Numeric value used as a lower bound for minor allele frequency.
  }
  \item{munge.log}{
    File log.
  }
}

\value{
   The function writes files in the ".sumstats" format, which can be used to estimate SNP
   heritability and genetic covariance using the ldsc() function. The function will also write a
   .log file that should be examined to ensure that column names are being interpret correctly.
}


\examples{

}

