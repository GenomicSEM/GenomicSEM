\name{sumstats}
\alias{sumstats}
\title{Allign summary statistics from univariate GWAS for a GWAS in GenomicSEM}
\description{
Function to process GWAS summary statistics files and prepare them for a GWAS in genomicSEM
}
\usage{
sumstats(filenames,reference,trait.names=NULL,se.logit=NULL,model=NULL,prev=NA,N=NULL,info.filter=.6,maf.filter=0.01,keep.indel=TRUE,parallel=TRUE,cores=NULL,log.prefix='sumstats', \dots)
}
\arguments{
  \item{filenames}{
    A list of file names, files are assumed to be located in the working directory unless an absolute path is provided.
  }
  \item{reference}{
    A reference file of SNPs to keep in your GWAS, one based on 1000 genomes phase 3 is provided.
  }
  \item{trait.names}{
    A vector of trait names which will be used as names for the munged filesi.
  }
  \item{se.logit}{
    A logical vector indicating whether the standard errors in each set of summary statistics is on the logit scale.
  }
  \item{model}{
    A character vector indicating the model underlying the summary statistics ('LIN','LOG','LPM','OLS'). See details.
  }
  \item{prev}{
    A vector of sample prevalences for dichotomous traits (used for 'LPM', may be set to NA for other models).
  }
  \item{N}{
    A vector of sample sizes (used for 'LPM' and 'OLS', may be set to NA for 'DEF' models).
  }
  \item{info.filter}{
    Numeric value used as a lower bound for imputation quality (INFO).
  }
  \item{maf.filter}{
    Numeric value used as a lower bound for minor allel frequency.
  }
  \item{keep.indel}{
    Whether to keep insertions/deletions or not.
  }
  \item{parallel}{
    Indicates whether the summary statistics files should be processed in parallel or serial fashion. Default is TRUE.
  }
  \item{cores}{
    Indicates how many cores to use when running in parallel. The default is NULL, in which case sumstats will use one
    fewer than the total number of cores available in the local environment.
  }
  \item{log.prefix}{
    Prefix for the log file(s).
  }
}

\value{
}

\details{
   The function parses multiple summary statistics files to extract association coefficients and their standard errors,
   and merge them, alongside the reference, in one data.frame to be input in GenomicSEM models that include single
   variants. The function ensures that the SNPs in each file are aligned to the same reference allele, attempts to
   filter strand issues and retains SNPs present in the reference file. The function can deal with GWAS of continous
   traits, dichotomous traits modeled by logistic regressions or dichotomous traits modeled by (misspecified) OLS
   regressions or mixed models. The function writes .log files that should be inspected to ensure that all column names
   were appropriately interpreted.
   'LIN' models are assumed for all summary statistics if no models are specified (model=NULL). For models other than
   'LIN', both effects and standard errors will be adjusted (LOG) or replaced by estimates based on p-value (OLS) or
   reference minor allele frequency, sample size and eventual sample prevalence (LPM). 'LPM' and 'OLS' specify summary
   statistics from a continuous or dichotomous trait obtained using OLS (or LMM). Note that the default model 'LIN'
   assumes the association coefficients to be standardized with respect to both predictors (genotypes) and outcome.
}


\examples{

}

