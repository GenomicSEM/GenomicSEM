\name{GSEM.process.covStruct}
\alias{GSEM.process.covStruct}
\title{build a convariance structure using LD score regression in R}
\description{
Function to run LD score regression (\link{https://github.com/bulik/ldsc}) to compute the genetic covariance between a series of traits based on genome wide summary statistics obtained from GWAS. The results generate by this function are necessary and sufficient to facilitate the fit of structural equation models (or other multivariate models) to the genetic covariance matrix.
}
\usage{
LDSC(traits,sample.prev,population.prev,ld,wld, sep_weights = FALSE, \dots)

}
\arguments{
   \item{traits}{a vector of strings which point to  LDSC munged files for trait you want to include in a Genomic SEM model.}
   \item{sample.prev}{ a vector of sample prevalences for dichotomous traits and \code{NA} for continous traits}
   \item{population.prev}{ a vector of population prevalences for dichotomous traits and \code{NA} for continous traits}
   \item{ld}{ string which contains the path to the folder in which the LD scores used in the analysis are located. Expects LD scores formated as required by the original LD score regression software. }
     \item{wld}{ string which contains the path to the folder in which the LD score weights used in the analysis are located. Expects LD scores formated as required by the original LD score regression software. }
       \item{sep_weights}{ ological which indicates wheter the weights are different form the LD scores used for the regression, defaults to FALSE. }
}

\value{
   The function returns a list with 5 named entries
   \item{S}{estimated genetic covariance matrix}
   \item{V}{varianve covariance matrix of the parameter estimates in S}
   \item{I}{matrix containing the "cross trait intercepts", or the error covariance between traits induced by overlap, in terms of subjects, between the GWASes on which the analyses are based}
   \item{N}{a vector contsaining the sample size (for the genetic variances) and the geometric mean of sample sizes (i.e. sqrt(N1,N2)) between two samples for the covariances }
  \item{m}{number of SNPs used to compute the LD scores with. }

}



\examples{

}
