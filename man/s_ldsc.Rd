\name{s_ldsc}
\alias{s_ldsc}
\title{Estimate genetic covariance matrices within functional annotations using multivariable Stratified LD score regression}
\description{
Function to run Stratified LD score regression.
}
\usage{
s_ldsc(traits,sample.prev=NULL,population.prev=NULL,ld,wld,frq,trait.names=NULL,n.blocks=200,ldsc.log=NULL,exclude_cont=TRUE, \dots)

}
\arguments{
   \item{traits}{A vector of file names which point to  LDSC munged files for trait you want to include.}
   \item{sample.prev}{A vector of sample prevalences for dichotomous traits and \code{NA} for continous traits. Default = NULL.}
   \item{population.prev}{A vector of population prevalences for dichotomous traits and \code{NA} for continous traits. Default = NULL.}
   \item{ld}{A folder (or folders) of partitioned LD scores used as the independent variable in S-LDSC.}
   \item{wld}{A folder of non-partitioned LD scores used as regression weights.}
   \item{frq}{A folder of allele frequency files.}
   \item{trait.names}{A character vector specifying how the traits should be named. These variable names can subsequently be used in later steps for model specification.}
   \item{n.blocks}{Number of blocks to use for the jacknive procedure which is used to estiamte entries in V, higher values will be optimal if you have a large number of variables and also slower, defaults to 200.}   
   \item{ldsc.log}{What to name the .log file if you want to overrride default to name file based on file names used as input.}   
   \item{exclude_cont}{Whether to exclude continuous annotations from S-LDSC estimation.}   
 }

\value{
   The function returns a list with 9 named entries
   \item{S}{The zero-order genetic covariance matrices for each annotation.}
   \item{V}{The zero-order sampling covariance matrices for each annotation.}
   \item{S_Tau}{The tau matrices for each annotation.}
   \item{V_Tau}{The tau sampling covariance matrices for each annotation.}
   \item{I}{matrix containing the "cross trait intercepts", or the error covariance between traits induced by overlap, in terms of subjects, between the GWASes on which the analyses are based}
   \item{N}{a vector contsaining the sample size (for the genetic variances) and the geometric mean of sample sizes (i.e. sqrt(N1,N2)) between two samples for the covariances }
   \item{m}{number of SNPs used to compute the LD scores with. }
   \item{Prop}{The proportional size of each annotation relative to the annotation containing all SNPs.}
   \item{Select}{A data.frame that codes flanking window and continuous annotations as 2 and all other annotations as 1. This is used by the `enrich` function to exclude the flanking window and continuous annotations from enrichment estimates.}
}



\examples{

}
