\name{ldsc}
\alias{ldsc}
\title{build a convariance structure using LD score regression in R}
\description{
Function to run LD score regression (\link{https://github.com/bulik/ldsc}) to compute the genetic covariance between a series of traits based on genome wide summary statistics obtained from GWAS. The results generate by this function are necessary and sufficient to facilitate the fit of structural equation models (or other multivariate models) to the genetic covariance matrix.
}
\usage{
ldsc(traits,sample.prev,population.prev,ld,wld,trait.names=NULL, sep_weights = FALSE,chr=22,n.blocks=200,ldsc.log=NULL,stand=FALSE,select=FALSE, \dots)

}
\arguments{
   \item{traits}{A vector of strings which point to  LDSC munged files for trait you want to include in a Genomic SEM model.}
   \item{sample.prev}{A vector of sample prevalences for dichotomous traits and \code{NA} for continous traits}
   \item{population.prev}{A vector of population prevalences for dichotomous traits and \code{NA} for continous traits}
   \item{ld}{String which contains the path to the folder in which the LD scores used in the analysis are located. Expects LD scores formated as required by the original LD score regression software. }
   \item{wld}{String which contains the path to the folder in which the LD score weights used in the analysis are located. Expects LD scores formated as required by the original LD score regression software. }
   \item{trait.names}{A character vector specifying how the traits should be named in the genetic covariance (S) matrix. These variable names can subsequently be used in later steps for model specification. If no value is provided, the function will automatically name the variables using the generic from of V1-VX.}
   \item{sep_weights}{Logical which indicates wheter the weights are different form the LD scores used for the regression, defaults to FALSE. }
   \item{chr}{number of chromosomes over which the LDSC weights are split, defalts to 22 (Human) but can be switched for other species }
   \item{n.blocks}{Number of blocks to use for the jacknive procedure which is used to estiamte entries in V, higher values will be optimal if you have a large number of variables and also slower, defaults to 200  }
   \item{ldsc.log}{How you want to name your .log file for ldsc. The default is NULL as the package will automatically name the log based in the file names unless a log file name is provided to this arugment. }
   \item{stand}{Whether you want the package to also output a genetic correlation and sampling correlation matrix. Default is FALSE.}
   \item{select}{Whether you want the package to estimate LDSC using only even or odd chromosomes by setting select to "ODD" and "EVEN" respectively. It can also be set to a set of numbers, such as c(1,3,10), to run ldsc on a specific chromosome or chromosomes.  Default is FALSE, in which case LDSC is estimated using all chromosomes.}             
   \item{chisq.max}{Maximum value of the squared Z statistics for a SNP that is considered in the LD-score regression. Default behaviour is to set to the maximum of 80 and N*0.001}                 
       }

\value{
   The function returns a list with 5 named entries
   \item{S}{estimated genetic covariance matrix}
   \item{V}{variance covariance matrix of the parameter estimates in S}
   \item{I}{matrix containing the "cross trait intercepts", or the error covariance between traits induced by overlap, in terms of subjects, between the GWASes on which the analyses are based}
   \item{N}{a vector contsaining the sample size (for the genetic variances) and the geometric mean of sample sizes (i.e. sqrt(N1,N2)) between two samples for the covariances }
   \item{m}{number of SNPs used to compute the LD scores with. }
}



\examples{

}
