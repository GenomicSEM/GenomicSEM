\name{hdl}
\alias{hdl}
\title{build a convariance structure using High defenition Likelihood estimation in R}
\description{
Function to run HDL  (\link{https://github.com/zhenin/HDL}) to compute the genetic covariance between a series of traits based on genome wide summary statistics obtained from GWAS. The results generate by this function are necessary and sufficient to facilitate the fit of structural equation models (or other multivariate models) to the genetic covariance matrix. HDL is more powerfull than LDSC but if the LD structure in the reference file mismatches the GWAS LD structure, LDSC seems to perfrom better, expescially for estiamtes of heritability. For medium samples (N > 50.000) with moderate SNP-h2 (snp h2 > 0.07) where the LD structure isnt similar we would recomend ldsc, especially for GWAS. If you have small GWAS ( N < 25.000) the extra power HDL provides is worth the downward bias in snp h2 estimates relative to ldsc. 
}
\usage{
ldsc(traits,sample.prev,population.prev,ld,wld,trait.names=NULL, sep_weights = FALSE, \dots)

}
\arguments{
   \item{traits}{A vector of strings which point to munged files for trait you want to include in a Genomic SEM model. the HDL function works with standard munged files}
   \item{sample.prev}{A vector of sample prevalences for dichotomous traits and \code{NA} for continous traits}
   \item{population.prev}{A vector of population prevalences for dichotomous traits and \code{NA} for continous traits}
   \item{LD.path}{String which contains the path to the folder in which the LD matrices used in the analysis are located. Expects LD matirices formated as required by the original HDL software. }
   \item{Nref]{Sample size of the reference file, default is 335265}
   \item{trait.names}{A character vector specifying how the traits should be named in the genetic covariance (S) matrix. These variable names can subsequently be used in later steps for model specification. If no value is provided, the function will automatically name the variables using the generic from of V1-VX.}
   \item{method}{sting, either "piecewise" whcih estimates the heritability or genetic covariance locally in chunks acros the genome and then sums these, or "jackknife"}
}

\value{
   The function returns a list with 5 named entries
   \item{S}{estimated genetic covariance matrix}
   \item{V}{variance covariance matrix of the parameter estimates in S}
   \item{I}{matrix containing the "cross trait intercepts", or the error covariance between traits induced by overlap, in terms of subjects, between the GWASes on which the analyses are based}

}



\examples{

}
