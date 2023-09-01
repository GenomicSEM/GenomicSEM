\name{userGWAS}
\alias{userGWAS}
\title{Create genetic covariance matrices for individual SNPs and estimate SNP effects for a user specified multivariate GWAS}
\description{
Function to obtain model estimates for a user-specified model across SNPs. 
}
\usage{
userGWAS(covstruc=NULL,SNPs=NULL,estimation="DWLS",model="",printwarn=TRUE,sub=FALSE,cores=NULL,toler=FALSE,SNPSE=FALSE,
         parallel=TRUE,GC="standard",MPI=FALSE,smooth_check=FALSE,TWAS=FALSE,std.lv=FALSE,fix_measurement=TRUE, \dots)

}
\arguments{
   \item{covstruc}{Output from Genomic SEM `ldsc` function}
   \item{SNPs}{Summary statistics file created using the 'sumstats' function}
   \item{estimation}{The estimation method to be used when running the factor model. The options are Diagonally Weighted Least Squares ("DWLS", this is the default) or Maximum Likelihood ("ML")}
   \item{model}{The user-specified model to use in model estimation using lavaan syntax. The SNP is referred to as 'SNP' in the model.}
   \item{printwarn}{Whether you want warnings and errors printed for each run. This can take up significant space across all SNPs, but the default is set to TRUE as these warnings may not be safe to ignore.} 
   \item{sub}{Whether you want to only output a piece of the model results (e.g., F1 ~ SNP). The argument takes a vector, as multiple pieces of the model result can be output.}
   \item{cores}{Indicates how many cores to use when running in parallel. The default is NULL, in which case sumstats will use 1 less than the total number of cores available in the local environment.}
   \item{toler}{The tolerance to use for matrix inversion.} 
   \item{SNPSE}{Whether the user wants to provide a different standard error (SE) of the SNP variance than the package default. The default is to use .0005 to reflect the fact that the SNP SE is assumed to be population fixed.}
   \item{parallel}{Whether the function should run using parallel or serial processing. Default = TRUE}
   \item{GC}{Level of Genomic Control (GC) you want the function to use. The default is 'standard' which adjusts the univariate GWAS standard errors by multiplying them by the square root of the univariate LDSC intercept. Additional options include 'conserv' which corrects standard errors using the univariate LDSC intercept, and 'none' which does not correct the standard errors.}
   \item{MPI}{Whether the function should use multi-node processing (i.e., MPI). Please note that this should only be used on a computing cluster on which the R package Rmpi is already installed.}
   \item{smooth_check}{Whether the function should save the consequent largest Z-statistic difference between the pre and post-smooth matrices.}
   \item{TWAS}{Whether the function is being used to estimate a multivariate TWAS using read_fusion output for the SNPs argument.} 
   \item{std.lv}{Optional argument to denote whether all latent variables are standardized using unit variance identification (default = FALSE)}
   \item{fix_measurement}{Optional argument to denote whether the measurement model should be fixed across all SNPs (default = TRUE)}
}

\value{
  The function outputs results from the multivariate GWAS. If the sub argument is used, it will output as many list objects as there are sub objects requested. 
  If the sub argument is FALSE (as is the package default), the function will ouput as many list objects as there are SNPs.
}



\examples{

}
