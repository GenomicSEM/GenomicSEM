# GenomicSEM

Code for structural equation modeling based on the summary statistics obtained from genome wide association studies (GWAS). The wiki pages on this github contain documentention on the use of the package. Untill explicitly stated otherwise the code on this github in an alpha version (0.0.1) and under active development. The code may thus produce undesired results on certain operating systems or when run concurrently with specific packadges or R versions. Feel free to raise issues if (or when...) the package produces undesired results, we will attempt to swiftly deal with known issues. 

There are types of models in Genomic SEM: ones that include effects of individual SNPs and ones that don't. The following goes through an example of each. 

# Models without SNP Effects

There are three primary steps necessary to running a structural model that does not include SNP effects. The first is to prepare the summary statistics for multivariate LD-Score regression. We note here that the user can use previous munge filed produced by the LD-Score regression packages, but that our package is also capable of preparing the summary statistics. In the process of preparing the summary statistics, it is also important to know that for case/control designs you should provide the total sample size and NOT the effective sample size. Both the original LDSC software and Genomic SEM correct for disproportionate numbers of cases/controls by using the input percentage of cases/controls for estimation of liability scale heritabilities and covariances. In order to use our package to prep the summary statistics [this will get filled in once the processing function is updated]


The second step is to run multivariate LD-Score regression to obtain the genetic covariance (i.e., S) matrix and corresponding sampling covariance matrix (i.e., V). This is achieved by running the GSEM.process.covStruct function. 

GSEM.process.covStruct(c("mood1920.sumstats.gz","misery1930.sumstats.gz", "irrit1940.sumstats.gz","hurt1950.sumstats.gz", "fedup1960.sumstats.gz", "nervous1970.sumstats.gz", "worry1980.sumstats.gz", "tense1990.sumstats.gz", "embarass2000.sumstats.gz", "nerves2010.sumstats.gz", "lonely2020.sumstats.gz", "guilt2030.sumstats.gz"), c(0.451,0.427,0.28,0.556,.406,.237,.568,.171,.478,.213,.177,.283), c(0.451,0.427,0.28,0.556,.406,.237,.568,.171,.478,.213,.177,.283),  "eur_w_ld_chr/", "eur_w_ld_chr/")

Things to know:
Enter the total sample size
Whatever reference file you used should name the reference and alternative allele as A1 and A2, respectively. 
