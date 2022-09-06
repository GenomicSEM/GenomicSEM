Log of version histories.

**0.0.1** Initial release  
Able to perfrom all analysis in preprint

**0.0.2** Second release.  
Small updates to existing function to make more user friendly and addition of 5 new functions (multiSNP, addGenes, read.fusion, userGWAS, multiGene). We note that some of the results presented on the wiki are slightly inconsistent with what is reported in the biorxiv preprint. This is due to using slightly different QC procedures. Overall, the primary difference is that we use the original LDSC munge function from Bulik-Sullivan et al. (2015) in the biorxiv paper and the output obtained from our own munge function in the wiki (our in-house munge function was still in development). We outline specific differences between these two munge functions in the wiki. The GWIS results are the most discrepant as we use the QC procedures from Nieuwboer et al. (2016) for comparative purposes in the biorxiv paper, but, for illustrative purposes, produce results using the default QC procedures for the Genomic SEM package on the wiki. Results produced from the previous version (0.0.1) of Genomic SEM (i.e. ldsc output) should still be compatible with the new version, but please e-mail us with any new errors or concerns.

**0.0.2b** Second release part b.  
Small updates to second release. Addition of function for parallel processing for userGWAS (userGWASpar). Both userGWAS and userGWASpar also now have their own wiki page to explain their usage. The munge and sumstats functions are updated to produce a .log file to clarify at what stages SNP are being removed in the QC pipeline. The sumstats and addSNPs fucntions are updated to allow for parallel processing. The userGWAS and commonfactorGWAS functions are updated to include additional information in the output (i.e., Z-statistics, p-values for the Q-index, model estimates and model chi-square and degrees of freedom for the Q-index and model chi-square).

**0.0.2c** Fix to ldsc function. March 10th, 2019   
We corrected a bug which impacts the weights used in ld-score regression and caused estimates produced to vary slightly. Substantive results are unlikely to be impacted, but we recommend reinstalling Genomic SEM if it was installed prior to this date.

**0.0.2d** Second release part d.  
Small updates to second release, including: more informative error messages, the addition of the fully standardized column for usermodel, and the ability to estimate models that do not include variables in your ldsc output. This last case can save a lot of time when you want to run subsets of certain models, but do not want to re-run ldsc to exclude a single variable. We also provided updates to any model estimating function (e.g., commonfactorGWAS, usermodel) to allow compatibility with an update to the lavaan package made on July 3rd, 2019. These updates were minimal, and should only affect point estimates and standard errors out to approximately the third or fourth decimal place.

**0.0.2e** Second release part e.  
Updates to second release that merge the addSNPs function, and parallel versions of the multivariate GWAS functions, into two functions: commonfactorGWAS and userGWAS. We combine these functions for two primary reasons. 1. It reduces the memory needs by creating the full S and V matrices one at a time, running the model, and discarding them, as opposed to creating all of them simultaneously. 2. It reduces the number of steps necessary to run a multivariate GWAS.

**0.0.2f** Second release part f. March 9th, 2020  
We discovered and corrected a bug that was introduced to the usermodel function on August 20th, 2019. For some model specifications, the bug had the potential to cause slight changes in point estimates and standard errors due to the sampling covariance (V) matrix not being correctly ordered. This bug did not affect the commonfactor, commonfactorGWAS, or userGWAS functions. All examples on the github wiki, and all examples that we have personally worked with, were not affected by this bug. In most cases, your results are likely to be unchanged. We recommend reinstalling to the most recent update of GenomicSEM, particularly if it was installed or updated between August 20th, 2019 and March 9, 2020.

**0.0.2g** Second release part g. March 27th, 2020  
A number of small code updates were made to make the package more user friendly. This included printing a .log file for the ldsc function, printing a warning if only one trait is provided to the ldsc function, additional checks for the sumstats function to ensure columns were correctly interpreted (e.g., ensuring the average ratio of estimate over standard error post-transformation is within a reasonable range), printing a note for the usermodel function that model fits are printed as NA for saturated models (i.e., df = 0), and printing a warning if mathematical operations (e.g., a '-' symbol) are used to name variables as this cannot be parsed by lavaan. An update to the commonfactorGWAS function ensures that the code will not raise an error and stop running due to models that fail to converge. Finally, the default SE of the SNP variance used in the sampling covariance matrix (V) has been raised to .0005, as opposed to the previous value of 1e-8. This value can still be toggled based on user preferences using the SNPSE argument. We found that the previous level of 1e-8 would often cause lavaan to raise benign warnings that took up unnecessary memory in the context of the GWAS functions. In addition, .0005 seems to be a more conservative and appropriate estimate of the SNP SE based on the European subset of the 1000 Genomes reference date. Our observation is that this update is highly unlikely to change Z-statistics produced from GWAS above the 4th decimal place, that the rank ordering of solutions before and after the update is maintained, and that qualitative conclusions remain the same.

**0.0.2h** Second release part h. May 15th, 2020  
We included a number of small code updates. This includes the ldsc function now printing a .log file and the ability to specify how this log file is named. We also now include the argument 'stand' to the ldsc function that allows the user to specify whether they would like to additionally output the genetic correlation matrix and its sampling covariance matrix. The message "SE could not be computed" has now been removed from the GWAS functions. This message was previously printed when lavaan printed NA for the standard error of one or more parameters, despite converging on a solution. Although Genomic SEM uses sandwich corrected standard errors calculated outside of lavaan, it was initially unclear whether this was diagnostic of problems that might affect matrices output from lavaan used for calculation of the sandwich corrected standard errors. Further testing revealed that an NA was often printed for SEs for parameters estimated very near 0, but that the matrices used for sandwich correction were unaffected.

Within the GWAS functions (userGWAS, commonfactorGWAS) we now include a GC argument that allows the user to specify the form of genomic control (GC) used to correct univariate GWAS standard errors (SEs). The new package default (GC = "stand") is to adjust the SEs by taking the product of the SEs and the square root of the univariate LDSC intercepts. The prior package default was to adjust the SEs by taking the product of the SEs and the univariate LDSC intercept, which can still be used by specifying GC = "conserv". The user also has the option to use uncorrected SEs by specifying GC = "none". A more detailed consideration of genomic control options within Genomic SEM can be found here: https://gist.github.com/MichelNivard/fcb22ab7401d7e9af95a5e1fd46ad5b7.

We also now offer the option to run the GWAS functions using MPI (i.e., multi-node processing). We note that MPI will not be useful for more simple models, as the computational overhead required to run MPI outweighs the run time benefits. For example, for the five-indicator p-factor model the model run for each SNP takes ~0.3 seconds, and in this case using MPI will generally slow down the run time. However, for more complicated models that may run on the order of 5-10 seconds per SNP, MPI may greatly speed up results and in some cases may help with memory issues. We note that using the MPI functionality requires that Rmpi already be installed on the computing cluster being used.

**0.0.3** Third release. September 23rd, 2020  
We introduced the third major update to Genomic SEM that includes the Stratified Genomic SEM extension. Stratified Genomic SEM can be used to examine enrichment of any parameter for a model estimated in Genomic SEM. This includes examining multivariate enrichment across traits by examining enrichment of factor variances. A separate help page is provided for Stratified Genomic SEM, which includes running the new enrich and s_ldsc functions.

**0.0.3b** Third release part b. March 8th, 2021  
A bug was fixed that would provide a downwardly biased estimate of CFI when using equality constraints across factor loadings for the usermodel function. In addition, a bug was fixed in the usermodel function that would remove any parameter estimates that were freely estimated at 0 (an unlikely situation in real data scenarios but one that arose when running simulations).

**0.0.3c** Third release part c. April 12th, 2021  
The usermodel, commonfactor, and userGWAS functions were updated with respect to how model chi-square is calculated. In the prior version, the functions would estimate a follow-up "residual model" that specified a fully saturated model reflecting the covariances and residual variances among the indicators after fixing the estimated effects from the specified. This produced the residual sampling covariance matrix necessary to estimate model chi-square. We have found that this residual sampling covariance matrix in fact reflects the observed sampling covariance matrix, such that an equivalent estimate of model chi-square can be obtained using the observed sampling covariance matrix paired with the residual genetic covariance matrix, calculated as the difference between the model implied and observed genetic covariance matrix. Through comparisons via simulation we confirmed that this produces the same model chi-square estimate and is also chi-square distributed. Using this updated version of calculating model chi-square has the advantage of not having to estimate the follow-up, residual model, which in the context of the userGWAS function can reduce run times by ~50%. In addition, these follow-up residual models may fail to converge or converge local estimates, whereas the new version of estimating model chi-square is guaranteed to produce a solution as long as the specified model converges. This update also includes a correction for ML estimation as it has come to our attention that lavaan employs a small sample size correction to the genetic covariance matrix based on the user provided sample size. We have found that this correction results in very slight deviations in model estimates based on the input sample size. The updated code has turned off this correction for ML estimation by setting sample.cov.rescale to FALSE in all cases.

**0.0.3d** Third release part d. April 30th, 2021  
The read_fusion, addGenes, and commonfactorGWAS function were updated and documented to be able to run multivariate TWAS, which we refer to as Transcriptome-wide Structural Equation Modeling (T-SEM). T-SEM takes univariate output from FUSION and combines it with the output from multivariable ld-score regression within the general Genomic SEM framework to allow for examining the effect of inferred gene expression within a multivariate framework. This includes the ability to examine the effect of gene expression on a common factor.

**0.0.3e** Third release part e. September 23rd, 2021  
The munge function and wiki were updated to reflect our recommendation to use the sum of effective sample sizes and a sample prevalence of 0.5 to produce an unbiased estimate of liability scale conversion. This is based on simulation results and derivations described in our recent preprint that can be found here: https://www.medrxiv.org/content/10.1101/2021.09.22.21263909v1

**0.0.3f** Third release part f. November 2nd, 2021    
Three updates of note were made to the package. The first was, per version 0.0.3e, to use the sum of effective sample size in the sumstats function when backing out logistic betas from Z-statistics for binary traits. The second was an additional update to the sumstats function to use MAF from the summary statistics file for OLS traits or Z-statistics only binary traits as MAF is used to back out a partially standardized beta and mismatch between the reference panel MAF and in-sample MAF may cause a minimal amount of bias. If in-sample MAF is not available, sumstats also now includes an additional 'betas' argument, which can be used to directly supply the column name of the betas in the summary stats if the continuous phenotype was known to be standardized prior to running the GWAS (thereby circumventing the need to use the reference panel MAF as the function will not perform any additional transformations). The third update was to ldsc and reflects a new package default to set the number of blocks used to produce block jackknife standard errors to 1 more than the number of non-redundant elements in the genetic covariance matrix when there are more than 18 traits being analyzed. This reflects the fact that ldsc typically uses 200 blocks, which yields a sampling covariance (V) matrix that produces inaccurate estimates of model chi-square when the number of independent elements in the genetic covariance matrix is > 200. This update is unlikely to qualitatively change any model results beyond model chi-square, but results for Q_SNP with > 18 traits prior to this updated will likely have produced hugely inflated estimates (e.g., half of the genome is identified as significant for Q_SNP) and results in this case should be re-examined.

**0.0.4** Fourth release. December 16th, 2021  
This update contains minimal changes for users (see Feature updates), but large changes for developers (see Code updates).

**Feature updates**  
- Functions where sanity checks have been added will now produce an error when provided wrong input before proceeding to the analyses. This prevents later termination as a result of incorrect input.
- The input to the `files` argument for `munge()` and `sumstats()` is changed to a vector, where previously both a list or vector were accepted.
  For now if a list is passed it is automatically changed to a vector and a deprecationwarning is printed, as this functionality will be removed in a future version.
- Parallel on Windows enabled for `userGWAS()`, `commonfactorGWAS()` and `sumstats()`
- Added parallel functionality for all OSes to `munge()` (use of `parallel` and `cores` arguments identical to other functions)

**v0.0.4 Code updates**
- Added PATCHNOTES.md to keep track of changes in each patch as replacement for Wiki version history
- NAMESPACE was changed to enable the use of private functions (not directly accessible to the end-user), these functions will always start with a period (.).
- Global code changes:
  - Started code restructuring toward a more function-based paradigm. This is a large change for developers, but users (should) not notice.
  - Large functions will gradually be chunked into smaller functions which can be more easily modified or replaced with more efficient versions. This will more easily enable use of a different backend (C++, Python, etc) in future versions.
  - Moving from `if (X == TRUE)` and `if (X == FALSE)` to `if (X)` and `if (!X)` for readability.
  - Minor changes to reduce code length where this did not impact readability or performance.
  - Numerous functions (such as LOG, tryCatch.W.E) have been changed to private, and are moved to utils.R to enable access from all functions in the package.
  - `cat(print(paste(...)))` are replaced by `.LOG` for readability.
  - Changed assignments with = to <- script-wide for clarity.
  - Added sanity checks to ensure user input is correct to prevent functions terminating halfway through analysis *in progress*
- Creation of matrices that are used across different GSEM functions such as Z_pre V_SNP and V_full have also been moved to utils.R. 
- Moved recognition of column names and its related checks and warnings to utils (as it is used in sumstats and munge). 
- `userGWAS()`, `sumstats()`, `commonfactorGWAS()` and `munge()` now have 1 main function each (*_main) which is used for both parallel and serial operation, to prevent copy-pasting.
  - `.userGWAS_main()` is the full analysis per SNP
  - `.commonfactorGWAS_main()` is the full analysis per SNP
  - `.sumstats_main()` is the entire parsing per file
  - `.munge_main()` is the entire parsing per file
- in `commonfactorGWAS()` specifically:
  - Moved stopCluster to `on.exit` straight after creation to ensure cluster closes when the script terminates early for whatever reason.
- in `userGWAS()` specifically: 
  - Some local variable names have been changed (such as k to n_phenotypes) to improve readability. 
  - Removed argument modelchi as it was unused.
  - Moved stopCluster to `on.exit` straight after creation to ensure cluster closes when the script terminates early for whatever reason.
- in `munge()` specifically:
  - added column.names argument allowing users to manually provide column names (e.g. column.names=list(SNP=mysnpcolumn))
  - added overwrite argument passed to gzip() to prevent halting on already existing files
  - Added arguments `parallel`, and `cores` to `munge()`
  - Added arguments `parallel`, and `cores` to `munge()` documentation
  - Enabled parallel munging of sumstats, some initial test results below
  - Parallel munging will store log files separately to prevent clutter
  - Removed N_Cases and N_Controls from recognized columns
- in `sumstats()` specificallly:
  - Changed parallel function used from `mclapply` to `foreach()` with PSOCK cluster to (a) be in line with other parallel functions, and (b) enable parallel functionality on Windows.
  - Changed recognizing column names such that for `munge()` any beta or logOR column is preferred over Z, for `sumstats` if linprob and OLS is set to FALSE, Z will never be interpreted as effect
  - Removed N_Cases and N_Controls from recognized columns
  - Changed warning for missing SE column to early breaking error in cases where it would cause rest of the code to fail.
- in `usermodel()` specifically:
  - Added breaking error when no LDSC trait name matches names in the model. The function would otherwise break with error `incorrect number of dimensions`
    
**Parallel munge test**  
**HomePC**: Windows 10, Ryzen7 3700X @ 3.60-4.4GHz, 48GB RAM, 970 EvoPlus 1TB  
**Server**: Linux Ubuntu, 2xEPYC 7H12 @ 2.25-3.3GHz, 512GB RAM, storage device unknown  
**Sumstats**: approximately 4M SNPs * 17 columns per file (12 neuroticism items from UKB)  

|N(files)|parallel|cores | HomePC<br>runtime (s) | Server<br>runtime (s) |
|--------|--------|------|-----------------------|-----------------------|
| 4      |  FALSE |  1   |    313                |   207                 |
| 4      |  TRUE  |  2   |    204                |   195                 |
| 4      |  TRUE  |  4   |    158                |   127                 |

|N(files)|parallel|cores | HomePC<br>runtime (s) | Server<br>runtime (s) |
|--------|--------|------|-----------------------|-----------------------|
| 8      |  FALSE |  1   |    454                |   401                 |
| 8      |  TRUE  |  2   |    354                |   330                 |
| 8      |  TRUE  |  4   |    245                |   199                 |
| 8      |  TRUE  |  8   |    206                |   135                 |

|N(files)|parallel|cores | HomePC<br>runtime (s) | Server<br>runtime (s) |
|--------|--------|------|-----------------------|-----------------------|
| 12     |  FALSE |  1   |   755                 |   772                 |
| 12     |  TRUE  |  2   |   505                 |   472                 |
| 12     |  TRUE  |  4   |   332                 |   273                 |
| 12     |  TRUE  |  8   |   282                 |   209                 |
| 12     |  TRUE  |  12  |   256                 |   147                 |

**Code update 04.02.2022**
- Fixed issue that would couse an early exit if SNPSE was a numeric
- Fixed issue that forced an SE column when linprob is set to TRUE
- Fixed issue that caused lack of column naames in the result of userGWAS 

**Code update 07.02.2022**
- Fixed 'unexpected end of input' error when installing package
- Fixed issue that caused parallel=TRUE in userGWAS to not return headers.

**0.0.5** Fifth release. February 23th, 2022
- Fixed an issue that may cause parallel clusters to not register and close properly.
- Changes: userGWAS and commonfactorGWAS:
  Now, before starting the main analysis (`*_main()`) both userGWAS and commonfcatorGWAS do a partial run (everything up to and including `lavaan::sem()`) of the first SNP
  The resulting object from this first run is passed to the subsequent per-SNP analyses, and from the object of the first run the following is pased to `lavaan::lavaan()` (the function underneath `lavaan::sem()`):
   - `basemodel@Options`: Many fit- and modeloptions like estimator, tolerances, etc.
   - `basemodel@ParTable`: Pre-formatted empty result table and some other paramters
   - `basemodel@Data`: Not the actual covariance matrix (unlike the name would suggest), but statistics about the data, such as N and 'moment' indicating a covariance matrix is used
   - `basemodel@Model`: The Lavaan-formatted model, instead of the model string
- This prevents lavaan from having to re-format all of this identically for each SNP (as all of the above is equal across SNPs). Significantly saving compute time.  
System: Windows 10, Ryzen7 3700X @ 3.60-4.4GHz, 48GB RAM, 970 EvoPlus 1TB  
UserGWAS of 100K SNPs from 12 summary statistics  
  
| parallel | cores | v0.0.4 runtime (s) | v0.0.5 runtime (s) |
|----------|-------|:------------------:|:------------------:|
| FALSE    | 1     |       27,302       |       25,979       |
| TRUE     | 2     |       11,659       |       10,115       |
| TRUE     | 4     |       6,199        |       5,505        |
| TRUE     | 8     |       4,123        |       3,478        |
| TRUE     | 12    |       3,549        |       2,863        |

- It seems likely that the greater saves in higher core counts are due to less memory IO (due to fewer operations) but I don't know how to test this definitively. Lower RAM usage in v0.0.5 does support this theory.

| parallel | cores | v0.0.4 max RAM (MB) | v0.0.5 max RAM (MB) |
|----------|-------|:-------------------:|:-------------------:|
| FALSE    | 1     |        1,407        |        1,398        |
| TRUE     | 2     |        2,832        |        2,447        |
| TRUE     | 4     |        3,497        |        2,837        |
| TRUE     | 8     |        4,764        |        3,615        |
| TRUE     | 12    |        6,103        |        4,680        |

Note: max RAM values are obtained running solely the 100K-userGWAS in an isolated test-environment, without any RStudio/IDE overhead or any other data. RAM usage in practical application is likely to be higher.

- Added **Parallel performance on Linux** to `README.md`.
- Testing shows that on Linux performance is best when OPENBLAS_NUM_THREADS is set to 1, and parallel is purely done through `cores` in GenomicSEM. Note the specific argument may vary depending on OS and R build.

System: Ubuntu 20.04, 2x EPYC 7H12 @ 2.6-3.3GHz, 2x256GB RAM
UserGWAS of 100K SNPs from 12 summary statistics, GenomicSEM v0.0.5

| cores | OPENBLAS unlimited | OPENBLAS limited |
|-------|:------------------:|:----------------:|
| 1     |       12,645       |      10,545      |
| 2     |       7,577        |      4,489       |
| 4     |       4,214        |      2,559       |
| 8     |       5,347        |      1,186       |
| 12    |       6,585        |       793        |
| 24    |       5,170        |       458        |

- Runtime improvements on Linux are largely comparable, and slightly more consistent than on Windows. Note the following results describe a **full** (1.8M SNPs) GWAS, not the 100K SNPs used in many previous tests.

System: Ubuntu 20.04, 2x EPYC 7H12 @ 2.6-3.3GHz, 2x256GB RAM  
userGWAS of **1.8M** SNPs from 12 summary statistics  

| cores | v0.0.4 runtime (s) | v0.0.5 runtime (s) |
|-------|:------------------:|:------------------:|
| 18    |       14,324       |       11,988       |
| 24    |       10,916       |       9,157        |
| 30    |       8,820        |       7,412        |
| 60    |       5,147        |       4,362        |
| 120   |       3,708        |       3,304        |

**Code update 21.03.2022**
- Fixed issue that caused CHR/BP/MAF merges to fail in userGWAS and commonfactorGWAS

**Code update 22.03.2022** `labeled 0.0.5b`
- Fixed issue that caused `infinite or missing values in x` errors when N(snps) < N(cores)

**Code update 24.03.2022** `labeled 0.0.5c`
- Fixed SNP order issue in commonfactorGWAS and userGWAS 