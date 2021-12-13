# GenomicSEM
**TO DO**
- Unify use of sanity checks across all user-accessible functions

**Changes for users**
- Functions where sanity checks have been added will now produce an error when provided wrong input before proceeding to the analyses. 
  This prevents later termination as a result of incorrect input.
- The input to the `files` argument for `munge()` and `sumstats()` is changed to a vector, where previously both a list or vector were accepted.
  For now if a list is passed it is automatically changed to a vector and a deprecationwarning is printed, as this functionality will be removed in a future version.
- Parallel on Windows enabled for `userGWAS()`
- Parallel on Windows enabled for `commonfactorGWAS()`
- Parallel on Windows enabled for `sumstats()`
- Added parallel functionality for all OSes to `munge()` (use of `parallel` and `cores` arguments identical to other functions)


**Code Update 13.12.2021**
- Fixed issue caused by typo in `.sumstats_main()`
- Changed warning for missing SE column in `sumstats()` to early breaking error as it would cause rest of the code to fail.
- Fixed issue caused by typo in `commonfactorGWAS()`
- Created commonfactorGWAS_main.R
- Moved main body of analysis code of `commonfactorGWAS()` to `.commonfactorGWAS_main()` in line with other functions
- Moved stopCluster in commonfactorGWAS to `on.exit` straight after creation to ensure cluster closes when the script terminates early for whatever reason.
- Enabled parallel for Windows for `commonfactorGWAS()` using PSOCK cluster, same as `userGWAS`.
- Global change: Moving from `if (X == TRUE)` and `if (X == FALSE)` to `if (X)` and `if (!X)` for readability.

**Code Update 10.12.2021**
- Fixed `files2 not found` issue in sumstats
- Changed recognizing column names such that for `munge()` any beta or logOR column is preferred over Z, for `sumstats` if linprob and OLS is set to FALSE, Z will never be interpreted as effect
- Removed N_Cases and N_Controls from recognized columns
- Added warn_for_missing argument to `get_renamed_colnames()` to modify which missing columns should produce a warning

**Code Update 09.12.2021**
- Fixed 'N should be numeric' issue when NULL or NA is passed.
- Fixed issue in `.check_one_of()`.
- Changed parallel function used in sumstats from `mclapply` to `foreach()` with PSOCK cluster to (a) be in line with other parallel functions, and (b) enable parallel functionality on Windows.
- Changed `.sumstats_main()` input to single file and single values to (a) be in line with other main functions and (b) ensure the loop only works with the required data and nothing extra
- Enforced vector as input for `files` argument in both `munge()` and `sumstats()`
- Enabled parallel for `sumstats()` on Windows
- Moved stopCluster in userGWAS to `on.exit` straight after creation to ensure cluster closes when the script terminates early for whatever reason.


**Code Update 08.12.2021**
- Added 'Changes for users' to README.md to keep track of changes users (may) notice
- Fixed `NA_provided not found` error
- Changed `sumstats()` to be compliant with changes made to munge
- Fixed `,,` issue as a result of transitioning to `.LOG()` that caused `ldsc()` to crash
- Moved body of sumtats munging code to `.munge_main()`
- Added arguments `parallel`, and `cores` to `munge()`
- Added arguments `parallel`, and `cores` to `munge()` documentation
- Enabled parallel munging of sumstats, some initial test results below
- Parallel munging will store log files separately to prevent clutter

Systems:  
HomePC: Windows 10, Ryzen7 3700X @ 3.60-4.4GHz, 48GB RAM, 970 EvoPlus 1TB  
Server: Linux Ubuntu, 2xEPYC 7H12 @ 2.25-3.3GHz, 512GB RAM, storage device unknown  
Sumstats: approximately 4M SNPs * 17 columns per file (12 neuroticism items from UKB)

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

  
**Code Update 07.12.2021**
- Fixed issue that caused multiple columns to be found in munge.R
- Merged changes (MAF=0 check) from main branch in sumstats.R
- Added basic sanity check functions to utils_sanitychecks.R
- Added basic sanity checks to `munge()` and `userGWAS()` inputs
- Sanity checks will be added to other functions in the future, to prevent functions terminating halfway through due to wrong input.
- Added missing info.filter to .sumstats_main()
- Unified use of .LOG across all scripts
- Added sanity checks to `commonfactorGWAS()`

**Code Update 06.12.2021**
- Fixed `object 'filenames' not found”` error when munging summary statistics

**Code Update 17.11.2021**
- Fixed '}' issue in utils.R
- Changed assignments with = to <- in munge.R

**Code Update 10.11.2021**
- Started code restructuring toward a more function-based paradigm. This is a large change for developers, but users (should) not notice.
- Large functions will gradually be chunked into smaller functions which can be more easily modified or replaced with more efficient versions. This will more easily enable use of a different backend (C++, Python, etc) in future versions.
- Enabled use of parallel on Windows in userGWAS.
- Global minor changes to reduce code length where this did not impact readability or performance.
- NAMESPACE was changed to enable the use of private functions (not accessible to the end-user), these functions will always start with a period (.).
- Numerous functions (such as LOG, tryCatch.W.E) have been changed to private, and are moved to utils.R to enable access from all functions in the package. 
- Creation of matrices that are used across different GSEM functions such as Z_pre V_SNP and V_full have also been moved to utils.R. 
- Additionally started replaing cat(print(paste(...))) commands with .LOG to unify its use. 
- Changed assignments with = to <- script-wide for clarity.
- Moved recognition of column names and its related checks and warnings to utils (as it is used in sumstats and munge). 
- userGWAS and sumstats now have 1 main function each (userGWAS_main and sumstats_main) which is used for both parallel and serial operation, to prevent copy-pasting.
  - userGWAS_main is the full analysis per SNP
  - sumstats_main is the entire parsing per file
- in userGWAS specifically: 
    - Some local variable names have been changed (such as k to n_phenotypes) to improve readability. 
    - Removed argument modelchi as it was unused.
- in munge specifically:
    - added column.names argument allowing users to manually provide column names (e.g. column.names=list(SNP=mysnpcolumn))
    - added overwrite argument passed to gzip() to prevent halting on already existing files


R-package which allows the user to fit structural equation models 
based on the summary statistics obtained from genome wide association studies (GWAS). Until explicitly stated otherwise the code on this github is an alpha version (now on version **0.0.3d**) and under active development. The code may thus produce undesired results on certain operating systems or when run concurrently with specific packages or R versions. Feel free to raise issues if (or when...) the package produces undesired results, we will attempt to swiftly deal with known issues. Please  **[visit the wiki](https://github.com/MichelNivard/GenomicSEM/wiki)** to get started, or **[check out the paper](https://www.nature.com/articles/s41562-019-0566-x)**. If you are having issues and not finding the answers anywhere on the wiki or FAQs page, we encourage you to post your question on the **[google group](https://groups.google.com/forum/#!forum/genomic-sem-users)**.

**Feature update**: GenomicSEM can now run [HDL](https://t.co/OBHihTb7rE?amp=1) a novel method for estimating heritability and genetic correlation that can in some cases outperform LDSC. See out tutorial [HERE](https://rpubs.com/MichelNivard/640145)

**Code Update**: The most recent code update is dated April 12th, 2021. This code update includes a statistically equivalent, but far more efficient, way of calculating model chi-square. In the context of a userGWAS model where model chi-square was requested this may decrease run-times by ~50%. We recommend reinstalling to the most recent update of GenomicSEM. However, please note that changes in Genomic SEM defaults are likely to produce slight changes in results relative to previous versions. For further details, see [version history](https://github.com/MichelNivard/GenomicSEM/wiki/Version-History)

**PGC worldwide lab meeting on genomicSEM**

Click below for a video which provides a very clear introduction to the method/paper:

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/ECwQS5UD3YM/0.jpg)](https://www.youtube.com/watch?v=ECwQS5UD3YM?t=3m36s)

**Contents of the wiki:**

Learn how to [install](https://github.com/MichelNivard/GenomicSEM/wiki/1.-Installing-GenomicSEM) `GenomicSEM`

Consider some of the [nuances of summary data, and know where to find summary data.](https://github.com/MichelNivard/GenomicSEM/wiki/2.-Important-resources-and-key-information)

Fit SEM models to GWAS summary data [without a SNP](https://github.com/MichelNivard/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects)

Run a GWAS where the [SNP is included](https://github.com/MichelNivard/GenomicSEM/wiki/4.-Common-Factor-GWAS) in the structural equation model.

Estimate [functional enrichment](https://github.com/GenomicSEM/GenomicSEM/wiki/6.-Stratified-Genomic-SEM) for any parameter in a Genomic SEM model (e.g., factor variances). 

Run multivariate TWAS using [`T-SEM`](https://github.com/GenomicSEM/GenomicSEM/wiki/7.-Transcriptome-wide-SEM-(T-SEM))

**Installation:**

We assume you are running R 3.4.1 or newer. We guarantee no backward or forward comparability. If something breaks please raise the issue on GitHub and we will try and fix it ASAP. 

First, you need to install the `devtools` package. You can do this from CRAN, launch R and then type

```[r]
install.packages("devtools")
```
Load the `devtools` package.

```[r]
library(devtools)
```

Now you are ready to install the latest version of `GenomicSEM`. Note that this will often raise 24 warnings about replacing previous imports; these warnings are safe to ignore.

```[r]
install_github("GenomicSEM/GenomicSEM")
```

That's it! You  are ready to start using `GenomicSEM` 

**License**

Copyright (C) 2018 Andrew Grotzinger, Mijke Rhemtulla, Hill F. Ip, Michel Nivard, & Elliot Tucker-Drob

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
