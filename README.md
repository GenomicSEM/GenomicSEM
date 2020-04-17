# GenomicSEM

R-package which allows the user to fit structural equation models 
based on the summary statistics obtained from genome wide association studies (GWAS). Until explicitly stated otherwise the code on this github is an alpha version (now on version **0.0.2**) and under active development. The code may thus produce undesired results on certain operating systems or when run concurrently with specific packages or R versions. Feel free to raise issues if (or when...) the package produces undesired results, we will attempt to swiftly deal with known issues. Please  **[visit the wiki](https://github.com/MichelNivard/GenomicSEM/wiki)** to get started, or **[check out the paper](https://www.nature.com/articles/s41562-019-0566-x)**. If you are having issues and not finding the answers anywhere on the wiki or FAQs page, we encourage you to post your question on the **[google group](https://groups.google.com/forum/#!forum/genomic-sem-users)**.

**Code Update**: On December 4th, 2019 we combined the addSNPs and multivariate GWAS functions and their parallel counterparts into a single function to reduce memory load and the number of steps in the analytic pipeline. Previous pipelines using addSNPs output can still be used, but the user will need to be sure to specify the correct arguments for the GWAS functions in the subsequent step. 

**CODE UPDATE 2** On January 8th, 2020 we corrected a bug in `sumstats`, which prepares summary statistics for GWAS, Prior to this fix, a small proportion of SNPs were omitted erroneously. In some GWAS of rigourously cleaned SNPs we found no SNPs affected at all, but from all SNP (i.e. analyses that retain low MAF and low INFO SNPs) more SNPs are affected (1-5%). Specifically the code removed SNPs with duplicate basepair positions, but did not consider chromosome when doing so. This bug did not affect Genomic SEM models without SNPs, and did not affect estimates for SNPs that were included in GWAS. The bug simply resulted in a small proportion of SNPs being arbitrarily removed from consideration.

**CODE UPDATE 3** On March 9th, 2020 we discovered and corrected a bug that was introduced to the usermodel function on August 20th, 2019. For some model specifications, the bug had the potential to cause slight changes in point estimates and standard errors due to the sampling covariance (V) matrix not being correctly ordered. This bug did not affect the commonfactor, commonfactorGWAS, or userGWAS functions. All examples on the github wiki, and all examples that we have personally worked with, were not affected by this bug. In most cases, your results are likely to be unchanged. We recommend reinstalling to the most recent update of GenomicSEM, particularly if it was installed or updated between August 20th, 2019 and March 9, 2020.

**PGC worldwide lab meeting on genomicSEM**

Click below for a video which provides a very clear introduction to the method/paper:

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/ECwQS5UD3YM/0.jpg)](https://www.youtube.com/watch?v=ECwQS5UD3YM?t=3m36s)

**Contents of the wiki:**

Learn how to [install](https://github.com/MichelNivard/GenomicSEM/wiki/1.-Installing-GenomicSEM) `GenomicSEM`

Consider some of the [nuances of summary data, and know where to find summary data.](https://github.com/MichelNivard/GenomicSEM/wiki/2.-Important-resources-and-key-information)

Fit SEM models to GWAS summary data [without a SNP](https://github.com/MichelNivard/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects).

Run a GWAS where the [SNP is included](https://github.com/MichelNivard/GenomicSEM/wiki/4.-Common-Factor-GWAS) in the structural equation model.

Combine [`GenomicSEM` and `OpenMX`](https://github.com/MichelNivard/GenomicSEM/wiki/6.-GenomicSEM-and-OpenMx).

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
install_github("MichelNivard/GenomicSEM")
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
