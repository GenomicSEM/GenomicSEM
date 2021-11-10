# GenomicSEM
**Code Update 10.11.2021**

- Started code restructuring toward a more function-based paradigm. This is a large change for developers, but users (should) not notice.
- Enabled use of parallel on Windows in userGWAS.
- Global minor changes to reduce code length where this did not impact readability or performance.
- Large functions will gradually be chunked into smaller functions which can be more easily modified or replaced with more efficient versions.
- This will more easily enable use of a different backend (C++, Python, etc) in future versions.
- NAMESPACE was changed to enable the use of private functions (not accessible to the end-user), these functions will always start with a period (.).
- Numerous functions (such as LOG, tryCatch.W.E) have been changed to private, and are moved to utils.R to enable access across the whole package. 
- Creation of matrices that are used across different GSEM functions such as Z_pre V_SNP and V_full have also been moved to utils.R. 
- Additionally started replaing cat(print(paste(...))) commands with .LOG to unify its use. 
  Changed assignments with = to <- script-wide for clarity.
- Moved recognition of column names and its related checks and warnings to utils (as it is used in sumstats and munge). 
- userGWAS and sumstats now have 1 main function each (userGWAS_main and sumstats_main) which is used for both parallel and serial operation, to prevent copy-pasting.
  - userGWAS_main is the full analysis per SNP
  - sumstats_main is the entire parsing per file
- Enabled user-provided column names in munge.
- in userGWAS specifically: 
    - Some variable names have been changed (such as k to n_phenotypes) to improve readability. 
    - Removed argument modelchi as it was unused.


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
