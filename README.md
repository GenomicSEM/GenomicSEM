# GenomicSEM

R-package which allows the user to fit structural equation models 
based on the summary statistics obtained from genome wide association studies (GWAS). Until explicitly stated otherwise the code on this github is an alpha version (now on version **0.0.2**) and under active development. The code may thus produce undesired results on certain operating systems or when run concurrently with specific packages or R versions. Feel free to raise issues if (or when...) the package produces undesired results, we will attempt to swiftly deal with known issues. Please  **[visit the wiki](https://github.com/MichelNivard/GenomicSEM/wiki)** to get started, or **[check out the paper](https://www.nature.com/articles/s41562-019-0566-x)**. If you are having issues and not finding the answers anywhere on the wiki or FAQs page, we encourage you to post your question on the **[google group](https://groups.google.com/forum/#!forum/genomic-sem-users)**.

**Feature update**: GenomicSEM can now run [HDL](https://t.co/OBHihTb7rE?amp=1) a novel method for estimating heritability and genetic correlation that can in some cases outperform LDSC. See out tutorial [HERE](https://rpubs.com/MichelNivard/640145)

**Code Update**: The most recent code update is dated May 15, 2020. We recommend reinstalling to the most recent update of GenomicSEM. However, please note that changes in Genomic SEM defaults are likely to produce slight changes in results relative to previous versions. For further details, see [version history](https://github.com/MichelNivard/GenomicSEM/wiki/Version-History)

**PGC worldwide lab meeting on genomicSEM**

Click below for a video which provides a very clear introduction to the method/paper:

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/ECwQS5UD3YM/0.jpg)](https://www.youtube.com/watch?v=ECwQS5UD3YM?t=3m36s)

**Contents of the wiki:**

Learn how to [install](https://github.com/MichelNivard/GenomicSEM/wiki/1.-Installing-GenomicSEM) `GenomicSEM`

Consider some of the [nuances of summary data, and know where to find summary data.](https://github.com/MichelNivard/GenomicSEM/wiki/2.-Important-resources-and-key-information)

Fit SEM models to GWAS summary data [without a SNP](https://github.com/MichelNivard/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects)

Run a GWAS where the [SNP is included](https://github.com/MichelNivard/GenomicSEM/wiki/4.-Common-Factor-GWAS) in the structural equation model.

Combine [`GenomicSEM` and `OpenMX`](https://github.com/MichelNivard/GenomicSEM/wiki/7.-GenomicSEM-and-OpenMx)

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
