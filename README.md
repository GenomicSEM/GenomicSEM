# GenomicSEM

R-package which allows the user to fit structural equation models 
based on the summary statistics obtained from genome wide association studies (GWAS). Until explicitly stated otherwise the code on this github in an alpha version (recently updated to version **0.0.2** but still very much in development) and under active development. The code may thus produce undesired results on certain operating systems or when run concurrently with specific packages or R versions. Feel free to raise issues if (or when...) the package produces undesired results, we will attempt to swiftly deal with known issues. Please  **[visit the wiki](https://github.com/MichelNivard/GenomicSEM/wiki)** to get started, or **[check out the preprint](https://www.biorxiv.org/content/early/2018/04/21/305029)**.

**PGC worldwide lab meeting on genomicSEM**

Click below for a video which provides a very clear introduction to the method/paper:

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/ECwQS5UD3YM/0.jpg)](https://www.youtube.com/watch?v=ECwQS5UD3YM?t=3m36s)

**Contents of the wiki:**

Learn how to [install](https://github.com/MichelNivard/GenomicSEM/wiki/1.-Installing-GenomicSEM) `GenomicSEM`

Consider some of the [nuances of summary data, and know where to find summary data.](https://github.com/MichelNivard/GenomicSEM/wiki/2.-Important-resources-and-key-information)

Fit SEM models to GWAS summary data [without a SNP](https://github.com/MichelNivard/GenomicSEM/wiki/3.-Models-without-SNP-effects))

Run a GWAS where the [SNP is included](https://github.com/MichelNivard/GenomicSEM/wiki/4.-Common-Factor-GWAS) in the structural equation model.

Combine [`GenomicSEM` and `OpenMX`](https://github.com/MichelNivard/GenomicSEM/wiki/5.-GenomicSEM-and-OpenMx)

**Instalation:**

We assume you are running R 3.4.1 or newer. We guarantee no backward or forward comparability. If something breaks please raise the issue on GitHub and we will try and fix it ASAP. 

First, you need to install the `devtools` package. You can do this from CRAN, launch R and then type

```[r]
install.packages("devtools")
```
Load the `devtools` package.

```[r]
library(devtools)
```

Now your ready to install the latest version of `GenomicSEM`

```[r]
install_github("MichelNivard/GenomicSEM")
```

That's it! You  are ready to start using `GenomicSEM` 

**License**

 Copyright (C) 2018 Andrew Grotzinger, Michel Nivard & Elliot Tucker-Drob.

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
