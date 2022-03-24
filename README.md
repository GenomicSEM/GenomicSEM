# GenomicSEM

R-package which allows the user to fit structural equation models 
based on the summary statistics obtained from genome wide association studies (GWAS). Until explicitly stated otherwise the code on this github is an alpha version (now on version **0.0.5c**) and under active development. The code may thus produce undesired results on certain operating systems or when run concurrently with specific packages or R versions. Feel free to raise issues if (or when...) the package produces undesired results, we will attempt to swiftly deal with known issues. Please  **[visit the wiki](https://github.com/MichelNivard/GenomicSEM/wiki)** to get started, or **[check out the paper](https://www.nature.com/articles/s41562-019-0566-x)**. If you are having issues and not finding the answers anywhere on the wiki or FAQs page, we encourage you to post your question on the **[google group](https://groups.google.com/forum/#!forum/genomic-sem-users)**.

**For Linux users** Please read [1.1 For Linux users in the wiki](https://github.com/GenomicSEM/GenomicSEM/wiki/1.1-For-Linux-users).
**v0.0.5c Patch notes**: This update contains several bug fixes.
**v0.0.5 Patch notes**: This update contains some bug fixes as well as optimizations making `userGWAS()` and `commonfactorGWAS()` 5-20% faster while using slightly less RAM, greater improvement at higher core counts in parallel (see [patchnotes](PATCHNOTES.md)).  
 

**v0.0.4 Patch notes**: This update is the start of a large code restructure and contains minimal changes for users (see below), but large changes for developers (see [patchnotes](PATCHNOTES.md)).  
**v0.0.4 Feature update**: `userGWAS()`, `commonfactorGWAS()` and `sumstats()` can now be run in parallel on Windows systems, additionally `munge()` can now be run in parallel as well (on both Linux/Mac and Windows). A number of functions now contain checks on input (ideally) preventing the function from stopping halway through analysis due to incorrect input. Input to `files` argument for `munge()` and `sumstats()` should now be a vector (a list will still work for now for backwards compatibility).

**v0.0.3d Feature update**: GenomicSEM can now run [HDL](https://t.co/OBHihTb7rE?amp=1) a novel method for estimating heritability and genetic correlation that can in some cases outperform LDSC. See out tutorial [HERE](https://rpubs.com/MichelNivard/640145)  
**v0.0.3d Code update**: The most recent code update is dated April 12th, 2021. This code update includes a statistically equivalent, but far more efficient, way of calculating model chi-square. In the context of a userGWAS model where model chi-square was requested this may decrease run-times by ~50%. We recommend reinstalling to the most recent update of GenomicSEM. However, please note that changes in Genomic SEM defaults are likely to produce slight changes in results relative to previous versions. For further details, see [patchnotes](PATCHNOTES.md)

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

**Parallel performance on Linux**

In some instances of R on Linux a parallel backend is automatically configured which by default uses the maximum number of cores (i.e. 1 thread per core).
It is unclear to me if this is base R or some other package that's doing this in the background. 
In GenomicSEM `parallel` and `foreach` are used to manage parallel operation which seem to work outside this scope. 
Combined, this causes the creation of far too many threads, precisely: `cores` argument * number of actual CPU cores. 
For example, a 16-core machine, with `cores=15` would spawn `16*15=240` R threads. 
This in turn causes CPU congestion, and causes a very significant drop in performance, especially in high core-count machines. In testing on a 256-core machine, 100K SNP userGWAS took 1,5 hours without limiting this backend, compared to under 10 minutes with limits. 

How to change this may depend on your Linux and/or R build, but as a catch-all you can use the following prior to running R:
```
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```
Note that this may change behavior of other programs or packages as well, so it is recommended to do this in a separate session, or from a separate script, e.g. create `RunGSEMAnalyses.sh`:
```
#!/bin/bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
/usr/lib/R/bin/exec/R --no-echo --no-restore --file=MyGSEMAnalysesRscript.R --args argument1 arument2
```
Best performance is achieved by leaving the values for these backends at `1` and maximizing the number of cores in GenomicSEM (only bound by CPU or RAM constraints).

Detailed test results can be found in [patchnotes](PATCHNOTES.md).


**License**

Copyright (C) 2018 Andrew Grotzinger, Matthijs van der Zee, Mijke Rhemtulla, Hill F. Ip, Michel Nivard, & Elliot Tucker-Drob

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
