# Part 1: Introduction

## Introduction

Welcome to the MERLIN (MEndelian Randomization for Linear INteraction)
package.

Traditional Mendelian Randomization (MR) relies on genome-wide
association study (GWAS) summary statistics to infer the average causal
effect of an exposure on an outcome. However, it is increasingly
recognized that causal effects may not be uniform across populations but
rather exhibit heterogeneity driven by environmental factors.

MERLIN is designed to estimate this causal effect heterogeneity. By
systematically integrating summary statistics from both GWAS and
Genome-Wide Interaction Studies (GWIS), MERLIN provides a robust
statistical framework to evaluate how environmental variables modulate
the causal pathway from exposure to outcome.

**Key Features**

- Summary-Data Driven: Only requires GWAS and GWIS summary statistics,
  bypassing the need for individual-level genotype data (which is often
  restricted due to privacy concerns).

- Robust Framework: Effectively accounts for sample overlap and linkage
  disequilibrium (LD) among instrumental variables.

- Comprehensive Output: Provides estimates, standard errors, and
  p-values for both the main causal effect and the interaction
  (heterogeneity) effect.

### Installation

The development version of MERLIN can be installed directly from GitHub
using the remotes package.

Note: MERLIN depends on the Rcpp and RcppArmadillo packages for
computational efficiency. This requires appropriate C++ compiler
configurations (Rtools for Windows, and Xcode command-line tools for Mac
OS/X).

To install the package, execute the following commands in your R
console:

    if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")

    remotes::install_github("shilab-ecnu/MERLIN")

Once installed, load the package:

``` r

library(MERLIN)
```

### Session information

``` r

sessionInfo()
```

``` text
R version 4.4.3 (2025-02-28)
Platform: aarch64-apple-darwin20
Running under: macOS 26.4.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib; LAPACK version 3.12.0

locale:
[1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] MERLIN_1.1.2

loaded via a namespace (and not attached):
 [1] R6_2.6.1            tidyselect_1.2.1    bit_4.6.0
 [4] tzdb_0.5.0          magrittr_2.0.4      glue_1.8.0
 [7] tibble_3.3.1        parallel_4.4.3      pkgconfig_2.0.3
[10] bit64_4.6.0-1       lifecycle_1.0.5     readr_2.2.0
[13] cli_3.6.5           vctrs_0.7.2         data.table_1.18.2.1
[16] compiler_4.4.3      tools_4.4.3         hms_1.1.4
[19] pillar_1.11.1       Rcpp_1.1.1          crayon_1.5.3
[22] rlang_1.1.7         vroom_1.7.0
```
