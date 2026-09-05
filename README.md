# MERLIN: MEndelian Randomization for Linear INteraction

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation](#installation)

## Overview 

`MERLIN` (Causal Heterogeneity using Summary Statistics) is a unified
Bayesian framework that jointly estimates average and heterogeneity
effects using summary data from genome-wide association studies (GWAS)
and Genome-Wide Interaction Studies (GWIS). 

The simulation and real data analyses code for reproduction can be found here: <a href="https://github.com/ydong-work/MERLIN_analysis/tree/main">REPRODUCTION</a>. 

## Documenation

Full documentation available here:
<https://shilab-ecnu.github.io/MERLIN/>

## System Requirements 

### Hardware requirements

`MERLIN` package requires only a standard computer with enough RAM to
support the in-memory operations.

### Software requirements

#### OS Requirements

This package is supported for *macOS*, *Windows* and *Linux*. The
package has been tested on the following systems:

-   macOS: Tahoe (26.3.1)

-   Windows: 11

-   Linux: Red Hat 9.5

#### Dependencies

```         
data.table
readr
Rcpp
TwoSampleMR (>= 0.7.9)
```

## Installation

Install the development version of *MERLIN* by use of the 'remotes'
package. Note that *MERLIN* depends on the 'Rcpp' and 'RcppArmadillo'
packages, which also require appropriate settings of Rtools and Xcode
for Windows and Mac OS/X, respectively.

To install this package, run the following command in R.

```         
install.packages("remotes")
remotes::install_github("shilab-ecnu/MERLIN")
```

-   This takes 4-7 minutes to install.


