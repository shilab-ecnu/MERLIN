---
editor_options: 
  markdown: 
    wrap: 72
---

# MERLIN: MEndelian Randomization for Linear INteraction

-   [Overview](#overview)
-   [Documentation](#documentation)
-   [System Requirements](#system-requirements)
-   [Installation](#installation)
-   [Demo](#demo)
    -   [Fit MERLIN using simulated
        data](#fit-merlin-using-simulated-data)
    -   [Fit MERLIN using real data](#fit-merlin-using-real-data)

# Overview {#overview}

`MERLIN` (Causal Heterogeneity using Summary Statistics) is a unified
Bayesian framework that jointly estimates average and heterogeneity
effects using summary data from genome-wide association studies (GWAS)
and Genome-Wide Interaction Studies (GWIS).

# Documenation

Full documentation available here:
<https://shilab-ecnu.github.io/MERLIN/>

# System Requirements {#system-requirements}

## Hardware requirements

`MERLIN` package requires only a standard computer with enough RAM to
support the in-memory operations.

## Software requirements

### OS Requirements

This package is supported for *macOS*, *Windows* and *Linux*. The
package has been tested on the following systems:

-   macOS: Tahoe (26.3.1)

-   Windows: 11

-   Linux: Red Hat 9.5

### Dependencies

```         
data.table
readr
Rcpp
TwoSampleMR (>= 0.7.9)
```

# Installation {#installation}

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

To update the package just run the
`remotes::install_github("shilab-ecnu/MERLIN")` command again.

# Demo {#demo}

This tutorial provides an introduction to the *MERLIN* package. R
package *MERLIN* implements MERLIN for causal heterogeneity using
summary statistics.

## Fit MERLIN using simulated data {#fit-merlin-using-simulated-data}

The simulation code in the paper is all in this hyperlink:
<a href="https://github.com/shilab-ecnu/MERLIN/tree/main/simulation">SIMULATION</a>.

We first generate the genotype data and the environmental variable:

```         
library(MERLIN)
library(mvtnorm)
library(MASS)
set.seed(2026)
```

```         
n_exp <- 80000
n_out <- 80000
m <- 100

# True Causal Parameters
b1 <- 0      
b4 <- 0.3    

# Generate Genotype matrix G
maf <- runif(m, 0.05, 0.5)
G <- matrix(rbinom((n_exp + n_out) * m, 2, rep(maf, each = n_exp + n_out)),
            nrow = n_exp + n_out, ncol = m)
G <- scale(G, center = TRUE, scale = FALSE)

# Generate Environmental variable E
E_exp <- sample(rep(c(-1, 1), each = n_exp / 2))
E_out <- sample(rep(c(-1, 1), each = n_out / 2))
E <- c(E_exp, E_out)
```

Now simulate the genetic effect sizes. The main genetic effects
$\gamma^{(G)}$ (we denote it as $\gamma_1$ in the tutorial code) and G×E
interaction effects $\gamma^{(GI)}$ (we denote it as $\gamma_3$ in the
tutorial code) are generated as correlated multivariate normal variables
with specified heritabilities.

```         
h_g1 <- 0.3    
h_g3 <- 0.1 
h_b <- 0.05     
cor_g1g3 <- 0.4  

sigma2g1 <- h_g1 / m
sigma2g3 <- h_g3 / m
sigma2b  <- h_b / m

cov_matrix <- matrix(c(sigma2g1, 
                      cor_g1g3 * sqrt(sigma2g1 * sigma2g3),
                      cor_g1g3 * sqrt(sigma2g1 * sigma2g3),
                      sigma2g3), ncol = 2)
                      
gamma1_3 <- rmvnorm(m, mean = c(0, 0), sigma = cov_matrix)
gamma_1x <- gamma1_3[, 1]
gamma_3x <- gamma1_3[, 2]
beta_2   <- rnorm(m, 0, sqrt(sigma2b))
```

Generate the exposure ($X$) and outcome ($Y$) variables with the genetic
effects defined.

```         
# Construct the GxE interaction term
GE <- G * E

# Define the correlation between residual errors
rhoxy <- 0.6  

# Calculate variances for the noise terms
var_noise_x <- 1 - h_g1 - h_g3
var_noise_y <- 1 - b1*b1 - h_b - b4*b4

# Construct the covariance matrix for the bivariate normal distribution
cov_xy <- rhoxy * sqrt(var_noise_x * var_noise_y)
sigma_noise <- matrix(c(var_noise_x, cov_xy,
                        cov_xy, var_noise_y), ncol = 2)

# Generate correlated noise for X and Y
noise_xy <- rmvnorm(n_exp + n_out, mean = c(0, 0), sigma = sigma_noise)
noise_x <- noise_xy[, 1]
noise_y <- noise_xy[, 2]

# Generate Exposure (X)
X <- G %*% gamma_1x + GE %*% gamma_3x + 0.1 * E + noise_x

# Generate Outcome (Y) 
Y <- X * b1 + G %*% beta_2 + X * E * b4 + 0.1 * E + noise_y

# Split into Exposure and Outcome datasets
exp_gwas <- X[1:n_exp]
exp_gwis <- X[1:n_exp]
exp_E    <- E[1:n_exp]

out_gwas <- Y[(n_exp + 1):(n_exp + n_out)]
out_gwis <- Y[(n_exp + 1):(n_exp + n_out)]
out_E    <- E[(n_exp + 1):(n_exp + n_out)]
```

We then conduct single-variant analysis to obtain the summary
statistics.

```         
get_sumstats <- function(G, pheno, interaction = FALSE, E = NULL) {
  betas <- numeric(ncol(G))
  ses <- numeric(ncol(G))
  
  for(i in 1:ncol(G)) {
    if(interaction) {
      # GWIS: Model with GxE interaction term
      model <- lm(pheno ~ G[, i] + E + G[, i]:E)
      betas[i] <- coef(model)[4]
      ses[i] <- summary(model)$coefficients[4, 2]
    } else {
      # GWAS: Standard additive genetic model
      model <- lm(pheno ~ G[, i])
      betas[i] <- coef(model)[2]
      ses[i] <- summary(model)$coefficients[2, 2]
    }
  }
  return(list(beta = betas, se = ses))
}

# Calculate Summary Statistics for Exposure Cohort
exp_gwas_sum <- get_sumstats(G[1:n_exp, ], exp_gwas)
exp_gwis_sum <- get_sumstats(G[1:n_exp, ], exp_gwis, interaction = TRUE, E = exp_E)

# Calculate Summary Statistics for Outcome Cohort
out_gwas_sum <- get_sumstats(G[(n_exp + 1):(n_exp + n_out), ], out_gwas)
out_gwis_sum <- get_sumstats(G[(n_exp + 1):(n_exp + n_out), ], out_gwis, 
                             interaction = TRUE, E = out_E)
```

We select genetic instruments using a p-value threshold. SNPs in either
the GWAS or GWIS analysis are included in the union set of instruments.

```         
p_threshold <- 5e-8;

pvals_gwas <- 2 * pnorm(-abs(exp_gwas_sum$beta / exp_gwas_sum$se))
iv_gwas <- which(pvals_gwas < p_threshold)

# Filter GWIS signals
pvals_gwis <- 2 * pnorm(-abs(exp_gwis_sum$beta / exp_gwis_sum$se))
iv_gwis <- which(pvals_gwis < p_threshold)

# Union of instruments and LD matrix
iv_union <- union(iv_gwas, iv_gwis)
R <- diag(length(iv_union))
```

Finally, we apply the MERLIN methods.

```         
gamma_hat <- exp_gwas_sum$beta[iv_union];
gamma3_hat <- exp_gwis_sum$beta[iv_union];
Gamma_hat <- out_gwas_sum$beta[iv_union];
Gamma3_hat <- out_gwis_sum$beta[iv_union];

se_gamma <- exp_gwas_sum$se[iv_union];
se_gamma3 <- exp_gwis_sum$se[iv_union];
se_Gamma <- out_gwas_sum$se[iv_union];
se_Gamma3 <- out_gwis_sum$se[iv_union];

rho_1 <- 0;
rho_2 <- 0;

res <- MERLIN(
  gammah1 = gamma_hat,
  gammah3 = gamma3_hat,
  Gammah1 = Gamma_hat,
  Gammah3 = Gamma3_hat,
  se1 = se_gamma,
  se2 = se_gamma3,
  se3 = se_Gamma,
  se4 = se_Gamma3,
  R = R,
  rho_1 = rho_1,
  rho_2 = rho_2,
  model = "standard",
  seed = 2026
)

str(res);

beta1_hat <- res$Beta1.hat;
se1_hat <- res$Beta1.se;
pval1 <- res$Beta1.pval;
beta4_hat <- res$Beta4.hat;
se4_hat <- res$Beta4.se;
pval4 <- res$Beta4.pval
```

beta1_hat, se1_hat, and pval1 are the estimated average causal effect,
corresponding standard error, and p-value of beta1_hat. beta4_hat,
se4_hat, and pval4 are the estimated heterogeneity causal effect,
corresponding standard error, and p-value of beta4_hat.

## Fit MERLIN using real data {#fit-merlin-using-real-data}

### The Testosterone-BD study with environmental factor sex

All the raw data for the real-data analyses in the replicated paper are
stored on
<a href="https://figshare.com/articles/dataset/Data_for_MERLIN/29910116">MERLIN
Dataset on Figshare</a>. Here, we take the dataset “The Testosterone–BD
study with the environmental factor sex” as an example.

Furthermore, we give an example to illustrate the implementation of
MERLIN for real data analysis. The following datasets
('Testosterone.GWAS.txt.gz', 'Testosterone.GWIS.txt.gz',
'BD.GWAS.txt.gz', 'BD.GWIS.txt.gz', 'g1000_eur.bed','g1000_eur.fam',
'g1000_eur.bim', 'all.bed') should be prepared. Download here:
<a href="https://figshare.com/articles/dataset/Data_for_MERLIN/29910116">MERLIN
Dataset on Figshare</a>

```         
expgwas <- "Testosterone.GWAS.txt.gz";
expgwis <- "Testosterone.GWIS.txt.gz";
outgwas <- "BD.GWAS.txt.gz";
outgwis <- "BD.GWIS.txt.gz";
stringname3 <- "g1000_eur";
block_file <- "all.bed";
```

'expgwas', 'expgwis', 'outgwas', and 'outgwis' are the dataset names for
exposure GWAS, exposure GWIS, outcome GWAS, and outcome GWIS,
respectively. Here, the environment variable is sex.

These four datasets must have the following format (note that it must be
tab-delimited): including columns as SNP, CHR, BP, A1, A2, BETA, SE, and
P.

**Example data format used for exposure and outcome summary statistics**

| SNP        | CHR |        BP | A1  | A2  |    BETA |     SE |      P |
|------------|----:|----------:|:---:|:---:|--------:|-------:|-------:|
| rs1000000  |  12 | 126890980 |  A  |  G  | -0.0023 | 0.0179 | 0.8969 |
| rs10000003 |   4 |  57561647 |  A  |  G  | -0.0099 | 0.0167 | 0.5535 |
| rs10000006 |   4 | 108826383 |  C  |  T  |  0.0118 | 0.0559 | 0.8329 |
| rs10000010 |   4 |  21618674 |  C  |  T  | -0.0113 | 0.0170 | 0.5070 |
| rs10000011 |   4 | 138223055 |  T  |  C  |  0.0246 | 0.0408 | 0.5470 |
| rs10000012 |   4 |   1357325 |  G  |  C  | -0.0128 | 0.0212 | 0.5455 |
| rs10000013 |   4 |  37225069 |  C  |  A  |  0.0062 | 0.0202 | 0.7600 |
| rs10000015 |   4 |  84143987 |  G  |  A  |  0.0366 | 0.0312 | 0.2404 |
| rs10000017 |   4 |  84778125 |  T  |  C  | -0.0151 | 0.0191 | 0.4292 |
| rs10000018 |   4 | 100458448 |  G  |  A  | -0.0123 | 0.0160 | 0.4410 |

If GWAS and GWIS data cannot be directly obtained, and the environmental
factor is a binary variable (e.g., sex), one can generate the required
GWAS and GWIS inputs for MERLIN by converting the sex-stratified summary
statistics (e.g., Testosterone.male.txt and Testosterone.female.txt) as
follows.

The GWAS summary statistics can be generated by meta-analyzing the male
and female data using inverse-variance weighting, as implemented in
METAL
(<a href="https://github.com/statgen/METAL">https://github.com/statgen/METAL</a>).
After installing the software, the analysis can be executed via the
command line (e.g., in Linux or other shell environments) using a
configuration file. A sample configuration file
'metal.config.Testosterone.txt' is available for download:
<a href="https://figshare.com/articles/dataset/Data_for_MERLIN/29910116">MERLIN
Dataset on Figshare</a>.

```         
metal metal.config.Testosterone.txt
```

The SNP effects and standard errors for GWIS summary statistics were
derived based on the following formula, assuming sex coded as Male=1,
Female=-1. Allele direction must be aligned before analyzing
sex-stratified data.

$$
\hat{b}_{gwis,j}=\frac{1}{2}(\hat{b}_{male,j}-{\hat{b}}_{female,j})
$$

$$
se(\hat{b}_{gwis,j})=\frac{1}{2}\sqrt{(se(\hat{b}_{male,j})^2+se(\hat{b}_{female,j})^2}
$$

'stringname3' is the name of the reference panel data. Here we use
samples from '1000 Genomes Project European panel', which is in plink
binary format. 'block_file' is used to partition the whole genome into
blocks.

The `matchpanel` function is used to match a GWAS/GWIS dataset with the
reference panel data, alongside initial data quality control. The output
includes a data frame (`$data`) and the corresponding storage path
(`$data_dir`).

```         
expgwas.match <- matchpanel(expgwas,stringname3)$data_dir;
expgwis.match <- matchpanel(expgwis,stringname3)$data_dir;
outgwas.match <- matchpanel(outgwas,stringname3)$data_dir;
outgwis.match <- matchpanel(outgwis,stringname3)$data_dir;
```

Having given that we have the formatted data, we can use the `ivselect`
function to screen the instrumental variables (IVs) and estimate the
correlations among those IVs. `plink_dir` specifies the local path to
the PLINK executable; if not provided, PLINK will be automatically
downloaded. `pval_cutoff_gwas` and `pval_cutoff_gwis` define the P-value
thresholds for the exposure GWAS and GWIS, respectively. `r2_cutoff` and
`kb_cutoff` are used in LD clumping to specify the $r^2$ threshold and
the physical distance (in kilobases) between SNPs. `maf_cutoff` sets the
threshold for minor allele frequency. `lam` denotes the shrinkage
parameter used in the regularization of the LD matrix. `CoreNum`
indicates the number of CPU cores to be used for parallel computation.
`intersect_mode` controls whether to merge GWAS and GWIS IVs using
intersection (default: union).

```         
plink_dir <- NULL;
pval_cutoff_gwas <- 5e-8;
pval_cutoff_gwis <- 5e-8;
r2_cutoff <- 0.5;
kb_cutoff <- 1024;
maf_cutoff <- 0.05;
lam <- 0.1;
coreNum <- 1;
intersect_mode <- FALSE;
```

```         
ivselect.res <- ivselect(expgwas.match, expgwis.match, outgwas.match, 
                         outgwis.match,
                         stringname3, block_file, plink_dir,
                         pval_cutoff_gwas, pval_cutoff_gwis, r2_cutoff, 
                         kb_cutoff, maf_cutoff, lam, coreNum, intersect_mode);
snp.causal <- ivselect.res$snp.causal;
gammah1 <- ivselect.res$gammah1;
gammah3 <- ivselect.res$gammah3;
Gammah1 <- ivselect.res$Gammah1;
Gammah3 <- ivselect.res$Gammah3;
se1 <- ivselect.res$se1;
se2 <- ivselect.res$se2;
se3 <- ivselect.res$se3;
se4 <- ivselect.res$se4;
R <- ivselect.res$R;
```

When the exposure and outcome samples are independent, the sample
correlation parameters `rho1` (for GWAS) and `rho2` (for GWIS) are set
to 0.

```         
rho1 <- 0; rho2 <- 0;
```

For data with sample overlap, $\rho_1$ and $\rho_2$ are estimated using
summary statistics among independent variants following
<a href="https://www.nature.com/articles/s41467-022-34164-1">Chen et al
(2022)</a>, we select independent SNPs using the truncated algorithm
($r^2$ threshold denoted by `ld_r2_thresh`). `pth` is the critical value
adapted to the truncated normal distribution in the estimation
procedure. `lambda` is the shrinkage turning parameter for LD estimator.

```         
ld_r2_thresh <- 0.001;
lambda <- 0.85;
pth <- 1.96;
RhoEst1 <- EstRhofun(expgwas, outgwas, stringname3, ld_r2_thresh, lambda, pth);
rho1 <- mean(RhoEst1$Rhores);
RhoEst2 <- EstRhofun(expgwis, outgwis, stringname3, ld_r2_thresh, lambda, pth);
rho2 <- mean(RhoEst2$Rhores);
```

Now we can fit MERLIN using the function *MERLIN*.

```         
res <- MERLIN(
  gammah1 = gammah1,
  gammah3 = gammah3,
  Gammah1 = Gammah1,
  Gammah3 = Gammah3,
  se1 = se1,
  se2 = se2,
  se3 = se3,
  se4 = se4,
  R = R,
  rho_1 = rho1,
  rho_2 = rho2,
  model = "standard",
  seed = 2026
)
str(res)
```

Check the convergence of the Gibbs sampler using `traceplot`.

```         
traceplot(res$Beta1res);
traceplot(res$Beta4res);
```

```         
MERLINbeta1 <- res$Beta1.hat;
MERLINse1 <- res$Beta1.se;
MERLINpvalue1 <- res$Beta1.pval;
MERLINbeta4 <- res$Beta4.hat;
MERLINse4 <- res$Beta4.se;
MERLINpvalue4 <- res$Beta4.pval;
```

```         
cat("The estimated average effect of Testosterone on BD: ", MERLINbeta1, 
    "\n with Standard error: ", MERLINse1, "and P-value: ", MERLINpvalue1);
```

```         
cat("The estimated heterogeneity effect of Testosterone on BD: ", MERLINbeta4, 
    "\n with Standard error: ", MERLINse4, "and P-value: ", MERLINpvalue4)
```

## Session information

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
