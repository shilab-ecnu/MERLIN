# MERLIN Functions

## Package index

This page provides a reference for the main user-facing functions,
helper functions, and lower-level C++ interfaces currently available in
**MERLIN**. It combines a compact package index with inline reference
entries, so users can read function descriptions, usage, arguments, and
return values on the same page.

### All functions

#### Main workflow

[`MERLIN()`](https://shilab-ecnu.github.io/MERLIN/reference/MERLIN-package.md)  
Run the MERLIN Bayesian Mendelian randomization model for average causal
effects and causal heterogeneity effects.

`ivselect()`  
Select instrumental variables, perform LD clumping, align summary
statistics, and estimate the LD correlation matrix.

`matchpanel()`  
Match a GWAS or GWIS summary-statistics file to a PLINK reference panel
and harmonise allele directions.

`EstRhofun()`  
Estimate sample-overlap correlation between two summary-statistics
files.

`traceplot()`  
Create a trace plot for saved MCMC samples.

#### Allele matching and summary-statistics helpers

`reverse()`  
Return the complementary nucleotide for a single allele.

`matchAellel()`  
Compare allele pairs and return sign factors for harmonising effect
directions.

`summaryQC()`  
Remove MHC-region SNPs and variants with extreme chi-square statistics.

`matchsnp()`  
Low-level C++ interface for matching SNPs between exposure, outcome, and
reference-panel data.

`matchscreen()`  
Low-level C++ interface for screening SNPs and matching them across
input files.

`ReadSNPinfo()`  
Low-level C++ helper for reading SNP information from a reference.

`Read_summarystat()`  
Low-level C++ helper for reading summary-statistics columns.

`select()`  
Low-level C++ helper for selecting character-vector elements by index.

`getLineNum()`  
Count the number of lines in a text file.

#### LD blocks and correlation matrices

`load_block_file()`  
Load genomic block definitions.

`test_blocks()`  
Assign variants to LD blocks and return block-level diagnostics.

`Cal_blockinf()`  
Calculate block membership information for variants.

`Cal_blockR()`  
Calculate block-wise LD information and independent variant indices.

`Cal_block_Rmatrix()`  
Calculate an LD correlation matrix for selected variants.

`Cal_block_Rvec()`  
Calculate a vectorised representation of block-wise LD correlations.

`IndepSummary()`  
Select approximately independent variants from summary statistics.

`LDclump()`  
Perform LD-based pruning or clumping on a correlation matrix.

`std_setdiff()`  
Low-level set-difference helper used during LD pruning.

#### MERLIN model kernels

`MRGEI_Gam3seo()`  
Run the standard MERLIN MCMC kernel.

`MRGEI_Gam3seo_addE2()`  
Run the continuous-environment MERLIN MCMC kernel.

`MRGEI_Gam3seo_binary()`  
Run the binary-environment MERLIN MCMC kernel.

`MRGEI_Gamseo_fixb1()`  
Run the `MO`(without outcome gwis) MERLIN MCMC kernel.

#### Matrix and numerical utilities

`comb()`  
Generate pairwise index combinations.

`Mat2Vec()`  
Convert the upper-triangular part of a matrix into a vector.

`Vec2Mat()`  
Reconstruct a symmetric matrix from a vectorised upper-triangular
representation.

`MatSum()`  
Calculate an element-wise vector summary used by truncated-normal
routines.

#### Normal and truncated-normal utilities

`phi()`  
Evaluate a standard normal cumulative distribution approximation.

`multiphi()`  
Apply `phi()` to a numeric vector.

`cdfNormal()`  
Evaluate a normal cumulative distribution function.

`MulticdfNormal()`  
Apply the standard normal cumulative distribution function to a numeric
vector.

`inverseNormal()`  
Evaluate a normal quantile function.

`MultiinverseNormal()`  
Apply the standard normal quantile function to a probability vector.

`truncEstfun()`  
Run a Gibbs sampler for bivariate truncated-normal correlation
estimation.

`testR()`  
Calculate a P-value for an estimated correlation coefficient.

## Detailed reference

### Main workflow

#### **`MERLIN()`**

Run the MERLIN Bayesian Mendelian randomization model. This wrapper
validates inputs and dispatches to a model-specific C++ MCMC routine.

##### Usage

``` r

MERLIN(
  gammah1,
  gammah3 = NULL,
  Gammah1,
  Gammah3 = NULL,
  se1 = NULL,
  se2 = NULL,
  se3 = NULL,
  se4 = NULL,
  R,
  rho_1 = NULL,
  rho_2 = NULL,
  model = c("standard", "continuous_E", "binary", "MO", "ME"),
  p1 = NULL,
  maxIter = 12000,
  burnin = 5000,
  thin = 10,
  seed = NULL
)
```

##### Arguments

| Argument | Description |
|----|----|
| `gammah1` | Numeric vector of exposure GWAS main-effect estimates. |
| `gammah3` | Numeric vector of exposure GWIS interaction-effect estimates. Required for `standard`, `continuous_E`, `binary`, and `MO`. |
| `Gammah1` | Numeric vector of outcome GWAS main-effect estimates. |
| `Gammah3` | Numeric vector of outcome GWIS interaction-effect estimates. Required for `standard`, `continuous_E`, `binary`, and `ME`. |
| `se1`, `se2`, `se3`, `se4` | Standard-error vectors corresponding to `gammah1`, `gammah3`, `Gammah1`, and `Gammah3`. |
| `R` | Numeric LD correlation matrix for the selected SNPs. |
| `rho_1` | Correlation parameter for the GWAS summary-statistics components. |
| `rho_2` | Correlation parameter for the GWIS summary-statistics components. Required for models using both GWIS terms. |
| `model` | MERLIN model variant. |
| `p1` | Proportion parameter required by the `binary` model. |
| `maxIter` | Total number of MCMC iterations. Must be divisible by `thin`. |
| `burnin` | Number of burn-in iterations. |
| `thin` | Thinning interval. |
| `seed` | Optional random seed for reproducibility. |

##### Value

A list returned by the selected C++ model kernel. The public workflow
typically uses elements such as `Beta1.hat`, `Beta1.se`, `Beta1.pval`,
`Beta4.hat`, `Beta4.se`, `Beta4.pval`, and saved MCMC samples such as
`Beta1res` and `Beta4res`.

#### **`ivselect()`**

Select instrumental variables and calculate the LD matrix used by
MERLIN.

##### Usage

``` r

ivselect(
  expgwas_dir,
  expgwis_dir = NULL,
  outgwas_dir,
  outgwis_dir = NULL,
  stringname3,
  block_file,
  plink_dir = NULL,
  pval_cutoff_gwas = 0.00000005,
  pval_cutoff_gwis = 0.00000005,
  r2_cutoff = 0.01,
  kb_cutoff = 1024,
  maf_cutoff = 0.05,
  lam = 0.1,
  coreNum = 1,
  intersect_mode = FALSE
)
```

##### Arguments

| Argument | Description |
|----|----|
| `expgwas_dir` | Path to exposure GWAS summary statistics. |
| `expgwis_dir` | Optional path to exposure GWIS summary statistics. |
| `outgwas_dir` | Path to outcome GWAS summary statistics. |
| `outgwis_dir` | Optional path to outcome GWIS summary statistics. |
| `stringname3` | Prefix of the PLINK reference panel files. |
| `block_file` | Genomic block definition file. |
| `plink_dir` | Path to the PLINK executable. If `NULL`, the function attempts to download PLINK via `bigsnpr::download_plink()`. |
| `pval_cutoff_gwas` | P-value threshold for exposure GWAS clumping. |
| `pval_cutoff_gwis` | P-value threshold for exposure GWIS clumping. |
| `r2_cutoff` | LD `r^2` threshold used by PLINK clumping. |
| `kb_cutoff` | Physical distance window in kilobases for PLINK clumping. |
| `maf_cutoff` | Minor allele frequency threshold for PLINK clumping. |
| `lam` | Regularisation parameter used when estimating the LD matrix. |
| `coreNum` | Number of CPU cores used in LD calculations. |
| `intersect_mode` | If `TRUE`, use the intersection of GWAS and GWIS IVs; otherwise use their union followed by LD clumping. |

##### Value

A list with selected SNPs, harmonised summary-statistics vectors,
standard-error vectors, and the LD matrix: `snp.causal`, `gammah1`,
`gammah3`, `Gammah1`, `Gammah3`, `se1`, `se2`, `se3`, `se4`, and `R`.

#### **`matchpanel()`**

Match a summary-statistics file to a reference panel and harmonise
effect directions.

##### Usage

``` r

matchpanel(ss, stringname3)
```

##### Arguments

| Argument | Description |
|----|----|
| `ss` | Path to a GWAS or GWIS summary-statistics file. The expected columns include `SNP`, `CHR`, `BP`, `A1`, `A2`, `BETA`, `SE`, and `P`. |
| `stringname3` | Prefix of the PLINK reference panel files. The function reads `paste0(stringname3, ".bim")`. |

##### Value

A list with `data`, the matched data frame, and `data_dir`, the path to
the written matched file. The output file is named by prefixing the
input file name with `match.`.

#### **`EstRhofun()`**

Estimate the correlation parameter induced by sample overlap.

##### Usage

``` r

EstRhofun(fileexposure, fileoutcome, stringname3, ld_r2_thresh, lam, pth, coreNum = 1)
```

##### Arguments

| Argument | Description |
|----|----|
| `fileexposure` | Path to exposure summary statistics. |
| `fileoutcome` | Path to outcome summary statistics. |
| `stringname3` | Prefix of the PLINK reference panel files. |
| `ld_r2_thresh` | LD threshold used when selecting approximately independent variants. |
| `lam` | Regularisation parameter for LD estimation. |
| `pth` | Truncation threshold or vector of thresholds used in the truncated-normal estimator. |
| `coreNum` | Number of CPU cores used in LD calculations. |

##### Value

A list containing `rhohat`, estimated correlation values; `pvalue`,
correlation test P-values; `pres`, the number of variants used for each
threshold; and `Rhores`, saved samples from the truncated-normal
correlation estimator.

#### **`traceplot()`**

Create a trace plot for a vector of saved MCMC samples.

##### Usage

``` r

traceplot(bhatpoint)
```

##### Arguments

| Argument | Description |
|----|----|
| `bhatpoint` | Numeric vector of saved MCMC samples, such as `res$Beta1res` or `res$Beta4res`. |

##### Value

A `ggplot` object showing sample index on the x-axis and sampled effect
values on the y-axis.

### Allele matching and summary-statistics helpers

#### **`reverse()`**

Return the complementary nucleotide for a single allele.

##### Usage

``` r

reverse(x)
```

##### Arguments

| Argument | Description |
|----|----|
| `x` | Allele string. Recognised values are `A`, `T`, `C`, and `G`. Other values are returned unchanged. |

##### Value

The complementary allele for recognised nucleotides.

#### **`matchAellel()`**

Compare allele pairs from summary statistics and the reference panel.

##### Usage

``` r

matchAellel(al_gwas, al_gtex)
```

##### Arguments

| Argument | Description |
|----|----|
| `al_gwas` | Two-column matrix or data frame of alleles from summary statistics. |
| `al_gtex` | Two-column matrix or data frame of alleles from the reference panel. |

##### Value

A numeric sign-factor vector. Values are `1` for matched alleles, `-1`
for swapped alleles, and `0` for unresolved mismatches.

#### **`summaryQC()`**

Remove variants from problematic genomic regions or with extreme
statistics.

##### Usage

``` r

summaryQC(mhcstart, mhcend, bh1, bh2, s12, s22, bp, chr, rsname,avbIndex, idx4panel, xbound, ybound)
```

##### Arguments

| Argument | Description |
|----|----|
| `mhcstart`, `mhcend` | Start and end base-pair positions defining the MHC exclusion region on chromosome 6. |
| `bh1`, `bh2` | Effect estimates for two traits or summary-statistics components. |
| `s12`, `s22` | Standard errors corresponding to `bh1` and `bh2`. |
| `bp`, `chr` | Base-pair positions and chromosomes for SNPs. |
| `rsname` | SNP identifiers. |
| `avbIndex` | Variant indices in the reference panel. |
| `idx4panel` | Indices used for panel-specific processing. |
| `xbound`, `ybound` | Chi-square thresholds for filtering exposure and outcome components. |

##### Value

A list of filtered effects, standard errors, positions, SNP names,
indices, and counts of removed SNPs: `pmhc`, `px`, and `py`.

#### **`matchsnp()`**

Low-level C++ interface for matching SNPs across exposure summary
statistics, outcome summary statistics, and a reference panel.

##### Usage

``` r

matchsnp(stringname1, stringname2, stringname3, matchExp)
```

##### Arguments

| Argument | Description |
|----|----|
| `stringname1` | Path to the first summary-statistics file, typically exposure. |
| `stringname2` | Path to the second summary-statistics file, typically outcome. |
| `stringname3` | Prefix of the PLINK reference panel files. |
| `matchExp` | Logical flag controlling matching orientation in the C++ routine. |

##### Value

A list containing matched effect estimates, standard errors,
chromosomes, positions, SNP identifiers, and reference-panel indices.

#### **`matchscreen()`**

Low-level C++ interface for P-value screening and SNP matching.

##### Usage

``` r

matchscreen(screenname, stringname1, stringname2, stringname3, pva_cutoff, matchExp = FALSE)
```

##### Arguments

| Argument | Description |
|----|----|
| `screenname` | Summary-statistics file used for initial P-value screening. |
| `stringname1`, `stringname2` | Additional summary-statistics files to match. |
| `stringname3` | Prefix of the PLINK reference panel files. |
| `pva_cutoff` | P-value threshold used for screening. |
| `matchExp` | Logical flag controlling matching orientation in the C++ routine. |

##### Value

A list of screened and matched summary-statistics components and
reference-panel indices.

#### **`ReadSNPinfo()`**

Read SNP information from a reference file into preallocated vectors.

##### Usage

``` r

ReadSNPinfo(stringname, A1, A2, rsname, chr, bp, morgan, N)
```

##### Arguments

| Argument     | Description                                                 |
|--------------|-------------------------------------------------------------|
| `stringname` | Path to the SNP information file.                           |
| `A1`, `A2`   | Preallocated vectors for allele codes.                      |
| `rsname`     | Preallocated vector for SNP identifiers.                    |
| `chr`, `bp`  | Preallocated vectors for chromosome and base-pair position. |
| `morgan`     | Preallocated vector for genetic-map positions.              |
| `N`          | Number of records to read.                                  |

##### Value

A list containing the populated SNP information vectors.

#### **`Read_summarystat()`**

Read summary-statistics columns into preallocated vectors.

##### Usage

``` r

Read_summarystat(stringname, SA1, SA2, rsname, betah, s2, pvalue, chr, bp, N)
```

##### Arguments

| Argument     | Description                                                 |
|--------------|-------------------------------------------------------------|
| `stringname` | Path to the summary-statistics file.                        |
| `SA1`, `SA2` | Preallocated allele-code vectors.                           |
| `rsname`     | Preallocated SNP identifier vector.                         |
| `betah`      | Preallocated vector for effect estimates.                   |
| `s2`         | Preallocated vector for standard errors.                    |
| `pvalue`     | Preallocated vector for P-values.                           |
| `chr`, `bp`  | Preallocated vectors for chromosome and base-pair position. |
| `N`          | Number of records to read.                                  |

##### Value

Called for its side effects on preallocated vectors. The R wrapper
returns invisibly.

#### **`select()`**

Select elements from a character vector by index.

##### Usage

``` r

select(vec_, idx_)
```

##### Arguments

| Argument | Description                |
|----------|----------------------------|
| `vec_`   | Character vector.          |
| `idx_`   | Numeric vector of indices. |

##### Value

A character vector containing selected elements.

#### **`getLineNum()`**

Count the number of lines in a file.

##### Usage

``` r

getLineNum(filename)
```

##### Arguments

| Argument   | Description       |
|------------|-------------------|
| `filename` | Path to the file. |

##### Value

Integer line count.

### LD blocks and correlation matrices

#### **`load_block_file()`**

Load genomic block definitions.

##### Usage

``` r

load_block_file(block_file)
```

##### Arguments

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `block_file` | Path to the genomic block definition file. |

##### Value

An object returned by the C++ loader containing block definitions.

#### **`test_blocks()`**

Assign variants to blocks and return diagnostic information.

##### Usage

``` r

test_blocks(bp, chr, block_file)
```

##### Arguments

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `bp`         | Integer vector of base-pair positions.     |
| `chr`        | Integer vector of chromosome numbers.      |
| `block_file` | Path to the genomic block definition file. |

##### Value

A list with block assignment and block-level diagnostic information.

#### **`Cal_blockinf()`**

Calculate block membership information for variants.

##### Usage

``` r

Cal_blockinf(bp, chr, block_file)
```

##### Arguments

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `bp`         | Integer vector of base-pair positions.     |
| `chr`        | Integer vector of chromosome numbers.      |
| `block_file` | Path to the genomic block definition file. |

##### Value

A list of block membership information.

#### **`Cal_blockR()`**

Calculate block-wise LD information and independent variant indices.

##### Usage

``` r

Cal_blockR(bp, chr, avbIndex, idx4panel, block_file, stringname3, ld_r2_thresh, coreNum, lam)
```

##### Arguments

| Argument       | Description                                            |
|----------------|--------------------------------------------------------|
| `bp`, `chr`    | Base-pair positions and chromosomes for selected SNPs. |
| `avbIndex`     | Zero-based reference-panel variant indices.            |
| `idx4panel`    | Zero-based indices for panel-specific processing.      |
| `block_file`   | Genomic block definition file.                         |
| `stringname3`  | Prefix of the PLINK reference panel files.             |
| `ld_r2_thresh` | LD threshold for selecting independent variants.       |
| `coreNum`      | Number of CPU cores.                                   |
| `lam`          | LD regularisation parameter.                           |

##### Value

A list containing block-wise LD results and independent variant indices.

`stringname3` \| Prefix of the PLINKreference panel files. \|  
`ld_r2_thresh` \| LD threshold for independent-SNP selection.\|  
`coreNum` \| Number of CPU cores. \|  
`lam` \| LD regularisation parameter. \|

##### Value

A list containing the LD matrix `R` and related block information.

#### **`Cal_block_Rvec()`**

Calculate a vectorised representation of block-wise LD correlations.

##### Usage

``` r

Cal_block_Rvec(bp, chr, avbIndex, block_file, stringname3, coreNum, lam)
```

##### Arguments

| Argument      | Description                                            |
|---------------|--------------------------------------------------------|
| `bp`, `chr`   | Base-pair positions and chromosomes for selected SNPs. |
| `avbIndex`    | Zero-based reference-panel variant indices.            |
| `block_file`  | Genomic block definition file.                         |
| `stringname3` | Prefix of the PLINK reference panel files.             |
| `coreNum`     | Number of CPU cores.                                   |
| `lam`         | LD regularisation parameter.                           |

##### Value

A list containing vectorised LD-correlation output and block
information.

#### **`IndepSummary()`**

Select approximately independent variants and return their summary
statistics.

##### Usage

``` r

IndepSummary(bp, chr, avbIndex, block_file, stringname3, bh1, bh2, se1, se2, coreNum, lam, ld_r2_thresh)
```

##### Arguments

| Argument       | Description                                             |
|----------------|---------------------------------------------------------|
| `bp`, `chr`    | Base-pair positions and chromosomes.                    |
| `avbIndex`     | Zero-based reference-panel variant indices.             |
| `block_file`   | Genomic block definition file.                          |
| `stringname3`  | Prefix of the PLINK reference panel files.              |
| `bh1`, `bh2`   | Effect estimates for two summary-statistics components. |
| `se1`, `se2`   | Standard errors corresponding to `bh1` and `bh2`.       |
| `coreNum`      | Number of CPU cores.                                    |
| `lam`          | LD regularisation parameter.                            |
| `ld_r2_thresh` | LD threshold for independent variant selection.         |

##### Value

A list containing independent-variant effects and standard errors,
including `bh1_ind`, `bh2_ind`, `se1_ind`, and `se2_ind`.

#### **`LDclump()`**

Perform LD pruning or clumping on a correlation matrix.

##### Usage

``` r

LDclump(R, ld_r2_thresh)
```

##### Arguments

| Argument       | Description                                           |
|----------------|-------------------------------------------------------|
| `R`            | LD correlation matrix.                                |
| `ld_r2_thresh` | Squared-correlation threshold used to prune variants. |

##### Value

An integer vector of retained variant indices.

#### **`std_setdiff()`**

Return the set difference between two integer vectors.

##### Usage

``` r

std_setdiff(x, y)
```

##### Arguments

| Argument | Description      |
|----------|------------------|
| `x`, `y` | Integer vectors. |

##### Value

Elements of `x` not present in `y`.

### MERLIN model kernels

The following functions are low-level C++ MCMC kernels. Most users
should call
[`MERLIN()`](https://shilab-ecnu.github.io/MERLIN/reference/MERLIN-package.md)
instead.

#### **`MRGEI_Gam3seo()`**

Run the standard MERLIN MCMC kernel.

##### Usage

``` r

MRGEI_Gam3seo(gammah1, gammah3, Gammah1, Gammah3, se1, se2, se3, se4, R, rho_1, rho_2, maxIter, burnin, thin)
```

##### Arguments

| Argument | Description |
|----|----|
| `gammah1`, `gammah3` | Exposure GWAS and GWIS effect estimates. |
| `Gammah1`, `Gammah3` | Outcome GWAS and GWIS effect estimates. |
| `se1`, `se2`, `se3`, `se4` | Corresponding standard errors. |
| `R` | LD correlation matrix. |
| `rho_1`, `rho_2` | Correlation parameters for GWAS and GWIS components. |
| `maxIter`, `burnin`, `thin` | MCMC control parameters. |

##### Value

A list of MERLIN estimates and saved MCMC samples.

#### **`MRGEI_Gam3seo_addE2()`**

Run the continuous-environment MERLIN MCMC kernel.

##### Usage

``` r

MRGEI_Gam3seo_addE2(gammah1, gammah3, Gammah1, Gammah3, se1, se2, se3, se4, R, rho_1, rho_2, maxIter, burnin, thin)
```

##### Arguments

Same as `MRGEI_Gam3seo()`.

##### Value

A list of MERLIN estimates and saved MCMC samples.

#### **`MRGEI_Gam3seo_binary()`**

Run the binary-environment MERLIN MCMC kernel.

##### Usage

``` r

MRGEI_Gam3seo_binary(gammah1, gammah3, Gammah1, Gammah3, se1, se2, se3, se4, R, rho_1, rho_2, p1, maxIter, burnin, thin)
```

##### Arguments

| Argument | Description |
|----|----|
| `gammah1`, `gammah3` | Exposure GWAS and GWIS effect estimates. |
| `Gammah1`, `Gammah3` | Outcome GWAS and GWIS effect estimates. |
| `se1`, `se2`, `se3`, `se4` | Corresponding standard errors. |
| `R` | LD correlation matrix. |
| `rho_1`, `rho_2` | Correlation parameters for GWAS and GWIS components. |
| `p1` | Binary-environment proportion parameter. |
| `maxIter`, `burnin`, `thin` | MCMC control parameters. |

##### Value

A list of MERLIN estimates and saved MCMC samples.

#### **`MRGEI_Gamseo_fixb1()`**

Run the `MO`(without outcome gwis) MERLIN MCMC kernel.

##### Usage

``` r

MRGEI_Gamseo_fixb1(gammah1, gammah3, Gammah1, se1, se2, se3, R, rho, b1, maxIter, burnin, thin)
```

##### Arguments

| Argument | Description |
|----|----|
| `gammah1`, `gammah3` | Exposure GWAS and GWIS effect estimates. |
| `Gammah1` | Outcome GWAS effect estimates. |
| `se1`, `se2`, `se3` | Corresponding standard errors. |
| `R` | LD correlation matrix. |
| `rho` | Correlation parameter. |
| `b1` | The main effect estimate (\beta_1), which can be either manually specified or obtained from an alternative MR method. When the “`MO`” mode is selected, the Egger method is used to estimate \beta_1. |
| `maxIter`, `burnin`, `thin` | MCMC control parameters. |

##### Value

A list of MERLIN estimates and saved MCMC samples.

#### **`Mat2Vec()`**

Convert the upper-triangular part of a matrix into a vector.

##### Usage

``` r

Mat2Vec(R)
```

##### Arguments

| Argument | Description    |
|----------|----------------|
| `R`      | Square matrix. |

##### Value

A vector containing the upper-triangular entries of `R`.

#### **`Vec2Mat()`**

Reconstruct a symmetric matrix from a vectorised upper-triangular
representation.

##### Usage

``` r

Vec2Mat(RV, p1)
```

##### Arguments

| Argument | Description                            |
|----------|----------------------------------------|
| `RV`     | Vectorised upper-triangular entries.   |
| `p1`     | Dimension of the target square matrix. |

##### Value

A symmetric matrix.

#### **`MatSum()`**

Calculate an element-wise vector summary used by truncated-normal
routines.

##### Usage

``` r

MatSum(y1, y2)
```

##### Arguments

| Argument   | Description      |
|------------|------------------|
| `y1`, `y2` | Numeric vectors. |

##### Value

A numeric vector returned by the C++ routine.

### Normal and truncated-normal utilities

#### **`phi()`**

Evaluate an approximation to the standard normal cumulative distribution
function.

##### Usage

``` r

phi(x)
```

##### Arguments

| Argument | Description    |
|----------|----------------|
| `x`      | Numeric value. |

##### Value

Approximate standard normal cumulative probability.

#### **`multiphi()`**

Apply `phi()` to each element of a numeric vector.

##### Usage

``` r

multiphi(x)
```

##### Arguments

| Argument | Description     |
|----------|-----------------|
| `x`      | Numeric vector. |

##### Value

Numeric vector of standard normal cumulative probabilities.

#### **`cdfNormal()`**

Evaluate a normal cumulative distribution function.

##### Usage

``` r

cdfNormal(x, mean, sd)
```

##### Arguments

| Argument | Description                |
|----------|----------------------------|
| `x`      | Numeric value.             |
| `mean`   | Normal mean.               |
| `sd`     | Normal standard deviation. |

##### Value

Normal cumulative probability.

#### **`MulticdfNormal()`**

Apply the standard normal cumulative distribution function to a vector.

##### Usage

``` r

MulticdfNormal(x)
```

##### Arguments

| Argument | Description     |
|----------|-----------------|
| `x`      | Numeric vector. |

##### Value

Numeric vector of standard normal cumulative probabilities.

#### **`inverseNormal()`**

Evaluate a normal quantile function.

##### Usage

``` r

inverseNormal(prob, mean, sd)
```

##### Arguments

| Argument | Description                |
|----------|----------------------------|
| `prob`   | Probability value.         |
| `mean`   | Normal mean.               |
| `sd`     | Normal standard deviation. |

##### Value

Normal quantile.

#### **`MultiinverseNormal()`**

Apply the standard normal quantile function to a probability vector.

##### Usage

``` r

MultiinverseNormal(prob)
```

##### Arguments

| Argument | Description         |
|----------|---------------------|
| `prob`   | Probability vector. |

##### Value

Numeric vector of standard normal quantiles.

#### **`truncEstfun()`**

Run a Gibbs sampler for estimating correlation from bivariate
truncated-normal data.

##### Usage

``` r

truncEstfun(a, b, x1, x2, maxIter, burnin, thin)
```

##### Arguments

| Argument   | Description                                         |
|------------|-----------------------------------------------------|
| `a`, `b`   | Lower and upper truncation bounds.                  |
| `x1`, `x2` | Numeric vectors of paired observations or Z-scores. |
| `maxIter`  | Total number of MCMC iterations.                    |
| `burnin`   | Number of burn-in iterations.                       |
| `thin`     | Thinning interval.                                  |

##### Value

A numeric vector of saved correlation samples.

#### **`testR()`**

Calculate a P-value for an estimated correlation coefficient.

##### Usage

``` r

testR(rho, n)
```

##### Arguments

| Argument | Description                                              |
|----------|----------------------------------------------------------|
| `rho`    | Estimated correlation coefficient.                       |
| `n`      | Number of variants or observations used in the estimate. |

##### Value

A P-value for testing the correlation estimate.

## Suggested public API

For website-facing documentation, the clearest public API is likely:

- [`MERLIN()`](https://shilab-ecnu.github.io/MERLIN/reference/MERLIN-package.md)
- `matchpanel()`
- `ivselect()`
- `EstRhofun()`
- `traceplot()`

The remaining functions are best documented as internal helpers or
low-level developer interfaces unless they are intentionally supported
for end users.
