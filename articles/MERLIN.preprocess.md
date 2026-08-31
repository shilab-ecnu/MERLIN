# Part 4: Data Preprocessing and Instrumental Variable Selection

## Overview of the Preprocessing Pipeline

Before fitting the MERLIN model on summary statistics, a rigorous data
preprocessing workflow is required. This workflow ensures that:

- Genetic Variants Alignment: Summary statistics from different studies
  (exposure and outcome) are perfectly aligned with the reference panel
  alleles.

- Instrumental Variable (IV) Selection: Valid genetic instruments are
  selected based on genome-wide significance thresholds from GWAS and/or
  GWIS.

- Linkage Disequilibrium (LD) Estimation: The correlation structure (LD
  matrix) among selected instruments is properly calculated.

- Sample Overlap Adjustment: Correlated estimation errors are accounted
  for if the exposure and outcome cohorts share overlapping
  participants.

## Step 1: Reference Panel Matching (`matchpanel`)

The matchpanel function cross-references your summary statistics with a
PLINK format reference panel (e.g., 1000 Genomes European panel). It
filters out SNPs with mismatched alleles, performing initial quality
control (QC).

**Parameter Specifications**

- `summary_file`: Path to the input GWAS/GWIS summary data
  (tab-delimited).
- `reference_panel`: Base name of the PLINK binary files (`.bed`,
  `.bim`, `.fam`).

**Execution**

Run the matching process for all four input datasets. The function
returns a list containing the filtered dataframe (`$data`) and the file
path to the newly saved matched dataset (`$data_dir`).

To illustrate the implementation of MERLIN, we applied the method to
investigate the causal effect of testosterone on Bipolar Disorder (BD),
using sex as the environmental variable to assess potential
gene–environment interaction effects.

The following datasets (‘Testosterone.GWAS.txt.gz’,
‘Testosterone.GWIS.txt.gz’, ‘BD.GWAS.txt.gz’, ‘BD.GWIS.txt.gz’,
‘g1000_eur.bed’,‘g1000_eur.fam’, ‘g1000_eur.bim’, ‘all.bed’) should be
prepared. Download here:
<https://figshare.com/articles/dataset/Data_for_MERLIN/29910116>.

``` r

# Define the names/paths of your raw files and reference panel
library(MERLIN)

expgwas     <- "Testosterone.GWAS.txt.gz"
expgwis     <- "Testosterone.GWIS.txt.gz"
outgwas     <- "BD.GWAS.txt.gz"
outgwis     <- "BD.GWIS.txt.gz"
stringname3 <- "g1000_eur"
block_file  <- "all.bed"
```

``` r

# Perform matching and extract the directories of the output matched data
expgwas.match <- matchpanel(expgwas, stringname3)$data_dir
expgwis.match <- matchpanel(expgwis, stringname3)$data_dir
outgwas.match <- matchpanel(outgwas, stringname3)$data_dir
outgwis.match <- matchpanel(outgwis, stringname3)$data_dir
```

The console output below shows the progress of processing the exposure
GWAS dataset. The remaining datasets are processed in the same manner
and are therefore omitted for brevity.

``` text
Starting processing of file: Testosterone.GWAS.txt.gz
Step 1/6: Reading reference panel (bim file)...
|--------------------------------------------------|
|==================================================|
Step 2/6: Reading GWAS summary statistics...
Rows: 9686364 Columns: 8── Column specification ─────────────────────────────────────────────────────────────────────────────────────────
Delimiter: " "
chr (3): SNP, A1, A2
dbl (5): CHR, BP, BETA, SE, P
ℹ Use `spec()` to retrieve the full column specification for this data.
ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
Step 3/6: Preprocessing GWAS data...
Step 4/6: Merging with reference panel...
Step 5/6: Matching allele directions...
Step 6/6: Saving matched data...
                                                                                                                      
Successfully processed 9686364 SNPs
Total processing time: 352.6 seconds
Output saved to: ./match.Testosterone.GWAS.txt.gz
```

## Step 2: Instrument Screening and LD Estimation (`ivselect`)

The `ivselect` function is the core preprocessing component. It selects
instrumental variables using Linkage Disequilibrium (LD) clumping and
extracts the effect sizes and standard errors for downstream analysis.

Detailed Parameter Interpretations

To customize the clumping and regularizing process, `ivselect` utilizes
the following parameters:

- `plink_dir`: Local path to the PLINK executable. If set to `NULL`, the
  package automatically downloads an appropriate version for your OS.

- `pval_cutoff_gwas / pval_cutoff_gwis`: Significance thresholds for
  extracting potential IVs from GWAS and GWIS respectively.

- `r2_cutoff`: The r^2 threshold for LD clumping. SNPs with correlations
  higher than this value within the specified physical distance will be
  clumped.

- `kb_cutoff`: The physical distance window (in kilobases) used for LD
  clumping.

- `maf_cutoff`: Minor Allele Frequency threshold for quality control.

- `lam`: Shrinkage turning parameter (\lambda) used for regularizing the
  estimated LD correlation matrix to guarantee positive definiteness.

- `coreNum`: Number of CPU cores allocated for parallel processing to
  speed up calculation over genomic blocks.

- `intersect_mode`: Logical flag. If `TRUE`, instruments must pass
  thresholds in both GWAS and GWIS (intersection). If `FALSE`,
  instruments passing thresholds in either GWAS or GWIS are included
  (union, default).

**Running `ivselect`**

``` r

# Define clumping and filtering thresholds
plink_dir        <- NULL
pval_cutoff_gwas <- 5e-8
pval_cutoff_gwis <- 5e-8
r2_cutoff        <- 0.5
kb_cutoff        <- 1024
maf_cutoff       <- 0.05
lam              <- 0.1
coreNum          <- 1
intersect_mode   <- FALSE
```

``` r

# Execute instrument screening
ivselect.res <- ivselect(expgwas.match, expgwis.match, outgwas.match, outgwis.match,
                         stringname3, block_file, plink_dir,
                         pval_cutoff_gwas, pval_cutoff_gwis, r2_cutoff, 
                         kb_cutoff, maf_cutoff, lam, coreNum, intersect_mode);
snp.causal <- ivselect.res$snp.causal;
gammah1    <- ivselect.res$gammah1;
gammah3    <- ivselect.res$gammah3;
Gammah1    <- ivselect.res$Gammah1;
Gammah3    <- ivselect.res$Gammah3;
se1        <- ivselect.res$se1;
se2        <- ivselect.res$se2;
se3        <- ivselect.res$se3;
se4        <- ivselect.res$se4;
R          <- ivselect.res$R;
```

The console output below shows the progress of IV selection.

``` text
Starting MERLIN IV selection process...
Step 1/7: Reading reference panel (bim file)...

Step 2/7: Reading GWAS summary statistics...

Step 3/7: Setting up PLINK...
  ✓ Using PLINK: user/plink/plink
  
Step 4/7: Running LD clumping for exposure GWAS...
PLINK v1.90b6.1 64-bit (28 May 2018)           www.cog-genomics.org/plink/1.9/
(C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ld_result1.log.
Options in effect:
  --bfile g1000_eur
  --clump match.Testosterone.GWAS.txt.gz
  --clump-field P
  --clump-kb 1024
  --clump-p1 5e-08
  --clump-r2 0.5
  --clump-snp-field SNP
  --maf 0.05
  --out ld_result1
......
--clump: 536 clumps formed from 9955 top variants.
Results written to ld_result1.clumped .

Step 5/7: Running LD clumping for exposure GWIS...
PLINK v1.90b6.1 64-bit (28 May 2018)           www.cog-genomics.org/plink/1.9/
(C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ld_result2.log.
Options in effect:
  --bfile g1000_eur
  --clump match.Testosterone.GWIS.txt.gz
  --clump-field P
  --clump-kb 1024
  --clump-p1 5e-08
  --clump-r2 0.5
  --clump-snp-field SNP
  --maf 0.05
  --out ld_result2
......
--clump: 393 clumps formed from 9182 top variants.
Results written to ld_result2.clumped.

Step 6/7: Selecting causal SNPs...
  - Mode: Union
PLINK v1.90b6.1 64-bit (28 May 2018)           www.cog-genomics.org/plink/1.9/
(C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ld_result3.log.
Options in effect:
  --bfile g1000_eur
  --clump union.snp.txt
  --clump-field P
  --clump-kb 1024
  --clump-p1 1
  --clump-r2 0.5
  --clump-snp-field SNP
  --maf 0.05
  --out ld_result3
......
--clump: 762 clumps formed from 880 top variants.
Results written to ld_result3.clumped.
Numbers of causal SNPs selected:342

Step 7/7: Calculating correlation matrix...

--------------------------------------------------
IV selection completed successfully!
Total SNPs selected: 342
Total processing time: 459.8 seconds
--------------------------------------------------
```

**Understanding the Extracted Outputs** The ivselect function returns a
structured list containing main components required by the MERLIN
estimation algorithm:

- `snp.causal`: Vector of RS IDs for the selected independent
  instrumental variables.

- `gammah1` / `gammah3`: Vector of estimated genetic main effects
  (\hat{\gamma}\_1) and interaction effects (\hat{\gamma}\_3) on the
  exposure.

- `Gammah1` / `Gammah3`: Vector of estimated genetic main effects
  (\hat{\Gamma}\_1) and interaction effects (\hat{\Gamma}\_3) on the
  outcome.

- `se1` / `se2` / `se3` / `se4`: Respective standard errors for
  `gammah1`, `gammah3`, `Gammah1`, and `Gammah3`.

- `R`: The estimated and regularized LD correlation matrix among the
  selected instruments.

## Step 3: Handling Sample Overlap (`EstRhofun`)

Causal effect estimation can be biased if the exposure dataset and
outcome dataset share overlapping participants. MERLIN explicitly
addresses this via sample correlation parameters `rho1` (for GWAS) and
`rho2` (for GWIS).

**Independent vs. Overlapping Designs**

- **Two-Sample MR (Independent Samples)**: If the exposure and outcome
  summary statistics come from completely separate cohorts, correlations
  due to sample overlap do not exist. You must set both parameters
  directly to zero:

``` r

rho1 <- 0
rho2 <- 0
```

- **One-Sample / Overlapping MR**: If the cohorts overlap, \rho_1 and
  \rho_2 are estimated using the `EstRhofun` function on summary
  statistics among independent variants following [Chen et al
  (2022)](https://www.nature.com/articles/s41467-022-34164-1).

**Parameter Specifications for `EstRhofun`**:

- ld_r2_thresh: Strict r^2 threshold to filter completely independent
  SNPs across the genome.

- lambad: Shrinkage tuning parameter for the LD estimator matrix during
  correlation calculation.

- pth: Critical value threshold adapted to the truncated normal
  distribution in the estimation routine.

``` r

# Set parameters for correlation estimation
ld_r2_thresh <- 0.001
lambad       <- 0.85
pth          <- 1.96

# Estimate correlation for GWAS summary statistics
RhoEst1 <- EstRhofun(expgwas, outgwas, stringname3, ld_r2_thresh, lambad, pth)
rho1    <- mean(RhoEst1$Rhores)

# Estimate correlation for GWIS summary statistics 
# (using exposure and outcome GWIS or relevant null datasets)
RhoEst2 <- EstRhofun(expgwis, outgwis, stringname3, ld_r2_thresh, lambad, pth)
rho2    <- mean(RhoEst2$Rhores)
```

With the matched data, screened instrumental variables, estimated LD
matrix `R`, and sample overlap parameters `rho1`/`rho2` successfully
prepared, the data pipeline is fully complete and ready for statistical
modeling.
