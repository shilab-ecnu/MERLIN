# Part 5: Model Fitting with a Discrete Binary Modifier (\`discrete_E\` and \`discrete_E_adj\`)

## Introduction to the Binary-Modifier Models

When the modifier E is binary, its mean and variance depend on the
proportion of individuals in the two groups. Consequently, GWAS and GWIS
summary statistics obtained from two cohorts may not be directly
comparable if the distribution of E differs between them. MERLIN
provides two model options for binary modifiers coded as E=-1 and E=1:

- **`discrete_E`** is used when the exposure and outcome cohorts have
  the same value of P(E=1). The common proportion does not need to be
  0.5.

- **`discrete_E_adj`** is used when P(E=1) differs between the exposure
  and outcome cohorts. This model harmonizes the exposure summary
  statistics to the modifier distribution in the outcome cohort before
  model fitting.

Both models require exposure GWAS, exposure GWIS, outcome GWAS, and
outcome GWIS summary statistics.

| Model option       | Exposure and outcome proportions | Additional arguments |
|:-------------------|:---------------------------------|:---------------------|
| `"discrete_E"`     | p\_{exp}=p\_{out}                | `p_common`           |
| `"discrete_E_adj"` | p\_{exp}\ne p\_{out}             | `p_exp`, `p_out`     |

Here, `p_common`, `p_exp`, and `p_out` always denote P(E=1) under the
-1/+1 coding. For example, if males are coded as 1, these quantities are
the proportions of males in the corresponding cohorts.

In this tutorial, a single simulated dataset with a balanced binary
modifier is used to generate G, X, and Y. The two scenarios differ only
in the indices used to construct the final exposure and outcome samples
for the marginal GWAS and GWIS regressions.

## Common Data Generation

We first generate all quantities shared by the two analyses. No
scenario-specific sampling is performed in this section.

``` r

library(MERLIN)
library(mvtnorm)
set.seed(2028)
```

### Common Parameters

Following the simulation code used in the study, `samp_num` specifies
the final sample size of each exposure or outcome analysis. The complete
genotype data contain enough individuals from both groups to construct
the required non-overlapping samples.

``` r

# Final sample size of each exposure or outcome analysis
samp_num <- 80000
m <- 100

# True causal parameters
b1 <- 0
b4 <- 0.3

# Variance parameters
h_g1 <- 0.3
h_g3 <- 0.1
h_b <- 0.05
cor_g1g3 <- 0

# Modifier proportions used in the two analyses
p_common <- 0.4
p_exp <- 0.4
p_out <- 0.6

# Instrument-selection threshold and sample correlations
p_threshold <- 5e-8
rho_1 <- 0
rho_2 <- 0
```

### Genotypes, Modifier, and Genetic Effects

The modifier is balanced, with half of the individuals coded as E=1 and
half as E=-1. The modifier is generated before X and Y and is used in
both phenotype models.

``` r

# Generate the complete genotype data
maf <- runif(m, 0.05, 0.5)
G <- matrix(
  rbinom(200000 * m, 2, rep(maf, each = 200000)),
  nrow = 200000,
  ncol = m
)
G <- scale(G, center = TRUE, scale = FALSE)

# Number of individuals available in the complete genotype data
Gcau_r <- nrow(G)

# Generate a balanced binary modifier
E <- sample(c(rep(1, Gcau_r / 2), rep(-1, Gcau_r / 2)))

stopifnot(mean(E == 1) == 0.5)

# Simulate genetic effects
sigma2g1 <- h_g1 / m
sigma2g3 <- h_g3 / m
sigma2b <- h_b / m

cov_matrix <- matrix(
  c(sigma2g1,
    cor_g1g3 * sqrt(sigma2g1 * sigma2g3),
    cor_g1g3 * sqrt(sigma2g1 * sigma2g3),
    sigma2g3),
  ncol = 2
)

gamma1_3 <- rmvnorm(m, mean = c(0, 0), sigma = cov_matrix)
gamma_1x <- gamma1_3[, 1]
gamma_3x <- gamma1_3[, 2]
beta_2 <- rnorm(m, 0, sqrt(sigma2b))
```

### Generating the Exposure and Outcome

The exposure and outcome are generated once for the entire source
population. Both subsequent analyses use these same values of G, E, X,
and Y.

``` r

GE <- G * E

# Correlation between residual errors represents unmeasured confounding
rhoxy <- 0.6
var_noise_x <- 1 - h_g1 - h_g3
var_noise_y <- 1 - b1 * b1 - h_b - b4 * b4
cov_xy <- rhoxy * sqrt(var_noise_x * var_noise_y)

sigma_noise <- matrix(
  c(var_noise_x, cov_xy,
    cov_xy, var_noise_y),
  ncol = 2
)

noise_xy <- rmvnorm(
  Gcau_r,
  mean = c(0, 0),
  sigma = sigma_noise
)

X <- G %*% gamma_1x + GE %*% gamma_3x + noise_xy[, 1]
Y <- X * b1 + G %*% beta_2 + X * E * b4 + noise_xy[, 2]

idx_E1 <- which(E == 1)
idx_Em1 <- which(E == -1)
```

### Function for Obtaining Summary Statistics

For a GWAS, the phenotype is regressed on each SNP. For a GWIS, the
coefficient and standard error of the G_j\times E term are extracted
from a model that also contains the main effects of G_j and E.

``` r

get_sumstats <- function(G, pheno, interaction = FALSE, E = NULL) {
  betas <- numeric(ncol(G))
  ses <- numeric(ncol(G))

  for (i in seq_len(ncol(G))) {
    if (interaction) {
      model <- lm(pheno ~ G[, i] + E + G[, i]:E)
      betas[i] <- coef(model)[4]
      ses[i] <- summary(model)$coefficients[4, 2]
    } else {
      model <- lm(pheno ~ G[, i])
      betas[i] <- coef(model)[2]
      ses[i] <- summary(model)$coefficients[2, 2]
    }
  }

  list(beta = betas, se = ses)
}
```

## Scenario 1: Matched Proportions (`discrete_E`)

For `discrete_E`, both final cohorts contain 32,000 individuals with E=1
and 48,000 with E=-1, giving P(E=1)=0.4 in both cohorts.

### Selecting the Exposure and Outcome Samples

The exposure and outcome indices are selected without replacement from
the complete data and are non-overlapping within this analysis.

``` r

n_exp_E1_common <- as.integer(samp_num * p_common)
n_exp_Em1_common <- samp_num - n_exp_E1_common
n_out_E1_common <- as.integer(samp_num * p_common)
n_out_Em1_common <- samp_num - n_out_E1_common

idx_exp_common <- c(
  idx_E1[seq_len(n_exp_E1_common)],
  idx_Em1[seq_len(n_exp_Em1_common)]
)

idx_out_common <- c(
  idx_E1[seq.int(n_exp_E1_common + 1,
                 n_exp_E1_common + n_out_E1_common)],
  idx_Em1[seq.int(n_exp_Em1_common + 1,
                  n_exp_Em1_common + n_out_Em1_common)]
)

stopifnot(length(idx_exp_common) == samp_num)
stopifnot(length(idx_out_common) == samp_num)
stopifnot(length(intersect(idx_exp_common, idx_out_common)) == 0)

exp_E_common <- E[idx_exp_common]
out_E_common <- E[idx_out_common]

c(
  exposure_n = length(idx_exp_common),
  outcome_n = length(idx_out_common),
  exposure_p = mean(exp_E_common == 1),
  outcome_p = mean(out_E_common == 1)
)
```

``` text
exposure_n  outcome_n exposure_p  outcome_p 
     80000      80000        0.4        0.4
```

### Obtaining Summary Statistics

The marginal regressions use the selected 80,000 exposure and 80,000
outcome individuals. The underlying G, E, X, and Y were generated in the
common section above.

``` r

exp_gwas_common <- get_sumstats(
  G[idx_exp_common, ], X[idx_exp_common]
)
exp_gwis_common <- get_sumstats(
  G[idx_exp_common, ], X[idx_exp_common],
  interaction = TRUE, E = E[idx_exp_common]
)
out_gwas_common <- get_sumstats(
  G[idx_out_common, ], Y[idx_out_common]
)
out_gwis_common <- get_sumstats(
  G[idx_out_common, ], Y[idx_out_common],
  interaction = TRUE, E = E[idx_out_common]
)
```

### Instrument Selection and Model Fitting

The SNPs are independent in this example. Instruments are selected from
the union of variants associated with the exposure in either the GWAS or
GWIS at P\<5\times10^{-8}.

``` r

pvals_gwas_common <- 2 * pnorm(
  -abs(exp_gwas_common$beta / exp_gwas_common$se)
)
pvals_gwis_common <- 2 * pnorm(
  -abs(exp_gwis_common$beta / exp_gwis_common$se)
)

iv_gwas_common <- which(pvals_gwas_common < p_threshold)
iv_gwis_common <- which(pvals_gwis_common < p_threshold)
iv_common <- union(iv_gwas_common, iv_gwis_common)
R_common <- diag(length(iv_common))

fit_discrete <- MERLIN(
  gammah1 = exp_gwas_common$beta[iv_common],
  gammah3 = exp_gwis_common$beta[iv_common],
  Gammah1 = out_gwas_common$beta[iv_common],
  Gammah3 = out_gwis_common$beta[iv_common],
  se1 = exp_gwas_common$se[iv_common],
  se2 = exp_gwis_common$se[iv_common],
  se3 = out_gwas_common$se[iv_common],
  se4 = out_gwis_common$se[iv_common],
  R = R_common,
  rho_1 = rho_1,
  rho_2 = rho_2,
  model = "discrete_E",
  p_common = p_common,
  seed = 2028
)

str(fit_discrete)
```

``` text
List of 6
 $ Beta1.hat : num -0.00379
 $ Beta1.se  : num 0.0166
 $ Beta1.pval: num 0.82
 $ Beta4.hat : num 0.307
 $ Beta4.se  : num 0.0101
 $ Beta4.pval: num 1.53e-203
```

The model automatically adjust the inputs using `p_common`. The returned
estimates are then transformed back to the original E=-1/+1 scale.

## Scenario 2: Different Proportions (`discrete_E_adj`)

The exposure sample again contains 32,000 individuals with E=1 and
48,000 with E=-1, whereas the outcome sample contains 48,000 individuals
with E=1 and 32,000 with E=-1. Thus, p\_{exp}=0.4 and p\_{out}=0.6.

### Selecting the Exposure and Outcome Samples

``` r

n_exp_E1_adj <- as.integer(samp_num * p_exp)
n_exp_Em1_adj <- samp_num - n_exp_E1_adj
n_out_E1_adj <- as.integer(samp_num * p_out)
n_out_Em1_adj <- samp_num - n_out_E1_adj

idx_exp_adj <- c(
  idx_E1[seq_len(n_exp_E1_adj)],
  idx_Em1[seq_len(n_exp_Em1_adj)]
)

idx_out_adj <- c(
  idx_E1[seq.int(n_exp_E1_adj + 1,
                 n_exp_E1_adj + n_out_E1_adj)],
  idx_Em1[seq.int(n_exp_Em1_adj + 1,
                  n_exp_Em1_adj + n_out_Em1_adj)]
)

stopifnot(length(idx_exp_adj) == samp_num)
stopifnot(length(idx_out_adj) == samp_num)
stopifnot(length(intersect(idx_exp_adj, idx_out_adj)) == 0)

exp_E_adj <- E[idx_exp_adj]
out_E_adj <- E[idx_out_adj]

c(
  exposure_n = length(idx_exp_adj),
  outcome_n = length(idx_out_adj),
  exposure_p = mean(exp_E_adj == 1),
  outcome_p = mean(out_E_adj == 1)
)
```

``` text
exposure_n  outcome_n exposure_p  outcome_p 
     80000      80000        0.4        0.6
```

### Obtaining Summary Statistics

The exposure and outcome summary statistics now differ from the first
scenario only because different individual indices are used.

``` r

exp_gwas_adj <- get_sumstats(
  G[idx_exp_adj, ], X[idx_exp_adj]
)
exp_gwis_adj <- get_sumstats(
  G[idx_exp_adj, ], X[idx_exp_adj],
  interaction = TRUE, E = E[idx_exp_adj]
)
out_gwas_adj <- get_sumstats(
  G[idx_out_adj, ], Y[idx_out_adj]
)
out_gwis_adj <- get_sumstats(
  G[idx_out_adj, ], Y[idx_out_adj],
  interaction = TRUE, E = E[idx_out_adj]
)
```

The GWAS and GWIS coefficients remain on the original E=-1/+1 scale. The
high-level
[`MERLIN()`](https://shilab-ecnu.github.io/MERLIN/reference/MERLIN-package.md)
function performs the required standardization and between-cohort
harmonization. The inputs should not be manually transformed before
model fitting.

### Instrument Selection and Model Fitting

``` r

pvals_gwas_adj <- 2 * pnorm(-abs(exp_gwas_adj$beta / exp_gwas_adj$se))
pvals_gwis_adj <- 2 * pnorm(-abs(exp_gwis_adj$beta / exp_gwis_adj$se))

iv_gwas_adj <- which(pvals_gwas_adj < p_threshold)
iv_gwis_adj <- which(pvals_gwis_adj < p_threshold)
iv_adj <- union(iv_gwas_adj, iv_gwis_adj)
R_adj <- diag(length(iv_adj))

fit_discrete_adj <- MERLIN(
  gammah1 = exp_gwas_adj$beta[iv_adj],
  gammah3 = exp_gwis_adj$beta[iv_adj],
  Gammah1 = out_gwas_adj$beta[iv_adj],
  Gammah3 = out_gwis_adj$beta[iv_adj],
  se1 = exp_gwas_adj$se[iv_adj],
  se2 = exp_gwis_adj$se[iv_adj],
  se3 = out_gwas_adj$se[iv_adj],
  se4 = out_gwis_adj$se[iv_adj],
  R = R_adj,
  rho_1 = rho_1,
  rho_2 = rho_2,
  model = "discrete_E_adj",
  p_exp = p_exp,
  p_out = p_out,
  seed = 2028
)

str(fit_discrete_adj)
```

``` text
List of 6
 $ Beta1.hat : num -0.0129
 $ Beta1.se  : num 0.0174
 $ Beta1.pval: num 0.458
 $ Beta4.hat : num 0.307
 $ Beta4.se  : num 0.0108
 $ Beta4.pval: num 6.6e-178
```

For this model, MERLIN first projects the exposure GWAS associations
from the exposure modifier mean to the outcome modifier mean. It then
expresses both GWIS inputs on the standardized modifier scale defined by
`p_out`. The final estimates are returned on the original E=-1/+1 scale.
These operations are implemented inside
[`MERLIN()`](https://shilab-ecnu.github.io/MERLIN/reference/MERLIN-package.md)
already and should not be repeated by the user.

The current implementation uses

\operatorname{Var}(\widehat{\gamma}\_{1,\mathrm{out}}) =
\operatorname{SE}(\widehat{\gamma}\_1)^2 +(\mu\_{out}-\mu\_{exp})^2
\operatorname{SE}(\widehat{\gamma}\_3)^2,

which treats the covariance between the exposure GWAS and exposure GWIS
association estimates as zero.

## Interpretation on the Original Binary-Modifier Scale

For both model options, the returned estimates correspond to the causal
effect model

\beta_1 + \beta_4 E,

with E coded as -1 or 1. The group-specific causal effects are therefore

\text{Causal effect for }E=-1=\beta_1-\beta_4,

and

\text{Causal effect for }E=1=\beta_1+\beta_4.

If `p_exp` and `p_out` are equal, `discrete_E` should be used with
`p_common` instead of `discrete_E_adj`.
