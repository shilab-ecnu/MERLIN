# Part 3: Model Fitting Using Simulated Data

## Introduction to the Simulation Framework

Before applying MERLIN to real-world summary statistics, it is highly
instructive to understand its behavior using simulated data. In this
section, we will simulate a scenario where an environmental factor $`E`$
interacts with genetic variants $`G`$ to influence both the exposure
$`X`$ and the outcome $`Y`$.

The underlying data generating mechanism can be defined as:
``` math
X=G\gamma_1+(G\times E)\gamma_3+0.1E+\epsilon_X,
```

``` math
Y=(\beta_1+\beta_4E)X+G\beta_2+0.1E+\epsilon_Y.
```
Where: - $`\beta_1`$: The average causal effect of $`X`$ on $`Y`$,
corresponding to $`\beta^{(A)}`$ in the manuscript.

- $`\beta_4`$: The causal interaction effect quantifying the extent to
  which the causal effect of $`X`$ on $`Y`$ is modified by $`E`$,
  corresponding to $`\beta^{(I)}`$ in the manuscript.

- $`\gamma_1`$: The main genetic effects on the exposure, corresponding
  to $`\boldsymbol{\gamma}^{(G)}`$ in the manuscript.

- $`\gamma_3`$: The genetic-by-environment interaction effects on the
  exposure, corresponding to $`\boldsymbol{\gamma}^{(GI)}`$ in the
  manuscript.

- $`\beta_2`$: The direct genetic effects on the outcome, representing
  horizontal pleiotropy and corresponding to
  $`{\boldsymbol{\beta}}^{(G)}`$ in the manuscript.

- $`\epsilon_X`$, $`\epsilon_Y`$: The residual error terms for the
  exposure and outcome, respectively. These capture the cumulative
  errors, crucially including the effects of unmeasured confounders,
  independent environmental noise, and measurement errors.

## Step-by-Step Simulation

First, we load the required packages and set a random seed for
reproducibility.

``` r

library(MERLIN)
library(mvtnorm)
library(MASS)
```

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:MERLIN':
    ## 
    ##     select

``` r

set.seed(2026)
```

**Generating Genotypes and Environmental Variables**

We simulate $`m = 100`$ independent SNPs for a total of 160,000
individuals, split equally into an exposure cohort (`n_exp = 80000`) and
an outcome cohort (`n_out = 80000`). The environmental variable $`E`$ is
simulated as a balanced binary modifier coded as $`-1`$ and $`1`$, with
equal proportions in both the exposure and outcome samples.

``` r

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

**Simulating Genetic Effects**

The main genetic effects ($`\gamma_1`$) and $`G \times E`$ interaction
effects ($`\gamma_3`$) are generated from a bivariate normal
distribution to allow for correlation (`cor_g1g3 = 0.4`) between them.
We control their variance using specified heritability parameters
(`h_g1` and `h_g3`).

``` r

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

**Generating Exposure and Outcome Phenotypes**

Using the mathematical model defined above, we generate the phenotypes
$`X`$ and $`Y`$, and then split the data into independent exposure and
outcome cohorts to mimic a two-sample Mendelian Randomization design.

``` r

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

## Obtaining Summary Statistics

In practice, researchers usually only have access to summary statistics.
We define a helper function `get_sumstats` to run single-variant linear
regressions and extract the effect sizes ($`\hat{\beta}`$) and standard
errors ($`SE`$).

``` r

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

## Instrument Selection and Model Fitting

**Selecting Instrumental Variables**

We select valid genetic instruments using a significance threshold
(e.g., $`P < 5e-8`$ for this toy example). SNPs that pass the threshold
in either the exposure GWAS or GWIS analysis are included in the union
set of instruments. Since SNPs are simulated independently here, the LD
correlation matrix $`R`$ is simply an identity matrix.

``` r

p_threshold <- 5e-8

# Filter GWAS signals
pvals_gwas <- 2 * pnorm(-abs(exp_gwas_sum$beta / exp_gwas_sum$se))
iv_gwas <- which(pvals_gwas < p_threshold)

# Filter GWIS signals
pvals_gwis <- 2 * pnorm(-abs(exp_gwis_sum$beta / exp_gwis_sum$se))
iv_gwis <- which(pvals_gwis < p_threshold)

# Union of instruments and LD matrix
iv_union <- union(iv_gwas, iv_gwis)
R <- diag(length(iv_union))
```

**Applying the MERLIN Method**

Finally, we format the inputs and execute the MERLIN function. Since we
simulated independent non-overlapping samples, the sample overlap
correlations (rho_1 and rho_2) are explicitly set to $`0`$.

``` r

# Extract values for the selected instruments
gamma_hat  <- exp_gwas_sum$beta[iv_union]
gamma3_hat <- exp_gwis_sum$beta[iv_union]
Gamma_hat  <- out_gwas_sum$beta[iv_union]
Gamma3_hat <- out_gwis_sum$beta[iv_union]

se_gamma  <- exp_gwas_sum$se[iv_union]
se_gamma3 <- exp_gwis_sum$se[iv_union]
se_Gamma  <- out_gwas_sum$se[iv_union]
se_Gamma3 <- out_gwis_sum$se[iv_union]

# Independent samples assumption
rho_1 <- 0
rho_2 <- 0

# Run standard MERLIN
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

str(res)
```

The console output below shows the progress and results of MERLIN.

``` text
Running MERLIN method (Model: standard)...

--------------------------------------------------
MERLIN Analysis Completed Successfully!
Total processing time: 0.9 seconds
--------------------------------------------------
List of 6
 $ Beta1.hat : num -0.00454
 $ Beta1.se  : num 0.0189
 $ Beta1.pval: num 0.81
 $ Beta4.hat : num 0.314
 $ Beta4.se  : num 0.0107
 $ Beta4.pval: num 3.93e-190
```

**Interpreting the Results**

The output `res` is a list containing the estimated parameters. We can
extract the average causal effect ($`\widehat{\beta}_1`$) and the
heterogeneity causal effect ($`\widehat{\beta}_4`$).

``` r

# Main causal effect estimates
beta1_hat <- res$Beta1.hat
se1_hat   <- res$Beta1.se
pval1     <- res$Beta1.pval

# Heterogeneity (interaction) causal effect estimates
beta4_hat <- res$Beta4.hat
se4_hat   <- res$Beta4.se
pval4     <- res$Beta4.pval

cat("Estimated Main Effect (b1):", round(beta1_hat, 4), "\n")
cat("Estimated Interaction Effect (b4):", round(beta4_hat, 4), "\n")
```

``` text
Estimated Main Effect (b1): -0.0045
Estimated Interaction Effect (b4): 0.3137
```

Crucially, recall that we deliberately introduced strong unmeasured
confounding into the simulation by correlating the residual error terms
of the exposure and outcome ($`\rho=0.6`$). Despite this confounding,
the estimated average causal effect was close to its true value of zero
($`\widehat{\beta}_1=-0.0045`$, $`P=0.810`$), whereas the estimated
interaction effect closely recovered its true value of 0.3
($`\widehat{\beta}_4=0.3137`$, $`P=3.93\times10^{-190}`$). These results
illustrate that MERLIN can recover both the average and
modifier-dependent causal effects in this simulated setting with a
balanced binary modifier.
