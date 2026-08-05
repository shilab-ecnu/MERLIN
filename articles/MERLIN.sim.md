# Part 3: Model Fitting Using Simulated Data

## Introduction to the Simulation Framework

Before applying MERLIN to real-world summary statistics, it is highly
instructive to understand its behavior using simulated data. In this
section, we will simulate a scenario where an environmental factor $`E`$
interacts with genetic variants $`G`$ to influence both the exposure
$`X`$ and the outcome $`Y`$.

The underlying data generating mechanism can be defined as:
``` math
 X = G\gamma_1 + (G \times E)\gamma_3 + \epsilon_X 
```
``` math
Y = (\beta_1 + \beta_4 E)X + \epsilon_Y
```
Where: - $`\gamma_1`$: Main genetic effects on the exposure.

- $`\gamma_3`$: Genetic-environmental interaction effects.

- $`\beta_1`$: The average (main) causal effect of $`X`$ on $`Y`$.

- $`\beta_4`$: The heterogeneity causal effect (the extent to which the
  causal effect is modified by the environment).

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

We simulate $`m = 1000`$ independent SNPs for a total of 160,000
individuals, split equally into an exposure cohort (`n_exp = 80000`) and
an outcome cohort (`n_out = 80000`). The environmental variable $`E`$ is
simulated from a standard normal distribution.

``` r

n_exp <- 80000
n_out <- 80000
m <- 1000

# True Causal Parameters
b1 <- 0      
b4 <- 0.3    

# Generate Genotype matrix G
maf <- runif(m, 0.05, 0.5)
G <- matrix(rbinom((n_exp + n_out) * m, 2, rep(maf, each = n_exp + n_out)),
            nrow = n_exp + n_out, ncol = m)
G <- scale(G, center = TRUE, scale = FALSE)

# Generate Environmental variable E
E_x <- rnorm(n_exp + n_out)
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
cor_g1g3 <- 0.4  

sigma2g1 <- h_g1 / m
sigma2g3 <- h_g3 / m

cov_matrix <- matrix(c(sigma2g1, 
                      cor_g1g3 * sqrt(sigma2g1 * sigma2g3),
                      cor_g1g3 * sqrt(sigma2g1 * sigma2g3),
                      sigma2g3), ncol = 2)
                      
gamma1_3 <- rmvnorm(m, mean = c(0, 0), sigma = cov_matrix)
gamma_1x <- gamma1_3[, 1]
gamma_3x <- gamma1_3[, 2]
```

**Generating Exposure and Outcome Phenotypes**

Using the mathematical model defined above, we generate the phenotypes
$`X`$ and $`Y`$, and then split the data into independent exposure and
outcome cohorts to mimic a two-sample Mendelian Randomization design.

``` r

# Construct the GxE interaction term
GE <- G * E_x

# Define the correlation between residual errors
rhoxy <- 0.6  

# Calculate variances for the noise terms
var_noise_x <- 1 - h_g1 - h_g3
var_noise_y <- 1

# Construct the covariance matrix for the bivariate normal distribution
cov_xy <- rhoxy * sqrt(var_noise_x * var_noise_y)
sigma_noise <- matrix(c(var_noise_x, cov_xy,
                        cov_xy, var_noise_y), ncol = 2)

# Generate correlated noise for X and Y
noise_xy <- rmvnorm(n_exp + n_out, mean = c(0, 0), sigma = sigma_noise)
noise_x <- noise_xy[, 1]
noise_y <- noise_xy[, 2]

# Generate Exposure (X)
X <- G %*% gamma_1x + GE %*% gamma_3x + noise_x

# Generate Outcome (Y) 
causal_effect <- b1 + b4 * E_x
Y <- X * causal_effect + noise_y

# Split into Exposure and Outcome datasets
exp_gwas <- X[1:n_exp]
exp_gwis <- X[1:n_exp]
exp_E    <- E_x[1:n_exp]

out_gwas <- Y[(n_exp + 1):(n_exp + n_out)]
out_gwis <- Y[(n_exp + 1):(n_exp + n_out)]
out_E    <- E_x[(n_exp + 1):(n_exp + n_out)]
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
(e.g., $`P < 0.01`$ for this toy example). SNPs that pass the threshold
in either the exposure GWAS or GWIS analysis are included in the union
set of instruments. Since SNPs are simulated independently here, the LD
correlation matrix $`R`$ is simply an identity matrix.

``` r

p_threshold <- 0.01

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

# Run MERLIN
res <- MERLIN(gamma_hat, gamma3_hat, Gamma_hat, Gamma3_hat,
              se_gamma, se_gamma3, se_Gamma, se_Gamma3, R, rho_1, rho_2)

str(res)
```

The console output below shows the progress and results of MERLIN.

``` text
Running MERLIN method (Model: standard)...

--------------------------------------------------
MERLIN Analysis Completed Successfully!
Total processing time: 408.1 seconds
--------------------------------------------------
List of 14
 $ Beta1.hat : num 0.0199
 $ Beta1.se  : num 0.0131
 $ Beta1.pval: num 0.129
 $ Beta4.hat : num 0.281
 $ Beta4.se  : num 0.0114
 $ Beta4.pval: num 6.91e-133
 $ gamma1    : num [1:571, 1] 0.0312 0.0153 0.0166 -0.0135 -0.0162 ...
 $ beta2     : num [1:571, 1] 0.00056 0.000281 0.000415 0.000604 0.005952 ...
 $ gamma3    : num [1:571, 1] 0.02122 -0.00824 -0.01623 -0.0089 0.00434 ...
 $ Beta1res  : num [1:1200, 1] 0.0147 0.0281 0.0283 0.012 0.0208 ...
 $ Beta4res  : num [1:1200, 1] 0.286 0.284 0.293 0.284 0.283 ...
 $ Sg12Res   : num [1:1200, 1] 0.000465 0.000461 0.000475 0.000451 0.000452 ...
 $ Sg22Res   : num [1:1200, 1] 2.20e-05 2.13e-05 1.97e-05 2.14e-05 2.16e-05 ...
 $ Sg32Res   : num [1:1200, 1] 0.000171 0.000184 0.000195 0.000167 0.000179 ...
```

**Interpreting the Results**

The output `res` is a list containing the estimated parameters. We can
extract the main causal effect ($`\hat{\beta}_1`$) and the heterogeneity
causal effect ($`\hat{\beta}_4`$).

``` r

# Main causal effect estimates
beta1_hat <- res$Beta1.hat
se1_hat   <- res$Beta1.se
pval1     <- res$Beta1.pval

# Heterogeneity (Interaction) causal effect estimates
beta4_hat <- res$Beta4.hat
se4_hat   <- res$Beta4.se
pval4     <- res$Beta4.pval

cat("Estimated Main Effect (b1):", round(beta1_hat, 4), "\n")
cat("Estimated Interaction Effect (b4):", round(beta4_hat, 4), "\n")
```

``` text
Estimated Main Effect (b1): 0.0199 
Estimated Interaction Effect (b4): 0.281
```

Crucially, recall that we deliberately introduced strong unmeasured
confounding into our simulation by correlating the residual error terms
of the exposure and outcome ($`\rho = 0.6`$). Despite this severe
confounding, you should observe that `beta1_hat` correctly remains
statistically non-significant (close to the true `b1 = 0`), while
`beta4_hat` is highly significant and accurately recovers the true
underlying interaction effect size of approximately 0.3. This
demonstrates MERLIN’s robust capability to estimate unbiased
environment-dependent causal effects even in the presence of strong
unmeasured confounders.
