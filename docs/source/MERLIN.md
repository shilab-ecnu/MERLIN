---
title: "MERLIN: MEndelian Randomization for Linear INteraction"
author: "\\small Yadong Yang, Minxi Bai, Jiacheng Miao, Jin Liu, Qiongshi Lv and Xingjie Shi"
date: "\\small \\today"
output:
  pdf_document: default
editor_options:
  markdown:
    wrap: 72
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
options(encoding = "UTF-8")
```



# Installation

```         
library(remotes)
install_github("shilab-ecnu/MERLIN")
```

Load the package using the following command:

```{r}
library(MERLIN)
```

# Fit MERLIN using simulated data
The repeated code in the paper is all in this hyperlink: <a href="https://github.com/shilab-ecnu/MERLIN/tree/main/simulation">SIMULATION</a>.

We first generate the genotype data and the environmental variable:

```{r, eval = FALSE}
library(mvtnorm)
library(MASS)
set.seed(2025)
```

```{r, eval = FALSE}
n_exp <- 80000;
n_out <- 80000;
m <- 1000;
h_g1 <- 0.3;    
h_g3 <- 0.1;   
b1 <- 0;   
b4 <- 0.3;       
cor_g1g3 <- 0.4; 
rhoxy <- 0.6;

maf <- runif(m, 0.05, 0.5);
G <- matrix(rbinom((n_exp + n_out) * m, 2, rep(maf, each = n_exp + n_out)),
            nrow = n_exp + n_out, ncol = m);
G <- scale(G, center = TRUE, scale = FALSE);
E_x <- rnorm(n_exp + n_out)
```

Now simulate the genetic effect sizes. The main genetic effects ($\gamma_1$) and G×E interaction effects ($\gamma_3$) are generated as correlated multivariate normal variables with specified heritabilities.

```{r, eval = FALSE}
sigma2g1 <- h_g1 / m;
sigma2g3 <- h_g3 / m;

cov_matrix <- matrix(c(sigma2g1, 
                      cor_g1g3 * sqrt(sigma2g1 * sigma2g3),
                      cor_g1g3 * sqrt(sigma2g1 * sigma2g3),
                      sigma2g3), ncol = 2);
                      
gamma1_3 <- rmvnorm(m, mean = c(0, 0), sigma = cov_matrix)
gamma_1x <- gamma1_3[, 1];
gamma_3x <- gamma1_3[, 2]
```

Generate the exposure ($X$) and outcome ($Y$) variables with the genetic effects defined.

```{r, eval = FALSE}
GE <- G * E_x;

noise_x <- rnorm(n_exp + n_out, sd = sqrt(1 - h_g1 - h_g3));
X <- G %*% gamma_1x + GE %*% gamma_3x + noise_x;

noise_y <- rnorm(n_exp + n_out, sd = sqrt(1 - b1^2 - b4^2));
Y <- X * b1 + GE %*% gamma_3x * b4 + noise_y;

exp_gwas <- X[1:n_exp];  
exp_gwis <- X[1:n_exp];   
exp_E <- E_x[1:n_exp];     
out_gwas <- Y[(n_exp + 1):(n_exp + n_out)];  
out_gwis <- Y[(n_exp + 1):(n_exp + n_out)];
out_E <- E_x[(n_exp + 1):(n_exp + n_out)]
```

We then conduct single-variant analysis to obtain the summary statistics.

```{r, eval = FALSE}
get_sumstats <- function(G, pheno, interaction = FALSE, E = NULL) {
  betas <- numeric(ncol(G));
  ses <- numeric(ncol(G));
  
  for(i in 1:ncol(G)) {
    if(interaction) {
      # Model with G×E interaction term
      model <- lm(pheno ~ G[, i] + E + G[, i]:E);
      betas[i] <- coef(model)[4];
      ses[i] <- summary(model)$coefficients[4, 2];
    } else {
      # Standard additive genetic model
      model <- lm(pheno ~ G[, i]);
      betas[i] <- coef(model)[2];
      ses[i] <- summary(model)$coefficients[2, 2]
    }
  }
  return(list(beta = betas, se = ses))
}

exp_gwas_sum <- get_sumstats(G[1:n_exp, ], exp_gwas);
exp_gwis_sum <- get_sumstats(G[1:n_exp, ], exp_gwis, interaction = TRUE, E = exp_E);

out_gwas_sum <- get_sumstats(G[(n_exp + 1):(n_exp + n_out), ], out_gwas);
out_gwis_sum <- get_sumstats(G[(n_exp + 1):(n_exp + n_out), ], out_gwis, 
                            interaction = TRUE, E = out_E)
```

We select genetic instruments using a p-value threshold. SNPs in either the GWAS or GWIS analysis are included in the union set of instruments.

```{r, eval = FALSE}
p_threshold <- 0.01;

pvals_gwas <- 2 * pnorm(-abs(exp_gwas_sum$beta / exp_gwas_sum$se));
iv_gwas <- which(pvals_gwas < p_threshold);

pvals_gwis <- 2 * pnorm(-abs(exp_gwis_sum$beta / exp_gwis_sum$se));
iv_gwis <- which(pvals_gwis < p_threshold);

iv_union <- union(iv_gwas, iv_gwas);
R <- diag(length(iv_union)) 
```

Finally, we apply the MERLIN methods. 

```{r, eval = FALSE}
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

res <- MERLIN(gamma_hat, gamma3_hat, Gamma_hat, Gamma3_hat,
            se_gamma, se_gamma3, se_Gamma, se_Gamma3, R, rho_1, rho_2);
str(res);

beta1_hat <- res$Beta1.hat;
se1_hat <- res$Beta1.se;
pval1 <- res$Beta1.pval;
beta4_hat <- res$Beta4.hat;
se4_hat <- res$Beta4.se;
pval4 <- res$Beta4.pval
```

beta1_hat, se1_hat, and pval1 are the estimated average causal effect, corresponding standard error, and p-value of beta1_hat. beta4_hat, se4_hat, and pval4 are the estimated heterogeneity causal effect, corresponding standard error, and p-value of beta4_hat.

# Real data
## The Testosterone-BD study with environmental factor sex

All the raw data for the real-data analyses in the replicated paper are stored on <a href="https://figshare.com/articles/dataset/Data_for_MERLIN/29910116">MERLIN Dataset on Figshare</a>. Here, we take the dataset “The Testosterone–BD study with the environmental factor sex” as an example.

Furthermore, we give an example to illustrate the implementation of MERLIN for real data analysis. The following datasets ('Testosterone.GWAS.txt.gz', 'Testosterone.GWIS.txt.gz', 'BD.GWAS.txt.gz', 'BD.GWIS.txt.gz', 'g1000_eur.bed','g1000_eur.fam', 'g1000_eur.bim', 'all.bed') should be prepared. Download here: <a href="https://figshare.com/articles/dataset/Data_for_MERLIN/29910116">MERLIN Dataset on Figshare</a>

```{r}
expgwas <- "Testosterone.GWAS.txt.gz";
expgwis <- "Testosterone.GWIS.txt.gz";
outgwas <- "BD.GWAS.txt.gz";
outgwis <- "BD.GWIS.txt.gz";
stringname3 <- "g1000_eur";
block_file <- "all.bed";
```

'expgwas', 'expgwis', 'outgwas', and 'outgwis' are the dataset names for exposure GWAS, exposure GWIS, outcome GWAS, and outcome GWIS, respectively. Here, the environment variable is sex.

These four datasets must have the following format (note that it must be tab-delimited): including columns as SNP, CHR, BP, A1, A2, BETA, SE, and P.

```{r, echo=FALSE}
example <- read.table("example.txt", header = TRUE)
```

```{r, echo=FALSE, results='asis'}
knitr::kable(example, digits = 6, caption = "\\label{exp}Data format used for exposure and outcome data.")
```


If GWAS and GWIS data cannot be directly obtained, and the environmental factor is a binary variable (e.g., sex), one can generate the required GWAS and GWIS inputs for MERLIN by converting the sex-stratified summary statistics (e.g., Testosterone.male.txt and Testosterone.female.txt) as follows.

The GWAS summary statistics can be generated by meta-analyzing the male and female data using inverse-variance weighting, as implemented in METAL (<a href="https://github.com/statgen/METAL">https://github.com/statgen/METAL</a>).  After installing the software, the analysis can be executed via the command line (e.g., in Linux or other shell environments) using a configuration file. A sample configuration file 'metal.config.Testosterone.txt' is available for download: <a href="https://figshare.com/articles/dataset/Data_for_MERLIN/29910116">MERLIN Dataset on Figshare</a>.


```         
metal metal.config.Testosterone.txt
```

The SNP effects and standard errors for GWIS summary statistics were derived based on the following formula, assuming sex coded as Male=1, Female=-1. Allele direction must be aligned before analyzing sex-stratified data.

```math
\hat{b}_{gwis,j}=\frac{1}{2}(\hat{b}_{male,j}-{\hat{b}}_{female,j})
```

```math
se(\hat{b}_{gwis,j})=\frac{1}{2}\sqrt{(se(\hat{b}_{male,j})^2+se(\hat{b}_{female,j})^2}
```

'stringname3' is the name of the reference panel data. Here we use samples from '1000 Genomes Project European panel', which is in plink binary format. 'block_file' is used to partition the whole genome into blocks.

The `matchpanel` function is used to match a GWAS/GWIS dataset with the reference panel data, alongside initial data quality control. The output includes a data frame (`$data`) and the corresponding storage path
(`$data_dir`).


```{r, eval = FALSE}
expgwas.match <- matchpanel(expgwas,stringname3)$data_dir;
expgwis.match <- matchpanel(expgwis,stringname3)$data_dir;
outgwas.match <- matchpanel(outgwas,stringname3)$data_dir;
outgwis.match <- matchpanel(outgwis,stringname3)$data_dir;
```

Having given that we have the formatted data, we can use the `ivselect` function to screen the instrumental variables (IVs) and estimate the correlations among those IVs. `plink_dir` specifies the local path to
the PLINK executable; if not provided, PLINK will be automatically downloaded. `pval_cutoff_gwas` and `pval_cutoff_gwis` define the P-value thresholds for the exposure GWAS and GWIS, respectively. `r2_cutoff` and
`kb_cutoff` are used in LD clumping to specify the $r^2$ threshold and the physical distance (in kilobases) between SNPs. `maf_cutoff` sets the threshold for minor allele frequency. `lam` denotes the shrinkage
parameter used in the regularization of the LD matrix. `CoreNum` indicates the number of CPU cores to be used for parallel computation. `intersect_mode` controls whether to merge GWAS and GWIS IVs using intersection (default: union).   

```{r, eval = FALSE}
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

```{r, eval = FALSE}
ivselect.res <- ivselect(expgwas.match, expgwis.match, outgwas.match, outgwis.match,
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

```{r, eval = FALSE}
rho1 <- 0; rho2 <- 0;
```

For overlap samples, Since $\rho_1$ and $\rho_2$ are estimated using summary statistics among independent variants, we select independent SNPs using the clumping algorithm ($r^2$ threshold denoted by `ld_r2_thresh`). `pth` is the critical value adapted to the truncated normal distribution in the estimation procedure. `lambda` is the shrinkage turning parameter for LD estimator.

```         
ld_r2_thresh <- 0.001;
lambad <- 0.85;
pth <- 1.96;
RhoEst1 <- EstRhofun(expgwas, outgwas, stringname3, ld_r2_thresh, lambad, pth);
rho1 <- mean(RhoEst1$Rhores);
RhoEst2 <- EstRhofun(expgwas, outgwas, stringname3, ld_r2_thresh, lambad, pth);
rho2 <- mean(RhoEst2$Rhores);
```

Now we can fit MERLIN using the function *MERLIN*. 

```{r, eval = FALSE}
res <- MERLIN(gammah1, gammah3, Gammah1, Gammah3,
             se1, se2, se3, se4, R, rho1, rho2)
str(res)
```

Check the convergence of the Gibbs sampler using `traceplot`.

```{r, eval = FALSE}
traceplot(res$Beta1res);
traceplot(res$Beta4res);
```

```{r, eval = FALSE}
MERLINbeta1 <- res$Beta1.hat;
MERLINse1 <- res$Beta1.se;
MERLINpvalue1 <- res$Beta1.pval;
MERLINbeta4 <- res$Beta4.hat;
MERLINse4 <- res$Beta4.se;
MERLINpvalue4 <- res$Beta4.pval;
```

```{r results="hide", eval = FALSE}
cat("The estimated main effect of Testosterone on BD: ", MERLINbeta1, 
    "\n with Standard error: ", MERLINse1, "and P-value: ", MERLINpvalue1);
```

```{r results="hide", eval = FALSE}
cat("The estimated interaction effect of Testosterone on BD: ", MERLINbeta4, 
    "\n with Standard error: ", MERLINse4, "and P-value: ", MERLINpvalue4)
```

## Figure 1

######Overlap plot
```{r, eval = FALSE}
library(ggplot2)  
library(dplyr)      

rm(list=ls())
data_bind <- read.table("data1.txt", header = TRUE)

colnames(data_bind) <- c("GEI.beta1", "GEI.beta1.pval", "GEI.beta4", "GEI.beta4.pval",  #MERLIN(p)
                         "GEI3.beta1", "GEI3.beta1.pval", "GEI3.beta4", "GEI3.beta4.pval",  #MERLIN
                         "LDP.beta1", "LDP.beta1.pval",  #LDP method
                         "raps(T).beta1", "raps(T).beta1.pval",  #MR.raps in twosamplemr
                         "IVW.beta1", "IVW.beta1.pval",  #MR.IVW in twosamplemr
                         "Egger(T).beta1", "Egger(T).beta1.pval",   #MR.Egger in twosamplemr
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
source("qqplot.r")
##qqplot
dat_method_b1 <- data_bind[, c("GEI3.beta1.pval", "LDP.beta1.pval", 
                               "raps(T).beta1.pval", "IVW.beta1.pval", "Egger(T).beta1.pval")]
colnames(dat_method_b1) <- c("MERLIN", "MR-LDP", "RAPS", "IVW", "MR-Egger")
col <- ncol(dat_method_b1)
dat_hg3hb2 <- data_bind[, c("cor_g1g3", "h_g3", "h_b2")]
dat_hg3hb2 <- do.call(rbind, lapply(1:col, function(i) dat_hg3hb2))
dat_boxplot_b1 <- data.frame(value = unlist(dat_method_b1), method = rep(names(dat_method_b1), each = nrow(dat_method_b1)))
dat_boxplot_b1$method <- factor(dat_boxplot_b1$method, 
                                levels = c("MERLIN", "MR-LDP", "RAPS", "IVW", "MR-Egger"))
dat_hg3hb2[, 1] <- factor(dat_hg3hb2[, 1])
dat_hg3hb2[, 2] <- factor(dat_hg3hb2[, 2])
dat_hg3hb2[, 3] <- factor(dat_hg3hb2[, 3])
dat <- cbind(dat_boxplot_b1, dat_hg3hb2)
qq_df <- dat %>%
  group_by(method, h_g3) %>%
  arrange(value, .by_group = TRUE) %>%
  mutate(
    n = n(),
    observed = -log10(value),
    expected = -log10(ppoints(n())),
    lower = -log10(qbeta(0.025, rank(value), n - rank(value) + 1)),
    upper = -log10(qbeta(0.975, rank(value), n - rank(value) + 1))
  )

qq_df$method <- factor(qq_df$method, levels = unique(qq_df$method))
updated_morandi_colors <- c("#F5BE8F", "#C1E0DB", "#CCD376", "#A28CC2", "#8498AB")
p <- ggplot(qq_df, aes(x = expected, y = observed, color = method)) +
  theme_bw()+ 
  geom_ribbon(
  aes(x = expected, ymin = lower, ymax = upper),
  fill = "gray90", alpha = 0.5, inherit.aes = TRUE)+
  geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dashed") +
  geom_point(size = 2, alpha = 1) +
  facet_grid(cols = vars(h_g3), 
                label = "label_parsed", 
                scales="fixed") +
  scale_color_manual(values = updated_morandi_colors) +
  ggtitle(" ") +
  theme(legend.text = element_text(size = 30, face = "bold"),
           legend.position = "none", 
           legend.key = element_blank() )+   
  labs(x = "Expected -log10(p)",
       y = "Observed -log10(p)",
       color = "Method") +                                                  
  theme(axis.text.x = element_text(size=30,face = "bold"),
        axis.text.y = element_text(size=30,face = "bold"),
        axis.title.x=element_text(size=30,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        strip.text.x = element_text(size = 30, face = "bold", colour = "black", angle = 0),
        strip.text.y = element_text(size = 30, face = "bold", colour = "black", angle = -90)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 15, face = "bold", color = "black"))+
  guides(fill=guide_legend(title=NULL)) +  
  theme(panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "grey", fill = NA, size = 0.5), 
        strip.background = element_rect(color = "grey", size = 0.3)) 
```




## Figure 2

```{r, eval = FALSE}

######Figure 1

library(ggplot2)    
library(dplyr)      

rm(list=ls())
data_bind <- read.table("data1.txt", header = TRUE)

colnames(data_bind) <- c("GEI.beta1", "GEI.beta1.pval", "GEI.beta4", "GEI.beta4.pval",  #MERLIN(p)
                         "GEI3.beta1", "GEI3.beta1.pval", "GEI3.beta4", "GEI3.beta4.pval",  #MERLIN
                         "LDP.beta1", "LDP.beta1.pval",  #LDP method
                         "raps(T).beta1", "raps(T).beta1.pval",  #MR.raps in twosamplemr
                         "IVW.beta1", "IVW.beta1.pval",  #MR.IVW in twosamplemr
                         "Egger(T).beta1", "Egger(T).beta1.pval",   #MR.Egger in twosamplemr
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
##Boxplot
dat_method_b1 <- data_bind[, c("GEI3.beta1", "LDP.beta1", "raps(T).beta1", "IVW.beta1", "Egger(T).beta1")]
colnames(dat_method_b1) <- c("MERLIN", "MR-LDP", "RAPS", "IVW", "MR-Egger")
col <- ncol(dat_method_b1)
dat_hg3hb2 <- data_bind[, c("cor_g1g3", "h_g3", "h_b2")]
dat_hg3hb2 <- do.call(rbind, lapply(1:col, function(i) dat_hg3hb2))
dat_boxplot_b1 <- data.frame(value = unlist(dat_method_b1), method = rep(names(dat_method_b1), each = nrow(dat_method_b1)))
dat_boxplot_b1$method <- factor(dat_boxplot_b1$method, 
                                levels = c("MERLIN", "MR-LDP", "RAPS", "IVW", "MR-Egger"))
dat_hg3hb2[, 1] <- factor(dat_hg3hb2[, 1])
dat_hg3hb2[, 2] <- factor(dat_hg3hb2[, 2])
dat_hg3hb2[, 3] <- factor(dat_hg3hb2[, 3])
dat <- cbind(dat_boxplot_b1, dat_hg3hb2)

updated_morandi_colors <- c("#F5BE8F", "#C1E0DB", "#CCD376", "#A28CC2", "#8498AB")
p <- ggplot(dat, aes(x = cor_g1g3, y = value, fill = method)) +
     geom_boxplot(position = position_dodge(width = 0.8)) +
     theme_bw()+
     facet_grid(rows = vars(h_b2), 
                cols = vars(h_g3), 
                label = "label_parsed", 
                scales="fixed") +
     scale_fill_manual(values = updated_morandi_colors) +
     geom_hline(aes(yintercept=b1), size = 1.3, linetype=5,col="red") +
     ggtitle(" ") +
     theme(legend.text = element_text(size = 30, face = "bold"),
           legend.position = "none")+                                                        
     theme(axis.text.x = element_text(size=30,face = "bold"),
           axis.text.y = element_text(size=30,face = "bold"),
           axis.title.x=element_text(size=30,face = "bold"),
           axis.title.y=element_text(size=30,face = "bold"),
           strip.text.x = element_text(size = 30, face = "bold", colour = "black", angle = 0),
           strip.text.y = element_text(size = 30, face = "bold", colour = "black", angle = -90)) +
     theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 15, face = "bold", color = "black"))+
     guides(fill=guide_legend(title=NULL)) +  
     theme(panel.spacing=unit(.05, "lines"),
           panel.border = element_rect(color = "grey", fill = NA, size = 0.5), 
           strip.background = element_rect(color = "grey", size = 0.3)) 
         
##Type 1 error
dat_method_b1_p <- data_bind[, c("GEI3.beta1.pval",  
                                 "LDP.beta1.pval", "raps(T).beta1.pval", 
                                 "IVW.beta1.pval", "Egger(T).beta1.pval")]
colnames(dat_method_b1_p) <- c("MERLIN", "MR-LDP", "RAPS", "IVW", "MR-Egger")
col <- ncol(dat_method_b1_p)
dat_hg3hb2 <- data_bind[, c("cor_g1g3", "h_g3", "h_b2")]
dat_hg3hb2 <- do.call(rbind, lapply(1:col, function(i) dat_hg3hb2))
dat_boxplot_b1_p <- data.frame(value = unlist(dat_method_b1_p), method = rep(names(dat_method_b1_p), each = nrow(dat_method_b1_p)))
dat_boxplot_b1_p$method <- factor(dat_boxplot_b1_p$method, 
                                levels = c("MERLIN", "MR-LDP", "RAPS", "IVW", "MR-Egger"))
dat_hg3hb2[, 1] <- factor(dat_hg3hb2[, 1])
dat_hg3hb2[, 2] <- factor(dat_hg3hb2[, 2])
dat_hg3hb2[, 3] <- factor(dat_hg3hb2[, 3])
dat <- cbind(dat_boxplot_b1_p, dat_hg3hb2)
proportion <- dat %>%
  group_by(method, cor_g1g3, h_g3, h_b2) %>%
  summarise(prop_beta_pval = sum(value < 0.05) / n(), .groups = "keep")

updated_morandi_colors <- c("#F5BE8F", "#C1E0DB", "#CCD376", "#A28CC2", "#8498AB")
p <- ggplot(proportion, aes(x = cor_g1g3, y = prop_beta_pval, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw()+
  facet_grid(rows = vars(h_b2), 
             cols = vars(h_g3), 
             label = "label_parsed",
             scales="free") +
  labs(x = "cor", y = "Type 1 error") +           
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = updated_morandi_colors) +
  geom_hline(yintercept = 0.05, size = 1.3, linetype = "dashed", color = "red") +
  ggtitle("  ") +
  theme(legend.text = element_text(size = 30, face = "bold"),
        legend.position = "none") +
  theme(axis.text.x = element_text(size=30,face = "bold"),
        axis.text.y = element_text(size=30,face = "bold"),
        axis.title.x=element_text(size=30,face = "bold"),
        axis.title.y=element_text(size=30,face = "bold"),
        strip.text.x = element_text(size = 30, face = "bold", colour = "black", angle = 0),
        strip.text.y = element_text(size = 30, face = "bold", colour = "black", angle = -90)) +
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 15, face = "bold", color = "black"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "grey", fill = NA, size = 0.5), 
        strip.background = element_rect(color = "grey", size = 0.3))  
        


##Power plot
data_bind <- read.table("data2.txt", header = TRUE)
colnames(data_bind) <- c("GEI.beta1", "GEI.beta1.pval", "GEI.beta4", "GEI.beta4.pval",  #MERLIN(p)
                         "GEI3.beta1", "GEI3.beta1.pval", "GEI3.beta4", "GEI3.beta4.pval",  #MERLIN
                         "LDP.beta1", "LDP.beta1.pval",  #LDP method
                         "raps(T).beta1", "raps(T).beta1.pval",  #MR.raps in twosamplemr
                         "IVW.beta1", "IVW.beta1.pval",  #MR.IVW in twosamplemr
                         "Egger(T).beta1", "Egger(T).beta1.pval",   #MR.Egger in twosamplemr
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
dat_method_b1 <- data_bind[, c("GEI3.beta1.pval",  
                               "LDP.beta1.pval", "raps(T).beta1.pval", 
                               "IVW.beta1.pval", "Egger(T).beta1.pval")]
colnames(dat_method_b1) <- c("MERLIN", "MR-LDP", "RAPS", "IVW", "MR-Egger")
col <- ncol(dat_method_b1)
dat_hg3hb2b1 <- data_bind[, c("beta1", "h_g3", "h_b2")]
dat_hg3hb2b1 <- do.call(rbind, lapply(1:col, function(i) dat_hg3hb2b1))
dat_boxplot_b1 <- data.frame(value = unlist(dat_method_b1), method = rep(names(dat_method_b1), each = nrow(dat_method_b1)))
dat_boxplot_b1$method <- factor(dat_boxplot_b1$method, 
                                levels = c("MERLIN", "MR-LDP", "RAPS", "IVW", "MR-Egger"))
dat_hg3hb2b1[, 1] <- factor(dat_hg3hb2b1[, 1],
                            levels = c("0", "0.05", "0.1", "0.15", "0.2"))
dat_hg3hb2b1[, 2] <- factor(dat_hg3hb2b1[, 2]) 
dat_hg3hb2b1[, 3] <- factor(dat_hg3hb2b1[, 3]) 
dat <- cbind(dat_boxplot_b1, dat_hg3hb2b1)
proportion <- dat %>%
  group_by(method, beta1, h_g3, h_b2) %>%
  summarise(prop_beta_pval = sum(value < 0.05) / n(), .groups = "keep")
      
h_b2_filtered <- paste0("h[2]^2==", 0.1)
proportion_filtered <- proportion[proportion$h_b2 == h_b2_filtered, ]
proportion_final <- subset(proportion_filtered, 
                            (!(method %in% c("MR-LDP", "RAPS", "IVW", "MR-Egger") & h_g3 %in% c("h[3]^2==0.15", "h[3]^2==0.3"))))
print(proportion_final)

updated_morandi_colors <- c("#C1E0DB", "#CCD376", "#A28CC2", "#8498AB")    
p <- ggplot(proportion_final, aes(x = beta1, y = prop_beta_pval, color = method, shape = h_g3, linetype = h_g3)) +
  theme_bw()+ 
  geom_line(data = proportion_final[proportion_final$method != "MERLIN", ], aes(group = method), linewidth = 1.3, alpha = 4) +  
  geom_point(data = proportion_final[proportion_final$method != "MERLIN", ], aes(group = method), size = 10, alpha = 0.5) +  
  geom_line(data = proportion_final[proportion_final$method == "MERLIN", ], aes(group = h_g3), color = "#F5BE8F", linewidth = 1.3, alpha = 4) +  
  geom_point(data = proportion_final[proportion_final$method == "MERLIN", ], aes(group = h_g3), color = "#F5BE8F", size = 10, alpha = 0.5) + 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +  
  labs(title = "", 
       x = expression(bold("β"["1"])), 
       y = "Power") +
  scale_color_manual(values = updated_morandi_colors) +
  theme(legend.text = element_text(size = 30, face = "bold"),
        legend.position = "none",
        legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=30,face = "bold"),
        axis.text.y = element_text(size=30,face = "bold"),
        axis.title.x = element_text(size=30,face = "bold"), 
        axis.title.y = element_text(size=30,face = "bold"), 
        strip.text.x = element_text(size = 30, face = "bold", colour = "black", angle = 0),
        strip.text.y = element_text(size = 30, face = "bold", colour = "black", angle = -90),
        plot.title = element_text(hjust = 0.5, vjust = -1, size = 15, face = "bold", color = "black")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1)) + 
  theme(panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "grey", fill = NA, size = 0.5), 
        strip.background = element_rect(color = "grey", size = 0.3)) 
        
```
## Figure 3

```{r, eval = FALSE}

######Figure 1

library(ggplot2)    
library(dplyr)  

##Boxplot
rm(list=ls())
data_bind <- read.table("data3.txt", header = TRUE)

colnames(data_bind) <- c("GEI.beta1", "GEI.beta1.pval", "GEI.beta4", "GEI.beta4.pval",  #MERLIN(p)
                         "GEI3.beta1", "GEI3.beta1.pval", "GEI3.beta4", "GEI3.beta4.pval",  #MERLIN
                         "LDP.beta1.female", "LDP.beta1.pval.female",  #LDP
                         "raps_T.beta1.female", "raps_T.beta1.pval.female",  #MR.raps in twosamplemr
                         "IVW.beta1.female", "IVW.beta1.pval.female", "IVW.beta1.se.female", #MR.IVW in twosamplemr
                         "Egger_T.beta1.female", "Egger_T.beta1.pval.female",  "Egger_T.beta1.se.female", #MR.Egger in twosamplemr
                         "n_female",
                         "LDP.beta1.male", "LDP.beta1.pval.male",  #LDP
                         "raps_T.beta1.male", "raps_T.beta1.pval.male",  #MR.raps in twosamplemr
                         "IVW.beta1.male", "IVW.beta1.pval.male", "IVW.beta1.se.male",  #MR.IVW in twosamplemr
                         "Egger_T.beta1.male", "Egger_T.beta1.pval.male", "Egger_T.beta1.se.male", #MR.Egger in twosamplemr
                         "n_male",
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
data_bind$sex_stratified_beta4_ivw <- 0.5 * (data_bind$IVW.beta1.male - data_bind$IVW.beta1.female)
data_bind$sex_stratified_beta4_ldp <- 0.5 * (data_bind$LDP.beta1.male - data_bind$LDP.beta1.female)
data_bind$sex_stratified_beta4_raps <- 0.5 * (data_bind$raps_T.beta1.male - data_bind$raps_T.beta1.female)
data_bind$sex_stratified_beta4_egger <- 0.5 * (data_bind$Egger_T.beta1.male - data_bind$Egger_T.beta1.female)

dat_method_b4 <- data.frame(data_bind[, c("GEI3.beta4", "sex_stratified_beta4_ldp", "sex_stratified_beta4_raps",
                                          "sex_stratified_beta4_ivw", "sex_stratified_beta4_egger")])
colnames(dat_method_b4) <- c("MERLIN", "Sex-stratified LDP", "Sex-stratified RAPS", "Sex-stratified IVW", "Sex-stratified Egger")
col <- ncol(dat_method_b4)
dat_hg3hb2 <- data_bind[, c("cor_g1g3", "h_g3", "h_b2")]
dat_hg3hb2 <- do.call(rbind, lapply(1:col, function(i) dat_hg3hb2))
dat_boxplot_b4 <- data.frame(value = unlist(dat_method_b4), method = rep(names(dat_method_b4), each = nrow(dat_method_b4)))
dat_boxplot_b4$method <- factor(dat_boxplot_b4$method, 
                                levels = c("MERLIN", "Sex-stratified LDP", "Sex-stratified RAPS", "Sex-stratified IVW", "Sex-stratified Egger"))
dat_hg3hb2[, 1] <- factor(dat_hg3hb2[, 1])
dat_hg3hb2[, 2] <- factor(dat_hg3hb2[, 2])
dat_hg3hb2[, 3] <- factor(dat_hg3hb2[, 3])
dat <- cbind(dat_boxplot_b4, dat_hg3hb2)
dat <- dat[dat$value!= 0, ]

updated_morandi_colors <- c("#F5BE8F", "#C1E0DB", "#CCD376", "#A28CC2", "#8498AB")
p <- ggplot(dat, aes(x = cor_g1g3, y = value, fill = method)) +
     geom_boxplot(position = position_dodge(width = 0.8)) +
     theme_bw()+
     facet_grid(rows = vars(h_b2), 
                cols = vars(h_g3), 
                label = "label_parsed", 
                scales = "fixed") +
     scale_fill_manual(values = updated_morandi_colors) +
     geom_hline(aes(yintercept=b4), size = 1.3, linetype=5,col="red") +
     ggtitle(" ") +
     theme(legend.text = element_text(size = 30, face = "bold"),
           legend.position = "none")+                                                        
     theme(axis.text.x = element_text(size=30,face = "bold"),
           axis.text.y = element_text(size=30,face = "bold"),
           axis.title.x=element_text(size=30,face = "bold"),
           axis.title.y=element_text(size=30,face = "bold"),
           strip.text.x = element_text(size = 30, face = "bold", colour = "black", angle = 0),
           strip.text.y = element_text(size = 30, face = "bold", colour = "black", angle = -90)) +
     theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 15, face = "bold", color = "black"))+
     guides(fill=guide_legend(title=NULL)) +  
     theme(panel.spacing=unit(.05, "lines"),
           panel.border = element_rect(color = "grey", fill = NA, size = 0.5), 
           strip.background = element_rect(color = "grey", size = 0.3)) +
     scale_y_continuous(limits = c(-0.5, 0.5))


#Power plot
rm(list=ls())
data_bind <- read.table("data4.txt", header = TRUE)
colnames(data_bind) <- c("GEI.beta1", "GEI.beta1.pval", "GEI.beta4", "GEI.beta4.pval",  #MERLIN(p)
                         "GEI3.beta1", "GEI3.beta1.pval", "GEI3.beta4", "GEI3.beta4.pval",  #MERLIN
                         "LDP.beta1.female", "LDP.beta1.pval.female",  #LDP
                         "raps_T.beta1.female", "raps_T.beta1.pval.female", "raps_T.beta1.se.female",  #MR.raps in twosamplemr
                         "IVW.beta1.female", "IVW.beta1.pval.female", "IVW.beta1.se.female", #MR.IVW in twosamplemr
                         "Egger_T.beta1.female", "Egger_T.beta1.pval.female",  "Egger_T.beta1.se.female",  #MR.Egger in twosamplemr
                         "n_female",
                         "LDP.beta1.male", "LDP.beta1.pval.male",  #LDP
                         "raps_T.beta1.male", "raps_T.beta1.pval.male", "raps_T.beta1.se.male",  #MR.raps in twosamplemr
                         "IVW.beta1.male", "IVW.beta1.pval.male", "IVW.beta1.se.male",  #MR.IVW in twosamplemr
                         "Egger_T.beta1.male", "Egger_T.beta1.pval.male", "Egger_T.beta1.se.male",  #MR.Egger in twosamplemr
                         "n_male",
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
p_male <- pmax(data$LDP.beta1.pval.male, 1e-300) 
ldp_male_z <- -qnorm(log(p_male / 2), lower.tail = TRUE, log.p = TRUE) * sign(data$LDP.beta1.male)
ldp_male_se <- abs(data$LDP.beta1.male / ldp_male_z)
data$LDP.beta1.se.male <- ldp_male_se

p_female <- pmax(data$LDP.beta1.pval.female, 1e-300)
ldp_female_z <- -qnorm(log(p_female / 2), lower.tail = TRUE, log.p = TRUE) * sign(data$LDP.beta1.female)
ldp_female_se <- abs(data$LDP.beta1.female / ldp_female_z)
data$LDP.beta1.se.female <- ldp_female_se
data$sex_z_ldp <- (data$LDP.beta1.male - data$LDP.beta1.female)/sqrt(data$LDP.beta1.se.male*data$LDP.beta1.se.male + data$LDP.beta1.se.female*data$LDP.beta1.se.female)
log_p_ldp <- pnorm(abs(data$sex_z_ldp), lower.tail = FALSE, log.p = TRUE)  # log(p) for upper tail
data$sex_p_ldp <- 2 * exp(log_p_ldp) 

data$sex_z_ivw <- (data$IVW.beta1.male - data$IVW.beta1.female)/sqrt(data$IVW.beta1.se.male*data$IVW.beta1.se.male + data$IVW.beta1.se.female*data$IVW.beta1.se.female)
data$sex_p_ivw <- 2*(1 - pnorm(abs(data$sex_z_ivw)))

data$sex_z_raps <- (data$raps_T.beta1.male - data$raps_T.beta1.female)/sqrt(data$raps_T.beta1.se.male*data$raps_T.beta1.se.male + data$raps_T.beta1.se.female*data$raps_T.beta1.se.female)
data$sex_p_raps <- 2*(1 - pnorm(abs(data$sex_z_raps)))

data$sex_z_egger <- (data$Egger_T.beta1.male - data$Egger_T.beta1.female)/sqrt(data$Egger_T.beta1.se.male*data$Egger_T.beta1.se.male + data$Egger_T.beta1.se.female*data$Egger_T.beta1.se.female)
data$sex_p_egger <- 2*(1 - pnorm(abs(data$sex_z_egger)))

dat_method_b4 <- data_bind[, c("GEI3.beta4.pval", "sex_p_ldp", "sex_p_raps", "sex_p_ivw", "sex_p_egger")]
colnames(dat_method_b4) <- c("MERLIN", "Sex-stratified LDP", "Sex-stratified RAPS", "Sex-stratified IVW", "Sex-stratified Egger")
col <- ncol(dat_method_b4)
dat_hg3hb2b4 <- data_bind[, c("beta4", "h_g3", "h_b2")]
dat_hg3hb2b4 <- do.call(rbind, lapply(1:col, function(i) dat_hg3hb2b4))
dat_boxplot_b4 <- data.frame(value = unlist(dat_method_b4), method = rep(names(dat_method_b4), each = nrow(dat_method_b4)))
dat_boxplot_b4$method <- factor(dat_boxplot_b4$method, 
                                levels = c("MERLIN", "Sex-stratified LDP", "Sex-stratified RAPS", "Sex-stratified IVW", "Sex-stratified Egger"))
dat_hg3hb2b4[, 1] <- factor(dat_hg3hb2b4[, 1],
                            levels = c("0", "0.05", "0.1", "0.15", "0.2"))
dat_hg3hb2b4[, 2] <- factor(dat_hg3hb2b4[, 2]) 
dat_hg3hb2b4[, 3] <- factor(dat_hg3hb2b4[, 3]) 
dat <- cbind(dat_boxplot_b4, dat_hg3hb2b4)
proportion <- dat %>%
  group_by(method, beta4, h_g3, h_b2) %>%
  summarise(prop_beta = sum(value < 0.05), .groups = "keep")
proportion$prop_beta_pval <- proportion$prop_beta/500

updated_morandi_colors <- c("#F5BE8F", "#C1E0DB", "#CCD376", "#A28CC2", "#8498AB")
n_methods <- length(levels(proportion$method))

linewidth_values <- c(1.3, 5, 3, 1, 1.3)   
alpha_values <- c(3, 0.8, 1.5, 3, 3)   
p <- ggplot(proportion, aes(x = beta4, y = prop_beta_pval, color = method)) +
     theme_bw()+ 
     geom_line(aes(group = method, linewidth = method, alpha = method)) +
     geom_point(aes(shape = method), size = 6, alpha = 0.5) +  
     geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +  # 添加y=0.05的红色虚线 
     facet_grid(rows = vars(h_b2), 
                cols = vars(h_g3), 
                label = "label_parsed",
                scales="free") +
     labs(title = "", 
          x = expression(bold("β"["4"])), 
          y = "Power") +
     scale_color_manual(values = updated_morandi_colors) +
     scale_linewidth_manual(values = linewidth_values) +
     scale_alpha_manual(values = alpha_values) +
     theme(legend.text = element_text(size = 30, face = "bold"),
           legend.position = "none",
           legend.title = element_blank()) +
    theme(axis.text.x = element_text(size=30,face = "bold"),
          axis.text.y = element_text(size=30,face = "bold"),
          axis.title.x = element_text(size=30,face = "bold"), 
          axis.title.y = element_text(size=30,face = "bold"), 
          strip.text.x = element_text(size = 30, face = "bold", colour = "black", angle = 0),
          strip.text.y = element_text(size = 30, face = "bold", colour = "black", angle = -90),
          plot.title = element_text(hjust = 0.5, vjust = -1, size = 15, face = "bold", color = "black")) +
     scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1)) + 
     guides(fill = FALSE) +
     theme(panel.spacing=unit(.05, "lines"),
           panel.border = element_rect(color = "grey", fill = NA, size = 0.5), 
           strip.background = element_rect(color = "grey", size = 0.3)) 


##Continue E boxplot
rm(list=ls())
data_bind <- read.table("data3.txt", header = TRUE)
colnames(data_bind) <- c("GEI.beta1", "GEI.beta1.pval", "GEI.beta4", "GEI.beta4.pval",
                         "GEI3.beta1", "GEI3.beta1.pval", "GEI3.beta4", "GEI3.beta4.pval",
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
dat_method_b4 <- data.frame(data_bind[, c("GEI3.beta4")])
colnames(dat_method_b4) <- c("MERLIN")
col <- ncol(dat_method_b4)
dat_hg3hb2 <- data_bind[, c("cor_g1g3", "h_g3", "h_b2")]
dat_hg3hb2 <- do.call(rbind, lapply(1:col, function(i) dat_hg3hb2))
dat_boxplot_b4 <- data.frame(value = unlist(dat_method_b4), method = rep(names(dat_method_b4), each = nrow(dat_method_b4)))
dat_boxplot_b4$method <- factor(dat_boxplot_b4$method, levels = c("MERLIN"))
dat_hg3hb2[, 1] <- factor(dat_hg3hb2[, 1])
dat_hg3hb2[, 2] <- factor(dat_hg3hb2[, 2])
dat_hg3hb2[, 3] <- factor(dat_hg3hb2[, 3])
dat <- cbind(dat_boxplot_b4, dat_hg3hb2)
dat <- dat[dat$value!= 0, ]

updated_morandi_colors <- c("#F5BE8F") 
p <- ggplot(dat, aes(x = cor_g1g3, y = value, fill = method)) +
     geom_boxplot(position = position_dodge(width = 0.8)) +
     theme_bw()+
     facet_grid(rows = vars(h_b2), 
                cols = vars(h_g3), 
                label = "label_parsed", 
                scales = "fixed") +
     scale_fill_manual(values = updated_morandi_colors) +
     geom_hline(aes(yintercept=b4), size = 1.3, linetype=5,col="red") +
     ggtitle(" ") +
     theme(legend.text = element_text(size = 30, face = "bold"),
           legend.position = "none")+                                                        
     theme(axis.text.x = element_text(size=30,face = "bold"),
           axis.text.y = element_text(size=30,face = "bold"),
           axis.title.x=element_text(size=30,face = "bold"),
           axis.title.y=element_text(size=30,face = "bold"),
           strip.text.x = element_text(size = 30, face = "bold", colour = "black", angle = 0),
           strip.text.y = element_text(size = 30, face = "bold", colour = "black", angle = -90)) +
     theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 15, face = "bold", color = "black"))+
     guides(fill=guide_legend(title=NULL)) +  
     theme(panel.spacing=unit(.05, "lines"),
           panel.border = element_rect(color = "grey", fill = NA, size = 0.5), 
           strip.background = element_rect(color = "grey", size = 0.3)) 

##continue E power plot
rm(list=ls())
data_bind <- read.table("data4.txt", header = TRUE)
colnames(data_bind) <- c("GEI.beta1", "GEI.beta1.pval", "GEI.beta4", "GEI.beta4.pval",
                         "GEI3.beta1", "GEI3.beta1.pval", "GEI3.beta4", "GEI3.beta4.pval",
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
dat_method_b4 <- data.frame(data_bind[, c("GEI3.beta4.pval")])
colnames(dat_method_b4) <- c("MERLIN")
col <- ncol(dat_method_b4)
dat_hg3hb2b4 <- data_bind[, c("beta4", "h_g3", "h_b2", "cor_g1g3")]
dat_hg3hb2b4 <- do.call(rbind, lapply(1:col, function(i) dat_hg3hb2b4))
dat_boxplot_b4 <- data.frame(value = unlist(dat_method_b4), method = rep(names(dat_method_b4), each = nrow(dat_method_b4)))
dat_boxplot_b4$method <- factor(dat_boxplot_b4$method, 
                                levels = c("MERLIN"))
dat_hg3hb2b4[, 1] <- factor(dat_hg3hb2b4[, 1], levels = c("0", "0.05", "0.1", "0.15", "0.2"))
dat_hg3hb2b4[, 2] <- factor(dat_hg3hb2b4[, 2]) 
dat_hg3hb2b4[, 3] <- factor(dat_hg3hb2b4[, 3]) 
dat_hg3hb2b4[, 4] <- factor(dat_hg3hb2b4[, 4], levels = c("0", "0.4", "0.8")) 
dat <- cbind(dat_boxplot_b4, dat_hg3hb2b4)
proportion <- dat %>%
  group_by(beta4, h_g3, h_b2, cor_g1g3) %>%
  summarise(prop_beta = sum(value < 0.05), .groups = "keep")
proportion$prop_beta_pval <- proportion$prop_beta/500

updated_morandi_colors <- c("#F5BE8F", "#F5BE8F", "#F5BE8F")  
p <- ggplot(proportion, aes(x = beta4, y = prop_beta_pval, color = cor_g1g3)) +
  theme_bw()+ 
  geom_line(aes(group = cor_g1g3), linewidth = 1.3, alpha = 4) + 
  geom_point(aes(shape = cor_g1g3), size = 6, alpha = 0.5) +  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +  
  facet_grid(rows = vars(h_b2), 
             cols = vars(h_g3), 
             label = "label_parsed",
             scales="fixed") +
  labs(title = "", 
       x = expression(bold("β"["4"])), 
       y = "Power") +
  scale_color_manual(values = updated_morandi_colors) +
  theme(legend.text = element_text(size = 30, face = "bold"),
        legend.position = "none",
        legend.title = element_blank()) +
  theme(axis.text.x = element_text(size=30,face = "bold"),
        axis.text.y = element_text(size=30,face = "bold"),
        axis.title.x = element_text(size=30,face = "bold"), 
        axis.title.y = element_text(size=30,face = "bold"), 
        strip.text.x = element_text(size = 30, face = "bold", colour = "black", angle = 0),
        strip.text.y = element_text(size = 30, face = "bold", colour = "black", angle = -90),
        plot.title = element_text(hjust = 0.5, vjust = -1, size = 15, face = "bold", color = "black")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1)) + 
  guides(fill = FALSE) +
  theme(panel.spacing=unit(.05, "lines"),
        panel.border = element_rect(color = "grey", fill = NA, size = 0.5), 
        strip.background = element_rect(color = "grey", size = 0.3)) 
```
