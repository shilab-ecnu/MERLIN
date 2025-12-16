###-----------------------------------------------
### This code generates simulated phenotypes.
###-----------------------------------------------

library(mvtnorm)
library(MASS)
library(bigsnpr)
library(bigstatsr)
library(bigreadr)
library(data.table)
setwd("/MERLIN")

## Load file 
cau_info <- read.table("/extractsnpID_info.txt") #causal snps
snpall <- fread("/ukb22828_c1_v3_set2_norepeat.bim") #all snps
sampleall <-  fread("/ukb22828_c1_v3_set2_norepeat.fam") #all samples
E_x <- read.table("/E_1-1.txt", header = T) #If we want to analyse continue case, change it to continue E
  
G_cau_tempbedfile <- "/G_cau_temp.bed" #causal genotypes
G_cau_rds <- snp_readBed(G_cau_tempbedfile, backingfile = tempfile())
Gcau <- snp_attach(G_cau_rds)
Gcau  <- Gcau$genotypes
Gcau <- Gcau[1:nrow(Gcau),1:ncol(Gcau)]#genotype matrix
Gcau[is.na(Gcau)] <- 0
str(Gcau)

## Set parameters 
b1 <- 0
b4 <- 0.3
h_g1 <- 0.3;  
h_g3 <- 0.1
h_b <- c(0)
cor_g1g3 <- 0.4

n <- 80000  #samples
m <- ncol(Gcau) #snps
rhoxy <- 0.6; #cor between noise

results_matrix <- matrix(, nrow = 1, ncol = 18)

tic <- proc.time()
pp <- 20242025
set.seed(pp)    
sigma2g1 <- h_g1/m
sigma2g3 <- h_g3/m
sigma2b  <- h_b/m

gamma1_3 <- rmvnorm(m, mean = c(0, 0), matrix(c(sigma2g1, cor_g1g3*sqrt(sigma2g1)*sqrt(sigma2g3), 
                     cor_g1g3*sqrt(sigma2g1)*sqrt(sigma2g3),sigma2g3), ncol=2))
gamma_1x <- c(gamma1_3[, 1])
gamma_3x <- c(gamma1_3[, 2])
     
beta_2   <- rnorm(m, 0, sqrt(sigma2b))

sgx <- 1 - h_g1 - h_g3
sgy <- 1 - b1*b1 - h_b - b4*b4

Gcau_r <- nrow(Gcau) 
Gmean_x <- matrix(rep(colMeans(Gcau), nrow(Gcau)), ncol=ncol(Gcau), byrow = T)
Gcau <- Gcau - Gmean_x
   
GE_x <- Gcau * E_x[,1]

noise_o <- rmvnorm(Gcau_r, mean = c(0, 0), matrix(c(sgx,rhoxy*sqrt(sgx)*sqrt(sgy), 
                     rhoxy*sqrt(sgx)*sqrt(sgy),sgy), ncol=2))
                     noise_x <- noise_o[, 1]
                     noise_y <- noise_o[, 2]  

sampleall <-  fread("data_use/ukb22828_c1_v3_set2_norepeat.fam")  
X <- Gcau %*% gamma_1x + GE_x %*% gamma_3x + noise_x
Y <- X * b1 + Gcau %*% beta_2 + X * E_x[,1] * b4 + noise_y
