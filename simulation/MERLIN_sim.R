

library(mvtnorm)
library(Rcpp)
library(MR.LDP)
library(TwoSampleMR)
library(PDSCE)
library(MASS)
library(bigsnpr)
library(bigstatsr)
library(bigreadr)
library(data.table)
setwd("¡«/MERLIN")
source("code/R/summary_lm.R")
sourceCpp("code/cpp/Gibbs_samp_Gamseo.cpp")
sourceCpp("code/cpp/Gibbs_samp_Gam3seo.cpp")
sourceCpp("code/cpp/truncated.cpp")

######load file
cau_info <- read.table("/extractsnpID_info.txt")
snpall <- fread("/ukb22828_c1_v3_set2_norepeat.bim")
sampleall <-  fread("/ukb22828_c1_v3_set2_norepeat.fam")
E_x <- read.table("/E_1-1.txt", header = T)
  
G_cau_tempbedfile <- "/G_cau_temp.bed"
G_cau_rds <- snp_readBed(G_cau_tempbedfile, backingfile = tempfile())
Gcau <- snp_attach(G_cau_rds)
Gcau  <- Gcau$genotypes
Gcau <- Gcau[1:nrow(Gcau),1:ncol(Gcau)]##genotype matrix
Gcau[is.na(Gcau)] <- 0
str(Gcau)

######set parameters 
pc <- 1
b1 <- 0
b4 <- 0.3
h_g1 <- 0.3;  
h_g3 <- 0.1
h_b <- c(0)
cor_g1g3 <- 0.4
intersect_choose <- 4

name_part <- paste0("No.", pc, "_b1=", b1, "_b4=", b4, "_hg3=", h_g3, "_hb=", h_b, "_corg1g3=", cor_g1g3, "_intersect", intersect_choose)
  
n <- 80000  #samples
m <- ncol(Gcau)
rhoxy <- 0.6;  

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

exp_gwas <- X[1:n, ]
exp_gwis <- X[1:n, ]
out_gwas <- Y[(2*n+1):(3*n), ]
out_gwis <- Y[(2*n+1):(3*n), ]
  
######plink to estimate effect size
exp_gwas_phy <- cbind(sampleall$V1[1:n], sampleall$V1[1:n], exp_gwas)
exp_gwis_phy <- cbind(sampleall$V1[1:n], sampleall$V1[1:n], exp_gwis) 
out_gwas_phy <- cbind(sampleall$V1[c((2*n+1):(3*n))], sampleall$V1[c((2*n+1):(3*n))], out_gwas) 
out_gwis_phy <- cbind(sampleall$V1[c((2*n+1):(3*n))], sampleall$V1[c((2*n+1):(3*n))], out_gwis)

exp_gwis_E <- data.frame(sampleall$V1[1:n], sampleall$V1[1:n], E_x[1:n,])
out_gwis_E <- data.frame(sampleall$V1[c((2*n+1):(3*n))], sampleall$V1[c((2*n+1):(3*n))], E_x[c((2*n+1):(3*n)),])

write.table(exp_gwas_phy, paste0("/plink_summarystats/", name_part, "_exp_gwas_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(exp_gwis_phy, paste0("/plink_summarystats/", name_part, "_exp_gwis_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(out_gwas_phy, paste0("/plink_summarystats/", name_part, "_out_gwas_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(out_gwis_phy, paste0("/plink_summarystats/", name_part, "_out_gwis_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(exp_gwis_E, paste0("/plink_summarystats/", name_part, "_exp_gwis_E.txt"), row.names = FALSE, col.names = FALSE)
write.table(out_gwis_E, paste0("/plink_summarystats/", name_part, "_out_gwis_E.txt"), row.names = FALSE, col.names = FALSE)

fit1cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --pheno /plink_summarystats/", name_part, "_exp_gwas_phy.txt --glm allow-no-covars --out /plink_summarystats/", name_part, "_exp_gwas_output")
system(fit1cmd)

fit3cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --linear interaction --pheno /plink_summarystats/", name_part, "_exp_gwis_phy.txt --no-psam-pheno --parameters 2,3 --covar /plink_summarystats/", name_part, "_exp_gwis_E.txt --out /plink_summarystats/", name_part, "_exp_gwis_output_plink2")
system(fit3cmd)

Fit1cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --pheno /plink_summarystats/", name_part, "_out_gwas_phy.txt --glm allow-no-covars --out /plink_summarystats/", name_part, "_out_gwas_output")
system(Fit1cmd)

Fit3cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --linear interaction --pheno /plink_summarystats/", name_part, "_out_gwis_phy.txt --no-psam-pheno --parameters 2,3 --covar /plink_summarystats/", name_part, "_out_gwis_E.txt --out /plink_summarystats/", name_part, "_out_gwis_output_plink2")
system(Fit3cmd)
 
######load plink results
fit1 <- fread(paste0("/plink_summarystats/", name_part, "_exp_gwas_output.PHENO1.glm.linear"))
fit3 <- fread(paste0("/plink_summarystats/", name_part, "_exp_gwis_output_plink2.PHENO1.glm.linear"))
Fit1 <- fread(paste0("/plink_summarystats/", name_part, "_out_gwas_output.PHENO1.glm.linear"))    
Fit3 <- fread(paste0("/plink_summarystats/", name_part, "_out_gwis_output_plink2.PHENO1.glm.linear")) 
  
fit3 <- fit3[fit3$TEST == "ADDxCOVAR1", ]
Fit3 <- Fit3[Fit3$TEST == "ADDxCOVAR1", ]

p_cutoff <- 1e-4

######choose iv
fit1_pcut <- fit1[fit1$P < 0.01, ]
fit3_pcut <- fit3[fit3$P < 0.01, ]
write.table(fit1_pcut, paste0("/plink_summarystats/", name_part, "_exp_gwas_output.PHENO1.glm.linear"), row.names = FALSE, col.names = TRUE, quote = FALSE) 
expgwas_LD_cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --clump-p1 5e-8 --clump-r2 0.01 --clump-kb 500 --clump /plink_summarystats/", name_part, "_exp_gwas_output.PHENO1.glm.linear --out /plink_summarystats/causal_snp/", name_part, "_expgwas_causal")
system(expgwas_LD_cmd)
  
write.table(fit3_pcut, paste0("/plink_summarystats/", name_part, "_exp_gwis_output_plink2_interaction.PHENO1.glm.linear"), row.names = FALSE, col.names = TRUE, quote = FALSE)
expgwis_LD_cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --clump-p1 5e-8 --clump-r2 0.01 --clump-kb 500 --clump /plink_summarystats/", name_part, "_exp_gwis_output_plink2_interaction.PHENO1.glm.linear --out /plink_summarystats/causal_snp/", name_part, "_expgwis_causal")
system(expgwis_LD_cmd)
  
expgwas_causal <- fread(paste0("/plink_summarystats/causal_snp/", name_part, "_expgwas_causal.clumped"))
causal_snp <- expgwas_causal$SNP
  
expgwis_causal <- fread(paste0("/plink_summarystats/causal_snp/", name_part, "_expgwis_causal.clumped"))
causal_snp_gwis <- expgwis_causal$SNP

causal_union <- union(causal_snp, causal_snp_gwis)
dup_snps <- intersect(causal_snp, causal_snp_gwis)
non_dup_snps <- setdiff(causal_union, dup_snps)
gwas_nodup <- intersect(causal_snp, non_dup_snps)
gwis_nodup <- intersect(causal_snp_gwis, non_dup_snps)

fit1_dup <- fit1[fit1$ID %in% dup_snps, ]
fit3_dup <- fit3[fit3$ID %in% dup_snps, ]
P_min <- pmin(fit1_dup$P, fit3_dup$P)
fit1_dup <- fit1[fit1$P %in% P_min, ]
fit3_dup <- fit3[fit3$P %in% P_min, ]
fit1_nodup <- fit1[fit1$ID %in% gwas_nodup, ]
fit3_nodup <- fit3[fit3$ID %in% gwis_nodup, ]
final_data <- rbind(fit1_dup, fit3_dup, fit1_nodup, fit3_nodup)
write.table(final_data, paste0("/plink_summarystats/", name_part, "_exp_unioncausal_output_plink2_interaction.PHENO1.glm.linear"), row.names = FALSE, col.names = TRUE, quote = FALSE)
union_LD_cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --clump-p1 1 --clump-r2 0.01 --clump-kb 500 --clump /plink_summarystats/", name_part, "_exp_unioncausal_output_plink2_interaction.PHENO1.glm.linear --out /plink_summarystats/causal_snp/", name_part, "_union_causal")
system(union_LD_cmd)
union_causal <- fread(paste0("/plink_summarystats/causal_snp/", name_part, "_union_causal.clumped"))
causal_snp_union <- union_causal$SNP
gwas_gwis_bothin <- causal_snp_union

#######MR analysis
cauidx <- match(gwas_gwis_bothin, fit1$ID)
R <- diag(length(cauidx))
cauidx_other <- match(causal_snp, fit1$ID) 
R_other <- diag(length(cauidx_other))
rhoe <- 0
rhoee <- 0

##MERLIN(p)
res1 <- MRGEI_Gamseo(fit1$BETA[cauidx], fit3$BETA[cauidx], Fit1$BETA[cauidx], fit1$SE[cauidx], fit3$SE[cauidx], Fit1$SE[cauidx], R, rhoe)
str(res1)
#MERLIN 
res2 <- MRGEI_Gam3seo(fit1$BETA[cauidx], fit3$BETA[cauidx], Fit1$BETA[cauidx], Fit3$BETA[cauidx], fit1$SE[cauidx], fit3$SE[cauidx], Fit1$SE[cauidx], Fit3$SE[cauidx], R, rhoe, rhoee)
str(res2)
##MR.LDP
gammah <- fit1$BETA[cauidx_other]
Gammah <- Fit1$BETA[cauidx_other]
segamma <- fit1$SE[cauidx_other]
seGamma <- Fit1$SE[cauidx_other]
gamma <- rep(0.01,length(cauidx_other))
alpha <- rep(0.01,length(cauidx_other))
sgga2 <- 0.01
sgal2 <- 0.01
maxIter = 10000;
diagnostics = FALSE
beta0 <- 0.1
epsStopLogLik <- 1e-6
out1 <- MRLDP_SimPXvb(gammah, Gammah, segamma, seGamma, gamma, alpha,  beta0, sgga2, sgal2, R_other,
                              0, epsStopLogLik, maxIter, model = 1);
out2 <- MRLDP_SimPXvb(gammah, Gammah, segamma, seGamma, gamma, alpha,  beta0, sgga2, sgal2, R_other,
                              1, epsStopLogLik, maxIter, model = 1);
tstatLD <- 2*(out1$tstat - out2$tstat)
pvalLD = pchisq(tstatLD, 1, lower.tail = F)

out3 <- MRLDP_SimPXvb(gammah, Gammah, segamma, seGamma, gamma, alpha,  beta0, sgga2, sgal2, R_other,
                           0, epsStopLogLik, maxIter, model = 2);
out4 <- MRLDP_SimPXvb(gammah, Gammah, segamma, seGamma, gamma, alpha,  beta0, sgga2, sgal2, R_other,
                           1, epsStopLogLik, maxIter, model = 2);
tstatLDP <- 2*(out3$tstat - out4$tstat)
pvalLDP = pchisq(tstatLDP, 1, lower.tail = F)
##MR.raps
raps.res <- mr_raps(gammah, Gammah, segamma, seGamma)
##ivw
ivw.res <- mr_ivw(gammah, Gammah, segamma, seGamma);
##mr.egger
egger.res <- mr_egger_regression(gammah, Gammah, segamma, seGamma)
  
results_matrix[pp, ] <- c(res1$Beta1.hat, res1$Beta1.pval, res1$Beta4.hat, res1$Beta4.pval,  #MRGEI_Gamseo
                        res2$Beta1.hat, res2$Beta1.pval, res2$Beta4.hat, res2$Beta4.pval,  #MRGEI_Gam3seo
                        out1$beta0, pvalLD,  #MR.LD
                        out3$beta0, pvalLDP,  #MR.LDP
                        raps.res$b, raps.res$pval,
                        ivw.res$b, ivw.res$pval,
                        egger.res$b, egger.res$pval
                        )
toc <- proc.time()
print((toc - tic)[3])

##remove files
setwd("/plink_summarystats/")
files <- list.files(full.names = TRUE)
files_to_remove <- files[grepl(name_part, files)]
if (length(files_to_remove) > 0) {
   file.remove(files_to_remove)
}
