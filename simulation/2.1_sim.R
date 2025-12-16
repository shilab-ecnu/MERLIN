###---------------------------------------------------------
### This code simulate MERLIN(p) and MERLIN in binary balance E
###---------------------------------------------------------

library(MERLIN)
library(MR.LDP)
library(TwoSampleMR)
setwd("/MERLIN")

n <- samp_num
exp_gwas <- X[1:n, ]
exp_gwis <- X[1:n, ]
out_gwas <- Y[(2*n+1):(3*n), ]
out_gwis <- Y[(2*n+1):(3*n), ]
  
## plink to estimate effect size 
exp_gwas_phy <- cbind(sampleall$V1[1:n], sampleall$V1[1:n], exp_gwas)
exp_gwis_phy <- cbind(sampleall$V1[1:n], sampleall$V1[1:n], exp_gwis) 
out_gwas_phy <- cbind(sampleall$V1[c((2*n+1):(3*n))], sampleall$V1[c((2*n+1):(3*n))], out_gwas) 
out_gwis_phy <- cbind(sampleall$V1[c((2*n+1):(3*n))], sampleall$V1[c((2*n+1):(3*n))], out_gwis)

exp_gwis_E <- data.frame(sampleall$V1[1:n], sampleall$V1[1:n], E_x[1:n,])
out_gwis_E <- data.frame(sampleall$V1[c((2*n+1):(3*n))], sampleall$V1[c((2*n+1):(3*n))], E_x[c((2*n+1):(3*n)),])

write.table(exp_gwas_phy, paste0("/plink_summarystats/exp_gwas_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(exp_gwis_phy, paste0("/plink_summarystats/exp_gwis_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(out_gwas_phy, paste0("/plink_summarystats/out_gwas_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(out_gwis_phy, paste0("/plink_summarystats/out_gwis_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(exp_gwis_E, paste0("/plink_summarystats/exp_gwis_E.txt"), row.names = FALSE, col.names = FALSE)
write.table(out_gwis_E, paste0("/plink_summarystats/out_gwis_E.txt"), row.names = FALSE, col.names = FALSE)

fit1cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --pheno /plink_summarystats/exp_gwas_phy.txt --glm allow-no-covars --out /plink_summarystats/exp_gwas_output")
system(fit1cmd)

fit3cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --linear interaction --pheno /plink_summarystats/exp_gwis_phy.txt --no-psam-pheno --parameters 2,3 --covar /plink_summarystats/exp_gwis_E.txt --out /plink_summarystats/exp_gwis_output_plink2")
system(fit3cmd)

Fit1cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --pheno /plink_summarystats/out_gwas_phy.txt --glm allow-no-covars --out /plink_summarystats/out_gwas_output")
system(Fit1cmd)

Fit3cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --linear interaction --pheno /plink_summarystats/out_gwis_phy.txt --no-psam-pheno --parameters 2,3 --covar /plink_summarystats/out_gwis_E.txt --out /plink_summarystats/out_gwis_output_plink2")
system(Fit3cmd)
 
## load plink results 
fit1 <- fread(paste0("/plink_summarystats/exp_gwas_output.PHENO1.glm.linear"))
fit3 <- fread(paste0("/plink_summarystats/exp_gwis_output_plink2.PHENO1.glm.linear"))
Fit1 <- fread(paste0("/plink_summarystats/out_gwas_output.PHENO1.glm.linear"))    
Fit3 <- fread(paste0("/plink_summarystats/out_gwis_output_plink2.PHENO1.glm.linear")) 
  
fit3 <- fit3[fit3$TEST == "ADDxCOVAR1", ]
Fit3 <- Fit3[Fit3$TEST == "ADDxCOVAR1", ]

p_cutoff <- 1e-4

## choose iv 
fit1_pcut <- fit1[fit1$P < 0.01, ]
fit3_pcut <- fit3[fit3$P < 0.01, ]
write.table(fit1_pcut, paste0("/plink_summarystats/exp_gwas_output.PHENO1.glm.linear"), row.names = FALSE, col.names = TRUE, quote = FALSE) 
expgwas_LD_cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --clump-p1 5e-8 --clump-r2 0.01 --clump-kb 500 --clump /plink_summarystats/exp_gwas_output.PHENO1.glm.linear --out /plink_summarystats/causal_snp/expgwas_causal")
system(expgwas_LD_cmd)
  
write.table(fit3_pcut, paste0("/plink_summarystats/exp_gwis_output_plink2_interaction.PHENO1.glm.linear"), row.names = FALSE, col.names = TRUE, quote = FALSE)
expgwis_LD_cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --clump-p1 5e-8 --clump-r2 0.01 --clump-kb 500 --clump /plink_summarystats/exp_gwis_output_plink2_interaction.PHENO1.glm.linear --out /plink_summarystats/causal_snp/expgwis_causal")
system(expgwis_LD_cmd)
  
expgwas_causal <- fread(paste0("/plink_summarystats/causal_snp/expgwas_causal.clumped"))
causal_snp <- expgwas_causal$SNP
  
expgwis_causal <- fread(paste0("/plink_summarystats/causal_snp/expgwis_causal.clumped"))
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
write.table(final_data, paste0("/plink_summarystats/exp_unioncausal_output_plink2_interaction.PHENO1.glm.linear"), row.names = FALSE, col.names = TRUE, quote = FALSE)
union_LD_cmd <- paste0("/plink2 --bfile /ukb22828_c1_v3_set2_norepeat --clump-p1 1 --clump-r2 0.01 --clump-kb 500 --clump /plink_summarystats/exp_unioncausal_output_plink2_interaction.PHENO1.glm.linear --out /plink_summarystats/causal_snp/union_causal")
system(union_LD_cmd)
union_causal <- fread(paste0("/plink_summarystats/causal_snp/union_causal.clumped"))
causal_snp_union <- union_causal$SNP
gwas_gwis_both <- causal_snp_union

## MR analysis 
cauidx <- match(gwas_gwis_both, fit1$ID)
R <- diag(length(cauidx))
cauidx_other <- match(causal_snp, fit1$ID) 
R_other <- diag(length(cauidx_other))
rhoe <- 0
rhoee <- 0

# MERLIN(p)
res1 <- MERLIN(fit1$BETA[cauidx], fit3$BETA[cauidx], Fit1$BETA[cauidx], fit1$SE[cauidx], fit3$SE[cauidx], Fit1$SE[cauidx], R, rhoe)
str(res1)
# MERLIN 
res2 <- MERLIN(fit1$BETA[cauidx], fit3$BETA[cauidx], Fit1$BETA[cauidx], Fit3$BETA[cauidx], fit1$SE[cauidx], fit3$SE[cauidx], Fit1$SE[cauidx], Fit3$SE[cauidx], R, rhoe, rhoee)
str(res2)
# MR.LDP
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
# raps
raps.res <- mr_raps(gammah, Gammah, segamma, seGamma)
# ivw
ivw.res <- mr_ivw(gammah, Gammah, segamma, seGamma);
# mr.egger
egger.res <- mr_egger_regression(gammah, Gammah, segamma, seGamma)
  
results_matrix <- c(res1$Beta1.hat, res1$Beta1.pval, res1$Beta4.hat, res1$Beta4.pval,  #MERLIN(p)
                        res2$Beta1.hat, res2$Beta1.pval, res2$Beta4.hat, res2$Beta4.pval,  #MERLIN
                        out1$beta0, pvalLD,  #MR.LD
                        out3$beta0, pvalLDP,  #MR.LDP
                        raps.res$b, raps.res$pval, #raps
                        ivw.res$b, ivw.res$pval, #ivw
                        egger.res$b, egger.res$pval #egger
                        )
toc <- proc.time()
print((toc - tic)[3])

## remove files 
setwd("/plink_summarystats/")
files <- list.files(full.names = TRUE)
files_to_remove <- files[grepl(files)]
if (length(files_to_remove) > 0) {
   file.remove(files_to_remove)
}
