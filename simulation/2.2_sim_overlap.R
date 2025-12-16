###---------------------------------------------
### This code simulate MERLIN under sample overlap
###---------------------------------------------

library(MERLIN)
setwd("/MERLIN")

d1 <- 20000 
d2 <- 20000
n <- samp_num
exp_gwas <- X[1:n, ]
exp_gwis <- X[1:n, ]
out_gwas <- Y[(2*n+1):(3*n), ]
out_gwis <- Y[(2*n+1):(3*n), ]
  
# simulate overlapping samples
# This is creating overlap samples between gamma1 & Gamma1 
o1 <- 1:d1
out_gwaso  <- Y[o1, ]
out_gwas <- c(out_gwaso, out_gwas[-o1])
     
# This is creating overlap samples between gamma3 & Gamma3
o2 <- 1:d2
out_gwiso  <- Y[(o2), ]
out_gwis <- c(out_gwiso, out_gwis[-o2])

exp_gwas_phy <- cbind(sampleall$V1[1:n], sampleall$V1[1:n], exp_gwas)
exp_gwis_phy <- cbind(sampleall$V1[1:n], sampleall$V1[1:n], exp_gwis) 
out_gwas_phy <- cbind(sampleall$V1[c(o1, (n+1+d1):(2*n))], sampleall$V1[c(o1, (n+1+d1):(2*n))], out_gwas) 
out_gwis_phy <- cbind(sampleall$V1[c(o2, (n+1+d1):(2*n))], sampleall$V1[c(o2, (n+1+d1):(2*n))], out_gwis)

exp_gwis_E <- data.frame(sampleall$V1[1:n], sampleall$V1[1:n], E_xx[1:n,])
out_gwis_E <- data.frame(sampleall$V1[c(o2, (n+1+d1):(2*n))], sampleall$V1[c(o2, (n+1+d1):(2*n))], E_xx[c(o2, (n+1+d1):(2*n)),])

write.table(exp_gwas_phy, paste0("/plink_summarystats/exp_gwas_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(exp_gwis_phy, paste0("/plink_summarystats/exp_gwis_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(out_gwas_phy, paste0("/plink_summarystats/out_gwas_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(out_gwis_phy, paste0("/plink_summarystats/out_gwis_phy.txt"), row.names = FALSE, col.names = FALSE)
write.table(exp_gwis_E, paste0("/plink_summarystats/exp_gwis_E.txt"), row.names = FALSE, col.names = FALSE)
write.table(out_gwis_E, paste0("/plink_summarystats/out_gwis_E.txt"), row.names = FALSE, col.names = FALSE)

## Plink
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
gwas_gwis_bothin <- causal_snp_union

## calculate rho caused by overlap 
# 1. independent snps 2. non causal snp
fit1u <- fit1[fit1$ID %in% ind_info$V2, ]
fit3u <- fit3[fit3$ID %in% ind_info$V2, ]
Fit1u <- Fit1[Fit1$ID %in% ind_info$V2, ]
Fit3u <- Fit3[Fit3$ID %in% ind_info$V2, ]

pth1 <- 1.96
x_gwas <- fit1u$BETA/fit1u$SE
x_gwis <- fit3u$BETA/fit3u$SE
y_gwas <- Fit1u$BETA/Fit1u$SE
y_gwis <- Fit3u$BETA/Fit3u$SE

indsnp_null_gwas <- intersect(which(abs(x_gwas) < pth1), which(abs(y_gwas) < pth1))
indsnp_null_gwis <- intersect(which(abs(x_gwis) < pth1), which(abs(y_gwis) < pth1))
indsnp_null <- intersect(indsnp_null_gwas, indsnp_null_gwis)
  
a = rep(-pth1, 2)
b = rep( pth1, 2)
rhohat1 <- EstRhofun(a, b, x_gwas[indsnp_null_gwas], y_gwas[indsnp_null_gwas], 4000, 1000, 10)
rhohat2 <- EstRhofun(a, b, x_gwis[indsnp_null_gwis], y_gwis[indsnp_null_gwis], 4000, 1000, 10)

rhoe <- mean(rhohat1)
rhoee <- mean(rhohat2)
## MR analysis 
cauidx <- match(gwas_gwis_bothin, fit1$ID)
R <- diag(length(cauidx))
cauidx_other <- match(causal_snp, fit1$ID) 
R_other <- diag(length(cauidx_other))

# MERLIN(p)
res1 <- MERLIN(fit1$BETA[cauidx], fit3$BETA[cauidx], Fit1$BETA[cauidx], fit1$SE[cauidx], fit3$SE[cauidx], Fit1$SE[cauidx], R, rhoe)
str(res1)
# MERLIN 
res2 <- MERLIN(fit1$BETA[cauidx], fit3$BETA[cauidx], Fit1$BETA[cauidx], Fit3$BETA[cauidx], fit1$SE[cauidx], fit3$SE[cauidx], Fit1$SE[cauidx], Fit3$SE[cauidx], R, rhoe, rhoee)
str(res2)

toc <- proc.time()
print((toc - tic)[3])

## remove files
setwd("/plink_summarystats/")
files <- list.files(full.names = TRUE)
files_to_remove <- files[grepl(files)]
if (length(files_to_remove) > 0) {
   file.remove(files_to_remove)
}