
## Plot Figure 3

library(ggplot2)    
library(dplyr)  

## Boxplot
rm(list=ls())
data_bind <- read.table("data3.txt", header = TRUE)

colnames(data_bind) <- c("MERLIN(p).beta1", "MERLIN(p).beta1.pval", "MERLIN(p).beta4", "MERLIN(p).beta4.pval",  #MERLIN(p)
                         "MERLIN.beta1", "MERLIN.beta1.pval", "MERLIN.beta4", "MERLIN.beta4.pval",  #MERLIN
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

dat_method_b4 <- data.frame(data_bind[, c("MERLIN.beta4", "sex_stratified_beta4_ldp", "sex_stratified_beta4_raps",
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


## Power plot
rm(list=ls())
data_bind <- read.table("data4.txt", header = TRUE)
colnames(data_bind) <- c("MERLIN(p).beta1", "MERLIN(p).beta1.pval", "MERLIN(p).beta4", "MERLIN(p).beta4.pval",  #MERLIN(p)
                         "MERLIN.beta1", "MERLIN.beta1.pval", "MERLIN.beta4", "MERLIN.beta4.pval",  #MERLIN
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

dat_method_b4 <- data_bind[, c("MERLIN.beta4.pval", "sex_p_ldp", "sex_p_raps", "sex_p_ivw", "sex_p_egger")]
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


## Continue E boxplot
rm(list=ls())
data_bind <- read.table("data3.txt", header = TRUE)
colnames(data_bind) <- c("MERLIN(p).beta1", "MERLIN(p).beta1.pval", "MERLIN(p).beta4", "MERLIN(p).beta4.pval",
                         "MERLIN.beta1", "MERLIN.beta1.pval", "MERLIN.beta4", "MERLIN.beta4.pval",
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
dat_method_b4 <- data.frame(data_bind[, c("MERLIN.beta4")])
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

## Continue E power plot
rm(list=ls())
data_bind <- read.table("data4.txt", header = TRUE)
colnames(data_bind) <- c("MERLIN(p).beta1", "MERLIN(p).beta1.pval", "MERLIN(p).beta4", "MERLIN(p).beta4.pval",
                         "MERLIN.beta1", "MERLIN.beta1.pval", "MERLIN.beta4", "MERLIN.beta4.pval",
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
dat_method_b4 <- data.frame(data_bind[, c("MERLIN.beta4.pval")])
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



