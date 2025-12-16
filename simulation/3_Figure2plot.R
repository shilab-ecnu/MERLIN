
###### Plot Figure 2
# beta1 means beta_A, beta4 means beta_I

library(ggplot2)    
library(dplyr)      

rm(list=ls())
data_bind <- read.table("data1.txt", header = TRUE)

colnames(data_bind) <- c("MERLIN(p).beta1", "MERLIN(p).beta1.pval", "MERLIN(p).beta4", "MERLIN(p).beta4.pval",  #MERLIN(p)
                         "MERLIN.beta1", "MERLIN.beta1.pval", "MERLIN.beta4", "MERLIN.beta4.pval",  #MERLIN
                         "LDP.beta1", "LDP.beta1.pval",  #LDP method
                         "raps(T).beta1", "raps(T).beta1.pval",  #MR.raps in twosamplemr
                         "IVW.beta1", "IVW.beta1.pval",  #MR.IVW in twosamplemr
                         "Egger(T).beta1", "Egger(T).beta1.pval",   #MR.Egger in twosamplemr
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
## Boxplot
dat_method_b1 <- data_bind[, c("MERLIN.beta1", "LDP.beta1", "raps(T).beta1", "IVW.beta1", "Egger(T).beta1")]
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
         
## Type 1 error
dat_method_b1_p <- data_bind[, c("MERLIN.beta1.pval",  
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
        


## Power plot
data_bind <- read.table("data2.txt", header = TRUE)
colnames(data_bind) <- c("MERLIN(p).beta1", "MERLIN(p).beta1.pval", "MERLIN(p).beta4", "MERLIN(p).beta4.pval",  #MERLIN(p)
                         "MERLIN.beta1", "MERLIN.beta1.pval", "MERLIN.beta4", "MERLIN.beta4.pval",  #MERLIN
                         "LDP.beta1", "LDP.beta1.pval",  #LDP method
                         "raps(T).beta1", "raps(T).beta1.pval",  #MR.raps in twosamplemr
                         "IVW.beta1", "IVW.beta1.pval",  #MR.IVW in twosamplemr
                         "Egger(T).beta1", "Egger(T).beta1.pval",   #MR.Egger in twosamplemr
                         "beta1", "beta4", "h_g3", "h_b2", "cor_g1g3")
dat_method_b1 <- data_bind[, c("MERLIN.beta1.pval",  
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
        
        