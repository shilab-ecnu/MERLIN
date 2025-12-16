
## Overlap plot

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
## QQplot
dat_method_b1 <- data_bind[, c("MERLIN.beta1.pval", "LDP.beta1.pval", 
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




