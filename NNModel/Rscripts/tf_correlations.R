library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggprism)
library(stringr)
library(dplyr)
theme_set(theme_pubclean())
library(EnvStats) # for adding sample size
# library(ggpomological) #theme
library(survival)
library(ggfortify)
library(ranger)
library("anytime")
library(chron)
library(ggsurvfit)
library(gtsummary)
library(coxphf)
library(survminer)

setwd("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts")
source("prepare_data.R")

try="_t23"

data <- get_data_2(try)

# pearl 


pearl <- data$prl$all_pearl

cor(pearl$ichorTF, pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(pearl$ctDNA_Burden, pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(pearl$VAF_pers, pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))


data_pearl_ichor <- ggscatter(pearl, x = "ichorTF", y = "Predictedlabel",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "ichorCNA tumor fraction", y = "ALFAssay tumor fraction")


data_pearl_fragle <- ggscatter(pearl, x = "ctDNA_Burden", y = "Predictedlabel",
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "Fragle tumor fraction", y = "ALFAssay tumor fraction")


data_pearl_OncoFOLLOW <- ggscatter(pearl, x = "VAF_pers", y = "Predictedlabel",
                               add = "reg.line",  # Add regressin line
                               add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                               conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "OncoFollowDNA VAF", y = "ALFAssay tumor fraction")


pearl_tf_correlations <- cowplot::plot_grid(data_pearl_ichor,   data_pearl_fragle, data_pearl_OncoFOLLOW,
                        labels = c("A", "B", "C"),
                        ncol = 3, nrow=1)

# title <- cowplot::ggdraw() + cowplot::draw_label("Correlation plot between ALFAssay and other assays/tools on Pearl data", fontface='bold')
# pearl_tf_correlations <- cowplot::plot_grid(title, pearl_tf_correlations, ncol=1, rel_heights=c(0.1, 1))



ggsave(paste0("plots", try, "/PearlTFCorrelations.pdf", sep=""),pearl_tf_correlations,width=17, height=8)
ggsave(paste0("plots", try, "/plots_progress/PearlTFCorrelations.pdf", sep=""),pearl_tf_correlations,width=17, height=8)

# neorhea
neorhea <- data$nrh$all_neorhea

cor(neorhea$ichorTF, neorhea$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(neorhea$ctDNA_Burden, neorhea$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(neorhea$VAF_perc, neorhea$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))


data_neorhea_ichor <- ggscatter(neorhea, x = "ichorTF", y = "Predictedlabel",
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.05, label.y = 0.5) + labs(x = "ichorCNA tumor fraction", y = "ALFAssay tumor fraction")


data_neorhea_fragle <- ggscatter(neorhea, x = "ctDNA_Burden", y = "Predictedlabel",
                               add = "reg.line",  # Add regressin line
                               add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                               conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.02, label.y = 0.50) + labs(x = "Fragle tumor fraction", y = "ALFAssay tumor fraction")


data_neorhea_Radar <- ggscatter(neorhea, x = "VAF_perc", y = "Predictedlabel",
                                   add = "reg.line",  # Add regressin line
                                   add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                                   conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "RaDaR VAF %", y = "ALFAssay tumor fraction")


neorhea_tf_correlations <- cowplot::plot_grid(data_neorhea_ichor,   data_neorhea_fragle, data_neorhea_Radar,
                                            labels = c("C", "D", "E"),
                                            ncol = 3, nrow=1)

title <- cowplot::ggdraw() + cowplot::draw_label("Correlation plot between ALFAssay and other assays/tools on Neorhea data", fontface='bold')
neorhea_tf_correlations <- cowplot::plot_grid(title, neorhea_tf_correlations, ncol=1, rel_heights=c(0.1, 1))


ggsave(paste0("plots", try, "/NeorheaTFCorrelations.pdf", sep=""),neorhea_tf_correlations,width=17, height=8)


# synergy

synergy <- data$snr$all_synergy

cor(synergy$ichorTF, synergy$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(synergy$ctDNA_Burden, synergy$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))


data_snr_ichor <- ggscatter(synergy, x = "ichorTF", y = "Predictedlabel",
                                add = "reg.line",  # Add regressin line
                                add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "ichorCNA tumor fraction", y = "ALFAssay tumor fraction") 


data_snr_fragle <- ggscatter(synergy, x = "ctDNA_Burden", y = "Predictedlabel",
                                 add = "reg.line",  # Add regressin line
                                 add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                                 conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "Fragle tumor fraction", y = "ALFAssay tumor fraction")



snr_tf_correlations <- cowplot::plot_grid(data_snr_ichor,   data_snr_fragle,
                                              labels = c("A", "B"),
                                              ncol = 2, nrow=1)

title <- cowplot::ggdraw() + cowplot::draw_label("Correlation plot between ALFAssay and other assays/tools on Synergy data", fontface='bold')
snr_tf_correlations <- cowplot::plot_grid(title, snr_tf_correlations, ncol=1, rel_heights=c(0.1, 1))

ggsave(paste0("plots", try, "/SynergyTFCorrelations.pdf",sep=""),snr_tf_correlations,width=13, height=8)


suppl_file_corr <- cowplot::plot_grid(snr_tf_correlations,   neorhea_tf_correlations,
                                 # labels = c("I", "II"),
                                 ncol = 1, nrow=2)

ggsave(paste0("plots", try, "/TFCorrelationsTraining.pdf", sep=""),suppl_file_corr,width=17, height=8)

