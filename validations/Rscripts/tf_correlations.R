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
source("survivalFunctionR.R")

data <- read.csv("/Users/alexandra/PhD/FragmentationPatterns/Data/MetaData/AllMetaData.csv", sep="\t")
data$status <- data$survivalStatus
data$BCFS <-  data$survivalTime

data$ichorPrediction[data$ichorTF>=0.03]<-"Positive"
data$ichorPrediction[data$ichorTF<0.03]<-"Negative" 
data$ctDNADetected[data$ctDNADetected==0]<-"Negative"
data$ctDNADetected[data$ctDNADetected==1]<-"Positive" 
data$ctDNA[data$ctDNADetected=="Fail"]<-"Fail" 

# predictions <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ctDNADetection_corrected2.csv", sep="\t")

predictions <- read.csv("~/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF.csv", sep="\t")

data_pred <- merge(data,predictions, by="PatientId" )
data_pred$Predictedlabel = data_pred$Predictedlabel/100

data_pred$ALFAssay[data_pred$Predictedlabel>=0.03]<-"Positive"
data_pred$ALFAssay[data_pred$Predictedlabel<0.03]<-"Negative" 


# pearl 
fragle_pred <-  read.csv("/Users/alexandra/PhD/PearlStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")

fragle_pred$PatientId <- str_replace(fragle_pred$Sample_ID, ".markDup.GCtagged", "")

data_pearl<-data_pred %>% filter(study=="Pearl")
data_pearl$VAF_pers <- data_pearl$VAF/100

data_pearl <- data_pearl %>% filter(!is.na(ichorTF))
data_pearl <- data_pearl %>% filter(!is.na(BCFS))
data_pearl <- data_pearl %>% filter(!is.na(ALFAssay))

data_pearl <- merge(fragle_pred, data_pearl, by="PatientId")

data_pearl$ctDNA_Burden_surv[data_pearl$ctDNA_Burden>=0.01]<-"HighctDNA" 
data_pearl$ctDNA_Burden_surv[data_pearl$ctDNA_Burden<0.01]<-"LowctDNA"


data_pearl_bl <- data_pearl %>% filter(timepoint=="BL")
data_pearl_d15 <- data_pearl %>% filter(timepoint=="C1")
data_pearl_pd <- data_pearl %>% filter(timepoint=="Surgery")

cor(data_pearl$ichorTF, data_pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(data_pearl$ctDNA_Burden, data_pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(data_pearl$VAF_pers, data_pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))


data_pearl_ichor <- ggscatter(data_pearl, x = "ichorTF", y = "Predictedlabel",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "ichorCNA tumor fraction", y = "ALFAssay tumor fraction")


data_pearl_fragle <- ggscatter(data_pearl, x = "ctDNA_Burden", y = "Predictedlabel",
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "Fragle tumor fraction", y = "ALFAssay tumor fraction")


data_pearl_OncoFOLLOW <- ggscatter(data_pearl, x = "VAF_pers", y = "Predictedlabel",
                               add = "reg.line",  # Add regressin line
                               add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                               conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "OncoFollowDNA tumor fraction", y = "ALFAssay tumor fraction")


pearl_tf_correlations <- cowplot::plot_grid(data_pearl_ichor,   data_pearl_fragle, data_pearl_OncoFOLLOW,
                        labels = c("A", "B", "C"),
                        ncol = 2, nrow=2)

title <- cowplot::ggdraw() + cowplot::draw_label("Correlation plot between ALFAssay and other assays/tools on Pearl data", fontface='bold')
pearl_tf_correlations <- cowplot::plot_grid(title, pearl_tf_correlations, ncol=1, rel_heights=c(0.1, 1))



ggsave("plots/PearlTFCorrelations.pdf",pearl_tf_correlations,width=11, height=8)


## neorhea
data_neorhea <- data_pred %>% filter(study=="Neorhea")
fragle_pred_neo <-  read.csv("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")

fragle_pred_neo$PatientId <- str_replace(fragle_pred_neo$Sample_ID, ".markDup.GCtagged", "")

data_neorhea <- merge(fragle_pred_neo, data_neorhea, by="PatientId")


data_neorhea <- data_neorhea %>% filter(!is.na(ctDNADetected))
data_neorhea <- data_neorhea %>% filter(!is.na(BCFS))
data_neorhea <- data_neorhea %>% filter(!is.na(ALFAssay))


data_neorhea$ctDNA_Burden_surv[data_neorhea$ctDNA_Burden>=0.1]<-"HighctDNA" 
data_neorhea$ctDNA_Burden_surv[data_neorhea$ctDNA_Burden<0.1]<-"LowctDNA"
data_neorhea$VAF_perc <- data_neorhea$VAF*100

data_neorhea_bl <- data_neorhea %>% filter(timepoint=="BL")
data_neorhea_c1d28 <- data_neorhea %>% filter(timepoint=="C1")
data_neorhea_surgery <- data_neorhea %>% filter(timepoint=="Surgery")


cor(data_neorhea$ichorTF, data_neorhea$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(data_neorhea$ctDNA_Burden, data_neorhea$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(data_neorhea$VAF_perc, data_neorhea$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))


data_neorhea_ichor <- ggscatter(data_neorhea, x = "ichorTF", y = "Predictedlabel",
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.05, label.y = 0.5) + labs(x = "ichorCNA tumor fraction", y = "ALFAssay tumor fraction")


data_neorhea_fragle <- ggscatter(data_neorhea, x = "ctDNA_Burden", y = "Predictedlabel",
                               add = "reg.line",  # Add regressin line
                               add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                               conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.02, label.y = 0.50) + labs(x = "Fragle tumor fraction", y = "ALFAssay tumor fraction")


data_neorhea_Radar <- ggscatter(data_neorhea, x = "VAF_perc", y = "Predictedlabel",
                                   add = "reg.line",  # Add regressin line
                                   add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                                   conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "RaDaR tumor fraction %", y = "ALFAssay tumor fraction")


neorhea_tf_correlations <- cowplot::plot_grid(data_neorhea_ichor,   data_neorhea_fragle, data_neorhea_Radar,
                                            labels = c("A", "B", "C"),
                                            ncol = 2, nrow=2)

title <- cowplot::ggdraw() + cowplot::draw_label("Correlation plot between ALFAssay and other assays/tools on Neorhea data", fontface='bold')
neorhea_tf_correlations <- cowplot::plot_grid(title, neorhea_tf_correlations, ncol=1, rel_heights=c(0.1, 1))


ggsave("plots/NeorheaTFCorrelations.pdf",neorhea_tf_correlations,width=11, height=8)


# snr
pred_snr <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_train_SNR_healthy_ichorTF.csv", sep="\t")
data_snr <- data %>% filter(study=="Synergy")
data_pred_snr <- merge(data_snr,pred_snr, by="PatientId")

data_pred_snr$Predictedlabel <- round(data_pred_snr$Predictedlabel/100,6)

data_pred_snr$ALFAssay[data_pred_snr$Predictedlabel>=0.03]<-"Positive"
data_pred_snr$ALFAssay[data_pred_snr$Predictedlabel<0.03]<-"Negative" 

fragle_pred_snr <-  read.csv("/Users/alexandra/PhD/SynergyStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")

fragle_pred_snr$PatientId <- str_replace(fragle_pred_snr$Sample_ID, ".markDup.GCtagged", "")

data_snr <- merge(fragle_pred_snr, data_pred_snr, by="PatientId")


data_snr$ctDNA_Burden_surv[data_snr$ctDNA_Burden>=0.01]<-"HighctDNA" 
data_snr$ctDNA_Burden_surv[data_snr$ctDNA_Burden<0.01]<-"LowctDNA"


data_snr_bl <- data_snr %>% filter(timepoint=="BL")
data_snr_c1d28 <- data_snr %>% filter(timepoint=="C1")
data_snr_surgery <- data_snr %>% filter(timepoint=="Surgery")


cor(data_snr$ichorTF, data_snr$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
cor(data_snr$ctDNA_Burden, data_snr$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))


data_snr_ichor <- ggscatter(data_snr, x = "ichorTF", y = "Predictedlabel",
                                add = "reg.line",  # Add regressin line
                                add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "ichorCNA tumor fraction", y = "ALFAssay tumor fraction") 


data_snr_fragle <- ggscatter(data_snr, x = "ctDNA_Burden", y = "Predictedlabel",
                                 add = "reg.line",  # Add regressin line
                                 add.params = list(color = "darkblue", fill = "lightgray"), # Customize reg. line
                                 conf.int = TRUE # Add confidence interval
) + stat_cor(method = "spearman", label.x = 0.1, label.y = 0.55) + labs(x = "Fragle tumor fraction", y = "ALFAssay tumor fraction")



snr_tf_correlations <- cowplot::plot_grid(data_snr_ichor,   data_snr_fragle,
                                              labels = c("A", "B"),
                                              ncol = 2, nrow=1)

title <- cowplot::ggdraw() + cowplot::draw_label("Correlation plot between ALFAssay and other assays/tools on Synergy data", fontface='bold')
snr_tf_correlations <- cowplot::plot_grid(title, snr_tf_correlations, ncol=1, rel_heights=c(0.1, 1))

ggsave("plots/SynergyTFCorrelations.pdf",snr_tf_correlations,width=13, height=8)



