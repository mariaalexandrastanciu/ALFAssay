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

predictions <- read.csv("~/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF_t2.csv", sep="\t")
# predictions <- read.csv("~/PhD/PyCharmProjects/ALFAssayNN/plots/results_validation_NRH_PRL_healthy_ichorTF03.csv", sep="\t")


data_pred <- merge(data,predictions, by="PatientId" )
data_pred$Predictedlabel = data_pred$Predictedlabel/100

data_pred$ALFAssay[data_pred$Predictedlabel>=0.03]<-"Positive"
data_pred$ALFAssay[data_pred$Predictedlabel<0.03]<-"Negative" 


# pearl ichorCNA data
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

ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "BCFS", "ichorPrediction","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "BCFS", "ichorPrediction","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "BCFS", "ichorPrediction","PD")

library(ggtext)
caption_text <- "KM on Pearl data using ichorCNA results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3


ggsave("plots_t2/PearlIchorSurvival.pdf", p3, width=12, height=9)



ggplot(data_pearl, aes(x = VAF_pers, y = ichorTF)) +
  geom_point()

ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "BCFS", "ALFAssay","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "BCFS", "ALFAssay","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "BCFS", "ALFAssay","PD")

library(ggtext)
caption_text <- "KM on Pearl data using ALFAssay results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3


ggsave("plots_t2/PearlALFAssayTFSurvival2.pdf", p3, width=12, height=9)


# with(data_pearl, CrossTable(ctDNADetected, ctDNA_Burden_surv))
# cor(data_pearl$ichorTF, data_pearl$ctDNA_Burden,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_pearl$ctDNA_Burden, data_pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_pearl$ctDNA_Burden, data_pearl$VAF,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_pearl$Predictedlabel, data_pearl$VAF,  use = "complete.obs" ,method = c( "spearman"))
# ggplot(data_pearl, aes(x = ctDNA_Burden, y = VAF)) +
#   geom_point()
# 
# ggplot(data_pearl, aes(x = Predictedlabel, y = VAF)) +
#   geom_point()


# data_pearl_bl <- data_pearl %>% filter(timepoint=="BL")
# data_pearl_d15 <- data_pearl %>% filter(timepoint=="C1")
# data_pearl_pd <- data_pearl %>% filter(timepoint=="Surgery")

ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "BCFS", "ctDNA_Burden_surv","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "BCFS", "ctDNA_Burden_surv","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "BCFS", "ctDNA_Burden_surv","PD")

library(ggtext)
caption_text <- "KM on Pearl data using Fragle results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3


ggsave("plots_t2/PearlFragleSurvival.pdf", p3, width=12, height=9)


## neorhea
data_neorhea <- data_pred %>% filter(study=="Neorhea")

data_neorhea <- data_neorhea %>% filter(!is.na(ctDNADetected))
data_neorhea <- data_neorhea %>% filter(!is.na(BCFS))
data_neorhea <- data_neorhea %>% filter(!is.na(ALFAssay))


data_neorhea_bl <- data_neorhea %>% filter(timepoint=="BL")
data_neorhea_c1d28 <- data_neorhea %>% filter(timepoint=="C1")
data_neorhea_surgery <- data_neorhea %>% filter(timepoint=="Surgery")


ggplot(data_neorhea, aes(x = VAF, y = Predictedlabel)) +
  geom_point()

ichor <- list()

ichor[[1]] <- survival_plot(data_neorhea_bl, "BCFS", "ctDNADetected","Baseline")
ichor[[2]] <- survival_plot(data_neorhea_c1d28, "BCFS", "ctDNADetected","C1D28")
ichor[[3]] <- survival_plot(data_neorhea_surgery, "BCFS", "ctDNADetected","Surgery")

library(ggtext)
caption_text <- "KM on Neorhea data using ctDNADetected results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3


ggsave("plots_t2/NeorheactDNADetectedSurvival.pdf", p3, width=12, height=9)


ichor <- list()

ichor[[1]] <- survival_plot(data_neorhea_bl, "BCFS", "ALFAssay","Baseline")
ichor[[2]] <- survival_plot(data_neorhea_c1d28, "BCFS", "ALFAssay","C1D28")
ichor[[3]] <- survival_plot(data_neorhea_surgery, "BCFS", "ALFAssay","Surgery")

library(ggtext)
caption_text <- "KM on Neorhea data using ALFAssay results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3


ggsave("plots_t2/NeorheaALFAssayTFSurvival2.pdf", p3, width=12, height=9)


### fragle neorhea

fragle_pred_neo <-  read.csv("/Users/alexandra/PhD/NeoRheaStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")

fragle_pred_neo$PatientId <- str_replace(fragle_pred_neo$Sample_ID, ".markDup.GCtagged", "")
# 
# data_neorhea<-data_neorhea %>% filter(study=="Neorhea")
# 
# data_neorhea <- data_neorhea %>% filter(!is.na(ichorTF))
# data_neorhea <- data_neorhea %>% filter(!is.na(BCFS))
# data_neorhea <- data_neorhea %>% filter(!is.na(ALFAssay))

data_neorhea <- merge(fragle_pred_neo, data_neorhea, by="PatientId")

data_neorhea$ctDNA_Burden_surv[data_neorhea$ctDNA_Burden>=0.1]<-"HighctDNA" 
data_neorhea$ctDNA_Burden_surv[data_neorhea$ctDNA_Burden<0.1]<-"LowctDNA"

# with(data_neorhea, CrossTable(ctDNADetected, ctDNA_Burden_surv))
# cor(data_neorhea$ichorTF, data_neorhea$ctDNA_Burden,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_neorhea$ctDNA_Burden, data_neorhea$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_neorhea$ctDNA_Burden, data_neorhea$VAF,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_neorhea$Predictedlabel, data_neorhea$VAF,  use = "complete.obs" ,method = c( "spearman"))
# ggplot(data_neorhea, aes(x = ctDNA_Burden, y = VAF)) +
#   geom_point()
# 
# ggplot(data_neorhea, aes(x = Predictedlabel, y = VAF)) +
#   geom_point()


data_neorhea_bl <- data_neorhea %>% filter(timepoint=="BL")
data_neorhea_c1d28 <- data_neorhea %>% filter(timepoint=="C1")
data_neorhea_surgery <- data_neorhea %>% filter(timepoint=="Surgery")

ichor <- list()

ichor[[1]] <- survival_plot(data_neorhea_bl, "BCFS", "ctDNA_Burden_surv","Baseline")
ichor[[2]] <- survival_plot(data_neorhea_c1d28, "BCFS", "ctDNA_Burden_surv","D15")
ichor[[3]] <- survival_plot(data_neorhea_surgery, "BCFS", "ctDNA_Burden_surv","PD")

library(ggtext)
caption_text <- "KM on Neorhea data using Fragle results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3


ggsave("plots_t2/NeorheaFragleSurvival.pdf", p3, width=12, height=9)




### snr predictions

pred_snr <- read.csv("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/plots/results_train_SNR_healthy_ichorTF_t2.csv", sep="\t")
data_snr <- data %>% filter(study=="Synergy")
data_pred_snr <- merge(data_snr,pred_snr, by="PatientId")

data_pred_snr$Predictedlabel <- round(data_pred_snr$Predictedlabel/100,6)

data_pred_snr$ALFAssay[data_pred_snr$Predictedlabel>=0.03]<-"Positive"
data_pred_snr$ALFAssay[data_pred_snr$Predictedlabel<0.03]<-"Negative" 

data_pred_snr_bl <- data_pred_snr %>% filter(timepoint=="BL")
data_pred_snr_c1d28 <- data_pred_snr %>% filter(timepoint=="C1")
data_pred_snr_surgery <- data_pred_snr %>% filter(timepoint=="Surgery")

ichor <- list()

ichor[[1]] <- survival_plot(data_pred_snr_bl, "BCFS", "ctDNADetected","W1")
ichor[[2]] <- survival_plot(data_pred_snr_c1d28, "BCFS", "ctDNADetected","W3")
ichor[[3]] <- survival_plot(data_pred_snr_surgery, "BCFS", "ctDNADetected","W13")

library(ggtext)
caption_text <- "KM on Synergy data using ichorCNA results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3


ggsave("plots_t2/SynergyichorCNASurvival.pdf", p3, width=12, height=9)


## alfassay
ichor <- list()

ichor[[1]] <- survival_plot(data_pred_snr_bl, "BCFS", "ALFAssay","W1")
ichor[[2]] <- survival_plot(data_pred_snr_c1d28, "BCFS", "ALFAssay","W3")
ichor[[3]] <- survival_plot(data_pred_snr_surgery, "BCFS", "ALFAssay","W13")

library(ggtext)
caption_text <- "KM on Synergy data using ALFAssay results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3


ggsave("plots_t2/SynergyALFAssayTFSurvival2.pdf", p3, width=12, height=9)

### fragle snr

fragle_pred_snr <-  read.csv("/Users/alexandra/PhD/SynergyStudy/FragmentationPatterns/fragle/Fragle.csv", sep=",")

fragle_pred_snr$PatientId <- str_replace(fragle_pred_snr$Sample_ID, ".markDup.GCtagged", "")

data_snr <- merge(fragle_pred_snr, data_pred_snr, by="PatientId")


data_snr$ctDNA_Burden_surv[data_snr$ctDNA_Burden>=0.01]<-"HighctDNA" 
data_snr$ctDNA_Burden_surv[data_snr$ctDNA_Burden<0.01]<-"LowctDNA"
# 
# with(data_snr, CrossTable(ctDNADetected, ctDNA_Burden_surv))
# with(data_snr, CrossTable(ctDNADetected, ALFAssay))



data_snr_bl <- data_snr %>% filter(timepoint=="BL")
data_snr_c1d28 <- data_snr %>% filter(timepoint=="C1")
data_snr_surgery <- data_snr %>% filter(timepoint=="Surgery")

ichor <- list()

ichor[[1]] <- survival_plot(data_snr_bl, "BCFS", "ctDNA_Burden_surv","W1")
ichor[[2]] <- survival_plot(data_snr_c1d28, "BCFS", "ctDNA_Burden_surv","W3")
ichor[[3]] <- survival_plot(data_snr_surgery, "BCFS", "ctDNA_Burden_surv","W13")

library(ggtext)
caption_text <- "KM on Neorhea data using Fragle results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                         # rremove("x.text"), 
                         labels = c("C", "D"),
                         ncol = 2)
p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
  labs(caption = caption_text)

p3 <- cowplot::plot_grid(p1, p2,
                         # rremove("x.text"), 
                         # labels = c("A", "B"),
                         nrow = 2)
# gplot <- survminer::arrange_ggsurvplots(splots_ctDNA, print = TRUE, ncol = 2, nrow = 2) # risk.table.height = 0.3


ggsave("plots_t2/SynergyFragleSurvival.pdf", p3, width=12, height=9)





# 
# cor(data_neorhea$ichorTF, data_neorhea$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_pearl$ichorTF, data_pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_pred_snr$ichorTF, data_pred_snr$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
# 
# 
# cor(data_neorhea$VAF, data_neorhea$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_pearl$VAF_pers, data_pearl$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
# cor(data_pred_snr$ichorTF, data_pred_snr$Predictedlabel,  use = "complete.obs" ,method = c( "spearman"))
# 
# 
# library(gmodels)
# data_neorhea <- data_neorhea %>% filter(!is.na(ctDNADetected))
# 
# with(data_neorhea, CrossTable(ctDNADetected, ALFAssay))
# 
# data_neorhea_bl %>% filter(ALFAssay=="Positive") %>% nrow
# 
# data_neorhea_c1d28 %>% filter(ALFAssay=="Positive") %>% nrow
# 
# data_neorhea_surgery %>% filter(ALFAssay=="Positive") %>% nrow
# 
# 
# with(data_pearl, CrossTable(ctDNADetected, ALFAssay))
