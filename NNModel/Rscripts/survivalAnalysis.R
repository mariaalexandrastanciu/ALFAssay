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
source("prepare_data.R")


try = "_t26"

data <- get_data_2(try)


# ichor pearl

data_pearl_bl <- data$prl$bl_pearl #data_pearl %>% filter(timepoint=="BL")
data_pearl_d15 <- data$prl$d15_pearl #data_pearl %>% filter(timepoint=="C1")
data_pearl_pd <- data$prl$pd_pearl # data_pearl %>% filter(timepoint=="Surgery")

ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "PFS", "ichorPrediction","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "PFS", "ichorPrediction","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "PFS", "ichorPrediction","PD")

library(ggtext)
# caption_text <- "KM on Pearl data using ichorCNA results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)
# 
# p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
#                          # rremove("x.text"), 
#                          labels = c("C", "D"),
#                          ncol = 2)
# p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 20, 5, 10)))+
#   labs(caption = caption_text)
# 
# p3 <- cowplot::plot_grid(p1, p2,
#                          # rremove("x.text"), 
#                          # labels = c("A", "B"),
#                          nrow = 2)


ggsave(paste0("plots",try, "/PearlIchorSurvival.pdf",sep=""), p1, width=12, height=4)


## alfassay pearl
ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "PFS", "ALFAssay","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "PFS", "ALFAssay","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "PFS", "ALFAssay","PD")

library(ggtext)
# caption_text <- "KM on Pearl data using ALFAssay results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)
final_plot <- plot_grid(
  ggdraw() + draw_label("KM on Pearl data using ALFAssay results", fontface = 'bold', x = 0.5, hjust = 0.5, size = 14),
  p1,
  ncol = 1,
  rel_heights = c(0.1, 1)  # Adjust title height
)

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


ggsave(paste0("plots", try, "/PearlALFAssayTFSurvival.pdf",sep=""), p1, width=12, height=4)
ggsave(paste0("plots", try, "/plots_progress/PearlALFAssayTFSurvival.pdf",sep=""), p1, width=12, height=4)

## pearl fragle

ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "PFS", "FragleDetection","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "PFS", "FragleDetection","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "PFS", "FragleDetection","PD")

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


ggsave(paste0("plots", try, "/PearlFragleSurvival.pdf", sep=""), p3, width=12, height=9)


## pearl oncofollow

ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "PFS", "ctDNADetected","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "PFS", "ctDNADetected","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "PFS", "ctDNADetected","PD")

library(ggtext)
caption_text <- "KM on Pearl data using OncoDNA results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)

final_plot <- plot_grid(
  ggdraw() + draw_label("KM on Pearl data using OncoDNA results", fontface = 'bold', x = 0.5, hjust = 0.5, size = 14),
  p1,
  ncol = 1,
  rel_heights = c(0.1, 1)  # Adjust title height
)


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


ggsave(paste0("plots", try, "/PearlOncoDNASurvival.pdf", sep=""), p1, width=12, height=4)
ggsave(paste0("plots", try, "/plots_progress/PearlOncoDNASurvival.pdf", sep=""), final_plot, width=12, height=4)



## alfassay and onco on pearl

pearl <- data$prl$all_pearl
pearl$OncoALFA[pearl$ctDNADetected=="Negative" & pearl$ALFAssay=="Negative"]<-"Negative" 
pearl$OncoALFA[pearl$ctDNADetected=="Negative" & pearl$ALFAssay=="Positive"]<-"NegativePositive"
pearl$OncoALFA[pearl$ctDNADetected=="Positive" & pearl$ALFAssay=="Negative"]<-"PositiveNegative" 
pearl$OncoALFA[pearl$ctDNADetected=="Positive" & pearl$ALFAssay=="Positive"]<-"Positive"

data_pearl_bl <- pearl %>% filter(timepoint=="BL")
data_pearl_d15 <- pearl%>% filter(timepoint=="D15")
data_pearl_pd <- pearl%>% filter(timepoint=="PD")

ichor <- list()

ichor[[1]] <- survival_plot(data_pearl_bl, "PFS", "OncoALFA","Baseline")
ichor[[2]] <- survival_plot(data_pearl_d15, "PFS", "OncoALFA","D15")
ichor[[3]] <- survival_plot(data_pearl_pd, "PFS", "OncoALFA","PD")

library(ggtext)
# caption_text <- "KM on Pearl data using ichorCNA results"

p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                         # rremove("x.text"), 
                         labels = c("A", "B"),
                         ncol = 2)


ggsave(paste0("plots",try, "/PearlOncoALFASurvival.pdf",sep=""), p1, width=12, height=4)


## radar neorhea 

data_neorhea_bl <- data$nrh$bl_neorhea #data_neorhea %>% filter(timepoint=="BL")
data_neorhea_c1d28 <- data$nrh$c1d28_neorhea #data_neorhea %>% filter(timepoint=="C1")
data_neorhea_surgery <- data$nrh$surhery_neorhea #data_neorhea %>% filter(timepoint=="Surgery")

ichor <- list()

ichor[[1]] <- survival_plot(data_neorhea_bl, "BCFS", "ctDNADetected","Baseline")
ichor[[2]] <- survival_plot(data_neorhea_c1d28, "BCFS", "ctDNADetected","C1D28")
ichor[[3]] <- survival_plot(data_neorhea_surgery, "BCFS", "ctDNADetected","Surgery")

library(ggtext)
caption_text <- "KM on Neorhea data using RaDaR results"

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


ggsave(paste0("plots", try, "/NeorheaRaDaRSurvival.pdf",sep=""), p3, width=12, height=9)


# aflassy neorhea

ichor <- list()


ichor[[1]] <- survival_plot(data_neorhea_bl, "BCFS", "ALFAssay","Baseline")
ichor[[2]] <- survival_plot(data_neorhea_c1d28, "BCFS", "ALFAssay","C1D28")
ichor[[3]] <- survival_plot(data_neorhea_surgery, "BCFS", "ALFAssay","Surgery")

library(ggtext)
# caption_text <- "KM on Neorhea data using ALFAssay results"

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


ggsave(paste0("plots", try, "/NeorheaALFAssayTFSurvival.pdf",sep=""), p3, width=12, height=9)
ggsave(paste0("plots", try, "/plots_progress/NeorheaALFAssayTFSurvival.pdf", sep=""), p3, width=12, height=9)



# fragle neorhea  - 
ichor <- list()

ichor[[1]] <- survival_plot(data_neorhea_bl, "BCFS", "FragleDetection","Baseline")
ichor[[2]] <- survival_plot(data_neorhea_c1d28, "BCFS", "FragleDetection","D15")
ichor[[3]] <- survival_plot(data_neorhea_surgery, "BCFS", "FragleDetection","PD")

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


ggsave(paste0("plots", try, "/NeorheaFragleSurvival.pdf",sep=""), p3, width=12, height=9)




### ichor synergy

data_pred_snr_bl <- data$snr$w1_snr # data_pred_snr %>% filter(timepoint=="BL")
data_pred_snr_c1d28 <- data$snr$w3_snr #    data_pred_snr %>% filter(timepoint=="C1")
data_pred_snr_surgery <- data$snr$w13_snr #data_pred_snr %>% filter(timepoint=="Surgery")

ichor <- list()

ichor[[1]] <- survival_plot(data_pred_snr_bl, "PFS", "ichorPrediction","W1")
ichor[[2]] <- survival_plot(data_pred_snr_c1d28, "PFS", "ichorPrediction","W3")
ichor[[3]] <- survival_plot(data_pred_snr_surgery, "PFS", "ichorPrediction","W13")

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


ggsave(paste0("plots", try, "/SynergyichorCNASurvival.pdf", sep=""), p3, width=12, height=9)


## alfassay synergy
ichor <- list()

ichor[[1]] <- survival_plot(data_pred_snr_bl, "PFS", "ALFAssay","W1")
ichor[[2]] <- survival_plot(data_pred_snr_c1d28, "PFS", "ALFAssay","W3")
ichor[[3]] <- survival_plot(data_pred_snr_surgery, "PFS", "ALFAssay","W13")

library(ggtext)
# caption_text <- "KM on Synergy data using ALFAssay results"

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


ggsave(paste0("plots", try, "/SynergyALFAssayTFSurvival.pdf", sep=""), p3, width=12, height=9)
ggsave(paste0("plots", try, "/plots_progress/SynergyALFAssayTFSurvival.pdf", sep=""), p3, width=12, height=9)

## fragle synergy

ichor <- list()

ichor[[1]] <- survival_plot(data_pred_snr_bl, "PFS", "FragleDetection","W1")
ichor[[2]] <- survival_plot(data_pred_snr_c1d28, "PFS", "FragleDetection","W3")
ichor[[3]] <- survival_plot(data_pred_snr_surgery, "PFS", "FragleDetection","W13")

library(ggtext)
caption_text <- "KM on Synergy data using Fragle results"

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


ggsave(paste0("plots", try, "/SynergyFragleSurvival.pdf", sep=""), p3, width=12, height=9)



