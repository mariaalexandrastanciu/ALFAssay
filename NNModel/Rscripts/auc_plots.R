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
library(pROC)

setwd("/Users/alexandra/PhD/PyCharmProjects/ALFAssayNN/Rscripts")


try = "_t10"

data <- get_data_2(try)


# ichor pearl
data_pearl <- data$prl$all_pearl

data_pearl_bl <- data$prl$bl_pearl #data_pearl %>% filter(timepoint=="BL")
data_pearl_d15 <- data$prl$d15_pearl #data_pearl %>% filter(timepoint=="C1")
data_pearl_pd <- data$prl$pd_pearl # data_pearl %>% filter(timepoint=="Surgery")

data_pearl$ctDNADetected
data_pearl$oncoDNADetection <- data_pearl$ctDNADetected
data_pearl$oncoDNADetection[data_pearl$oncoDNADetection=="Negative"]<-0
data_pearl$oncoDNADetection[data_pearl$oncoDNADetection=="Positive"]<-1
data_pearl$oncoDNADetection <- as.numeric(data_pearl$oncoDNADetection)
data_pearl$Label <- as.numeric(data_pearl$Label)

controls <- data$ctrl %>% filter(study=="Healthy")
controls$oncoDNADetection <- controls$ctDNADetected
controls$oncoDNADetection[controls$oncoDNADetection=="Negative"]<-0
controls$oncoDNADetection[controls$oncoDNADetection=="Positive"]<-1
controls$oncoDNADetection <- as.numeric(controls$oncoDNADetection)

truth <- c(controls$Label,data_pearl$Label)

response_oncodna <- c(controls$Label,data_pearl$Label)
predicted_oncodna <- c(controls$oncoDNADetection,data_pearl$oncoDNADetection)

roc_obj_onco_dna <- roc(response_oncodna, predicted_oncodna)
auc(roc_obj)

df_oncodna <- data.frame(response_oncodna, predicted_oncodna)
df_oncodna <- na.omit(df_oncodna)
##plot
# pdf(paste0("plots", try, "/plots_progress/oncodna_roc.pdf", sep=""), width = 800, height = 800)
response_oncodna <- factor(df_oncodna$response_oncodna, levels = c(0, 1), labels = c("Negative", "Positive"))

# Create ROCit object
roc_oncodna <- rocit(score = df_oncodna$predicted_oncodna, class = df_oncodna$response_oncodna)


# Plot ROC curve with shaded CI
plot(
  roc_oncodna,
  legend = FALSE,
  YIndex = FALSE,
  line.width = 2,
  col = "#D95F02"  # use a different color from previous one
)


# Add title
title("OncoDNA classification ROC curve on Pearl", cex.main = 1.5)

# Add AUC and CI to plot
auc_val <- roc_oncodna$AUC
auc_text <- sprintf("AUC = %.3f", auc_val)
text(x = 0.4, y = 0.1, labels = auc_text, cex = 1.2, font = 2, col = "#333333", adj = c(0, 0))

dev.off()



data_pearl$ALFAssay[data_pearl$ALFAssay=="Negative"]<-0
data_pearl$ALFAssay[data_pearl$ALFAssay=="Positive"]<-1
data_pearl$ALFAssay <- as.numeric(data_pearl$ALFAssay)

controls$ALFAssay[controls$ALFAssay=="Negative"]<-0
controls$ALFAssay[controls$ALFAssay=="Positive"]<-1
controls$ALFAssay <- as.numeric(controls$ALFAssay)

response_alfassay <- c(controls$Label,data_pearl$Label)
predicted_alfassay <- c(controls$ALFAssay,data_pearl$ALFAssay)

roc_obj_alfassay <- roc(response_alfassay, predicted_alfassay)
auc(roc_obj_alfassay)


## aflassay roc curve
df_aflassay <- data.frame(response_alfassay, predicted_alfassay)
df_aflassay <- na.omit(df_aflassay)
##plot
# pdf(paste0("plots", try, "/plots_progress/alfassay_roc.pdf", sep=""), width = 800, height = 800)
response_aflassay <- factor(df_aflassay$response_alfassay, levels = c(0, 1), labels = c("Negative", "Positive"))

# Create ROCit object
roc_aflassay <- rocit(score = df_aflassay$predicted_alfassay, class = df_aflassay$response_alfassay)


# Plot ROC curve with shaded CI
plot(
  roc_aflassay,
  legend = FALSE,
  YIndex = FALSE,
  line.width = 2,
  col = "#D95F02"  # use a different color from previous one
)


# Add title
title("ALFAssay classification ROC curve on Pearl", cex.main = 1.5)

# Add AUC and CI to plot
auc_val <- roc_aflassay$AUC
auc_text <- sprintf("AUC = %.3f", auc_val)
text(x = 0.4, y = 0.1, labels = auc_text, cex = 1.2, font = 2, col = "#333333", adj = c(0, 0))

dev.off()



## combined
patientID <-  c(controls$PatientId,data_pearl$PatientId)
df <- data.frame(patientID,predicted_alfassay, predicted_oncodna, truth)

# Save to CSV
write.csv(df, "~/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t10/plots_progress/predictions_comparison.csv", row.names = FALSE)



### classifier
# Load required libraries
library(caret)
library(randomForest)
library(pROC)
library(ggplot2)
library(dplyr)

# Read data
df <- read.csv("~/PhD/PyCharmProjects/ALFAssayNN/Rscripts/plots_t10/plots_progress/predictions_comparison.csv")
df <- data.frame(patientID,predicted_alfassay, predicted_oncodna, truth)
df <- na.omit(df)
df$rowID <- seq_len(nrow(df))
# Convert 'truth' to factor with proper levels: Positive = 1, Negative = 0
df$truth <- factor(df$truth, levels = c(0, 1), labels = c("Negative", "Positive"))

# Set seed for reproducibility
set.seed(123)

# Define cross-validation settings
cv_control <- trainControl(
  method = "cv",                   # k-fold CV
  number = 10,                     # 10 folds
  classProbs = TRUE,               # Needed for AUC
  summaryFunction = twoClassSummary, # For ROC-based evaluation
  savePredictions = "final",      # Save out-of-fold predictions
  sampling = "up"                 # Optional: Handle class imbalance (can use "down", "smote", or remove)
)

# Train Random Forest model using cross-validation, optimizing for ROC
rf_model <- train(
  truth ~ predicted_alfassay + predicted_oncodna,
  data = df,
  method = "rf",
  metric = "ROC",
  trControl = cv_control,
  tuneLength = 3
)

# Extract predictions and compute ROC curve
preds <- rf_model$pred
roc_obj <- roc(preds$obs, preds$Positive, levels = c("Negative", "Positive"), direction = "<")

auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

# CI for ROC curve itself
ci_obj <- ci.se(roc_obj, specificities = seq(0, 1, l = 25))

# Plot ROC with CI and flipped axes
plot(
  roc_obj,
  col = "#1c61b6",
  cex.main =1,
  lwd = 2,
  main = "Random Forest ROC Curve with Confidence Interval",
  xlim = c(1, 0),       # Reverse x-axis: 1 → 0
  ylim = c(0, 1),       # Standard y-axis: 0 → 1
  asp = 1,
  legacy.axes = TRUE    # Prevents automatic axis extension
)

# Add shaded confidence interval
plot(ci_obj, type = "shape", col = rgb(0.1, 0.4, 0.9, alpha = 0.2), add = TRUE)

# Add legend with AUC and CI
legend_text <- sprintf("AUC = %.3f [%.3f–%.3f]", auc_val, auc_ci[1], auc_ci[3])
legend("bottomright", legend = legend_text, bty = "n", cex = 1)



## ggplot
library(ggplot2)

# Ensure correct AUC and CI (you should already have this)
auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj)

# Format AUC label
auc_label <- sprintf("AUC = %.3f [%.3f–%.3f]", auc_val, auc_ci[1], auc_ci[3])

# Build ROC dataframe (properly sorted and unique)
roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
) %>%
  distinct() %>%
  arrange(desc(specificity))  # Maintain descending order for correct X-axis

# CI data
ci_vals <- ci.se(roc_obj, specificities = seq(0, 1, l = 100))
ci_df <- data.frame(
  specificity = as.numeric(rownames(ci_vals)),
  lower = ci_vals[, 1],
  upper = ci_vals[, 3]
)

# Plot with fixed axis, labels, and smooth continuous line
ggplot() +
  geom_ribbon(data = ci_df, aes(x = specificity, ymin = lower, ymax = upper),
              fill = "skyblue", alpha = 0.3) +
  geom_line(data = roc_df, aes(x = specificity, y = sensitivity),
            color = "#1c61b6", size = 1.2) +
  coord_cartesian(xlim = c(1, 0), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = "Random Forest ROC Curve",
    x = "Specificity",
    y = "Sensitivity"
  ) +
  annotate(
    "text",
    x = 0.35, y = 0.1,
    label = auc_label,
    hjust = 0,
    size = 5,
    fontface = "bold",
    color = "#333333"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

### rocit

# Load required library
library(ROCit)
# pdf(paste0("plots", try, "/plots_progress/oncodna_alfassay_roc.pdf", sep=""), width = 800, height = 800)
# Extract predicted probabilities and true labels from your trained model
# 'Positive' corresponds to class 1
preds <- rf_model$pred
# If you used cross-validation, there may be multiple rows per sample (from tuning); filter best tune
best_mtry <- rf_model$bestTune$mtry
preds <- preds[preds$mtry == best_mtry, ]

# Create ROCit object using your model's probabilities and ground truth
ROCit_obj <- rocit(score = preds$Positive, class = preds$obs)

# Compute confidence interval
ci_obj <- ciROC(ROCit_obj)

# Plot ROC curve with shaded CI
plot(
  ROCit_obj,
  legend = FALSE,
  YIndex = FALSE,
  line.width = 2,
  col = "darkblue"
)

# Add CI ribbon to the plot
plot(ci_obj, add = TRUE, col = rgb(0.1, 0.4, 0.9, alpha = 0.9))
title("OncoDNA and ALFAAssay classification ROC curve", cex.main = 1.5)
# Extract and format AUC + CI
auc_val <- ROCit_obj$AUC
ci_vals <- ci.auc(roc_obj)
auc_text <- sprintf("AUC = %.3f [%.3f–%.3f]", auc_val, ci_vals[1], ci_vals[2])

# Add AUC label to plot (bottom right)
text(
  x = 0.4, y = 0.1,
  labels = auc_text,
  cex = 1.2,
  font = 2,
  col = "black",
  adj = c(0, 0)
)

  dev.off()
  
  best_mtry <- rf_model$bestTune$mtry
  preds <- rf_model$pred %>% filter(mtry == best_mtry)
  
  # Add back row IDs to match with patient IDs
  # CARET drops extra columns, so we need to merge based on row names
  # Use row index in df as "rowID" and match it to preds$rowIndex
  df_ids <- df %>% select(rowID, patientID)
  
  # Merge to associate predictions with patient ID
  preds_with_ids <- preds %>%
    left_join(df_ids, by = c("rowIndex" = "rowID")) %>%
    select(patientID, Positive)
  
  # Rename output for clarity
  final_predictions_df <- preds_with_ids %>% dplyr::rename(PatientId = patientID, PredictedProbability = Positive)
  
  # View result
  head(final_predictions_df)
  
  data_pearl_rf <- merge(data_pearl, final_predictions_df, by ="PatientId")
  data_pearl_rf$PredictedProbability[data_pearl_rf$PredictedProbability>0.5]<-1
  # data_pearl_rf$PredictedProbability <- factor(data_pearl_rf$PredictedProbability , levels = c(0, 1), labels = c("Negative", "Positive"))
  data_pearl_rf$PredictedProbability[data_pearl_rf$PredictedProbability==0]<-"Negative"
  data_pearl_rf$PredictedProbability[data_pearl_rf$PredictedProbability==1]<-"Positive" 
  data_pearl_bl <- data_pearl_rf %>% filter(timepoint=="BL")
  data_pearl_d15 <- data_pearl_rf %>% filter(timepoint=="D15")
  data_pearl_pd <-  data_pearl_rf %>% filter(timepoint=="PD")
  
  ichor <- list()
  
  ichor[[1]] <- survival_plot(data_pearl_bl, "PFS", "PredictedProbability","Baseline")
  ichor[[2]] <- survival_plot(data_pearl_d15, "PFS", "PredictedProbability","D15")
  ichor[[3]] <- survival_plot(data_pearl_pd, "PFS", "PredictedProbability","PD")
  
  library(ggtext)
  caption_text <- "KM on Pearl data using ALFAssayOncoDNA classifier results"
  
  p1 <- cowplot::plot_grid(ichor[[1]]$plot , ichor[[2]]$plot,
                           # rremove("x.text"), 
                           labels = c("A", "B"),
                           ncol = 2)
  
  p2 <- cowplot::plot_grid(ichor[[3]]$plot, 
                           # rremove("x.text"), 
                           labels = c("C", "D"),
                           ncol = 2)
  
  
  # p2  <- p2 + theme(plot.caption = element_textbox_simple(padding = margin(10, 10, 5, 10)))+
  #   labs(caption = caption_text)
  
  p3 <- cowplot::plot_grid(p1, p2,
                           # rremove("x.text"), 
                           # labels = c("A", "B"),
                           nrow = 2)

  ggsave(paste0("plots",try, "/plots_progress//PearlOncoDNAALFAssaySurvival.pdf",sep=""), p3, width=12, height=9)
  