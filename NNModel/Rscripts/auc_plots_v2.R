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


try = "_t24"

data <- get_data_2(try)


# ichor pearl
data_pearl <- data$prl$all_pearl
# 
# data_pearl_bl <- data$prl$bl_pearl #data_pearl %>% filter(timepoint=="BL")
# data_pearl_d15 <- data$prl$d15_pearl #data_pearl %>% filter(timepoint=="C1")
# data_pearl_pd <- data$prl$pd_pearl # data_pearl %>% filter(timepoint=="Surgery")

data_pearl$ctDNADetected
data_pearl$oncoDNADetection <- data_pearl$ctDNADetected
data_pearl$oncoDNADetection[data_pearl$oncoDNADetection=="Negative"]<-0
data_pearl$oncoDNADetection[data_pearl$oncoDNADetection=="Positive"]<-1
data_pearl$oncoDNADetection <- as.numeric(data_pearl$oncoDNADetection)

model_classif_data <- data_pearl %>% select(PatientId,oncoDNADetection,ichorTF,Predictedlabel)
model_classif_data <- na.omit(model_classif_data)
model_classif_data %<>% dplyr::rename(#"Patient"="Patient ID" ,
  "ALFAssay" = "Predictedlabel" )

write.table( model_classif_data, "~/PhD/PyCharmProjects/ALFAssay/validations/data/classification_data.csv", quote = FALSE, row.names = FALSE)


model_classif_data_ichorTF <- data_pearl %>% select(PatientId,VAF,ichorTF,Predictedlabel)
model_classif_data_ichorTF$VAF <- model_classif_data_ichorTF$VAF/100
model_classif_data_ichorTF$ichorTF[model_classif_data_ichorTF$ichorTF<0.03]<-0
model_classif_data_ichorTF$ichorTF[model_classif_data_ichorTF$ichorTF>=0.03]<-1
model_classif_data_ichorTF <- na.omit(model_classif_data_ichorTF)
model_classif_data_ichorTF %<>% dplyr::rename(#"Patient"="Patient ID" ,
  "ALFAssay" = "Predictedlabel" , "oncoDNADetection" = "VAF")

write.table( model_classif_data_ichorTF, "~/PhD/PyCharmProjects/ALFAssay/validations/data/classification_data_ichorTF.csv", quote = FALSE, row.names = FALSE)


# 
# 
#   # Define cross-validation settings
#   cv_control <- trainControl(
#     method = "cv",                   # k-fold CV
#     number = 10,                     # 10 folds
#     classProbs = TRUE,               # Needed for AUC
#     summaryFunction = twoClassSummary, # For ROC-based evaluation
#     savePredictions = "final",      # Save out-of-fold predictions
#     sampling = "up"                 # Optional: Handle class imbalance (can use "down", "smote", or remove)
#   )
#   
#   # Train Random Forest model using cross-validation, optimizing for ROC
#   rf_model <- train(
#     truth ~ Predictedlabel,
#     data = df,
#     method = "rf",
#     metric = "ROC",
#     trControl = cv_control,
#     tuneLength = 3
#   )
#   
  # # Extract predictions and compute ROC curve
  # preds <- rf_model$pred
  # roc_obj <- roc(preds$obs, preds$Positive, levels = c("Negative", "Positive"), direction = "<")
  # 
  # roc_onco_truth <- roc(data_pearl$oncoDNADetection, data_pearl$ALFAssay)
  # auc(roc_obj)
  # 
  # 
  
  ############################################################
  ## 1.  Setup
  ############################################################
  # Install once if needed:
  # install.packages(c("tidymodels", "ranger"))
  
  library(tidymodels)      # tidymodels meta-package
  library(ranger)          # fast random-forest engine
  
  set.seed(42)             # reproducibility
  
  # ------------------------------------------------------------------
  # df must contain:
  #   - truth : 0/1 (or TRUE / FALSE)  → turn it into a factor
  #   - predictedLabel : numeric
  # ------------------------------------------------------------------
  data_pearl$truth <- data_pearl$oncoDNADetection
  df <- data_pearl %>% select(truth, Predictedlabel)
  df <- na.omit(df)
  df <- df %>% 
    mutate(truth = factor(truth, levels = c(0, 1)))  # make it a factor
  
  
  ############################################################
  ## 2.  One outer train / test split (80 / 20, stratified)
  ############################################################
  split_obj  <- initial_split(df, prop = 0.8, strata = truth)
  train_data <- training(split_obj)
  test_data  <- testing(split_obj)
  
  
  ############################################################
  ## 3.  k-fold CV on training data (here: 5 folds, 3 repeats)
  ############################################################
  folds <- vfold_cv(train_data, v = 5, repeats = 3, strata = truth)
  
  rf_spec_big <- rand_forest(
    mtry  = 1,
    trees = 5000,      # << more trees → finer probabilities
    min_n = 5
  ) %>%
    set_mode("classification") %>%
    set_engine("ranger", importance = "impurity")
  
  rf_big <- workflow() %>%
    add_formula(truth ~ Predictedlabel) %>%
    add_model(rf_spec_big) %>%
    fit(data = train_data)
  
  test_preds <- predict(rf_big, test_data, type = "prob") %>%
    bind_cols(test_data) %>%
    mutate(truth = forcats::fct_relevel(truth, "1"))
  
  ############################################################
  ## 1.  “Exact” ROC  + 95 % CI for AUC & Curve
  ############################################################
  library(pROC)
  roc_obj <- roc(truth ~ .pred_1,
                 data   = test_preds,
                 levels = c("1", "0"),
                 quiet  = TRUE)
  
  auc_ci <- ci.auc(roc_obj, conf.level = 0.95, boot.n = 2000)
  auc_txt <- sprintf("AUC = %.3f (95%% CI %.3f–%.3f)",
                     auc(roc_obj), auc_ci[1], auc_ci[3])
  
  ci_curve <- ci.se(roc_obj,
                    specificities = seq(0, 1, by = 0.01),
                    boot.n        = 1000,
                    conf.level    = 0.95)
  
  spec_grid <- as.numeric(rownames(ci_curve))
  ci_df <- data.frame(
    fpr      = 1 - spec_grid,
    sens_low = ci_curve[, 1],
    sens_mid = ci_curve[, 2],
    sens_hi  = ci_curve[, 3]
  )
  
  roc_df <- data.frame( fpr = (1 - roc_obj$specificities), tpr = roc_obj$sensitivities)
  

  ############################################################
  ## 3.  Plot exact CI ribbon + smoothed curve
  ############################################################
  library(ggplot2)
  
  ggplot(roc_df, aes(x = fpr, y = tpr)) +geom_abline(linetype = "dashed") 
  
  test <- ggplot() +
    # Exact step ROC
    geom_step(
      data = roc_df,
      aes(x = fpr, y = tpr),
      direction = "hv",
      linewidth = 1.2
    ) +
    geom_abline(linetype = "dashed") +
    coord_equal() +
    labs(
      title    = "Exact ROC curve with 95 % confidence band",
      subtitle = auc_txt,
      x        = "False Positive Rate (1 – Specificity)",
      y        = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal(base_size = 13)
  ggsave(paste0("plots",try, "/plots_progress//test.pdf",sep=""), test, width=12, height=9)
  
  
  
