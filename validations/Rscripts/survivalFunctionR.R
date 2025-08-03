survival_plot = function(dataset, outcome, eval_var, timepoint_var){
  
  # all_data_timepoint <- all_data %>% filter(timepoint==timepoint_var)
  lst <- list(outcome=dataset[[outcome]],status=dataset$status,eval_var=dataset[[eval_var]])
  surv_table <- as.data.frame(do.call(cbind,lst))
  surv_table$outcome <- as.numeric(surv_table$outcome)
  surv_table$status <- as.numeric(surv_table$status)
  print(surv_table)
  surv_table %<>% filter(!is.na(outcome))
  
  eval_labels = unique(surv_table$eval_var)
  
  surv_dif_p_value <- paste0( "p = ", round((survdiff(formula = Surv(outcome,status) ~ eval_var, data = surv_table))$p,2))
  print(surv_dif_p_value)
  
  linelistsurv_fit <-  surv_fit(Surv(outcome, status) ~  eval_var, data = surv_table)
  names(linelistsurv_fit$strata) <- gsub("eval_var=", "", names(linelistsurv_fit$strata))
  print(linelistsurv_fit)
  
  plot <- survminer::ggsurvplot(
    linelistsurv_fit, 
    data = surv_table,          # again specify the data used to fit linelistsurv_fit_sex 
    conf.int = FALSE,              # do not show confidence interval of KM estimates
    surv.scale = "percent",        # present probabilities in the y axis in %
    break.time.by = 500,            # present the time axis with an increment of 10 days
    xlab = "Follow-up days",
    ylab = outcome,
    pval = T,                      # print p-value of Log-rank test 
    pval.coord = c(100,.50),        # print p-value at these plot coordinates
    risk.table = T,                # print the risk table at bottom 
    legend.title = timepoint_var,       # legend characteristics
    # legend.labs = linelistsurv_fit$strata,
    font.legend = 10, 
    palette = "Dark2",             # specify color palette 
    surv.median.line = "hv",       # draw horizontal and vertical lines to the median survivals
    ggtheme = theme_light() ,       # simplify plot background
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.y.text = FALSE
  )
  return(plot)
}