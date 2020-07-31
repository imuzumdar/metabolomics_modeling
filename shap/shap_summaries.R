library(SHAPforxgboost)
library(rmarkdown)

final_xgb_obj = final_model$finalModel
y_data = xgb_wide_full_85_imp_final_data$treatment

# creating X_data and removing response
X_data = xgb_wide_full_85_imp_final_data
X_data$treatment = NULL

# creating shap values
shap_values <- shap.values(xgb_model = final_xgb_obj, X_train = X_data)
shap_values$mean_shap_score

# show that rowSum is the output 
shap_data <- shap_values$shap_score
shap_data[, BIAS := shap_values$BIAS0]
pred_mod <- predict(final_xgb_obj, as.matrix(X_data), ntreelimit = 25)
shap_data[, `:=`(rowSum = round(rowSums(shap_data),6), pred_mod = round(pred_mod,6))]
shap_val_table = rmarkdown::paged_table(shap_data[1:14,])


# To prepare the long-format data:
shap_long <- shap.prep(xgb_model = final_xgb_obj, X_train = X_data)

#shap plot
shap.plot.summary(shap_long)

# dilluted points
shap.plot.summary(shap_long, x_bound  = 1.2, dilute = 10)

# examples of partial dependence plots
g1 <- shap.plot.dependence(data_long = shap_long, x = 'metabolite_262.1', y = 'metabolite_42.1', 
                           color_feature = 'Column_WV') + 
  ggtitle("(A) SHAP values of 42.1 vs. measurements of 262.1")

g1


# force plot
plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score)
shap.plot.force_plot(plot_data)

# without zoom
shap.plot.force_plot(plot_data, zoom_in = FALSE)

