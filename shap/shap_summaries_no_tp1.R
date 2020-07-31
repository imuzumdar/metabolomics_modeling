# same as other shap file but for new model. Could easily automate the contents of this folder
# to work for any xgb model

library(SHAPforxgboost)
library(rmarkdown)

bt_model = best_treatment_model$finalModel #(bt = best treatment)
y_data_bt = no_tp1_reduced_data$treatment

# creating X_data and removing response
X_data_bt = no_tp1_reduced_data
X_data_bt$treatment = NULL

X_data_bt_matrix = as.matrix(X_data_bt)

# creating shap values
shap_values_bt <- shap.values(xgb_model = bt_model, X_train = X_data_bt_matrix)
shap_values_bt$mean_shap_score

# show that rowSum is the output 
shap_data_bt <- shap_values_bt$shap_score
shap_data_bt[, BIAS := shap_values_bt$BIAS0]
pred_mod_bt <- predict(bt_model, X_data_bt_matrix, ntreelimit = 50)
shap_data_bt[, `:=`(rowSum = round(rowSums(shap_data_bt),6), pred_mod_bt = round(pred_mod_bt,6))]
shap_val_table_bt = rmarkdown::paged_table(shap_data_bt[1:14,])


# To prepare the long-format data:
shap_long_bt <- shap.prep(xgb_model = bt_model, X_train = X_data_bt_matrix)

#shap plot
shap.plot.summary(shap_long_bt)

# examples of partial dependence plots
g1_bt <- shap.plot.dependence(data_long = shap_long_bt, x = 'metabolite_86.2', 
                              y = 'metabolite_104.2', 
                           color_feature = 'Column_WV') + 
  ggtitle("(A) SHAP values of 104.2 vs. measurements of 86.2")

g1_bt


# force plot
plot_data_bt <- shap.prep.stack.data(shap_contrib = shap_values_bt$shap_score)
shap.plot.force_plot(plot_data_bt)

# without zoom
shap.plot.force_plot(plot_data_bt, zoom_in = FALSE)

