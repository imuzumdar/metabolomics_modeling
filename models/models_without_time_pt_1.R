# recall line 156 from models.R (below), which subsets variables based on ANOVA feature selection

# wide_subset_metabolites = c(wide_subset_metabolites_1, wide_subset_metabolites_2,
#                             wide_subset_metabolites_3, wide_subset_metabolites_4)

# we now construct a new subset of columns using only time periods 2, 3, 4 (NO TIME PERIOD 1)
# note: all data is based on 85 percent imputation of the full data set and is wide data.
# We exclude these facts from variable names for ease of naming.
new_wide_subset = c(wide_subset_metabolites_2, wide_subset_metabolites_3, 
                    wide_subset_metabolites_4)

new_wide_subset_columns = c("treatment", new_wide_subset)
no_tp1_data = wide_full_85_imp_data[, new_wide_subset_columns] #no_tp1 = no time period 1

# specify train control function
repeated_train = trainControl(method = "repeatedcv", number = 7, repeats = 3, verboseIter = TRUE)

# tuneGrid
tune_grid_no_tp1 = expand.grid(
                                nrounds = seq(from = 25, to = 200, by = 25),
                                eta = c(0.05, 0.1, 0.3, 0.4),
                                max_depth = c(1, 2, 3, 4, 5, 6),
                                gamma = 0,
                                colsample_bytree = .8,
                                min_child_weight = 1,
                                subsample = .75
                              )

# now we build first xgb model (NO TUNE GRID)
xgb_no_tp1 = train(treatment ~., data = no_tp1_data, method = "xgbTree",
                   trControl = repeated_train, preProcess = c("center", "scale"))

print(xgb_no_tp1) # 81 percent accuracy, sad drop :(
trellis.par.set(caretTheme())
plot(xgb_no_tp1, highlight = TRUE)
ggplot(xgb_no_tp1)
plot(xgb_no_tp1, metric = "Kappa", plotType = "level", scales = list(x = list(rot = 90)))
plot(varImp(xgb_no_tp1))

# now we build second xgb model (YES TUNE GRID)
xgb_no_tp1.2 = train(treatment ~., data = no_tp1_data, method = "xgbTree", 
                   tuneGrid = tune_grid_no_tp1,
                   trControl = repeated_train, preProcess = c("center", "scale"))

print(xgb_no_tp1.2) # 76 percent accuracy, even sadder drop :(
trellis.par.set(caretTheme())
plot(xgb_no_tp1.2, highlight = TRUE)
ggplot(xgb_no_tp1.2)
plot(xgb_no_tp1.2, metric = "Kappa", plotType = "level", scales = list(x = list(rot = 90)))
plot(varImp(xgb_no_tp1.2))

# model reduction 
varImp_features_xgb_no_tp1 = varImp(xgb_no_tp1)$importance
imp_features_xgb_no_tp1 = 
  rownames(varImp_features_xgb_no_tp1)[which(varImp_features_xgb_no_tp1$Overall > 0)]

no_tp1_columns = c("treatment", imp_features_xgb_no_tp1)
no_tp1_reduced_data = no_tp1_data[, no_tp1_columns]

# tuning nrounds and max_depth
tg_reduced_no_tp1 <-  expand.grid(nrounds = c(25, 50, 75, 100, 125, 150, 175, 200), 
                                  max_depth = c(1, 2, 3, 4, 5, 6), 
                                  eta = .4,
                                  gamma = 0,
                                  min_child_weight = 1,
                                  subsample = .75, 
                                  colsample_bytree = .8)

xgb_reduced_no_tp1 = train(treatment ~., data = no_tp1_reduced_data, method = "xgbTree", 
                           tuneGrid = tg_reduced_no_tp1,
                           trControl = repeated_train, preProcess = c("center", "scale"))

print(xgb_reduced_no_tp1)
plot(xgb_reduced_no_tp1, highlight = TRUE)
ggplot(xgb_reduced_no_tp1)
plot(xgb_reduced_no_tp1, metric = "Kappa", plotType = "level", scales = list(x = list(rot = 90)))
plot(varImp(xgb_reduced_no_tp1))

# FINAL model reduction 
varImp_final_no_tp1 = varImp(xgb_reduced_no_tp1)$importance
imp_features_final_no_tp1 = rownames(varImp_final_no_tp1)[which(varImp_final_no_tp1$Overall > 0)]

final_no_tp1_columns = c("treatment", imp_features_final_no_tp1)
final_no_tp1_data = no_tp1_data[, final_no_tp1_columns]

xgb_final_no_tp1 = train(treatment ~., data = final_no_tp1_data, method = "xgbTree", 
                         tuneGrid = tg_reduced_no_tp1, # same parameter grid as prev model
                         trControl = repeated_train, preProcess = c("center", "scale"))

# model plots
print(xgb_final_no_tp1)
ggplot(xgb_final_no_tp1)
plot(xgb_final_no_tp1, metric = "Kappa", plotType = "level", scales = list(x = list(rot = 90)))
plot(varImp(xgb_final_no_tp1))

# removing 233.2 and training on same parameter grid
fin_no_tp1 = final_no_tp1_data
fin_no_tp1$metabolite_233.2 = NULL

fin_xgb_grid <-  expand.grid(nrounds = c(25, 50, 75, 100, 125, 150, 175, 200), 
                             max_depth = 1, 
                             eta = .4,
                             gamma = 0,
                             min_child_weight = 1,
                             subsample = .75, 
                             colsample_bytree = .8)

xgb_fin_no_tp1 = train(treatment ~., data = fin_no_tp1, method = "xgbTree", 
                         tuneGrid = fin_xgb_grid,
                         trControl = repeated_train, preProcess = c("center", "scale"))

# 88 percent accuracy!!!!
print(xgb_fin_no_tp1)
ggplot(xgb_fin_no_tp1)
plot(xgb_fin_no_tp1, metric = "Kappa", plotType = "level", scales = list(x = list(rot = 90)))
plot(varImp(xgb_fin_no_tp1))


# best treatment model
best_treatment_model = xgb_reduced_no_tp1


