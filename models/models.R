library(xgboost)
library(caret)
library(caTools)
library(e1071)
library(ElemStatLearn)

# simple preprocessing to convert time to factor variable
long_full_85_imp_data$time = as.factor(long_full_85_imp_data$time)
long_full_90_imp_data$time = as.factor(long_full_90_imp_data$time)



# defining function for creating grouped k-folds
group_cv = function(x, k = length(unique(x))) {
  dat = data.frame(index = seq(along = x), group = x)
  groups = data.frame(group = unique(dat$group))
  group_folds = createFolds(groups$group, returnTrain = TRUE, k = k)
  group_folds = lapply(group_folds, function(x, y) y[x,,drop = FALSE], y = groups)
  dat_folds = lapply(group_folds, function(x, y) merge(x, y), y = dat)
  lapply(dat_folds, function(x) sort(x$index))
}

folds = group_cv(long_full_85_imp_data$id)

# RANDOM FOREST WITH LOOCV
train.control.loocv = trainControl(method = "LOOCV")
train.control.groupedKfolds = trainControl(method = "cv", index = folds)

ncol_85 = ncol(long_full_85_imp_data)
ncol_90 = ncol(long_full_90_imp_data)

rf_params = expand.grid(mtry = c(2))

rf_85_mfa = train(time ~., data = long_full_85_imp_data[c(2:ncol_85)], method = "rf", 
                  trControl = train.control.groupedKfolds, tuneLength = 15)

# finding two predictors in mtry
rf_85_mfa_2_pred = train(time ~., data = long_full_85_imp_data[c(2:ncol_85)], method = "rf", 
                  tuneGrid = rf_params, trControl = train.control.groupedKfolds)

print(rf_85_mfa)
predictors(rf_85_mfa)
plot(rf_85_mfa)

# 2 pred
print(rf_85_mfa_2_pred)
predictors(rf_85_mfa_2_pred)
plot(rf_85_mfa_2_pred)
plot(varImp(rf_85_mfa_2_pred), top = 2)

# var importance for random forest model
rf_varImp = varImp(rf_85_mfa)
rf_feature_importance = rf_varImp$importance
plot(rf_varImp, top = 30)

rf_most_important_features = rf_feature_importance[order(-rf_feature_importance$Overall),]

### reducing number of predictors
#preprocessing
rf_feature_importance$metabolites = rownames(rf_feature_importance)
rf_feature_imp_ordered = rf_feature_importance[order(-rf_feature_importance$Overall),]

# 30 metabolites plus treatment
treatment_var = which(rf_feature_imp_ordered$metabolites == "treatment1")
most_imp_features = rf_feature_imp_ordered$metabolites[c(1:30)]

reduced_data_rf_85 = long_full_85_imp_data[,c(most_imp_features, "time", "treatment")]

#reduced model
rf_85_mfa_reduced = train(time ~., data = reduced_data_rf_85, method = "rf", 
                  trControl = train.control.groupedKfolds, tuneLength = 15)

print(rf_85_mfa_reduced) # 100 percent accuracy, wow!!!!
predictors(rf_85_mfa)
plot(rf_85_mfa_reduced)

## reducing to only two predictors
two_most_imp_features = rf_feature_imp_ordered$metabolites[c(1:2)]

reduced_data_rf_85_two_pred = long_full_85_imp_data[,c(two_most_imp_features, "time", 
                                                       "treatment")]

# parsimonious model
rf_85_mfa_reduced_two_pred = train(time ~., data = reduced_data_rf_85_two_pred, method = "rf", 
                                   trControl = train.control.groupedKfolds)

print(rf_85_mfa_reduced_two_pred) # 100 percent accuracy, wow!!!!
predictors(rf_85_mfa_reduced_two_pred)
plot(rf_85_mfa_reduced_two_pred)

plot(varImp(rf_85_mfa_reduced_two_pred))

## reducing to four predictors and removing lidocaine
four_most_imp_features_no_lidocaine = rf_feature_imp_ordered$metabolites[c(2:5)]

reduced_data_rf_85_four_pred_minus_lidocaine = 
  long_full_85_imp_data[,c(four_most_imp_features_no_lidocaine, "time", "treatment")]

# parsimonious model
rf_85_mfa_reduced_four_pred_minus_lidocaine = 
  train(time ~., data = reduced_data_rf_85_four_pred_minus_lidocaine, method = "rf", 
                                   trControl = train.control.groupedKfolds)

print(rf_85_mfa_reduced_four_pred_minus_lidocaine) # 100 percent accuracy, wow!!!!
predictors(rf_85_mfa_reduced_four_pred_minus_lidocaine)
plot(rf_85_mfa_reduced_four_pred_minus_lidocaine)

plot(varImp(rf_85_mfa_reduced_four_pred_minus_lidocaine))

# four metabolites in this model
metabolites[c(546,558,205,252)]

# SVM MODEL FOR PREDICTING GROUP
train.control.groupedKfolds = trainControl(method = "cv", index = folds)

svm_Linear_85 = train(treatment ~., 
                    data = wide_full_85_imp_data[c(2:ncol(wide_full_85_imp_data))], 
                    method = "svmLinear2",
                    trControl = trainControl(method = "LOOCV"),
                    preProcess = c("center", "scale"),
                    tuneLength = 15)

print(svm_Linear_85)

svm_poly_85 = train(treatment ~., 
                    data = wide_full_85_imp_data[c(2:ncol(wide_full_85_imp_data))], 
                    method = "svmPoly",
                    trControl = trainControl(method = "LOOCV"),
                    preProcess = c("center", "scale"))

print(svm_poly_85)

### concerned about over-fitting, using anova_tests_treatment results to identify metabolites
  # significant in determining treatment vs control group (including sig interactions).

sig_metabolites_treatment = which(anova_tests_treatment$treatment_significance == "yes")
sig_metabolites_interaction = which(anova_tests_treatment$interaction_significance == "yes")

sig_metabolites_treatment_and_interaction = unique(c(sig_metabolites_treatment,
                                                     sig_metabolites_interaction))

treatment_metabolites = paste("metabolite", sig_metabolites_treatment_and_interaction, sep = "_")

# significant and has more than 85 percent of values
intersect_columns = which(treatment_metabolites %in% colnames(long_full_85_imp_data))

subset_columns = c("treatment", "time", treatment_metabolites[intersect_columns])
subset_long_full_85_imp_data = long_full_85_imp_data[, subset_columns]

# doing it for wide data
wide_subset_metabolites_1 = paste(treatment_metabolites[intersect_columns], "1", sep = ".")
wide_subset_metabolites_2 = paste(treatment_metabolites[intersect_columns], "2", sep = ".")
wide_subset_metabolites_3 = paste(treatment_metabolites[intersect_columns], "3", sep = ".")
wide_subset_metabolites_4 = paste(treatment_metabolites[intersect_columns], "4", sep = ".")

wide_subset_metabolites = c(wide_subset_metabolites_1, wide_subset_metabolites_2,
                            wide_subset_metabolites_3, wide_subset_metabolites_4)

wide_subset_columns = c("treatment", wide_subset_metabolites)
subset_wide_full_85_imp_data = wide_full_85_imp_data[, wide_subset_columns]

# saving data to load for python
wide_data = write.xlsx(subset_wide_full_85_imp_data, "wide_data_xgb.xlsx")
long_data = write.xlsx(subset_long_full_85_imp_data, "long_data_xgb.xlsx")
                                                       
# building new SVM
svm_Linear_85_reduced = train(treatment ~., data = subset_long_full_85_imp_data, 
                              method = "svmLinear2",
                              trControl = train.control.groupedKfolds,
                              preProcess = c("center", "scale"))

print(svm_Linear_85_reduced)

svm_poly_85_reduced = train(treatment ~., data = subset_long_full_85_imp_data, 
                              method = "svmPoly",
                              trControl = train.control.groupedKfolds,
                              preProcess = c("center", "scale"))

print(svm_poly_85_reduced)
plot(varImp(svm_poly_85_reduced))

# trying using wide format to yield usable kappa statistic
svm_Linear_85_reduced_wide = train(treatment ~., data = subset_wide_full_85_imp_data, 
                              method = "svmLinear2",
                              trControl = trainControl(method = "repeatedcv", number = 7, 
                                                      repeats = 3),
                              preProcess = c("center", "scale"))

print(svm_Linear_85_reduced_wide)

svm_poly_85_reduced_wide = train(treatment ~., data = subset_wide_full_85_imp_data, 
                            method = "svmPoly",
                            trControl = trainControl(method = "repeatedcv", number = 7, 
                                                     repeats = 3),
                            preProcess = c("center", "scale"))

print(svm_poly_85_reduced_wide) # 85.7% accuracy with .714 kappa!!!! HUGE IMPROVEMENT

# useful metabolites
treatment_metabolites[intersect_columns]
metabolites[c(11, 42, 104, 233, 262, 270, 275, 339, 387, 551, 29, 86, 212, 245, 341, 343, 411, 
              450, 509)]

# TRY RANDOM FOREST MODEL TO PREDICT TREATMENT
rf_85_reduced_wide = train(treatment ~., data = subset_wide_full_85_imp_data, 
                           method = "rf", 
                           trControl = trainControl(method = "repeatedcv", number = 7, 
                                                    repeats = 3), 
                           tuneLength = 15)

rf_85_reduced = train(treatment ~., data = subset_long_full_85_imp_data, 
                      method = "rf", 
                      trControl = train.control.groupedKfolds, 
                      tuneLength = 15)

print(rf_85_reduced_wide)
plot(varImp(rf_85_reduced_wide))

print(rf_85_reduced)
plot(varImp(rf_85_reduced))

# XGBOOST MODEL TO PREDICT TREATMENT
xgb_85_reduced_wide = train(treatment ~., data = subset_wide_full_85_imp_data, 
                                 method = "xgbTree",
                                 trControl = trainControl(method = "repeatedcv", number = 7, 
                                                          repeats = 3, verboseIter = TRUE),
                                 preProcess = c("center", "scale"))

print(xgb_85_reduced_wide) # 80 percent accurate, decent but could be better. Moving to python 
                           # jupyter notebook to try and calculate shapley values

trellis.par.set(caretTheme())
plot(xgb_85_reduced_wide, highlight = TRUE)
ggplot(xgb_85_reduced_wide)
plot(xgb_85_reduced_wide, metric = "Kappa", plotType = "level",scales = list(x = list(rot = 90)))

plot(varImp(xgb_85_reduced_wide))

# only metabolite 343, 262, 86, 42, 104, 387

# further subsetting data for calculations
new_mb = paste("metabolite", c(343, 262, 86, 42, 104, 387), sep = "_")
new_wide_subset_metabolites_1 = paste(new_mb, "1", sep = ".")
new_wide_subset_metabolites_2 = paste(new_mb, "2", sep = ".")
new_wide_subset_metabolites_3 = paste(new_mb, "3", sep = ".")
new_wide_subset_metabolites_4 = paste(new_mb, "4", sep = ".")

new_wide_subset_metabolites = c(new_wide_subset_metabolites_1, new_wide_subset_metabolites_2,
                                new_wide_subset_metabolites_3, new_wide_subset_metabolites_4)

new_wide_subset_columns = c("treatment", new_wide_subset_metabolites)
xgb_wide_full_85_imp_data = subset_wide_full_85_imp_data[, new_wide_subset_columns]

xgb_85_new_reduced_wide = train(treatment ~., data = xgb_wide_full_85_imp_data, 
                                method = "xgbTree",
                                trControl = trainControl(method = "repeatedcv", number = 7, 
                                                     repeats = 3, verboseIter = TRUE),
                                preProcess = c("center", "scale"))

print(xgb_85_new_reduced_wide)

plot(varImp(xgb_85_new_reduced_wide))

# last effort for model reduction 
varImp_features_xgb = varImp(xgb_85_new_reduced_wide)$importance
imp_features_xgb = rownames(varImp_features_xgb)[which(varImp_features_xgb$Overall > 0)]

# matching params from both previous xgb models. for time, assuming those are population params
xgb_grid <-  expand.grid(nrounds = 50, 
                        max_depth = 2, 
                        eta = .4,
                        gamma = 0,
                        min_child_weight = 1,
                        subsample = .75,
                        colsample_bytree = .8)

# tuning nrounds and max_depth
xgb_grid2 <-  expand.grid(nrounds = c(25, 50, 75, 100, 125, 150, 175, 200), 
                         max_depth = c(1, 2, 3, 4, 5, 6), 
                         eta = .4,
                         gamma = 0,
                         min_child_weight = 1,
                         subsample = .75,
                         colsample_bytree = .8)

final_wide_subset_columns = c("treatment", imp_features_xgb)
xgb_wide_full_85_imp_final_data = subset_wide_full_85_imp_data[, final_wide_subset_columns]

xgb_85_very_reduced_wide = train(treatment ~., data = xgb_wide_full_85_imp_final_data, 
                                 method = "xgbTree",
                                 trControl = trainControl(method = "repeatedcv", number = 7, 
                                                          repeats = 3, verboseIter = TRUE),
                                 preProcess = c("center", "scale"),
                                 tuneGrid = xgb_grid)

xgb_85_very_reduced_wide2 = train(treatment ~., data = xgb_wide_full_85_imp_final_data, 
                                 method = "xgbTree",
                                 trControl = trainControl(method = "cv", number = 7, 
                                                          verboseIter = TRUE),
                                 preProcess = c("center", "scale"),
                                 tuneGrid = xgb_grid2)

repeated_train = trainControl(method = "repeatedcv", number = 7, repeats = 3, verboseIter = TRUE)

final_model = train(treatment ~., data = xgb_wide_full_85_imp_final_data, 
                    method = "xgbTree",
                    trControl = repeated_train,
                    preProcess = c("center", "scale"),
                    tuneGrid = xgb_grid2)

# model 1 and 2 output
print(xgb_85_very_reduced_wide)
print(xgb_85_very_reduced_wide2) # 92 percent accuracy!!!! nice.

# model 2 plots
trellis.par.set(caretTheme())
plot(xgb_85_very_reduced_wide2, highlight = TRUE)
ggplot(xgb_85_very_reduced_wide2)
plot(xgb_85_very_reduced_wide2, metric = "Kappa", plotType = "level",
     scales = list(x = list(rot = 90)))

plot(varImp(xgb_85_very_reduced_wide2))

# final model check
print(final_model)
plot(final_model, highlight = TRUE)
ggplot(final_model)
plot(final_model, metric = "Kappa", plotType = "level",
     scales = list(x = list(rot = 90)))
plot(varImp(final_model))

