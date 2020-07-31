# metabolomics_modeling
This project seeks to identify metabolites which predict ischemic/reperfusive events and determine if humanin treatment significantly alters blood metabolomics. The dataset contains information on 736 metabolomic measurements over four times points in 14 animal subjects. Seven of these subjects belong to a humanin treatment group; the rest belong to a control group.

There are two main phases to this research:

1. Predict whether subjects are in treatment or control group based on metabolomic data.
2. Predict stage of ischemic events using metabolite measurements across time series data.

Some cool outcomes of this project:

1. Automated all feature engineering and parametric test scripts
2. Identified 5 metabolites which predict ischmeic events extremely well, leading to potential +50k cost savings
3. Developed ANOVA-based feature selection method which improved base model accuracy from 35 to 85.7 percent (useful tool to prevent overfitting on small datasets!)

## Treatment vs. Control Group Analysis

To understand the relationship between the humanin treatment and various metabolite measurements, we built gradient boosted tree (GBT) models to predict whether a given pig was in the treatment or control group, based on their metabolite concentrations at four time points. In building these models, we first performed ANOVA based feature selection, a novel feature engineering technique (currently evaluating its benefits to other datasets and model instances) which selects variables based on variational analyses. For each metabolite, we fit an ANOVA model that considered the main effects of the treatment and time period variables, as well as their interaction. If the treatment variable was significant as either a main effect or an interaction, it was included in the GBT model. Otherwise, it was discarded. This procedure reduced the dimensionality of the dataset from 736 metabolites to 19.

Following feature selection, we developed GBT models using the 19 previously selected metabolites to predict whether pigs were in the treatment or control group. We performed repeated cross validation with 7 folds and 3 repeated draws on a grid of 8 hyperparameters to tune the mode and optimize its accuracy. On each iteration of model fitting, we created feature importance plots and eliminated variables with zero mode importance. We continued this process until all variables were satisfied.

In the end, we found a model using four metabolites which predicted the treatment variable with 85.7% accuracy. For our model, we also calculated several explainability metrics using shap scores, a machine learning interpretability technique derived from game theory, and created plots for them. Shap scores determine the specific contribution that each metabolites makes to each of the models predictions, explaining more thoroughly how the model is making its predictions.

## Ischemic Event Analysis

Before performing any analysis, we had to deal with large amounts of missing information in the provided dataset. With the assumption that this data was missing at randome (MAR), we initially tried to use kNN imputation to estimate the missing values in the dataset. However, each record was missing at least 28 metabolite measurements (out of 736). As a result, we did not have any complete records with which kNN imputation could reference in filling the missing data. Instead, we chose to use MFA multiple-imputation, a relatively new method for data imputation. We first eliminated all metabolite variables which were missing more than 10 percent of their values, since imputing these values would be unreliable. Then, with the remaining metabolites, we imputed 10 datasets, fit models, and pooled the results over each dataset.

To predict the stage of ischemic events in the subjects, we trained a series of random forest models. We used feature importance metrics to reduce our models and perform model selection. On each iteration of model selection, we conducted 13-fold custom-indexed cross validation (to account for the repeated measures) to measure the accuracy of each model and then eliminated all variables with zero feature importance. Eventually, we obtained a model using only five predictors (out of 736) which yielded 100% accuracy on each test fold of the data. We do recognize that our model may be susceptible to "regression to the mean" on a larger dataset. However, given the overall findings from the analysis, we do believe that our model may generalize well on larger datasets. 

