# packages to load
library(rstatix)
library(lsmeans)
library(mice)
library(mitml)
library(DMwR)
library(lmer)
library(lmerTest)

#simple preprocessing
metabolites = colnames(pm_data)
colnames(pm_data)[1:736] = paste("metabolite_", 1:736, sep = "")

# SAMPLE test to inform analysis
dummy_data = data.frame(pm_data$metabolite_1, pm_data$time_period, pm_data$id)
colnames(dummy_data) = c("metabolite_1", "time_period", "id")

formula_anova = paste("metabolites_2", "~", "time_period + Error(id/time_period)")

dummy_test = anova_test(formula=as.formula('metabolite_1 ~ time_period + Error(id/time_period)'), 
                        data = dummy_data)

dummy = aov(formula = as.formula('metabolite_1 ~ time_period + Error(id/time_period)'), 
                   data = pm_data)

lsmeans(dummy, pairwise ~ treatment, adjust="tukey") 

get_anova_table(dummy_test)

# sample code
#pwc <- anxiety %>% group_by(time) %>% 
# pairwise_t_test(score ~ group, p.adjust.method = "tukey")

dummy_data %>% pairwise_t_test(metabolite_1 ~ time_period, paired = TRUE, 
                               p.adjust.method = "tukey")

# going to have to use combination of anova_test and aov function

### ANOVA TESTS FOR EACH METABOLITE

# create some null vectors and data frame
metabolite_name = rep(NA, 736)
metabolite_pathway = rep(NA, 736)
metabolite_subpathway = rep(NA, 736)
f_statistic = rep(NA, 736)
p_value = rep(NA, 736)
significance = rep(NA, 736)
effect_size = rep(NA, 736)
# pairwise_comparison = rep(NA, 736)
period_1_and_2_estimated_difference = rep(NA, 736)
period_1_and_2_p_value = rep(NA, 736)
period_1_and_2_significance = rep(NA, 736)
period_1_and_4_estimated_difference = rep(NA, 736)
period_1_and_4_p_value = rep(NA, 736)
period_1_and_4_significance = rep(NA, 736)
period_2_and_3_estimated_difference = rep(NA, 736)
period_2_and_3_p_value = rep(NA, 736)
period_2_and_3_significance = rep(NA, 736)
period_2_and_4_estimated_difference = rep(NA, 736)
period_2_and_4_p_value = rep(NA, 736)
period_2_and_4_significance = rep(NA, 736)

anova_tests = data.frame(metabolite_name, metabolite_pathway, metabolite_subpathway, 
                        f_statistic, p_value, significance, effect_size,
                        period_1_and_2_estimated_difference,
                        period_1_and_2_p_value,
                        period_1_and_2_significance,
                        period_1_and_4_estimated_difference,
                        period_1_and_4_p_value,
                        period_1_and_4_significance,
                        period_2_and_3_estimated_difference,
                        period_2_and_3_p_value,
                        period_2_and_3_significance,
                        period_2_and_4_estimated_difference,
                        period_2_and_4_p_value,
                        period_2_and_4_significance)

anova_tests_no_imp = anova_tests

# loop to conduct all tests
for (num in 1:736) {
  anova_tests_no_imp = anova_analyzer_no_imp(pm_data, metabolite_data, num, anova_tests_no_imp)
}

#save output to excel wb
rm_anova_excel = write.xlsx(anova_tests_no_imp, "anova_results.xlsx")

#save output to excel wb with all comparisons
rm_anova_excel_comparisons = write.xlsx(anova_tests_no_imp, "anova_results_all_comparisons.xlsx")

sig_metabolites = which(anova_tests_no_imp$significance == "yes")
anova_tests_no_imp_sig_metabolites = anova_tests_no_imp[sig_metabolites, ]

rm_anova_excel_significant_metabolies = 
  write.xlsx(anova_tests_no_imp_sig_metabolites, "anova_results.xlsx", sheetName = "Sheet2",
             append = TRUE)

##### WORK FOR IMPUTATION USING MICE ########

# checking percent of missing values in each metabolite
pMiss <- function(x){sum(is.na(x))/length(x)*100}
missing_more_than_10_percent = apply(pm_data,2,pMiss)[which(apply(pm_data,2,pMiss) > 10)]
missing_more_than_14_percent = apply(pm_data,2,pMiss)[which(apply(pm_data,2,pMiss) > 14)]
length(missing_more_than_10_percent) #number of metabolites missing more than 10 percent
length(missing_more_than_14_percent) #number of metabolites missing more than 14 percent

missing_data_pattern = md.pattern(pm_data)

missing_more_than_10_percent_names = names(missing_more_than_10_percent)
pm_data_90_percent_or_more = 
  pm_data[, -which(names(pm_data) %in% missing_more_than_10_percent_names)]

# perform multiple imputation on pm_data (metabolite values for each pig)
# start by identifying grouping variable
initial <- mice(pm_data_90_percent_or_more, maxit=0)
pred = initial$predictorMatrix
pred[,"id"]<- -2 # set ID as class variable for 2l.norm
pred["id",]<- -2 # may not be necessary

#Work on method
meth = initial$method
meth[which(meth == "pmm")] = "2l.norm"

# must convert id to integer beforehand
pm_data_90_percent_or_more$id = as.numeric(pm_data_90_percent_or_more$id)
pm_data_90_percent_or_more$time_period = as.numeric(pm_data_90_percent_or_more$time_period)

imputation.df<-mice(data = pm_data_90_percent_or_more,m = 5, method = meth, 
                    predictorMatrix = pred)

# convert id back to factor
pm_data_90_percent_or_more$id = as.factor(pm_data_90_percent_or_more$id)
pm_data_90_percent_or_more$time_period = as.factor(pm_data_90_percent_or_more$time_period)

anova_tests_no_imp = anova_tests

# missing metabolites
missing_metabolites = which(colSums(is.na(pm_data)) > 0)

## TRYING PCA IMPUTATION
imp_data_90 = imputePCA(pm_data_90_percent_or_more[, 1:592], ncp = 4, method = "Regularized", 
                        coeff.ridge=1, threshold = 1e-05, nb.init = 2,maxiter = 1000)
imp_data_90_table = imp_data_90$completeObs

class(imp_data_90_table)
imp_data_frame_90 = as.data.frame(imp_data_90_table)

# adding back time period variables
imp_data_frame_90$id = pm_data$id
imp_data_frame_90$time_period = pm_data$time_period

### INCLUDING TREATMENT GROUP

# to adjust SS calculations: options(contrasts=c("contr.sum","contr.poly"))

# SAMPLE test to inform analysis
dummy_data = data.frame(full_pm_data$metabolite_2, full_pm_data$time_period, full_pm_data$id,
                        full_pm_data$treatment)
colnames(dummy_data) = c("metabolite_2", "time_period", "id", "treatment")

formula_anova = paste("metabolite_2", "~", "time_period*treatment + (1|id:treatment)")
formula_anova_test = paste("metabolite_2", "~", "time_period*treatment + Error(id/time_period)")

dummy_test = anova_test(formula = as.formula(formula_anova_test), 
                        data = dummy_data)

dummy = lmer(formula = as.formula(formula_anova_test), 
             data = full_pm_data)

dummy_test = anova_test(data = full_pm_data, dv = metabolite_2, 
                        wid = id,between = treatment, within = time_period)

lsmeans(dummy, pairwise ~ time_period, adjust="tukey") 

get_anova_table(dummy_test)

# edit old stuff
### ANOVA TESTS FOR EACH METABOLITE

anova_tests_treatment = data.frame(metabolite_name, metabolite_pathway, metabolite_subpathway)

# create some null vectors and data frame
metabolite_name = rep(NA, 736)
metabolite_pathway = rep(NA, 736)
metabolite_subpathway = rep(NA, 736)

# adding new columns
anova_tests_treatment$interaction_f_statistic = rep(NA, 736)
anova_tests_treatment$interaction_p_value = rep(NA, 736)
anova_tests_treatment$interaction_significance = rep(NA, 736)
anova_tests_treatment$interaction_effect_size = rep(NA, 736)
anova_tests_treatment$treatment_f_statistic = rep(NA, 736)
anova_tests_treatment$treatment_p_value = rep(NA, 736)
anova_tests_treatment$treatment_significance = rep(NA, 736)
anova_tests_treatment$treatment_effect_size = rep(NA, 736)
anova_tests_treatment$time_period_f_statistic = rep(NA, 736)
anova_tests_treatment$time_period_p_value = rep(NA, 736)
anova_tests_treatment$time_period_significance = rep(NA, 736)
anova_tests_treatment$time_period_effect_size = rep(NA, 736)

# pairwise_comparisons time period
anova_tests_treatment$period_1_and_2_estimated_difference = rep(NA, 736)
anova_tests_treatment$period_1_and_2_p_value = rep(NA, 736)
anova_tests_treatment$period_1_and_2_significance = rep(NA, 736)
anova_tests_treatment$period_1_and_4_estimated_difference = rep(NA, 736)
anova_tests_treatment$period_1_and_4_p_value = rep(NA, 736)
anova_tests_treatment$period_1_and_4_significance = rep(NA, 736)
anova_tests_treatment$period_2_and_3_estimated_difference = rep(NA, 736)
anova_tests_treatment$period_2_and_3_p_value = rep(NA, 736)
anova_tests_treatment$period_2_and_3_significance = rep(NA, 736)
anova_tests_treatment$period_2_and_4_estimated_difference = rep(NA, 736)
anova_tests_treatment$period_2_and_4_p_value = rep(NA, 736)
anova_tests_treatment$period_2_and_4_significance = rep(NA, 736)

anova_tests_treatment$control_treatment_estimated_difference = rep(NA, 736)
anova_tests_treatment$control_treatment_p_value = rep(NA, 736)
anova_tests_treatment$control_treatment_significance = rep(NA, 736)

# loop to conduct all tests
for (num in 1:736) {
  anova_tests_treatment = anova_analyzer_treatment_group(full_pm_data, metabolite_data, num, 
                                                          anova_tests_treatment)
}

#save output to excel wb
rm_anova_excel = write.xlsx(anova_tests_treatment, "anova_results_treatment.xlsx")




