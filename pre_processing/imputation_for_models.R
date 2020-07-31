library(devtools)
library(reshape2)
library(tidyr)
library(missMDA)
library(VIM)
library(naniar)

#simple preprocessing
metabolites = colnames(full_pm_data)[1:736]
colnames(full_pm_data)[1:736] = paste("metabolite_", 1:736, sep = "")

# checking percent of missing values in each metabolite
pMiss <- function(x){sum(is.na(x))/length(x)*100}
full_missing_more_than_10_pct = apply(full_pm_data,2,pMiss)[which(apply(full_pm_data,2,pMiss) > 10)]
full_missing_more_than_15_pct = apply(full_pm_data,2,pMiss)[which(apply(full_pm_data,2,pMiss) > 15)]

#number of metabolites missing more than 10 percent
full_len_10 = length(full_missing_more_than_10_pct)
#number of metabolites missing more than 15 percent
full_len_15 = length(full_missing_more_than_15_pct)

# getting rid of metabolites with more than 10 percent missing observations
full_missing_more_than_10_pct_names = names(full_missing_more_than_10_pct)
full_pm_data_90_percent_or_more = 
  full_pm_data[, -which(names(full_pm_data) %in% full_missing_more_than_10_pct_names)]

# getting rid of metabolites with more than 15 percent missing observations
full_missing_more_than_15_pct_names = names(full_missing_more_than_15_pct)
full_pm_data_85_percent_or_more = 
  full_pm_data[, -which(names(full_pm_data) %in% full_missing_more_than_15_pct_names)]

# number of metabolites in data frames
nm_85_pct = ncol(full_pm_data_85_percent_or_more) - 3
num_90_pct = ncol(full_pm_data_90_percent_or_more) - 3

# converting from long to wide format for imputation
wide_full_85 = reshape(full_pm_data_85_percent_or_more, 
                       v.names = colnames(full_pm_data_85_percent_or_more)[1:nm_85_pct], 
                       timevar = "time_period",
                       idvar = "id", direction = "wide",
                      )

wide_full_90 = reshape(full_pm_data_90_percent_or_more, 
                       v.names = colnames(full_pm_data_90_percent_or_more)[1:num_90_pct], 
                       timevar = "time_period",
                       idvar = "id", direction = "wide",
                      )

# visualize missing data
gg_miss_var(wide_full_85)
matrixplot(wide_full_85, sortby = 2)

gg_miss_var(wide_full_90)
matrixplot(wide_full_90, sortby = 2)

matrixplot(full_pm_data, sortby = 2)

# MFA imputation instead of PCA imputation to generalize to different classes

# 85 percent data
wide_full_85_imp = imputeMFA(wide_full_85[c(3:ncol(wide_full_85))], ncp = 4,
                             group = rep(nm_85_pct, 4), method = "Regularized")

# convert to data frame
wide_full_85_imp_data = as.data.frame(wide_full_85_imp$completeObs)
wide_full_85_imp_data = cbind(wide_full_85[c(1:2)], wide_full_85_imp_data)

# 90 percent data
wide_full_90_imp = imputeMFA(wide_full_90[c(3:ncol(wide_full_90))], ncp = 4,
                             group = rep(num_90_pct, 4), method = "Regularized")

# convert to data frame
wide_full_90_imp_data = as.data.frame(wide_full_90_imp$completeObs)
wide_full_90_imp_data = cbind(wide_full_90[c(1:2)], wide_full_90_imp_data)

### MULTIPLE IMPUTATION USING FAMD
wide_full_85_MIFAMD = MIFAMD(wide_full_85, ncp = 4, method = "Regularized", nboot = 10)
wide_full_90_MIFAMD = MIFAMD(wide_full_90, ncp = 4, method = "Regularized", nboot = 10)

wide_full_85_MIFAMD_data = wide_full_85_MIFAMD$res.MI
wide_full_90_MIFAMD_data = wide_full_90_MIFAMD$res.MI

# creating varying vector to convert to long format


# CONVERTING DATA FROM WIDE TO LONG FORMAT FOR MODELING
long_full_85_imp_data = reshape(data=wide_full_85_imp_data, idvar="id",
                                varying = colnames(wide_full_85_imp_data)[3:2494],
                                new.row.names = 1:56,
                                direction="long")

long_full_90_imp_data = reshape(data=wide_full_90_imp_data, idvar="id",
                                varying = colnames(wide_full_90_imp_data)[3:2386],
                                new.row.names = 1:56,
                                direction="long")

# doing it now for MI data
# define reshape function

reshape_wide_to_long = function(df) {

    new_data = reshape(data=df, idvar="id", varying = colnames(df)[3:ncol(df)],
                       new.row.names = 1:56, direction="long") 
    return(new_data)
}

long_full_85_MIFAMD_data = lapply(wide_full_85_MIFAMD_data, FUN = reshape_wide_to_long)
long_full_90_MIFAMD_data = lapply(wide_full_90_MIFAMD_data, FUN = reshape_wide_to_long)






