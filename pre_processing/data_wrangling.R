# parameters for jvm
options(java.parameters = "-Xmx8000m")
options(java.parameters = "-Xmx8g" )

# packages to load
library(xlsx)
library(plyr)
library(dplyr)
library(ggplot2)
library(nnet)

### CONTROL ###
# reading in pig data
pig_metabolite_data = read.xlsx("pig_data.xlsx", sheetName = "OrigScale", startRow = 12, 
                     colIndex = c(2, 43:70), header = FALSE)

# transpose data
pig_metabolite_data = t(pig_metabolite_data)
pig_metabolite_data = as.data.frame(pig_metabolite_data)

colnames(pig_metabolite_data) = pig_metabolite_data[1,]

# remove first row with metabolite names
pig_metabolite_data = pig_metabolite_data[-1,]

# convert columns to numeric as opposed to character
pm_data <- mutate_all(pig_metabolite_data, function(x) as.numeric(x))

# adding time period for classification
pm_data$time_period = as.factor(rep(c(1,2,3,4), 7))

pig_ids = c(rep(6204, 4), rep(6217, 4), rep(6192, 4), rep(6254, 4), rep(6256, 4), rep(6193, 4), 
            rep(6660, 4))

pm_data$id = pig_ids
pm_data$id = as.factor(pm_data$id)

### METABOLITE INFO ###
# reading in data about each metabolite
metabolite_data = read.xlsx("pig_data.xlsx", sheetName = "OrigScale", startRow = 11, 
                                colIndex = c(1:9), header = TRUE)

### TREATMENT GROUP ###
# reading in pig data
pig_metabolite_data_treatment = read.xlsx("pig_data.xlsx", sheetName = "OrigScale", startRow = 12, 
                                colIndex = c(2, 15:42), header = FALSE)

# transpose data
pig_metabolite_data_treatment = t(pig_metabolite_data_treatment)
pig_metabolite_data_treatment = as.data.frame(pig_metabolite_data_treatment)

colnames(pig_metabolite_data_treatment) = pig_metabolite_data_treatment[1,]

# remove first row with metabolite names
pig_metabolite_data_treatment = pig_metabolite_data_treatment[-1,]

# convert columns to numeric as opposed to character
pm_data_treatment <- mutate_all(pig_metabolite_data_treatment, function(x) as.numeric(x))

# adding time period for classification
pm_data_treatment$time_period = as.factor(rep(c(1,2,3,4), 7))

pig_ids_treatment = c(rep(6260, 4), rep(6187, 4), rep(6190, 4), rep(6253, 4), rep(6209, 4), 
                      rep(6652, 4), rep(6699, 4))

pm_data_treatment$id = pig_ids_treatment
pm_data_treatment$id = as.factor(pm_data_treatment$id)

### FULL DATA SET ###
full_pm_data <- rbind(pm_data, pm_data_treatment)
full_pm_data$treatment = c(rep(0, 28), rep(1, 28))
full_pm_data$treatment = as.factor(full_pm_data$treatment)



