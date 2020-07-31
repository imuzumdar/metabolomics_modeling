anova_analyzer_no_imp <- function(pig_data, mb_data, number, anova_data){
  
  # forming dummy data set
  col_name = paste("metabolite", number, sep = "_")
  test_data = data.frame(pig_data[[col_name]], pig_data[["time_period"]], 
                          pig_data[["id"]])
  colnames(test_data) = c(col_name, "time_period", "id")
  
  #formula string
  formula_string = paste(col_name,"~", "time_period + Error(id/time_period)")
  
  # updating data table metabolite info
  anova_data$metabolite_name[number] = mb_data$BIOCHEMICAL[number]
  anova_data$metabolite_pathway[number] = mb_data$SUPER_PATHWAY[number]
  anova_data$metabolite_subpathway[number] = mb_data$SUB_PATHWAY[number]
  
  # checking if anova is even possible
  output = tryCatch(
    {
      # conducting repeated measures anova
      anova = anova_test(formula = as.formula(formula_string), data = test_data)
      anova_table = get_anova_table(anova)
      anova_for_pairwise_comp = aov(formula = as.formula(formula_string), data = test_data)
      
      # pairwise comparisons
      # pairwise_formula_string = paste("pairwise", "~", "time_period")
      pairwise_table = summary(lsmeans(anova_for_pairwise_comp, pairwise ~ time_period, 
                                       adjust="tukey")$contrasts)
      
      # determining significance of metabolite
      sig = NA
      if (anova_table$p <= .05) {
        sig = "yes"
      } else {
        sig = "no"
      }
      
      # determining sig of difference in metabolite measure from period 1 to 2
      sig_12 = NA
      if (pairwise_table$p.value[1] <= .05) {
        sig_12 = "yes"
      } else {
        sig_12 = "no"
      }
      
      # determining sig of difference in metabolite measure from period 1 to 4
      sig_14 = NA
      if (pairwise_table$p.value[3] <= .05) {
        sig_14 = "yes"
      } else {
        sig_14 = "no"
      }
      
      # determining sig of difference in metabolite measure from period 2 to 3
      sig_23 = NA
      if (pairwise_table$p.value[4] <= .05) {
        sig_23 = "yes"
      } else {
        sig_23 = "no"
      }
      
      # determining sig of difference in metabolite measure from period 2 to 4
      sig_24 = NA
      if (pairwise_table$p.value[5] <= .05) {
        sig_24 = "yes"
      } else {
        sig_24 = "no"
      }
      
      # updating statistics from testing
      anova_data$f_statistic[number] = anova_table$F
      anova_data$p_value[number] = anova_table$p
      anova_data$significance[number] = sig
      anova_data$effect_size[number] = anova_table$ges
      anova_data$period_1_and_2_estimated_difference[number] = pairwise_table$estimate[1]
      anova_data$period_1_and_2_p_value[number] = pairwise_table$p.value[1]
      anova_data$period_1_and_2_significance[number] = sig_12
      anova_data$period_1_and_4_estimated_difference[number] = pairwise_table$estimate[3]
      anova_data$period_1_and_4_p_value[number] = pairwise_table$p.value[3]
      anova_data$period_1_and_4_significance[number] = sig_14
      anova_data$period_2_and_3_estimated_difference[number] = pairwise_table$estimate[4]
      anova_data$period_2_and_3_p_value[number] = pairwise_table$p.value[4]
      anova_data$period_2_and_3_significance[number] = sig_23
      anova_data$period_2_and_4_estimated_difference[number] = pairwise_table$estimate[5]
      anova_data$period_2_and_4_p_value[number] = pairwise_table$p.value[5]
      anova_data$period_2_and_4_significance[number] = sig_24
      
      return(anova_data)
    },
    error=function(cond) {
      message(paste("metabolite does not have enough variables:", col_name))
      # message(cond)
      # Choose a return value in case of error
      return(NA)
    }
  )
  
  # checking if error went through or not
  if (is.na(output)) {
    return(anova_data)
  }
  else {
    return(output)
  }
}
