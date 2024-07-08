# neg binomial test of wb ~ features
glmnb_total_worm_on_feature <- function(col_numb, df_of_asvs, df_colnames, total_worm_covariate){
  set.seed(8675309)
  print(col_numb)
  summary_of_glm <- try({summary(glm.nb(total_worm_covariate~df_of_asvs[,col_numb]))}, silent=T)
  return(c(summary_of_glm$coefficients[2,1], summary_of_glm$coefficients[2,4],df_colnames[col_numb]))
}

# Function to perform the Adonis 2 test and extract results, explicitly for NAE testing with ABX interaction
extract_adonis_results <- function(data, column_name, bc_dist) {
  # Construct the formula for the Adonis 2 test
  formula <- as.formula(paste0("bc_dist~", column_name, "*drug"))
  
  # Perform the Adonis 2 test
  adonis_result <- adonis2(formula, data = data, permutations = 10000)
  
  # Create a data frame with the column name and p-values
  return(tibble(column = column_name,
                p_value = adonis_result$`Pr(>F)`[1],  # Main effect p-value
                interaction_p_value = adonis_result$`Pr(>F)`[3])) # Interaction p-value
}

# Function to perform correlation test between a metabolite and ASVs, mediation shares a few similarities with approach implemented here
# (https://github.com/fjw536/AnorexiaGutMicrobiome/blob/main/Mediation_codes.r)
parallel_calculate_med_and_cors <- function(metab, metab_name, ASV_set, dependent_var) {
  cor_values <- sapply(names(ASV_set), function(asv_name) {
    # Exclude the current ASV being tested from control variables
    control_vars <- ASV_set[, -which(names(ASV_set) == asv_name), drop = FALSE]
    
    # Check if the control_vars matrix is empty
    if (ncol(control_vars) == 0) {
      warning("Control variables matrix is empty.")
      return(NULL)
    }
    
    # Convert control_vars to numeric matrix
    control_vars <- as.matrix(as.data.frame(control_vars))
    
    # Check if control_vars is numeric
    if (!is.numeric(control_vars)) {
      stop("Control variables matrix must be numeric.")
    }
    
    # Calculate simple correlation
    cor_test <- cor.test(metab, ASV_set[[asv_name]])
    
    # Calculate partial correlation
    res_pcor <- nptest::np.cor.test(x = metab, y = ASV_set[[asv_name]], 
                                    z = control_vars, partial = TRUE, R = 1000)
    # Initiate mediation
    initial_df <- data.frame(cbind(metab, ASV_set[[asv_name]], dependent_var)) 
    names(initial_df) <- c("X", "M", "Y") 
    combined_df <- initial_df %>%
      mutate(XX = ntile(X, 4))
    
    mediator_on_X <- lm(M ~ XX, data = combined_df) 
    
    whole_model <- glm.nb(Y ~ XX + M, control=glm.control(maxit = 5000), 
                          data = combined_df)
    set.seed(1001)
    med_results <- try({(mediate(mediator_on_X, whole_model, treat = "XX", mediator = "M",
                                 control.value=1, treat.value=4, boot = TRUE, sims = 10))})
    
    if(class(med_results)!="try-error"){
      results_summary <- summary(med_results)
      
      # Create list with results from both tests
      list(metab_name = metab_name, asv_name = asv_name, 
           simple_estimate = cor_test$estimate, simple_p_value = cor_test$p.value,
           partial_estimate = res_pcor$estimate, partial_p_value = res_pcor$p.value,
           ACME_est = results_summary$d.avg,
           ACME_pval = results_summary$d.avg.p,
           ADE_est = results_summary$z.avg,
           ADE_pval = results_summary$z.avg.p,
           prop_mediated_est = results_summary$n.avg, 
           prop_mediated_pval = results_summary$n.avg.p)
    }
    else{
      list(metab_name = metab_name, asv_name = asv_name, 
           simple_estimate = cor_test$estimate, simple_p_value = cor_test$p.value,
           partial_estimate = res_pcor$estimate, partial_p_value = res_pcor$p.value,
           ACME_est = NA,
           ACME_pval = NA,
           ADE_est = NA,
           ADE_pval = NA,
           prop_mediated_est = NA, 
           prop_mediated_pval = NA)
    }
  })
  return(cor_values)
}