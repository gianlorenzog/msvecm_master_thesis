#--- General Likelihood Ratio (LR) Test
# H0: Restricted model is sufficient.
lr_test <- function(ll_unrestricted, ll_restricted, num_params_unrestricted, num_params_restricted) {
  # Calculate LR statistic
  LR <- 2 * (ll_unrestricted - ll_restricted)
  
  # Degrees of freedom
  df <- num_params_unrestricted - num_params_restricted
  
  # Handle extremely large LR statistic values gracefully
  if (LR > 1e6) {
    p_value <- 0  # Practically zero for extreme values
  } else {
    # Compute the p-value
    p_value <- 1 - pchisq(LR, df)
  }
  
  # Determine significance stars with more robust levels
  significance <- ""
  if (p_value <= 0.001) {
    significance <- "****"
  } else if (p_value <= 0.01) {
    significance <- "***"
  } else if (p_value <= 0.05) {
    significance <- "**"
  } else if (p_value <= 0.1) {
    significance <- "*"
  }
  
  # Format the p-value
  p_value_formatted <- ifelse(p_value < 2.2e-16, "< 2.2e-16", format(p_value, digits = 4))
  
  # Return results with a clean summary
  return(list(
    LR_statistic = LR,
    degrees_of_freedom = df,
    P_Value = p_value,
    Formatted_P_Value = p_value_formatted,
    Significance = significance
  ))
}

# Custom function to format and display LR test results
print_lr_test <- function(test_results) {
  cat("Likelihood Ratio (LR) Test Results:\n")
  cat(sprintf("  LR Statistic: %.3f\n", test_results$LR_statistic))
  cat(sprintf("  Degrees of Freedom: %d\n", test_results$degrees_of_freedom))
  cat(sprintf("  P-Value: %s %s\n", test_results$Formatted_P_Value, test_results$Significance))
}
