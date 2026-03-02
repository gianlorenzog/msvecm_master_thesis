# Hansen (1992) Likelihood Ratio Test with formatted output
hansen_likelihood_ratio_test <- function(ll_unrestricted, ll_restricted, num_params) {
  # Calculate the test statistic
  lr_stat <- 2 * (ll_unrestricted - ll_restricted)
  
  # Handle extremely large LR statistic values gracefully
  if (lr_stat > 1e6) {
    p_value <- 0  # Practically zero for extreme values
  } else {
    # Compute the p-value
    df <- num_params
    p_value <- 1 - pchisq(lr_stat, df)
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
    LR_Statistic = lr_stat,
    Degrees_of_Freedom = num_params,
    P_Value = p_value,
    Formatted_P_Value = p_value_formatted,
    Significance = significance
  ))
}

# Custom function to format and display Hansen LR test results
print_test_results <- function(test_results) {
  cat("Hansen Likelihood Ratio Test Results:\n")
  cat(sprintf("  LR Statistic: %.3f\n", test_results$LR_Statistic))
  cat(sprintf("  Degrees of Freedom: %d\n", test_results$Degrees_of_Freedom))
  cat(sprintf("  P-Value: %s %s\n", test_results$Formatted_P_Value, test_results$Significance))
}
