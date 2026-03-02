# Function to calculate AIC and BIC
calculate_aic_bic <- function(log_likelihood, num_params, num_obs) {
  aic <- -2 * log_likelihood + 2 * num_params
  bic <- -2 * log_likelihood + log(num_obs) * num_params
  return(list(AIC = aic, BIC = bic))
}