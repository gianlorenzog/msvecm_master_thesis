# Test Granger Causality: Wald test - MSVECM Output
#
# we test whether the lagged coeff of the causal variable (x) are jointly zero in each regime. 
# The following function is built to work with the MSVECM script (manual implementation)
#
# NOTE: revised function, integrating noise injection and handling for zero standard errors.
# Introduced small random noise during the bootstrap process to ensure variability: noise_level
#
# Recall that MSVECM has as default "const", so we replicated the same here
#
wald_test <- function(Pi, Gamma, regime, variables, lag, std_errors_gamma, std_errors_pi, include_Pi = FALSE, lambda = 1e-6) {
  # Validate regime
  if (regime > length(Pi)) {
    stop("Invalid regime number. Please choose a valid regime.")
  }
  
  # Validate lag
  num_lags <- dim(Gamma[[regime]])[3]
  if (lag > num_lags) {
    stop("Invalid lag number. Please choose a valid lag.")
  }
  
  # Validate variables
  num_vars <- nrow(Pi[[regime]])
  if (any(variables > num_vars) || any(variables < 1)) {
    stop("Invalid variable indices. Please choose valid variables.")
  }
  
  # Extract Gamma coefficients for the specified lag and regime
  gamma_hat <- c(
    Gamma[[regime]][variables[1], variables[1], lag],  # Row 2, Col 2
    Gamma[[regime]][variables[2], variables[2], lag]   # Row 1, Col 1
  )
  std_errors_gamma_subset <- c(
    std_errors_gamma[variables[1]],  # Row 2
    std_errors_gamma[variables[2]]   # Row 1
  )
  
  # Optionally include Pi coefficients
  if (include_Pi) {
    pi_hat <- c(
      Pi[[regime]][variables[1], variables[1]],  # Row 2, Col 2
      Pi[[regime]][variables[2], variables[2]]   # Row 1, Col 1
    )
    std_errors_pi_subset <- c(
      std_errors_pi[variables[1]],  # Row 2
      std_errors_pi[variables[2]]   # Row 1
    )
    
    # Combine Gamma and Pi coefficients and SEs
    gamma_hat <- c(gamma_hat, pi_hat)
    std_errors <- c(std_errors_gamma_subset, std_errors_pi_subset)
  } else {
    std_errors <- std_errors_gamma_subset
  }
  
  # Validate that lengths match
  if (length(gamma_hat) != length(std_errors)) {
    stop("Length of gamma_hat and std_errors does not match. Check indices or inputs.")
  }
  
  # Replace zero standard errors with a small positive value for numerical stability
  std_errors[std_errors == 0] <- 1e-6
  
  # Construct variance-covariance matrix
  var_gamma <- diag(std_errors^2)
  var_gamma_regularized <- var_gamma + lambda * diag(nrow(var_gamma))
  
  # Compute the Wald test statistic
  wald_stat <- t(gamma_hat) %*% solve(var_gamma_regularized) %*% gamma_hat
  wald_stat <- as.numeric(wald_stat)
  
  # Degrees of freedom
  df <- length(gamma_hat)
  
  # Critical value from chi-squared distribution (5% significance level)
  chi_critical <- qchisq(0.95, df)
  
  # Decision
  decision <- if (wald_stat > chi_critical) {
    paste(
      "Reject the null hypothesis: Granger causality exists in regime", regime, 
      "for lag", lag, "and levels (if included)."
    )
  } else {
    paste(
      "Fail to reject the null hypothesis: No Granger causality in regime", regime, 
      "for lag", lag, "and levels (if included)."
    )
  }
  
  # Return results
  list(
    regime = regime,
    variables = variables,
    lag = lag,
    include_Pi = include_Pi,
    wald_stat = wald_stat,
    chi_critical = chi_critical,
    decision = decision,
    std_errors = std_errors
  )
}



