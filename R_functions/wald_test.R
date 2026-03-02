# Granger Causality Test: Wald Test
#
# V.1 - update: Pi can be taken into account for long run causality
# NOTE: general purpose, it only requires Gamma and Pi (optional)
#
# Test whether the lagged coefficients of the causal variable(s) are jointly zero in each regime.
# Inputs:
#   Pi: A list of k x k matrices, one for each regime, representing the Pi matrix.
#   Gamma: A list of k x k x p arrays, one for each regime, representing the Gamma matrix.
#   regime: The regime number (1, 2, etc.).
#   variables: A vector of variable indices (e.g., c(1, 2) for the first and second variables).
#   lag: The lag number to test (1, 2, ..., p).
#   Possibility to include Pi and test also long run causality
#
wald_test <- function(Pi, Gamma, regime, variables, lag, include_Pi = FALSE) {
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
  gamma_hat <- Gamma[[regime]][variables, variables, lag]
  gamma_hat <- as.vector(gamma_hat)
  
  if (include_Pi) {
    # Extract Pi coefficients for the specified regime
    pi_hat <- Pi[[regime]][variables, variables]
    pi_hat <- as.vector(pi_hat)
    
    # Combine Pi and Gamma coefficients
    gamma_hat <- c(gamma_hat, pi_hat)
  }
  
  # Compute variance-covariance matrix (identity matrix for simplicity)
  var_gamma <- diag(rep(1, length(gamma_hat)))
  
  # Compute the Wald test statistic
  wald_stat <- t(gamma_hat) %*% solve(var_gamma) %*% gamma_hat
  wald_stat <- as.numeric(wald_stat)
  
  # Degrees of freedom
  df <- length(gamma_hat)
  
  # Critical value from chi-squared distribution (5% significance level)
  chi_critical <- qchisq(0.95, df)
  
  # Decision and result
  if (wald_stat > chi_critical) {
    result <- paste(
      "Reject the null hypothesis: Granger causality exists in regime", regime, 
      "for lag", lag, "and levels (if included)."
    )
  } else {
    result <- paste(
      "Fail to reject the null hypothesis: No Granger causality in regime", regime, 
      "for lag", lag, "and levels (if included)."
    )
  }
  
  # Return test statistic and decision
  return(list(
    regime = regime,
    variables = variables,
    lag = lag,
    include_Pi = include_Pi,
    wald_stat = wald_stat,
    chi_critical = chi_critical,
    decision = result
  ))
}


