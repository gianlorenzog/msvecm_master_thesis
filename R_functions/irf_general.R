calculate_irf_switching <- function(Pi_list, Gamma_list, intercept_list, shock_vector,
                                           horizon, transition_matrix, initial_regime) {
  # Number of variables
  K <- length(shock_vector)
  
  # Extract regime-specific parameters
  #Pi_list <- Pi_list
  #Gamma_list <- Gamma_list
  #intercept_list <- intercept_list
  
  # Initialize IRF matrix
  irf <- matrix(0, nrow = horizon + 1, ncol = K)
  irf[1, ] <- shock_vector  # Apply initial shock
  
  # Start with the initial regime
  current_regime <- initial_regime
  
  for (h in 2:(horizon + 1)) {
    irf_prev <- matrix(irf[h - 1, ], ncol = 1)  # Convert to column vector
    
    # Get regime-specific parameters for the current regime
    Pi <- Pi_list[[current_regime]]
    Gamma <- Gamma_list[[current_regime]]
    intercept <- intercept_list[[current_regime]]
    lag_order <- dim(Gamma)[3]  # Number of lags
    
    # Long-run effect
    long_run_effect <- Pi %*% irf_prev
    
    # Short-run effects over lag_order
    short_run_effect <- Reduce(`+`, lapply(1:lag_order, function(p) {
      if ((h - p) > 0) {
        Gamma[, , p] %*% matrix(irf[h - p, ], ncol = 1)
      } else {
        matrix(0, nrow = K, ncol = 1)  # Return zero matrix for invalid indices
      }
    }))
    
    # Combine effects
    combined_effect <- long_run_effect + short_run_effect
    if (!is.null(intercept)) {
      combined_effect <- combined_effect + matrix(intercept, ncol = 1)  # Add intercept
    }
    
    # Store the IRF for this step
    irf[h, ] <- t(combined_effect)
    
    # Transition to the next regime based on transition probabilities
    current_regime <- sample(1:length(Pi_list), 1, prob = transition_matrix[current_regime, ])
  }
  
  return(irf)
}