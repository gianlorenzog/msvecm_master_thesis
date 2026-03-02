# IRF computation for MS-VECM
#
# Previous version IRF_msvecm_v1 and IRF_MSwM didn't consider regime changes over time (MS)
# but fixed  S_t to a single regime throughout the IRF horizon, essentially treating it as a standard VECM.
#
# V.2 - Second Version (with lapply() + Reduce())
#
# Slightly more verbose but clearer logic, explicitly handles dimension mismatch with lapply() and Reduce().
# Slower due to lapply() and Reduce(), but safer for edge cases, as Reduce() handles accumulation explicitly.
#
#
# - consider Multi-Variable Shock to reflect how the system responds to the simultaneous shocks.
# - shock size should be inserted in the vector as follows: shock_variable <- c(1, 0, 0.5, -1)
# - single shock --> just fill one element (shock_variable <- c(1, 0, 0, 0))
#
#
calculate_irf_switching <- function(msvecm_model, regime, shock_vector, horizon = 10, num_simulations = 1000) {
  # Extract model parameters for the chosen regime
  Pi <- msvecm_model$vecm_params[[regime]]$Pi
  Gamma <- msvecm_model$vecm_params[[regime]]$Gamma
  intercept <- msvecm_model$vecm_params[[regime]]$intercept
  K <- nrow(Pi)
  lag_order <- dim(Gamma)[3]
  
  # Initialize matrices for IRFs
  irf <- matrix(0, nrow = horizon + 1, ncol = K)
  irf[1, ] <- shock_vector  # Apply initial shock
  
  # Simulate responses
  for (h in 2:(horizon + 1)) {
    # Ensure proper shape of previous IRF step
    irf_prev <- matrix(irf[h - 1, ], ncol = 1)  # Convert to column vector
    
    # Long-run effect
    long_run_effect <- Pi %*% irf_prev
    
    # Short-run effect (loop over lags)
    short_run_effect <- Reduce(`+`, lapply(1:lag_order, function(p) {
      if ((h - p) > 0) {
        Gamma[, , p] %*% matrix(irf[h - p, ], ncol = 1)
      } else {
        matrix(0, nrow = K, ncol = 1)  # Return a zero matrix for invalid indices
      }
    }))
    
    # Combine effects
    combined_effect <- long_run_effect + short_run_effect
    if (!is.null(intercept)) {
      combined_effect <- combined_effect + matrix(intercept, ncol = 1)  # Add intercept
    }
    
    # Store in IRF matrix (as a row vector)
    irf[h, ] <- t(combined_effect)
  }
  
  return(irf)
}

