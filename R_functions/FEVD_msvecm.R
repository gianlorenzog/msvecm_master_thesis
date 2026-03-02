#  FEVD with IRFs calculation integrated  - uses irf_ms_v.2
#
# Update 1: dynamically switches both regimes (via calculate_irf_switching) and 
# Sigma based on the current regime at each horizon step.
#
# Update 2: load irf_ms_v.2 directly inside the script
#
# Other features:
# - Normaled the FEVD across shocks ensures that the proportions sum to 1 for each variable and horizon
# - Modified FEVD computation to ensure non-negative contributions
# - integrated function for IRF matrix generation uses "gen_irf_matrix_v2" algorithm
#
#
source(file = "R_functions/irf_msvecm_v2.R")
compute_fevd_irf <- function(msvecm_model, regime, horizon, num_simulations, transition_matrix) {
  # Generate IRF Matrix
  k <- ncol(msvecm_model$vecm_params[[regime]]$Sigma)  # Number of variables
  irf_matrix <- array(0, dim = c(horizon + 1, k, k))  # Initialize IRF matrix
  
  for (shock in 1:k) {
    shock_vector <- rep(0, k)  # Create shock vector
    shock_vector[shock] <- 1   # Apply shock to variable `shock`
    
    # Compute IRFs for this shock, considering regime-switching
    irf_result <- calculate_irf_switching(msvecm_model, regime, shock_vector, horizon, transition_matrix)
    irf_matrix[, , shock] <- irf_result
  }
  
  # Compute FEVD
  fevd <- array(0, dim = c(horizon + 1, k, k))  # FEVD array
  
  for (var in 1:k) {  # Loop over dependent variables
    for (shock in 1:k) {  # Loop over shocks
      # Numerator: Contribution of shocks to `shock` for `var`
      numerator <- sapply(1:(horizon + 1), function(h) {
        sum(sapply(1:h, function(l) {
          irf_l <- irf_matrix[l, , ]  # IRF matrix at lag `l`
          
          # Dynamically determine regime and use its Sigma
          current_regime <- sample(1:length(msvecm_model$vecm_params), 1, prob = transition_matrix[regime, ])
          Sigma_dynamic <- msvecm_model$vecm_params[[current_regime]]$Sigma
          
          contribution <- as.numeric(irf_l[var, ] %*% Sigma_dynamic %*% matrix(irf_l[, shock], ncol = 1))
          max(contribution, 0)  # Ensure non-negative contributions
        }))
      })
      
      # Denominator: Total variance of forecast errors for variable `var`
      denominator <- sapply(1:(horizon + 1), function(h) {
        sum(sapply(1:h, function(l) {
          irf_l <- irf_matrix[l, , ]
          
          # Dynamically determine regime and use its Sigma
          current_regime <- sample(1:length(msvecm_model$vecm_params), 1, prob = transition_matrix[regime, ])
          Sigma_dynamic <- msvecm_model$vecm_params[[current_regime]]$Sigma
          
          as.numeric(irf_l[var, ] %*% Sigma_dynamic %*% matrix(irf_l[var, ], ncol = 1))
        }))
      })
      
      # Compute FEVD and handle division by zero
      fevd[, var, shock] <- ifelse(denominator > 0, numerator / denominator, 0)
    }
  }
  
  # Normalize FEVD across shocks for each variable and horizon
  for (h in 1:(horizon + 1)) {
    for (var in 1:k) {
      total <- sum(fevd[h, var, ])
      if (total > 0) {
        fevd[h, var, ] <- fevd[h, var, ] / total  # Normalize to ensure sum = 1
      }
    }
  }
  
  return(fevd)
}
