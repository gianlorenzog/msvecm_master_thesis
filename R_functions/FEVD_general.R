# FEVD using MS IRFs - general purpose
#
# Update 1: dynamically determines the current regime and switch Sigma
# based on the current regime at each horizon step.
#
# NOTE: the script doesn't compute irf and gen the irf matrix internally like
# the other tailored for the manual MSVECM model 
#
# Other features:
# - Normaled the FEVD across shocks ensures that the proportions sum to 1 for each variable and horizon
# - Modified FEVD computation to ensure non-negative contributions
#
#
compute_fevd_matrix <- function(irf_matrix, Sigma_list, horizon, transition_matrix, initial_regime) {
  k <- ncol(Sigma_list[[1]])  # Number of variables (assumes all Sigma matrices have the same dimensions)
  fevd <- array(0, dim = c(horizon + 1, k, k))  # FEVD array
  
  for (var in 1:k) {  # Loop over dependent variables
    for (shock in 1:k) {  # Loop over shocks
      # Numerator: Contribution of shocks to `shock` for `var`
      numerator <- sapply(1:(horizon + 1), function(h) {
        # Initialize regime tracking for this horizon step
        current_regime <- initial_regime
        
        sum(sapply(1:h, function(l) {
          irf_l <- irf_matrix[l, , ]  # IRF matrix at lag `l`
          
          # Update regime after the first step based on transition probabilities
          if (l > 1) {
            current_regime <- sample(1:length(Sigma_list), 1, prob = transition_matrix[current_regime, ])
          }
          
          # Use the Sigma matrix for the current regime
          Sigma_dynamic <- Sigma_list[[current_regime]]
          
          # Compute contribution of shocks
          contribution <- as.numeric(irf_l[var, ] %*% Sigma_dynamic %*% matrix(irf_l[, shock], ncol = 1))
          max(contribution, 0)  # Ensure non-negative contributions
        }))
      })
      
      # Denominator: Total variance of forecast errors for variable `var`
      denominator <- sapply(1:(horizon + 1), function(h) {
        # Initialize regime tracking for this horizon step
        current_regime <- initial_regime
        
        sum(sapply(1:h, function(l) {
          irf_l <- irf_matrix[l, , ]
          
          # Update regime after the first step based on transition probabilities
          if (l > 1) {
            current_regime <- sample(1:length(Sigma_list), 1, prob = transition_matrix[current_regime, ])
          }
          
          # Use the Sigma matrix for the current regime
          Sigma_dynamic <- Sigma_list[[current_regime]]
          
          # Compute total variance
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


