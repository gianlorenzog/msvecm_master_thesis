# Generate IRF Matrix (MS-VECM)
#
# This script is based on irf_ms_v.2
#
source(file = "R_functions/irf_msvecm_v2.R")
generate_irf_matrix <- function(msvecm_model, regime, horizon, num_simulations) {
  k <- ncol(msvecm_model$vecm_params[[regime]]$Sigma)  # Number of variables
  irf_matrix <- array(0, dim = c(horizon + 1, k, k))  # Initialize IRF matrix
  
  for (shock in 1:k) {
    shock_vector <- rep(0, k)  # Create shock vector
    shock_vector[shock] <- 1   # Apply shock to variable `shock`
    
    # Compute IRFs for this shock
    irf_result <- calculate_irf_switching(msvecm_model, regime, shock_vector, horizon, num_simulations)
    irf_matrix[, , shock] <- irf_result
  }
  
  return(irf_matrix)
}