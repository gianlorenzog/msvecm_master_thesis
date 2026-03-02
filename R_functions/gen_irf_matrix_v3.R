# Generate IRF Matrix (MS-VECM)
#
# This script is based on irf_ms_v.3
#
source(file = "R_functions/irf_msvecm_v3.R")
gen_irf_matrix<- function(msvecm_model, transition_matrix, horizon, initial_regime = 1) {
  # Number of variables
  k <- ncol(msvecm_model$vecm_params[[1]]$Pi)  # Assume same dimensions for all regimes
  
  # Initialize the IRF matrix (horizon + 1 x variables x variables)
  irf_matrix <- array(0, dim = c(horizon + 1, k, k))
  
  for (shock in 1:k) {
    # Create shock vector for variable `shock`
    shock_vector <- rep(0, k)
    shock_vector[shock] <- 1  # Apply unit shock to the `shock` variable
    
    # Compute the regime-switching IRF for this shock
    irf_result <- calculate_irf_switching(
      msvecm_model = msvecm_model,
      shock_vector = shock_vector,
      transition_matrix = transition_matrix,
      horizon = horizon,
      initial_regime = initial_regime
    )
    
    # Store the result in the IRF matrix
    irf_matrix[, , shock] <- irf_result
  }
  
  return(irf_matrix)
}