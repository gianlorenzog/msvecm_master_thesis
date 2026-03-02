# CFEVD
# General purpose function: only requires an FEVD matrix as input, making it reusable across models
# It does not depend on how the FEVD is computed, as long as the input array is correctly formatted.
#
compute_cfevd <- function(fevd_result) {
  # Dimensions of the FEVD array
  horizons <- dim(fevd_result)[1] - 1  # Horizon count (h)
  variables <- dim(fevd_result)[2]    # Number of variables (k)
  shocks <- dim(fevd_result)[3]       # Number of shocks (k)
  
  # Initialize CFEVD array
  cfevd <- array(0, dim = c(horizons + 1, variables, shocks))
  
  # Compute cumulative FEVD for each variable-shock pair
  for (var in 1:variables) {
    for (shock in 1:shocks) {
      cfevd[, var, shock] <- cumsum(fevd_result[, var, shock])
    }
  }
  
  return(cfevd)
}
