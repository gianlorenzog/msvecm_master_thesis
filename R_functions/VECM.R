library(MASS)  # For generalized inverse
library(mvtnorm)  # For multivariate normal density
# v.4
# ADDED: possibility to include and estimate the intercept terms
# 
#
fit_vecm <- function(data, lag_order = 1, rank = 1, include = c("const", "none"), max_iter = 100, tol = 1e-6) {
  # Match include argument
  include <- match.arg(include, c("const", "none"))
  
  # Standardize data
  data <- scale(data)
  T <- nrow(data)
  K <- ncol(data)
  diff_data <- diff(data)
  
  # Create lagged terms
  lagged_data <- embed(data, lag_order + 1)
  lagged_data <- lagged_data[, -(1:K)]  # Remove the current observation
  lagged_data <- matrix(lagged_data, ncol = K * lag_order, byrow = FALSE)  # Ensure proper shape
  diff_data <- diff_data[lag_order:(T - 1), , drop = FALSE]  # Align differences with lagged terms
  
  # Initialize VECM parameters
  set.seed(123)
  Pi <- matrix(runif(K^2, -0.5, 0.5), nrow = K, ncol = K)  # Cointegration matrix
  Gamma <- array(runif(K^2 * lag_order, -0.5, 0.5), dim = c(K, K, lag_order))  # Short-term matrices
  Sigma <- diag(K)  # Covariance matrix
  intercept <- if (include == "const") rep(0, K) else NULL  # Initialize intercept if included
  
  # Constrain the rank of Pi
  constrain_rank <- function(Pi, rank) {
    svd_decomp <- svd(Pi)
    u <- svd_decomp$u
    v <- svd_decomp$v
    d <- svd_decomp$d
    d[(rank + 1):length(d)] <- 0  # Keep only the largest 'rank' singular values
    # Reconstruct Pi with reduced rank
    return(u[, 1:rank, drop = FALSE] %*% diag(d[1:rank], nrow = rank, ncol = rank) %*% t(v[, 1:rank, drop = FALSE]))
  }
  Pi <- constrain_rank(Pi, rank)
  
  # Decompose Pi into Alpha and Beta
  decompose_pi <- function(Pi, rank) {
    svd_decomp <- svd(Pi)
    alpha <- svd_decomp$u[, 1:rank, drop = FALSE]  # Loading matrix (adjustment coefficients)
    beta <- svd_decomp$v[, 1:rank, drop = FALSE]  # Cointegration matrix
    return(list(alpha = alpha, beta = beta))
  }
  decomposition <- decompose_pi(Pi, rank)
  alpha <- decomposition$alpha
  beta <- decomposition$beta
  
  log_likelihood <- NULL
  fitted_values <- matrix(0, nrow = T - lag_order, ncol = K)
  residuals <- matrix(0, nrow = T - lag_order, ncol = K)
  prev_log_likelihoods <- c()  # Store last few log-likelihoods to detect oscillation
  
  for (iter in 1:max_iter) {
    # Compute residuals and fitted values
    for (t in 1:(T - lag_order)) {
      # Cointegration terms (long-run dynamics): Use only the first K columns of lagged_data
      lagged_terms <- Pi %*% matrix(lagged_data[t, 1:K], ncol = 1)
      
      # Short-run dynamics (Gamma matrices)
      short_run_terms <- rowSums(sapply(1:lag_order, function(p) {
        if (t - p + 1 > 0) {  # Ensure indices are valid
          Gamma[, , p] %*% matrix(diff_data[t - p + 1, ], ncol = 1)
        } else {
          matrix(0, nrow = K, ncol = 1)  # Return zero matrix for invalid indices
        }
      }))
      
      # Intercept term
      intercept_term <- if (!is.null(intercept)) intercept else rep(0, K)
      
      # Fitted values and residuals
      fitted_values[t, ] <- lagged_terms + short_run_terms + intercept_term
      residuals[t, ] <- diff_data[t, ] - fitted_values[t, ]
    }
    
    # Compute likelihood
    densities <- apply(residuals, 1, function(res) {
      dmvnorm(res, mean = rep(0, K), sigma = Sigma)
    })
    new_log_likelihood <- sum(log(pmax(densities, 1e-12)))
    prev_log_likelihoods <- c(prev_log_likelihoods, new_log_likelihood)
    
    # Convergence check
    if (length(prev_log_likelihoods) > 5) {
      diff_log <- abs(diff(prev_log_likelihoods[(length(prev_log_likelihoods) - 4):length(prev_log_likelihoods)]))
      if (all(diff_log < tol)) {
        cat("Convergence reached at iteration:", iter, "\n")
        break
      }
    }
    log_likelihood <- new_log_likelihood
    cat("Iteration:", iter, "Log-Likelihood:", new_log_likelihood, "\n")
    
    # Update intercept only if include = "const"
    if (include == "const") {
      intercept <- colMeans(diff_data - (fitted_values - intercept))
    }
    
    # Constrain Pi to the specified rank
    Pi <- constrain_rank(Pi, rank)
  }
  
  # Calculate the number of estimated parameters
  num_params_Pi <- rank * (K + rank)
  num_params_Gamma <- K^2 * lag_order
  num_params_Sigma <- K * (K + 1) / 2  # Covariance matrix parameters
  num_params_intercept <- if (include == "const") K else 0
  total_num_params <- num_params_Pi + num_params_Gamma + num_params_Sigma + num_params_intercept
  
  list(
    Pi = Pi,
    Gamma = Gamma,
    Sigma = Sigma,
    intercept = intercept,
    alpha = alpha,
    beta = beta,
    log_likelihood = log_likelihood,
    fitted_values = fitted_values,
    residuals = residuals,
    num_params = total_num_params  # Store total number of parameters
  )
}

