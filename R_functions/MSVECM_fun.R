library(MASS)
library(mvtnorm)
# v.9
# Code works
# Added some Safeguards in Backward Pass, Debug checks, Normalized next_probs to Avoid Bias
# TO DO: Fix smoothed probabilities
#
markov_switching_vecm <- function(data, num_states = 2, lag_order = 1, rank = 1, include = c("const", "none"), max_iter = 100, tol = 1e-6, debug = FALSE) {
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
  
  # Define the decompose_pi function to get Alpha and Beta from Pi
  decompose_pi <- function(Pi, rank) {
    svd_decomp <- svd(Pi)
    alpha <- svd_decomp$u[, 1:rank, drop = FALSE]  # Loading matrix (adjustment coefficients)
    beta <- svd_decomp$v[, 1:rank, drop = FALSE]  # Cointegration matrix
    return(list(alpha = alpha, beta = beta))
  }
  
  # Initialize Markov Switching parameters
  set.seed(123)
  pi <- rep(1 / num_states, num_states)  # Initial state probabilities
  trans_probs <- matrix(0.9 / (num_states - 1), nrow = num_states, ncol = num_states) +
    diag(0.1, num_states)  # High self-transition probabilities
  trans_probs <- trans_probs / rowSums(trans_probs)
  
  # Initialize state-dependent VECM parameters
  vecm_params <- list()
  model_specific <- list()  # To store Alpha, Beta, and intercepts for each regime
  
  for (s in 1:num_states) {
    Pi <- matrix(runif(K^2, -0.5, 0.5), nrow = K, ncol = K)  # Random initialization of Pi
    Gamma <- array(runif(K^2 * lag_order, -0.5, 0.5), dim = c(K, K, lag_order))
    intercept <- if (include == "const") rep(0, K) else NULL  # Proper initialization
    Sigma <- diag(K)
    
    # Decompose Pi into Alpha and Beta
    decomposition <- decompose_pi(Pi, rank)
    alpha <- decomposition$alpha
    beta <- decomposition$beta
    Pi <- alpha %*% t(beta)  # Reconstruct Pi to ensure proper rank
    
    vecm_params[[s]] <- list(
      Pi = Pi,
      Gamma = Gamma,
      intercept = intercept,
      Sigma = Sigma
    )
    model_specific[[s]] <- list(
      alpha = alpha,
      beta = beta,
      intercept = intercept
    )
  }
  
  log_likelihood <- NULL
  forward_probs <- matrix(0, T - lag_order, num_states)  # Forward probabilities
  smoothed_probs <- matrix(1 / num_states, T - lag_order, num_states)  # Initialize to uniform probabilities
  fitted_values <- array(0, dim = c(T - lag_order, K, num_states))  # Fitted values
  residuals <- array(0, dim = c(T - lag_order, K, num_states))  # Residuals
  
  # EM Algorithm
  for (iter in 1:max_iter) {
    # E-step: Compute residuals and forward probabilities
    for (s in 1:num_states) {
      for (t in 1:(T - lag_order)) {
        # Long-run dynamics
        lagged_terms <- vecm_params[[s]]$Pi %*% matrix(lagged_data[t, 1:K], ncol = 1)
        # Short-run dynamics
        short_run_terms <- rowSums(sapply(1:lag_order, function(p) {
          Gamma <- vecm_params[[s]]$Gamma
          if (t - p + 1 > 0) {
            Gamma[, , p] %*% matrix(diff_data[t - p + 1, ], ncol = 1)
          } else {
            matrix(0, nrow = K, ncol = 1)  # Return zero matrix for invalid indices
          }
        }))
        # Reshape intercept to ensure compatibility, only if included
        intercept <- if (include == "const") matrix(vecm_params[[s]]$intercept, nrow = K, ncol = 1) else matrix(0, nrow = K, ncol = 1)
        # Fitted values
        fitted_values[t, , s] <- (lagged_terms + short_run_terms + intercept)[, 1]
        # Residuals
        residuals[t, , s] <- diff_data[t, ] - fitted_values[t, , s]
      }
      # Compute densities for each state
      densities <- apply(residuals[, , s], 1, function(res) {
        dmvnorm(res, mean = rep(0, K), sigma = vecm_params[[s]]$Sigma)
      })
      forward_probs[, s] <- pmax(densities, 1e-12)  # Clamp probabilities to avoid underflow
    }
    
    # Normalize forward probabilities
    forward_probs <- forward_probs / rowSums(forward_probs)
    if (any(abs(rowSums(forward_probs) - 1) > 1e-6)) {
      cat("Warning: Forward probabilities are not properly normalized.\n")
    }
    gamma <- forward_probs  # Filtered state probabilities
    
    # M-step: Update intercepts based on smoothed probabilities, only if included
    if (include == "const") {
      for (s in 1:num_states) {
        numerator <- colSums(sapply(1:(T - lag_order), function(t) {
          gamma[t, s] * (diff_data[t, ] - (fitted_values[t, , s] - vecm_params[[s]]$intercept))
        }))
        denominator <- sum(gamma[, s])
        vecm_params[[s]]$intercept <- numerator / denominator
        vecm_params[[s]]$intercept <- vecm_params[[s]]$intercept[1:K]
        model_specific[[s]]$intercept <- vecm_params[[s]]$intercept
      }
    } else {
      for (s in 1:num_states) {
        vecm_params[[s]]$intercept <- NULL
        model_specific[[s]]$intercept <- NULL
      }
    }
    # Update Sigma
    for (s in 1:num_states) {
      weighted_residuals <- array(0, dim = c(K, K))  # Initialize matrix to accumulate weighted outer products
      for (t in 1:(T - lag_order)) {
        res <- matrix(residuals[t, , s], ncol = 1)  # Ensure residuals are treated as column vectors
        weighted_residuals <- weighted_residuals + gamma[t, s] * (res %*% t(res))  # Weighted outer product
      }
      #vecm_params[[s]]$Sigma <- weighted_residuals / sum(gamma[, s])  # Normalize by the sum of weights
      #vecm_params[[s]]$Sigma <- weighted_residuals / sum(gamma[, s]) + diag(1e-6, K) # Regularization for Sigma
      vecm_params[[s]]$Sigma <- diag(diag(weighted_residuals / sum(gamma[, s])) + 1e-6)  # Simplifies Sigma
    }
    
    # Update Transition Probabilities
    for (i in 1:num_states) {
      for (j in 1:num_states) {
        numerator <- sum(smoothed_probs[1:(T - lag_order - 1), i] * forward_probs[2:(T - lag_order), j])
        denominator <- sum(smoothed_probs[1:(T - lag_order - 1), i])
        trans_probs[i, j] <- ifelse(denominator > 0, numerator / denominator, 0)
      }
    }
    trans_probs <- trans_probs / rowSums(trans_probs)  # Normalize rows
    
    # Backward pass for smoothed probabilities
    #smoothed_probs[T - lag_order, ] <- gamma[T - lag_order, ]
    smoothed_probs[T - lag_order, ] <- gamma[T - lag_order, ] / sum(gamma[T - lag_order, ])  # Robust initialization
    for (t in (T - lag_order - 1):1) {
      next_probs <- smoothed_probs[t + 1, ] %*% trans_probs
      
      # Safeguard against numerical underflow
      max_val <- max(next_probs)
      if (max_val > 0) {
        next_probs <- next_probs / max_val
      }
      
      if (any(is.na(next_probs)) || sum(next_probs) == 0) {
        cat("Warning: NA or zero sum in next_probs at time", t, "\n") # Print or debug intermediate values
        smoothed_probs[t, ] <- gamma[t, ]
      } else {
        smoothed_probs[t, ] <- gamma[t, ] * next_probs
        #smoothed_probs[t, ] <- smoothed_probs[t, ] / sum(smoothed_probs[t, ])
        epsilon <- 1e-10
        smoothed_probs[t, ] <- smoothed_probs[t, ] / (sum(smoothed_probs[t, ]) + epsilon)
        
        # Additional check: Ensure probabilities sum to 1
        if (abs(sum(smoothed_probs[t, ]) - 1) > 0.01) {
          cat("Warning: Smoothed probabilities at time", t, "do not sum to 1. Sum:",
              sum(smoothed_probs[t, ]), "\n")
        }
      }
      
      # Debugging specific time points (controlled by the debug parameter)
      if (debug && t %% 50 == 0) {  # Print debug info every 50 steps if debug is TRUE
        cat("Time:", t, "\nNext Probs:", next_probs, "\nSmoothed Probs:", smoothed_probs[t, ], "\n")
      }
    }

    
    # Compute log-likelihood
    new_log_likelihood <- sum(log(pmax(rowSums(forward_probs), 1e-12)))
    if (!is.null(log_likelihood) && abs(new_log_likelihood - log_likelihood) < tol) {
      break
    }
    log_likelihood <- new_log_likelihood
    cat("Iteration:", iter, "Log-Likelihood:", new_log_likelihood, "\n")
  }
  
  # Calculate number of parameters
  n_params_trans <- num_states * (num_states - 1)
  n_params_pi <- num_states * (rank * (K + K - rank))
  n_params_gamma <- num_states * (lag_order * K^2)
  n_params_sigma <- num_states * (K * (K + 1) / 2)
  n_params_intercept <- if (include == "const") num_states * K else 0
  total_params <- n_params_trans + n_params_pi + n_params_gamma + n_params_sigma + n_params_intercept
  
  list(
    transition_probs = trans_probs,
    vecm_params = vecm_params,
    filtered_probs = gamma,
    smoothed_probs = smoothed_probs,
    log_likelihood = log_likelihood,
    fitted_values = fitted_values,
    residuals = residuals,
    model_specific = model_specific,  # Contains Alpha, Beta, and intercepts for each regime
    n_params = total_params  # Number of estimated parameters
  )
}
