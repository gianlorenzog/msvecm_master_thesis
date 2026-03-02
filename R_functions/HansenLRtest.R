hansen_lr_msvecm <- function(
    y, 
    B = 499,
    seed = 123,
    lag_order = 1,
    rank = 1,
    include = c("const","none"),
    num_states = 2,
    verbose = TRUE
) {
  include <- match.arg(include)
  
  if (!is.null(seed)) set.seed(seed)
  
  Y <- as.matrix(y)
  Tobs <- nrow(Y)
  
  #------------------------------------------------------------
  # 1. Fit H0: linear VECM (no switching)
  #------------------------------------------------------------
  fit_H0 <- fit_linear_vecm(Y, lag_order, rank, include)
  
  # Gaussian log-likelihood
  E0 <- fit_H0$residuals
  Sigma0 <- crossprod(E0) / nrow(E0)
  k <- ncol(E0)
  T_eff <- nrow(E0)
  ll0 <- -0.5 * T_eff * (k * log(2*pi) + log(det(Sigma0)) + k)
  
  #------------------------------------------------------------
  # 2. Fit H1: your MS-VECM
  #------------------------------------------------------------
  fit_H1 <- markov_switching_vecm(
    data = Y,
    num_states = num_states,
    lag_order = lag_order,
    rank = rank,
    include = include,
    max_iter = 100,
    tol = 1e-6
  )
  
  ll1 <- fit_H1$log_likelihood
  LR_obs <- 2 * (ll1 - ll0)
  
  if (verbose) {
    cat("H0 LogLik:", ll0, "\n")
    cat("H1 LogLik:", ll1, "\n")
    cat("Observed LR:", LR_obs, "\n")
  }
  
  #------------------------------------------------------------
  # 3. Bootstrap under H0
  #------------------------------------------------------------
  LR_boot <- numeric(B)
  
  for (b in 1:B) {
    if (verbose && b %% max(1, B %/% 10) == 0) cat("Bootstrap", b, "of", B, "\n")
    
    # Simulate bootstrap sample under H0
    Yb <- simulate_linear_vecm(
      fit_H0,
      T = Tobs,
      lag_order = lag_order,
      include = include
    )
    
    # Refit H0
    f0b <- fit_linear_vecm(Yb, lag_order, rank, include)
    Eb <- f0b$residuals
    Sigma0b <- crossprod(Eb) / nrow(Eb)
    ll0b <- -0.5 * nrow(Eb) * (k*log(2*pi) + log(det(Sigma0b)) + k)
    
    # Refit H1 (MS-VECM)
    f1b <- try(
      markov_switching_vecm(
        data = Yb,
        num_states = num_states,
        lag_order = lag_order,
        rank = rank,
        include = include,
        max_iter = 100,
        tol = 1e-6
      ),
      silent = TRUE
    )
    
    if (inherits(f1b, "try-error")) {
      LR_boot[b] <- NA
    } else {
      LR_boot[b] <- 2 * (f1b$log_likelihood - ll0b)
    }
  }
  
  LR_boot <- LR_boot[is.finite(LR_boot)]
  
  pval <- mean(LR_boot >= LR_obs)
  crit <- quantile(LR_boot, c(.90,.95,.99))
  
  list(
    LR_obs = LR_obs,
    p_value = pval,
    crit = crit,
    LR_boot = LR_boot,
    logLik = c(H0 = ll0, H1 = ll1),
    fit_H0 = fit_H0,
    fit_H1 = fit_H1
  )
}
