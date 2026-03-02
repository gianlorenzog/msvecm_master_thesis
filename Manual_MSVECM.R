# Libraries
library(tseries)

# setwd
rm(list=ls())
getwd()
setwd("/Users/gianlorenzo/Desktop/MSVECM_Rcode")

# import enviroment
load("obj.rData")
# df.xts contains the log of the TS at levels
# lev is df.xts converted as data.frame
diff_xts = diff(df.xts)
data = as.matrix(lev)
#--------------------------------------
###### MANUAL MS-VECM ######
#--------------------------------------
# Simulated Data
# set.seed(123)
# simdata <- matrix(rnorm(400), ncol = 4)


#--- Run a VECM Model without Regime Shifts
source(file = "R_functions/VECM.R")
basic_VECM <- fit_vecm(data, lag_order = 2, rank = 2, include = "const")
print(basic_VECM)


#--- Run a MSVECM
source(file = "R_functions/MSVECM_fun.R")
msVECM <- markov_switching_vecm(data, num_states = 2, lag_order = 2, rank = 2, include = "const")


filtered_probs <- msVECM$filtered_probs
smoothed_probs <- msVECM$smoothed_probs


#---- Plot Filtered/Smoothed State Probabilities (Aggregate by regime)
par(mfrow = c(1, 1))   # par(mfrow = c(2, 1))  
# Plot filtered probabilities
matplot(filtered_probs, type = "l", lty = 1, col = c("blue", "red"),
        ylab = "Filtered Probability", xlab = "Time",
        main = "Filtered State Probabilities for Markov Switching VECM")
legend("topright", legend = c("State 1", "State 2"), col = c("blue", "red"), lty = 1)
#
# Plot smoothed probabilities
matplot(smoothed_probs, type = "l", lty = 1, col = c("blue", "red"),
        ylab = "Smoothed Probability", xlab = "Time",
        main = "Smoothed State Probabilities for Markov Switching VECM")
legend("topright", legend = c("State 1", "State 2"), col = c("blue", "red"), lty = 1)
#----


#---- Plot Filtered/Smoothed State Probabilities (Separate Regime)
source(file = "R_functions/plot_prob.R")
plot_prob(y= smoothed_probs, prob = "smoothed")
plot_prob(y= filtered_probs, prob = "filtered")

#---- Plot smoothed and filtered probabilities (Aggregate prob) with regime splitting (2 graphs)
source(file = "R_functions/plot_ag_prob.R")
plot_ag_prob(smoothed_probs, filtered_probs, time)


plot(smoothed_probs[, 1], type = "l", col = "red", ylim = c(0, 1),
     ylab = "Probabilities", xlab = "Time",
     main = "Smoothed vs Filtered Probabilities - Regime 1")
polygon(c(1:length(filtered_probs[, 1]), rev(1:length(filtered_probs[, 1]))),
        c(filtered_probs[, 1], rep(0, length(filtered_probs[, 1]))),
        col = rgb(0.5, 0.5, 0.5, alpha = 0.5), border = NA)
lines(filtered_probs[, 1], col = "grey", lwd = 1)
   
   
plot(smoothed_probs[, 1], type = "l", col = "red", ylim = c(0, 1),
     ylab = "Probabilities", xlab = "Time",
     main = "Smoothed vs Filtered Probabilities - Regime 1")
polygon(c(1:length(filtered_probs[, 1]), rev(1:length(filtered_probs[, 1]))),
        c(filtered_probs[, 1], rep(0, length(filtered_probs[, 1]))),
        col = rgb(0.8, 0.8, 0.8, alpha = 0.5), border = NA)
lines(filtered_probs[, 1], col = "black", lwd = 1)


#---- Plot a variable with probabilities 
titles <- c("Eurostoxx50", "GDP", "HICP", "M1")

# Loop and plot it with the respective title
par(mfrow = c(4, 1), mar = c(4, 4, 2, 1))  # Adjust layout for 4 variables
for (i in 1:dim(data)[2]) {
  plot(data[, i], type = "l", col = "blue", ylab = "Value", xlab = "Time",
       main = paste("Variable", titles[i], "with Filtered and Smoothed Probabilities"))
  
  # Rescale probabilities to match the scale of the variable
  rescaled_filtered <- filtered_probs[, 1] * (max(data[, i]) - min(data[, i])) + min(data[, i])
  rescaled_smoothed <- smoothed_probs[, 1] * (max(data[, i]) - min(data[, i])) + min(data[, i])
  
  # Add filtered and smoothed probabilities to the plot
  lines(rescaled_filtered, col = "red", lty = 2, lwd = 2)
  lines(rescaled_smoothed, col = "green", lty = 3, lwd = 2)
  
  legend("topright", legend = c("Variable", "Filtered Prob", "Smoothed Prob"),
         col = c("blue", "red", "green"), lty = c(1, 2, 3))
}
#----


#---- Plot all variables with regime bands (both regimes in the same variable's graph)
regime <- max.col(smoothed_probs)  # Get the most probable regime at each time point
titles <- c("Eurostoxx50", "GDP", "HICP", "M1")  # Variable titles

par(mfrow = c(4, 1), mar = c(4, 4, 2, 1))  # Adjust layout for 4 variables
for (j in 1:dim(data)[2]) {  # Use `j` for the outer loop to avoid overwriting
  plot(data[, j], type = "l", col = "black", ylab = "Value", xlab = "Time",
       main = paste("Variable", titles[j], "with Regime Bands"))
  
  # Add regime bands
  for (t in 1:length(data[, j])) {  # Use `t` for the inner loop
    rect(t - 0.5, min(data[, j]), t + 0.5, max(data[, j]),
         col = ifelse(regime[t] == 1, rgb(1, 0, 0, alpha = 0.2),  # Red for regime 1
                      rgb(0, 0, 1, alpha = 0.2)), border = NA)  # Blue for regime 2
  }
  lines(data[, j], col = "black")  # Add the variable line
}
#----


#---- Plot a single variable with separate regime bands
source(file = "R_functions/plot_bands_onevar.R")

# Plot the SX5E variable with significant regime bands
plot_bands_onevar(data[, "SX5E"], smoothed_probs, threshold = 0.5,
                  main = "SX5E with Significant Regime Bands")
#---- 


#---- Plot the histogram and density for this state and variable (8 Graphs)  -->   VERIFICARE I RISULTATI
# (8 graphs - every variables in each states)
plot_vecm_state_dependent_distributions <- function(data, smoothed_probs, threshold = 0.7, bins = 30) {
  # Get the names of the variables
  var_names <- colnames(data)
  num_vars <- ncol(data)
  
  # Set up layout for each variable across two states
  par(mfrow = c(num_vars, 2), mar = c(4, 4, 2, 2))
  
  for (var_index in 1:num_vars) {
    variable <- data[, var_index]
    
    for (state in 1:2) {
      # Identify time points where the smoothed probability for the current state is above the threshold
      state_indices <- which(smoothed_probs[, state] > threshold)
      
      # Extract values of the current variable for this state
      state_values <- variable[state_indices]
      
      # Plot the histogram and density for this state and variable
      hist(
        state_values, breaks = bins, probability = TRUE, main = paste(var_names[var_index], "- State", state),
        xlab = var_names[var_index], col = ifelse(state == 1, "blue", "red"), border = "white"
      )
      
      # Overlay density curve
      lines(density(state_values), col = "black", lwd = 2)
    }
  }
  
  # Reset layout
  par(mfrow = c(1, 1))
}
#
plot_vecm_state_dependent_distributions(data, smoothed_probs)
#---- 



#--------------------------------------
###### POST-ESTIMATION TESTS ######
#--------------------------------------

#--- STORE VALUES FOR TESTS
# Restricted Model- parameters do not change across regimes.
# Unrestricted Model - parameters vary with regimes.
#
# Log-likelihoods (from models)
ll_restricted <- basic_VECM$log_likelihood  # Log-likelihood from single-regime VECM
ll_unrestricted <- msVECM$log_likelihood # Log-likelihood from MS-VECM

# Number of parameters
num_params_rs <- basic_VECM$num_params # Restricted model parameters
num_params_un <- msVECM$n_params # Unrestricted model parameters

num_obs <- 328             # Number of observations (not directly used in LR test)
#---



#--- General Likelihood Ratio (LR) Test
# H0: Restricted model is sufficient.
source(file = "R_functions/LRtest.R")

# Perform the test
general_LR <- lr_test(ll_unrestricted, ll_restricted, num_params_un, num_params_rs)
print(general_LR)
#---


#--- Hansen (1992) Likelihood Ratio Test
# Tests parameter stability across regimes or structural breaks.
# H0: No regime shifts (parameters constant).
source(file = "R_functions/LRtest_hansen1992.R")
num_params <- num_params_un
# Perform the test
hansenLR <- hansen_likelihood_ratio_test(ll_unrestricted, ll_restricted, num_params)
print(hansenLR)
#---


#---  AIC and BIC
# Function to calculate AIC and BIC
source(file = "R_functions/calculate_aic_bic.R")

calculate_aic_bic(ll_unrestricted, num_params_un, 328) # AIC/BIC - MS-VECM
calculate_aic_bic(ll_restricted, num_params_rs, 328) # AIC/BIC - VECM
#---- 




#--- Normality Test of the Residuals
# Extract Markov Switching Model Residuals
# res1 = msVECM$residuals[,,1][,1] # regime 1 residuals Eurostoxx50
titles <- c("Eurostoxx50", "GDP", "HICP", "M1")  # Variable titles
res <- list() # Initialize the outer list
for (i in 1:2) {
  res[[i]] <- list() # Initialize the inner list for each i
  for (j in 1:4) {
    res[[i]][[j]] <- msVECM$residuals[,,i][,j]
  }
  names(res[[i]]) <- titles # Assign names to the inner list
}
names(res) <- c("res_rg1", "res_rg2") # Assign names to the outer list


# Jarque Bera test      H_0 : Normality. If p-value < 0.05 we reject H_0
#regime 1
for (i in 1:4){
  print(jarque.bera.test(res$res_rg1[[i]]))
  cat("regime 1 - residuals:", titles[i], "\n\n")
}
#regime 2
for (i in 1:4){
  print(jarque.bera.test(res$res_rg2[[i]]))
  cat("regime 2 - residuals:", titles[i], "\n\n")
}

## QQ-Norm Plots for Residuals
par(mfrow = c(2, 4))  # Layout for 8 plots (2 rows, 4 columns)

# Regime 1
for (i in 1:4) {
  qqnorm(res$res_rg1[[i]], main = paste(titles[i], "- Reg 1 Res"))
  qqline(res$res_rg1[[i]], col = 2, distribution = qnorm)
}

# Regime 2
for (i in 1:4) {
  qqnorm(res$res_rg2[[i]], main = paste(titles[i], "- Reg 2 Res"))
  qqline(res$res_rg2[[i]], col = 2, distribution = qnorm)
}


# skweness, kurtosis, mean, median
#regime 1
for (i in 1:4){
  cat("regime 1 - residuals:", titles[i], "\n")
  cat("Skewness:", skewness(res$res_rg1[[i]], na.rm = TRUE), "\n")
  cat("Kurtosis:", kurtosis(res$res_rg1[[i]], na.rm = TRUE), "\n")
  cat("Mean:", mean(res$res_rg1[[i]], na.rm = TRUE), "\n")
  cat("Median:", median(res$res_rg1[[i]], na.rm = TRUE), "\n\n")
}
#regime 2
for (i in 1:4){
  cat("regime 2 - residuals:", titles[i], "\n")
  cat("Skewness:", skewness(res$res_rg2[[i]], na.rm = TRUE), "\n")
  cat("Kurtosis:", kurtosis(res$res_rg2[[i]], na.rm = TRUE), "\n")
  cat("Mean:", mean(res$res_rg2[[i]], na.rm = TRUE), "\n")
  cat("Median:", median(res$res_rg2[[i]], na.rm = TRUE), "\n\n")
}

# density distribution of the returns - regime 1
par(mfrow = c(1, 1))
for (i in 1:4) {
  data <- data.frame(Returns = as.numeric(res$res_rg1[[i]]))
  print(
    ggplot(data, aes(x = Returns)) +
      geom_density(fill = 'darkorange', color = 'black', alpha = .6) +
      theme_minimal() +
      stat_function(fun = dnorm, color = 'blue', size = 1,
                    args = list(mean = mean(data$Returns, na.rm = TRUE),
                                sd = sd(data$Returns, na.rm = TRUE))) +
      labs(y = "Density", x = "Returns", title =  paste("residuals", titles[i]))
  )
  rm(data)
}

# regime 2
par(mfrow = c(1, 1))
for (i in 1:4) {
  data <- data.frame(Returns = as.numeric(res$res_rg2[[i]]))
  print(
    ggplot(data, aes(x = Returns)) +
      geom_density(fill = 'darkorange', color = 'black', alpha = .6) +
      theme_minimal() +
      stat_function(fun = dnorm, color = 'blue', size = 1,
                    args = list(mean = mean(data$Returns, na.rm = TRUE),
                                sd = sd(data$Returns, na.rm = TRUE))) +
      labs(y = "Density", x = "Returns", title =  paste("residuals", titles[i]))
  )
  rm(data)
}

# ACF and PACF of the residuals
res1 = as.data.frame(res$res_rg1)  # res regime 1
res2 = as.data.frame(res$res_rg2)  # res regime 2
res = cbind(res1, res2)
colnames(res) <- c("Eurostoxx50 reg1", "GDP reg1", "HICP reg1", "M1 reg1",
                   "Eurostoxx50 reg2", "GDP reg2", "HICP reg2", "M1 reg2")

par(mfrow=c(4,4))
for (i in 1:ncol(res)) { # Loop over each column of residuals
  acf(coredata(res[, i]), main = paste("ACF of", colnames(res)[i])) # ACF plot
  pacf(coredata(res[, i]), main = paste("PACF of", colnames(res)[i])) # PACF plot
}
#---


#--- Test Granger Causality - Wald test
# test whether the lagged coefficients of the causal variable (x) are jointly zero in each regime. 
source(file = "R_functions/wald_test_msvecm.R")

# STEP 1: Extract Coeff and SE
# Extract Coeffiecients
Gamma_list <- c()
Gamma_list$Gamma_regime1 <- msVECM$vecm_params[[1]]$Gamma
Gamma_list$Gamma_regime2 <- msVECM$vecm_params[[2]]$Gamma

Pi_list <- c()
Pi_list$Pi_regime1 <- msVECM$vecm_params[[1]]$Pi
Pi_list$Pi_regime2 <- msVECM$vecm_params[[2]]$Pi

#--- Extract Standard Errors from coeff
# Compute standard errors for Gamma - REGIME 1
bootstrap_gamma <- replicate(500, {
  boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
  boot_model <- markov_switching_vecm(boot_data, include = "const")
  boot_model$vecm_params[[1]]$Gamma  # regime 1
})
se_gamma_1 <- apply(bootstrap_gamma, 1, sd)

# Compute standard errors for Pi - REGIME 1
bootstrap_pi <- replicate(500, {
  boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
  boot_model <- markov_switching_vecm(boot_data, include = "const")
  boot_model$vecm_params[[1]]$Pi # regime 1
})
se_pi_1 <- apply(bootstrap_pi, 1, sd)

# Compute standard errors for Gamma - REGIME 2
bootstrap_gamma <- replicate(500, {
  boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
  boot_model <- markov_switching_vecm(boot_data, include = "const")
  boot_model$vecm_params[[2]]$Gamma  # regime 2
})
se_gamma_2 <- apply(bootstrap_gamma, 1, sd)

# Compute standard errors for Pi - REGIME 2
bootstrap_pi <- replicate(500, {
  boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
  boot_model <- markov_switching_vecm(boot_data, include = "const")
  boot_model$vecm_params[[2]]$Pi  # regime 2
})
se_pi_2 <- apply(bootstrap_pi, 1, sd)
#--- 

# STEP 2: Perform Wald test
# GDP, HICP, and M1 Granger-causing SX5E (and vice versa):
wald_test(Pi_list, Gamma_list, 1, c(2,1), 2, se_gamma_1, se_pi_1, TRUE) #GDP -SX5E
wald_test(Pi_list, Gamma_list, 1, c(3,1), 2, se_gamma_1, se_pi_1, TRUE) #HICP-SX5E
wald_test(Pi_list, Gamma_list, 1, c(4,1), 2, se_gamma_1, se_pi_1, TRUE) #M1 - SX5E

wald_test(Pi_list, Gamma_list, 1, c(2,1), 2, se_gamma_2, se_pi_2, TRUE) #SX5E -GDP
wald_test(Pi_list, Gamma_list, 2, c(3,1), 2, se_gamma_2, se_pi_2, TRUE) #SX5E-HICP
wald_test(Pi_list, Gamma_list, 2, c(4,1), 2, se_gamma_2, se_pi_2, TRUE) #SX5E - M1


# GDP Granger-causing HICP and M1 (and vice versa):
wald_test(Pi_list, Gamma_list, 1, c(3,4), 2, se_gamma_1, se_pi_1, TRUE) #HICP - M1
wald_test(Pi_list, Gamma_list, 1, c(4,3), 2, se_gamma_1, se_pi_1, TRUE) #M1 - HICP

wald_test(Pi_list, Gamma_list, 2, c(3,4), 2, se_gamma_2, se_pi_2, TRUE) #HICP - M1
wald_test(Pi_list, Gamma_list, 2, c(4,3), 2, se_gamma_2, se_pi_2, TRUE) #M1 - HICP

# M1 Granger-causing GDP (and vice versa):
wald_test(Pi_list, Gamma_list, 1, c(3,2), 2, se_gamma_1, se_pi_1, TRUE) #HICP - GDP
wald_test(Pi_list, Gamma_list, 1, c(2,3), 2, se_gamma_1, se_pi_1, TRUE) #GDP - HICP

wald_test(Pi_list, Gamma_list, 2, c(3,2), 2, se_gamma_2, se_pi_2, TRUE) #HICP - GDP
wald_test(Pi_list, Gamma_list, 2, c(2,3), 2, se_gamma_2, se_pi_2, TRUE) #GDP - HICP


#---
# Visualize Regime Characteristics:
# Plot the variances (diagonal elements of Sigma) across regimes to interpret differences:
diag_sigma_1 <- diag(msVECM$vecm_params[[1]]$Sigma)
diag_sigma_2 <- diag(msVECM$vecm_params[[2]]$Sigma)
plot(1:length(diag_sigma_1), diag_sigma_1, type = "b", col = "blue", ylim = range(diag_sigma_1, diag_sigma_2))
points(1:length(diag_sigma_2), diag_sigma_2, type = "b", col = "red")
legend("topright", legend = c("Regime 1", "Regime 2"), col = c("blue", "red"), lty = 1)


#------------------------------------------------------------
# Impulse Response Function
#------------------------------------------------------------
#
#--- IRF for Regime Switching
shock_vector <- c(0, 1, 1, 1)  # Multi-variable shock
horizon <- 10
num_simulations <- 1000
regime =1


source(file = "R_functions/irf_msvecm_v2.R")
irf <- calculate_irf_switching(msVECM, regime, shock_vector, 
                               horizon, num_simulations)

# Plot IRF
matplot(seq_len(nrow(irf)) - 1, irf, type = "l", lty = 1,
        xlab = "Horizon", ylab = "Response", main = "Switching IRF")
legend("topright", legend = paste("Variable", 1:ncol(irf)),
       col = 1:ncol(irf), bty = "n", lty = 1)

#------------------------
# CIRF
calc_cirf <- function(irf, horizon = 10){
  # Initialize CIRF matrix
  cirf <- matrix(0, nrow = horizon + 1, ncol = ncol(irf))
  
  # Compute cumulative IRF by summing IRF over horizons for each variable
  for (h in 1:(horizon + 1)) {
    cirf[h, ] <- colSums(irf[1:h, , drop = FALSE])
  }
  
  return(cirf)
}

cirf_result  <- calc_cirf(irf, horizon)

# easier alternative
cirf_result <- apply(irf, 2, cumsum)


# Plot CIRF
par(mfrow = c(1, 1))
matplot(seq_len(nrow(cirf_result)) - 1, cirf_result, type = "l", lty = 1,
        xlab = "Horizon", ylab = "Cumulative Response", main = "Cumulative Switching IRF")
legend("topright", legend = paste("Variable", 1:ncol(cirf_result)),
       col = 1:ncol(cirf_result), bty = "n", lty = 1)

#
#------------------------------------------------------------
# Forecast Error Variance Decomposition
#------------------------------------------------------------
#
#--- FEVD - general approach

# Generate IRF matrix v2 --> uses irf_msvecm_v2
source(file = "R_functions/gen_irf_matrix_v2.R")
irf_matrix <- generate_irf_matrix(msVECM, regime = 1, horizon = horizon, 
                                  num_simulations = num_simulations)

# FEVD - general
Sigma <- list(msVECM$vecm_params[[1]]$Sigma, msVECM$vecm_params[[2]]$Sigma)

source(file = "R_functions/FEVD_general.R") # works with both v2 and v3
fevd_result <- compute_fevd_matrix(irf_matrix, Sigma_list = Sigma, horizon, 
                                   transition_matrix = msVECM$transition_probs,
                                   initial_regime = 1)


# Verify FEVD sums
fevd_sums <- apply(fevd_result, c(1, 2), sum)
print(fevd_sums)  # Should now be exactly 1


## plot multiple variables for comparison
# Create a 2x2 plotting layout
par(mfrow = c(2, 2))  # 2 rows, 2 columns of plots

# Plot for Variable 1
matplot(0:10, fevd_result[, 1, ], type = "l", lty = 1, xlab = "Horizon", ylab = "FEVD",
        main = "FEVD of Euro Stoxx 50")

# Plot for Variable 2
matplot(0:10, fevd_result[, 2, ], type = "l", lty = 1, xlab = "Horizon", ylab = "FEVD",
        main = "FEVD of GDP")

# Plot for Variable 3
matplot(0:10, fevd_result[, 3, ], type = "l", lty = 1, xlab = "Horizon", ylab = "FEVD",
        main = "FEVD of HICP")

# Plot for Variable 4
matplot(0:10, fevd_result[, 4, ], type = "l", lty = 1, xlab = "Horizon", ylab = "FEVD",
        main = "FEVD of M1")

par(mfrow = c(1, 1))
plot.new()
legend("topright", legend = paste("Shock in Variable", 1:4), col = 1:4, bty = "n", lty = 1)
#--- 

# Reset layout to single plot
par(mfrow = c(1, 1))



#--- CFEVD
source(file = "R_functions/CFEVD.R")

# Compute CFEVD from the FEVD results
cfevd_result <- compute_cfevd(fevd_result)


## plot multiple variables for comparison
# Create a 2x2 layout for CFEVD plots
par(mfrow = c(2, 2))  # 2 rows, 2 columns of plots

# Plot for Variable 1
matplot(0:horizon, cfevd_result[, 1, ], type = "l", lty = 1, xlab = "Horizon", ylab = "CFEVD",
        main = "CFEVD of Euro Stoxx 50")

# Plot for Variable 2
matplot(0:horizon, cfevd_result[, 2, ], type = "l", lty = 1, xlab = "Horizon", ylab = "CFEVD",
        main = "CFEVD of GDP")

# Plot for Variable 3
matplot(0:horizon, cfevd_result[, 3, ], type = "l", lty = 1, xlab = "Horizon", ylab = "CFEVD",
        main = "CFEVD of HICP")

# Plot for Variable 4
matplot(0:horizon, cfevd_result[, 4, ], type = "l", lty = 1, xlab = "Horizon", ylab = "CFEVD",
        main = "CFEVD of M1")

# Reset layout to single plot
par(mfrow = c(1, 1))
plot.new()
legend("topright", legend = paste("Shock in Variable", 1:4), col = 1:4, bty = "n", lty = 1)
##---------------------------
# END
##---------------------------






#---------------------------
# Alternatives for IRF and FEVD
#---------------------------

#--- Alternative - IRF for Regime Switching
shock_vector <- c(0, 1, 1, 1)  # Multi-variable shock
horizon <- 10
num_simulations <- 1000
regime =1

source(file = "R_functions/irf_msvecm_v3.R")
irf <- calculate_irf_switching(msVECM, shock_vector,
                               msVECM$transition_probs,
                               horizon = 10, initial_regime = 1)

##---------------------------


#--- Alternative 1 - FEVD - Integrated approach
# irf and irf_matrix computed inside the function
#
#--- FEVD - integrated function - uses "irf_msvecm_v2"
source(file = "R_functions/FEVD_msvecm.R")

# Parameters
horizon <- 10
num_simulations <- 1000
regime <- 1
transition_matrix <- msVECM$transition_probs

# Compute FEVD with integrated IRF computation
fevd_result <- compute_fevd_irf(msvecm_model = msVECM, regime = regime, 
                                horizon = horizon, 
                                num_simulations = num_simulations, 
                                transition_matrix = transition_matrix)
##---------------------------


#--- Alternative 2 - FEVD
#
# Generate IRF matrix v3 --> uses irf_msvecm_v3
# irf_msvecm is more stochastic
source(file = "R_functions/gen_irf_matrix_v3.R")

irf_matrix <- gen_irf_matrix(
  msvecm_model = msVECM,
  transition_matrix = msVECM$transition_probs,
  horizon = 10,
  initial_regime = 2
)

# FEVD - general
Sigma <- list(msVECM$vecm_params[[1]]$Sigma, msVECM$vecm_params[[2]]$Sigma)

source(file = "R_functions/FEVD_general.R") # works with both v2 and v3
fevd_result <- compute_fevd_matrix(irf_matrix, Sigma_list = Sigma, horizon, 
                                   transition_matrix = msVECM$transition_probs,
                                   initial_regime = 1)
##---------------------------

