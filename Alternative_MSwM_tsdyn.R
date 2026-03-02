# Libraries
library(ggplot2)
library(xts)
# library(e1071) # skewness
library(fBasics) # skewness

library(MSwM) # MS
library(vars) # VAR/VECM
library(tsDyn) # VECM
library(dplyr)    # to use %>% -- loads also "last"
library(tseries)

library(magrittr) # to use %>%
library(purrr)
library(tidyverse)
library(data.table)

library(mvtnorm) # dmvnorm

# library(stargazer)
# library(broom)
# library(xtable)

# library(highfrequency) # makeReturns
###

# setwd
rm(list=ls())
getwd()
setwd("/Users/gianlorenzo/Desktop/TESI")

# import enviroment
load("obj.rData")
# df.xts contains the log of the TS at levels
# lev is df.xts converted as data.frame
diff_xts = diff(df.xts)
data = as.matrix(lev)
#--------------------------------------
# DESCRIPTIVE STATS & PLOTS
#--------------------------------------
titles <- c("Eurostoxx50", "GDP", "HICP", "M1")
#
# Loop and plot it with the respective title
par(mfrow=c(4,1))
for (i in 1:dim(df.xts)[2]) {
  plot(as.zoo(df.xts[, i]),
     col = "steelblue",
     lwd = 2,
     ylab = "Logarithm",
     xlab = "Date",
     main = titles[i])
}


# QQNORM
par(mfrow=c(2,2))
for (i in 1:dim(df.xts)[2]) {
  qqnorm(df.xts[,i], main = titles[i])
  qqline(df.xts[,i], col = 2, distribution= qnorm)
}

# skweness, kurtosis, mean, median
for (i in 1:dim(df.xts)[2]) {
  cat("Column:", colnames(df.xts)[i], "\n")
  cat("Skewness:", skewness(df.xts[, i], na.rm = TRUE), "\n")
  cat("Kurtosis:", kurtosis(df.xts[, i], na.rm = TRUE), "\n\n")
  cat("Mean:", mean(df.xts[, i], na.rm = TRUE), "\n")
  cat("Median:", median(df.xts[, i], na.rm = TRUE), "\n")
  cat("SD:", sd(df.xts[, i], na.rm = TRUE), "\n\n")
}


#density distribution of the log values
par(mfrow = c(1, 1))
for (i in 1:dim(df.xts)[2]) {
  data <- data.frame(Returns = as.numeric(df.xts[, i]))
  print(
    ggplot(data, aes(x = Returns)) +
      geom_density(fill = 'darkorange', color = 'black', alpha = .6) +
      theme_minimal() +
      stat_function(fun = dnorm, color = 'blue', size = 1,
                    args = list(mean = mean(data$Returns, na.rm = TRUE),
                                sd = sd(data$Returns, na.rm = TRUE))) +
      labs(y = "Density", x = "Log Values", title = titles[i])
  )
  rm(data)
}

#density distribution of the log differences
par(mfrow = c(1, 1))
for (i in 1:dim(df.xts)[2]) {
  data <- data.frame(Returns = as.numeric(diff_xts[, i]))
  print(
    ggplot(data, aes(x = Returns)) +
      geom_density(fill = 'darkorange', color = 'black', alpha = .6) +
      theme_minimal() +
      stat_function(fun = dnorm, color = 'blue', size = 1,
                    args = list(mean = mean(data$Returns, na.rm = TRUE),
                                sd = sd(data$Returns, na.rm = TRUE))) +
      labs(y = "Density", x = "Log Differences", title = titles[i])
  )
  rm(data)
}

#
#--------------------------------------
### TESTS ###
#--------------------------------------
#
# Normality Test
# Jarque Bera test      H_0 : Normality. If p-value < 0.05 we reject H_0
for (i in 1:dim(df.xts)[2]) {
  cat("jarque bera test for column:", colnames(df.xts)[i], "\n")
  print(jarque.bera.test(df.xts[,i]))
  cat("\n")  # Add a blank line for readability
}

# ACF
par(mfrow=c(2,2))
for (i in 1:dim(df.xts)[2]) {
  acf(coredata(df.xts[,i]), main = titles[i])
}
# PACF
par(mfrow=c(2,2))
for (i in 1:dim(df.xts)[2]) {
  pacf(coredata(df.xts[,i]), main = titles[i])
}

# UNIT ROOT
for (i in 1:dim(df.xts)[2]) {
  cat("ADF Test for column:", colnames(df.xts)[i], "\n")
  print(adf.test(coredata(df.xts[, i])))
  cat("\n")  # Add a blank line for readability
}
# H_0 = no stationarity.

# Alternative
# ur.df(df.xts$SX5E, type = "trend", lags = 6, selectlags = "AIC")

# optimal number of lags (tsdyn)
lags.select(lev, lag.max = 10, include = "trend", 
            fitMeasure ="LL", sameSample = TRUE
)

# Johansen cointegration test
johansen_test <- ca.jo(lev, type = "trace", K = 2, 
                       ecdet = "const", spec = "longrun")
summary(johansen_test)
#
#--------------------------------------
###### MODEL (MSwM + tsdyn) ######
#--------------------------------------
# ca.jorls --> Fit a VECM model (OLS)
# vecm_model <- cajorls(johansen_test, r = 2)
#
# Fit the VECM model using tsDyn --> MLE
vecm_model_mle <- VECM(lev, lag = 2, r = 2, include = "const", estim = "ML") # const
summary(vecm_model_mle)

# Test of the cointegrating rank
rank.test(vecm_model_mle, "trace", 1)
rank.test(vecm_model_mle, "trace", 2)

# Extract the data used in the VECM model
vecm_data <- vecm_model_mle$model

# rename col with bad syntax
colnames(vecm_data)[8:11] <- c('SX5Elag', 'GDPlag','HICPlag', 'M1lag') 
colnames(vecm_data)[12:15] <- c('SX5Elag2', 'GDPlag2','HICPlag2', 'M1lag2') 
# colnames(vecm_data)[16:19] <- c('SX5Elag3', 'GDPlag3','HICPlag3', 'M1lag3') 
# colnames(vecm_data)[20:23] <- c('SX5Elag4', 'GDPlag4','HICPlag4', 'M1lag4') 
# colnames(vecm_data)[24:27] <- c('SX5Elag5', 'GDPlag5','HICPlag5', 'M1lag5') 


# Construct a formula manually
dep_var <- colnames(vecm_data)[1]
# ind_vars <- colnames(vecm_data)[-1] # all --> doesn't work
ind_vars <- colnames(vecm_data)[cbind(-1,-7)] # variables + lags
# ind_vars <- colnames(vecm_data)[2:6] # variables + ECT 1 & 2
# ind_vars <- colnames(vecm_data)[5:7] # ECT + Intercept/Trend
# ind_vars <- colnames(vecm_data)[5:6] # ECT

formula <- as.formula(paste(dep_var, "~", paste(ind_vars, collapse = " + ")))
any(is.na(vecm_data)) 
vecm_data2 <- na.omit(vecm_data)

# Fit the lm model using the constructed formula and cleaned data
lm_model <- lm(formula, data = vecm_data2)

# Determine the correct length of sw = num. coeff + 1 (for the intercept)
coef(lm_model)
num_coefficients <- length(coef(lm_model)) 

# Define sw
sw <- rep(TRUE, num_coefficients + 1)

# Fit the Markov Switching model with 2 regimes
ms_model <- msmFit(lm_model, k = 2, sw = sw)
summary(ms_model)

# graphs with regime probabilities:
# par(mfrow=c(3,1)) 
# plotProb(ms_model,which=1) # filtered/smoothed prob
# plotProb(ms_model,which=2) # regime 1 prob & shaded area
# plotProb(ms_model,which=3) # regime 2 prob & shaded area


par(mfrow=c(3,1))
plotReg(ms_model, regime = 1)
plotReg(ms_model, regime = 2)


# Store smoothed and filtered probabilities
smoothed_probs <- ms_model@Fit@smoProb  # Matrix 325 x 2
filtered_probs <- ms_model@Fit@filtProb


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
plot_prob(y= filtered_probs, prob = "filtered")
plot_prob(y= smoothed_probs, prob = "smoothed")


#---- Plot smoothed and filtered probabilities (Aggregate prob) with regime splitting (2 graphs)
source(file = "R_functions/plot_ag_prob.R")
plot_ag_prob(smoothed_probs, filtered_probs, time)


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
ll_restricted <- logLik(vecm_model_mle)  # Log-likelihood for the restricted model
ll_unrestricted <- ms_model@Fit@logLikel # Log-likelihood for the unrestricted model
sum(log(ms_model@Fit@margLik))           # Check log-likelihood unrestricted

num_params_rs <- vecm_model_mle$nparB    # vecm_model_mle$npar
num_params_un <- 32                      # Number of restrictions
# check num parameteres unrestricted model (msmFit)
dim(ms_model@Coef)[1] * dim(ms_model@Coef)[2] + ncol(ms_model@transMat) + length(ms_model@std) 

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
res1 = msmResid(ms_model, regime = 1)  # ms_model@Fit@error[,1]
res2 = msmResid(ms_model, regime = 2)  # ms_model@Fit@error[,2]
res = cbind(res1, res2)
#
# Jarque Bera test      H_0 : Normality. If p-value < 0.05 we reject H_0
for (i in 1:dim(res)[2]){
  print(jarque.bera.test(res[,i]))
}

# QQNORM
par(mfrow=c(1,2))
for (i in 1:dim(res)[2]){
  qqnorm(res[,i]); qqline(res[,i], col = 2, distribution= qnorm)
}

# skweness, kurtosis, mean, median
for (i in 1:dim(res)[2]){
  cat("Skewness:", skewness(res[,i], na.rm = TRUE), "\n")
  cat("Kurtosis:", kurtosis(res[,i], na.rm = TRUE), "\n")
  cat("Mean:", mean(res[,i], na.rm = TRUE), "\n")
  cat("Median:", median(res[,i], na.rm = TRUE), "\n\n")
}

#density distribution of the residuals
par(mfrow = c(1, 1))
for (i in 1:2) {
  data <- data.frame(Returns = as.numeric(res[, i]))
  print(
    ggplot(data, aes(x = Returns)) +
      geom_density(fill = 'darkorange', color = 'black', alpha = .6) +
      theme_minimal() +
      stat_function(fun = dnorm, color = 'blue', size = 1,
                    args = list(mean = mean(data$Returns, na.rm = TRUE),
                                sd = sd(data$Returns, na.rm = TRUE))) +
      labs(y = "Density", x = "residuals", title = paste("regime", i))
  )
  rm(data)
}


# ACF and PACF of the residuals
par(mfrow=c(2,2))
for (i in 1:ncol(res)) { # Loop over each column of residuals
  acf(coredata(res[, i]), main = paste("ACF of", colnames(res)[i])) # ACF plot
  pacf(coredata(res[, i]), main = paste("PACF of", colnames(res)[i])) # PACF plot
}
#---


#--- Test Granger Causality
# test whether the lagged coefficients of the causal variable (x) are jointly zero in each regime. 
# we will use the Wald test
source(file = "R_functions/wald_test.R")

# Extract Coefficients matrices: Pi, Gamma, intercept
source(file = "R_functions/extract_coeff_mswm.R")
Coeff <- extract_coeff()


# GDP, HICP, and M1 Granger-causing SX5E (and vice versa):
wald_test(Coeff$Pi_list, Coeff$Gamma_list, 1, c(2,1), 2, TRUE) # GDP - SX5E
wald_test(Coeff$Pi_list, Coeff$Gamma_list, 1, c(3,1), 2, TRUE) # HICP - SX5E
wald_test(Coeff$Pi_list, Coeff$Gamma_list, 1, c(4,1), 2, TRUE) # M1 - SX5E

wald_test(Coeff$Pi_list, Coeff$Gamma_list, 2, c(2,1), 2, TRUE) # SX5E - GDP
wald_test(Coeff$Pi_list, Coeff$Gamma_list, 2, c(3,1), 2, TRUE) # SX5E - HICP
wald_test(Coeff$Pi_list, Coeff$Gamma_list, 2, c(4,1), 2, TRUE) # SX5E - M1


# GDP Granger-causing HICP and M1 (and vice versa):
wald_test(Coeff$Pi_list, Coeff$Gamma_list, 1, c(3,4), 2, TRUE) # HICP - M1
wald_test(Coeff$Pi_list, Coeff$Gamma_list, 1, c(4,3), 2, TRUE) # M1 - HICP

# M1 Granger-causing GDP (and vice versa):
wald_test(Coeff$Pi_list, Coeff$Gamma_list, 2, c(3,2), 2, TRUE) # HICP - GDP
wald_test(Coeff$Pi_list, Coeff$Gamma_list, 2, c(2,3), 2, TRUE) # GDP - HICP


#------------------------------------------------------------
# Impulse Response Function
#------------------------------------------------------------
#
# Extract Coefficients matrices: Pi, Gamma, intercept
source(file = "R_functions/extract_coeff_mswm.R")
Coeff <- extract_coeff()


# IRF
source(file = "R_functions/irf_general.R")
shock_vector <- c(1, 0, 0, 0)  # Shock to the first variable
horizon = 10
transition_matrix = ms_model@transMat

irf <- calculate_irf_switching(Coeff$Pi_list, Coeff$Gamma_list, 
                               intercept_list = Coeff$int_list, 
                               shock_vector, horizon, transition_matrix,
                               initial_regime = 1) 


# Plot IRF
matplot(seq_len(nrow(irf)) - 1, irf, type = "l", lty = 1,
        xlab = "Horizon", ylab = "Response", main = "Switching IRF")
legend("topright", legend = paste("Variable", 1:ncol(irf)),
       col = 1:ncol(irf), bty = "n", lty = 1)

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

cirf_result <- calc_cirf(irf, horizon)

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
#--- FEVD
source(file = "R_functions/FEVD_general.R")

# Ensure Sigma is positive definite
Sigma <- list()
Sigma$Sigma_1 <- cov(Coeff$error_list$err_reg1)
Sigma$Sigma_2 <- cov(Coeff$error_list$err_reg2)


# source(file = "R_functions/gen_irf_matrix_mswm.R")
generate_irf_matrix <- function(horizon) {
  k <- ncol(Coeff$Pi_list[[1]])   # Number of variables
  irf_matrix <- array(0, dim = c(horizon + 1, k, k))  # Initialize IRF matrix
  
  for (shock in 1:k) {
    shock_vector <- rep(0, k)  # Create shock vector
    shock_vector[shock] <- 1   # Apply shock to variable `shock`
    
    # Compute IRFs for this shock
    irf_result <- calculate_irf_switching(Coeff$Pi_list, Coeff$Gamma_list, 
                                          intercept_list = Coeff$int_list, 
                                          shock_vector, horizon, transition_matrix,
                                          initial_regime = 1) 
    irf_matrix[, , shock] <- irf_result
  }
  
  return(irf_matrix)
}


# Generate IRF matrix
irf_matrix <- generate_irf_matrix(horizon = 10)

# Normalize IRF values to ensure stability
# irf_matrix <- irf / max(abs(irf))


horizon = 10
transition_matrix = ms_model@transMat
initial_regime = 1

# Compute FEVD using the updated function
fevd_result <- compute_fevd_matrix(irf_matrix, Sigma, horizon,
                                   transition_matrix, initial_regime)

# Verify FEVD sums across shocks for each variable and horizon
fevd_sums <- apply(fevd_result, c(1, 2), sum)
print(fevd_sums)  # Should now be exactly 1

# Plot FEVD for variable 1 due to shock in variable 2
plot(0:horizon, fevd_result[, 1, 2], type = "l", xlab = "Horizon", ylab = "FEVD",
     main = "FEVD of Variable 1 due to Shock in Variable 2 (Normalized)")


## plot multiple variables for comparison
# Create a 2x2 plotting layout
par(mfrow = c(2, 3))  # 2 rows, 2 columns of plots

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

plot.new()
legend("topright", legend = paste("Shock in Variable", 1:4), col = 1:4, bty = "n", lty = 1)

# Reset layout to single plot
par(mfrow = c(1, 1))



#--- CFEVD
source(file = "R_functions/CFEVD.R")

# Compute CFEVD from the FEVD results
cfevd_result <- compute_cfevd(fevd_result)


## plot multiple variables for comparison
# Create a 2x2 layout for CFEVD plots
par(mfrow = c(2, 3))  # 2 rows, 2 columns of plots

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
#par(mfrow = c(1, 1))
plot.new()
legend("topright", legend = paste("Shock in Variable", 1:4), col = 1:4, bty = "n", lty = 1)
#---




