# Extract Coefficients from ms_model (MSwM package)
# 
# procedure: fit 4 times with msmFit and store the results, binding them together
# 
extract_coeff <- function(data = vecm_data) {
  
  # Dataset
  dataset <- data[, -7]
  any(is.na(dataset)) 
  dataset <- na.omit(dataset)
  
  # Initialize Gamma list
  Gamma_list <- list(
    Gamma_regime1 = array(NA, dim = c(4, 4, 2)),
    Gamma_regime2 = array(NA, dim = c(4, 4, 2))
  )
  
  # Initialize Pi list
  Pi_list <- list(
    Pi_regime1 = array(NA, dim = c(4, 4)),
    Pi_regime2 = array(NA, dim = c(4, 4))
  )
  
  # Initialize intercept list
  int_list <- list(
    int_regime1 = array(NA, dim = c(4, 1)),
    int_regime2 = array(NA, dim = c(4, 1))
  )
  
  # Initialize error list
  error_list <- list(
    err_reg1 = matrix(NA, nrow = nrow(dataset), ncol = 4),
    err_reg2 = matrix(NA, nrow = nrow(dataset), ncol = 4)
  )
  
  
  # Fit model and extract coefficients
  for (i in 1:4) {
    dep_var <- colnames(dataset)[i]
    ind_vars <- setdiff(colnames(dataset), dep_var)
    formula <- as.formula(paste(dep_var , "~", paste(ind_vars, collapse = " + ")))
    lm_run <- lm(formula, data = dataset)
    ms_fit <- msmFit(object = lm_run, k = 2, sw = rep(TRUE, ncol(dataset) + 1))
    
    # Extract Gamma row vector for each regime
    Gamma_row_regime1 <- ms_fit@Coef[1, -1:-6]  # Regime 1 coefficients
    Gamma_row_regime2 <- ms_fit@Coef[2, -1:-6]  # Regime 2 coefficients
    
    # Extract Pi row vector for each regime
    Pi_row_regime1 <- ms_fit@Coef[1, 1:4]  # Regime 1 coefficients
    Pi_row_regime2 <- ms_fit@Coef[2, 1:4]  # Regime 2 coefficients
    
    # Extract intercept row vector for each regime
    int_row_regime1 <- ms_fit@Coef[1, 1]  # Regime 1 coefficients
    int_row_regime2 <- ms_fit@Coef[2, 1]  # Regime 2 coefficients
    
    # Extract error row vector for each regime
    err_row_regime1 <- ms_model@Fit@error[,1]
    err_row_regime2 <- ms_model@Fit@error[,2]
    
    # Append to the Gamma list
    Gamma_list$Gamma_regime1[, , 1][i, ] <- matrix(as.numeric(Gamma_row_regime1[1:4]))
    Gamma_list$Gamma_regime1[, , 2][i, ] <- matrix(as.numeric(Gamma_row_regime1[5:8]))
    
    Gamma_list$Gamma_regime2[, , 1][i, ] <- matrix(as.numeric(Gamma_row_regime2[1:4]))
    Gamma_list$Gamma_regime2[, , 2][i, ] <- matrix(as.numeric(Gamma_row_regime2[5:8]))
    
    # Append to the Pi list
    Pi_list$Pi_regime1[i, ] <- as.numeric(Pi_row_regime1[1:4])
    Pi_list$Pi_regime2[i, ] <- as.numeric(Pi_row_regime2[1:4])
    
    # Append to the intercept list
    int_list$int_regime1[i, 1] <- as.numeric(int_row_regime1)
    int_list$int_regime2[i, 1] <- as.numeric(int_row_regime2)
    
    # Append to the error list
    error_list$err_reg1[, i] <- as.numeric(err_row_regime1)  # Add error for regime 1
    error_list$err_reg2[, i] <- as.numeric(err_row_regime2)  # Add error for regime 2
  }
  
  # Return both lists as a single list
  return(list(Pi_list = Pi_list, Gamma_list = Gamma_list, int_list = int_list, 
              error_list = error_list))
}










