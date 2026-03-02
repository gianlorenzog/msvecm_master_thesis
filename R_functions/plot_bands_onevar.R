#---- Plot a single variable with separate regime bands
plot_bands_onevar <- function(variable, state_probs, time = NULL, threshold = 0.7, main = "Significant Regime Bands") {
  # Ensure lengths match by trimming the variable to match state_probs
  variable <- variable[-(1:(length(variable) - nrow(state_probs)))]  # Trim the variable to match state_probs length
  
  # If time is not provided, create a default time sequence
  if (is.null(time)) {
    time <- 1:nrow(state_probs)
  } else {
    # Trim time to match state_probs
    time <- time[-(1:(length(time) - nrow(state_probs)))]
  }
  
  # Set up plot layout for two regimes
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))  # Two rows for two states
  
  # Regime colors
  regime_colors <- c("lightgreen", "lightgray")
  
  for (state in 1:2) {
    # Identify time indices where the state probability is above the threshold
    significant_indices <- which(state_probs[, state] > threshold)
    
    # Split into consecutive segments for each significant regime period
    segments <- split(significant_indices, cumsum(c(1, diff(significant_indices) != 1)))
    
    # Create plot for the variable with significant bands for this regime
    plot(time, variable, type = "l", col = "black", lwd = 2, 
         ylab = "SX5E", xlab = "Time", main = paste("Regime", state))
    
    # Draw significant bands
    for (segment in segments) {
      if (length(segment) > 1) {  # Only draw if there are multiple consecutive points
        polygon(
          x = c(time[segment], rev(time[segment])),
          y = c(rep(min(variable, na.rm = TRUE), length(segment)),
                rev(rep(max(variable, na.rm = TRUE), length(segment)))),
          col = regime_colors[state], border = NA
        )
      }
    }
    
    # Redraw the variable line on top of the bands
    lines(time, variable, col = "black", lwd = 2)
  }
  
  # Reset layout
  par(mfrow = c(1, 1))
}
