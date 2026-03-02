#
# Function to plot smoothed and filtered probabilities (Aggregate prob) with regime splitting (2 graphs)
#
plot_ag_prob <- function(smoothed_probs_input, filtered_probs_input, time = NULL) {
  # Create local copies to avoid modifying the original objects
  smoothed_probs <- smoothed_probs_input
  filtered_probs <- filtered_probs_input
  
  # Ensure smoothed_probs and filtered_probs have the same number of rows
  if (nrow(smoothed_probs) != nrow(filtered_probs)) {
    if (nrow(smoothed_probs) > nrow(filtered_probs)) {
      smoothed_probs <- smoothed_probs[2:nrow(smoothed_probs), ]  # Drop the first row of smoothed_probs
    } else {
      filtered_probs <- filtered_probs[2:nrow(filtered_probs), ]  # Drop the first row of filtered_probs
    }
  }
  
  # Initialize 'time' inside the function
  if (is.null(time) || length(time) != nrow(smoothed_probs)) {
    time <- 1:nrow(smoothed_probs)  # Set 'time' as a sequence matching the rows
  }
  
  # Ensure time length matches the number of rows
  if (length(time) != nrow(smoothed_probs)) {
    stop("The length of 'time' must match the number of rows in smoothed_probs.")
  }
  
  # Set up plot layout for two regimes
  par(mfrow = c(2, 1), mar = c(5, 4, 4, 4) + 0.1)  # Two rows for two states with extra space for dual axes
  
  for (state in 1:2) {
    # Plot smoothed probabilities on primary y-axis (left)
    plot(
      time, smoothed_probs[, state], type = "l", col = "red", lwd = 2,
      ylim = c(0, 1), ylab = "Smoothed Probabilities", xlab = "Time",
      main = paste("Regime", state)
    )
    
    # Add filtered probabilities on secondary y-axis (right)
    par(new = TRUE)
    plot(
      time, filtered_probs[, state], type = "h", col = "grey", lwd = 1,
      ylim = c(0, 1), axes = FALSE, xlab = "", ylab = ""
    )
    axis(4)  # Add right-side axis for filtered probabilities
    mtext("Filtered Probabilities", side = 4, line = 3)
  }
  
  # Reset layout
  par(mfrow = c(1, 1))
}

