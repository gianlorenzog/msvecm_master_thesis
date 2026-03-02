plot_prob <- function(y, prob = c("filtered", "smoothed")) {
  # Match the argument to ensure valid input
  prob <- match.arg(prob)
  
  # Define colors based on probability type
  state_colors <- if (prob == "filtered") c("red", "green") else c("blue", "purple")
  
  # Ensure the input `y` is a matrix or data frame of probabilities
  if (!is.matrix(y) && !is.data.frame(y)) {
    stop("Input 'y' must be a matrix or data frame of probabilities.")
  }
  
  # Set up plot layout
  par(mfrow = c(ncol(y), 1))  # Subplots for each state
  for (s in 1:ncol(y)) {
    plot(y[, s], type = "l", col = state_colors[s %% length(state_colors) + 1], lwd = 2,
         ylab = paste(prob, "probabilities"),  # Dynamic y-axis label
         xlab = "Time", main = paste("State", s, "Probability"))
    abline(h = 0.5, lty = 2, col = "gray")  # Reference line for regime threshold
  }
  par(mfrow = c(1, 1))  # Reset layout
}

