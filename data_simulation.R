library(MASS)  # Per la simulazione multivariata
library(xts)   # Per oggetti time series estesi

# Parametri
set.seed(123)
n_obs <- 500
n_regimes <- 2

# Date fittizie (giornaliere, solo per demo – puoi usare mensili ecc.)
start_date <- as.Date("2000-01-01")
dates <- seq(start_date, by = "days", length.out = n_obs)

# Parametri dei regimi
mu_list <- list(
  c(0.002, 0.001, 0.001, 0.0005),   # Regime 1
  c(-0.001, -0.0005, 0.002, 0.001)  # Regime 2
)

sigma_list <- list(
  matrix(c(0.0005, 0, 0, 0,
           0, 0.0001, 0, 0,
           0, 0, 0.0002, 0,
           0, 0, 0, 0.0001), 4),
  
  matrix(c(0.001, 0, 0, 0,
           0, 0.0002, 0, 0,
           0, 0, 0.0004, 0,
           0, 0, 0, 0.0002), 4)
)

# Matrice di transizione markoviana
P <- matrix(c(0.95, 0.05,
              0.10, 0.90), nrow=2, byrow=TRUE)

# Funzione per simulare stati markoviani
simulate_states <- function(P, n) {
  states <- numeric(n)
  states[1] <- sample(1:nrow(P), 1)
  for (t in 2:n) {
    states[t] <- sample(1:nrow(P), 1, prob = P[states[t-1], ])
  }
  return(states)
}

# Simula gli stati latenti
states <- simulate_states(P, n_obs)

# Simulazione dati condizionata agli stati
Y <- matrix(NA, nrow = n_obs, ncol = 4)
colnames(Y) <- c("StockIndex", "GDP", "CPI", "M1")

for (t in 1:n_obs) {
  mu <- mu_list[[states[t]]]
  sigma <- sigma_list[[states[t]]]
  Y[t, ] <- mvrnorm(1, mu = mu, Sigma = sigma)
}


######### FINE

# Trasforma in livelli cumulando i rendimenti
# Y_levels <- apply(Y, 2, cumsum)

# Crea oggetto xts
# Y_xts <- xts(Y_levels, order.by = dates)
# colnames(Y_xts) <- c("StockIndex", "GDP", "CPI", "M1")

# Plot delle serie simulate
plot.zoo(Y_xts, screens = 1, col = 1:4, lty = 1,
         xlab = "Time", ylab = "Simulated Levels",
         main = "Serie Storiche Simulate con Markov Switching")
legend("topleft", legend = colnames(Y_xts), col = 1:4, lty = 1)
